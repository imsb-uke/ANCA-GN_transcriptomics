import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
from matplotlib.colorbar import Colorbar
import matplotlib.colors as mcolors
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import seaborn as sns
import scanpy as sc
from scipy.sparse import csr_matrix
import scvi
import celltypist
from celltypist import models
import muon as mu
import collections
from tqdm import tqdm
import matplotlib.cm as cm

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


def scale(adata, layer="log1p"):
    adata.X = adata.layers[layer].copy()
    sc.pp.scale(adata)
    adata.layers["scaled"] = adata.X.copy()


def log_normalize(adata, target_sum=1e4):
    # normalize counts
    restore_raw_counts(adata)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    # log transform
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()


def protein_clr(adata):
    restore_raw_counts(adata)
    mu.prot.pp.clr(adata, inplace=True)
    adata.layers["clr"] = adata.X.copy()


def filter_genes_and_cells(adata, min_cells=10, min_genes=100):
    if min_genes is not None:
        print(f"Total number of cells: {adata.n_obs}")
        sc.pp.filter_cells(adata, min_genes=min_genes)
        print(f"Number of cells after gene filter: {adata.n_obs}")
    if min_cells is not None:
        print(f"Total number of genes: {adata.n_vars}")
        # Min 20 cells - filters out 0 count genes
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(f"Number of genes after cell filter: {adata.n_vars}")


def restore_raw_counts(adata):
    try:
        adata.X = adata.layers["counts"].copy()
    except KeyError:
        print("No raw counts layer found, assuming counts are already raw")
        adata.layers["counts"] = adata.X.copy()


def scale(adata, **kwargs):
    sc.pp.scale(adata, **kwargs)
    adata.layers["scaled"] = adata.X


# def preprocess(adata, min_cells=20):
#     filter_genes_and_cells(adata, min_cells=min_cells)
#     log_normalize(adata)
#     sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="patient")
# scale(adata)


def adapt_marker_genes_to_data(adata, marker_genes):
    if type(marker_genes) == list:
        return find_duplicates(adata.var_names, marker_genes)

    marker_genes_in_data = dict()
    for ct, markers in marker_genes.items():
        markers_found = find_duplicates(adata.var_names, markers)
        if len(markers_found) > 0:
            marker_genes_in_data[ct] = markers_found
    return marker_genes_in_data


def find_duplicates(list1, list2):
    return list(set(list1).intersection(list2))


def convert_to_sparse(X):
    X = sc._utils.as_sparse(X, copy=True)
    return X


def load_yaml(file_path):
    with open(file_path, "r") as f:
        data = yaml.safe_load(f)
    return data


def get_adata_stats(adata, sample_col="sample"):
    # check number of cells per sample and number of expressed genes per cell and number cell expressing a gene
    if sample_col is not None:
        min_cells_per_sample = adata.obs[sample_col].value_counts().min()
        max_cells_per_sample = adata.obs[sample_col].value_counts().max()
    else:
        min_cells_per_sample = None
        max_cells_per_sample = None
    # number expressed genes
    min_n_expressed_genes_per_cell = (adata.X > 0).sum(axis=1).min()
    max_n_expressed_genes_per_cell = (adata.X > 0).sum(axis=1).max()
    # number of cells expressing a gene
    min_n_cells_expressing_gene = (adata.X > 0).sum(axis=0).min()
    output = (
        f"min cells per sample: {min_cells_per_sample},\n"
        f"max cells per sample: {max_cells_per_sample},\n"
        f"min n expressed genes per cell: {min_n_expressed_genes_per_cell},\n"
        f"max n expressed genes per cell: {max_n_expressed_genes_per_cell},\n"
        f"min n cells expressing gene: {min_n_cells_expressing_gene}\n"
    )
    print(output)


def integrate_with_totalvi(
    rna_data,
    protein_data,
    batch_key="sample",
    n_top_genes=2000,
    categorical_covariate_keys=None,
    continuous_covariate_keys=None,
    save_path=None,
    prefix=None,
    empirical_protein_background_prior=False,
    **train_kwargs,
):
    # make sure we dont change the original data
    rna_data = rna_data.copy()
    protein_data = protein_data.copy()

    print("Preparing RNA data")
    rna_data.X = rna_data.layers["counts"].copy()
    log_normalize(rna_data)

    print("Selecting highly variable genes")
    sc.pp.highly_variable_genes(
        rna_data,
        batch_key=batch_key,
        # flavor="seurat_v3",
        flavor="seurat",
        # flavor="cell_ranger",
        n_top_genes=n_top_genes,
        # span=0.75,
        subset=False,
        # layer="counts",
        layer="log1p",
    )
    rna_subset = rna_data[:, rna_data.var["highly_variable"]]
    print("Stats after filtering highly variable genes")
    get_adata_stats(rna_subset, sample_col=batch_key)

    rna_subset.layers["counts"] = csr_matrix(rna_subset.layers["counts"]).toarray()

    rna_subset.obsm["protein_counts"] = pd.DataFrame(
        protein_data.layers["counts"].toarray(),
        index=rna_subset.obs.index,
        columns=protein_data.var_names,
    )
    # check number of cell without protein expression
    n_cells_without_protein = (
        (rna_subset.obsm["protein_counts"]).sum(axis=1) == 0
    ).sum()
    print(f"Number of cells without protein expression: {n_cells_without_protein}")

    print("Setting up anndata object for TOTALVI")
    scvi.model.TOTALVI.setup_anndata(
        rna_subset,
        batch_key=batch_key,
        protein_expression_obsm_key="protein_counts",
        layer="counts",
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
    )

    print("Creating model")
    model = scvi.model.TOTALVI(
        rna_subset,
        latent_distribution="normal",
        n_layers_decoder=2,
        empirical_protein_background_prior=empirical_protein_background_prior,
    )

    print("Starting training")
    model.train(**train_kwargs)

    # save the trained model
    if save_path is not None:
        if prefix is None:
            prefix = ""
        model.save(dir_path=save_path, save_anndata=True, overwrite=True, prefix=prefix)
    return model, rna_subset


def run_celltypist(rna_data, variant=None):
    models.download_models(
        force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"]
    )
    model_low = models.Model.load(model="Immune_All_Low.pkl")
    model_high = models.Model.load(model="Immune_All_High.pkl")
    predictions_high = celltypist.annotate(
        rna_data, model=model_high, majority_voting=True
    )
    predictions_high_adata = predictions_high.to_adata()
    rna_data.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
        rna_data.obs.index, "majority_voting"
    ]
    rna_data.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
        rna_data.obs.index, "conf_score"
    ]
    predictions_low = celltypist.annotate(
        rna_data, model=model_low, majority_voting=True
    )
    predictions_low_adata = predictions_low.to_adata()
    rna_data.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
        rna_data.obs.index, "majority_voting"
    ]
    rna_data.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
        rna_data.obs.index, "conf_score"
    ]


def get_marker_results_for_marker_list(marker_results, marker_list):
    results_subset = marker_results[marker_results.names.isin(marker_list)]
    results_subset = results_subset.sort_values(by="names")
    return results_subset


def find_overlap(list_1, list_2):
    return list(set(list_1).intersection(list_2))


def check_markers_in_rna_and_protein(rna_data, protein_data, markers: list):
    rna_markers = find_overlap(rna_data.var_names, markers)
    # print out all markers that are not in the rna data
    print(f"Markers in RNA data: {rna_markers}")
    print(f"Markers not in RNA data: {set(markers) - set(rna_markers)}")
    protein_markers = find_overlap(protein_data.var_names, markers)
    # print out all markers that are not in the protein data
    print(f"Markers in protein data: {protein_markers}")
    print(f"Markers not in protein data: {set(markers) - set(protein_markers)}")
    return rna_markers, protein_markers


def get_genes_with_highest_log_fold_changes(de_results, n_genes=10, clusters=None):
    if clusters is None:
        clusters = de_results["group"].unique()
    result_dict = {}
    try:
        clusters = sorted(clusters, key=lambda x: int(x))
    except:
        clusters = sorted(clusters)
    for cluster in clusters:
        subdf = de_results[de_results.group == cluster].sort_values(
            by="logfoldchanges", ascending=False
        )
        relevant_genes = subdf.head(n_genes)
        print(f"{cluster}:", ", ".join(relevant_genes.names.to_list()))
        result_dict[cluster] = relevant_genes.names.to_list()
    return result_dict


def extract_clusters_per_gene(marker_results):
    pass


def run_de_analysis(adata, labeling, data_type="rna", n_genes=10, min_expression=0.1):
    adata.uns["log1p"] = {}
    adata.uns["log1p"]["base"] = None
    if data_type == "protein" or data_type == "cite":
        layer = "clr"
    else:
        layer = [key for key in adata.layers.keys() if "log" in key][0]
    print(f"Running DE analysis for {labeling}")
    sc.tl.rank_genes_groups(
        adata,
        groupby=labeling,
        method="wilcoxon",
        key_added=f"dea_{labeling}",
        layer=layer,
        use_raw=False,
        pts=True,
    )

    all_marker_results = sc.get.rank_genes_groups_df(adata, None, key=f"dea_{labeling}")
    filtered_marker_results, best_markers = filter_de_markers(
        adata,
        labeling,
        all_marker_results,
        min_expression=min_expression,
        n_genes=n_genes,
    )
    return all_marker_results, filtered_marker_results, best_markers


def filter_de_markers(
    adata, cluster_key, all_marker_results, min_expression=0.1, n_genes=10
):
    clusters = adata.obs[cluster_key].unique().tolist()
    filtered_marker_results = all_marker_results[
        (all_marker_results["pvals_adj"] < 0.05)
        & (all_marker_results["logfoldchanges"] > 0.5)
        & (all_marker_results["pct_nz_group"] > min_expression)
    ]
    best_markers = get_genes_with_highest_log_fold_changes(
        filtered_marker_results, n_genes=n_genes, clusters=clusters
    )
    try:
        best_markers = collections.OrderedDict(
            sorted(best_markers.items(), key=lambda x: int(x[0]))
        )
        best_markers = collections.OrderedDict(
            sorted(best_markers.items(), key=lambda x: int(x[0]))
        )
    except:
        # best_markers = collections.OrderedDict(sorted(best_markers.items()))
        best_markers = collections.OrderedDict(best_markers.items())
    return filtered_marker_results, best_markers


def get_top_n_markers(top_markers, n=5):
    top_n_markers = collections.OrderedDict(
        (key, val[0:n]) for key, val in top_markers.items() if len(val) > 0
    )
    return top_n_markers


def delete_leiden_clustering_results(adata, only_uns=False):
    if not only_uns:
        leiden_cols = [col for col in adata.obs.columns if "leiden" in col]
        adata.obs = adata.obs.drop(columns=leiden_cols)

    leiden_keys = [key for key in adata.uns.keys() if "leiden" in key]
    for key in leiden_keys:
        del adata.uns[key]

    return adata


def merge_annotations(
    adata, base_annot_col, extra_annot_col, new_col_name="merged_celltypes"
):
    # use all string values that cannot be converted to float or integer in the extra col
    # use these annotations instead of the base annotations
    adata.obs[new_col_name] = adata.obs[base_annot_col]
    adata.obs.loc[
        adata.obs[extra_annot_col].apply(lambda x: not is_number(x)), new_col_name
    ] = adata.obs[extra_annot_col]
    return adata


def is_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def subset_and_cluster(
    adata,
    cluster_col,
    cluster_name,
    representation="totalvi",
    resolutions=[0.5, 0.6, 0.75, 1.0],
    show=True,
):
    adata_sub = adata[adata.obs[cluster_col] == cluster_name, :].copy()
    sc.pp.neighbors(adata_sub, use_rep=f"X_{representation}")
    sc.tl.umap(adata_sub)
    sc.tl.leiden(adata_sub)
    delete_leiden_clustering_results(adata_sub)
    for resolution in resolutions:
        sc.tl.leiden(adata_sub, key_added=f"leiden_{resolution}", resolution=resolution)

    if show:
        colors = [
            key
            for key in adata_sub.obs.keys()
            if "leiden" in key and representation in key
        ]
        fig = sc.pl.umap(
            adata_sub,
            color=colors,
            legend_loc="on data",
            # save=f"_leiden_{representation}.png",
            return_fig=True,
        )
    return adata_sub


def get_cluster_names(adata, cluster_col):
    try:
        for idx in sorted(adata.obs[cluster_col].unique(), key=lambda x: int(x)):
            print(f'"{idx}": "",')
    except:
        for idx in sorted(adata.obs[cluster_col].unique()):
            print(f'"{idx}": "",')


def rank_genes_with_threshold(adata, de_results, cluster_key, threshold=0.2):
    # order genes by logfold changes per cluster

    clusters = adata.obs[cluster_key].unique().tolist()
    de_clusters = de_results["group"].unique().tolist()
    assert set(de_clusters).issubset(
        set(clusters)
    ), "Clusters in adata and de_results do not match"

    result_dict = {}
    for cluster in tqdm(de_clusters):
        subdf = de_results[de_results.group == cluster].sort_values(
            by="logfoldchanges", ascending=False
        )
        # calculate fraction of cells that express each gene in the cluster
        sub_adata = adata[adata.obs[cluster_key] == cluster, :]
        genes = subdf.names.to_list()
        logfolds = subdf.logfoldchanges.to_list()
        n_expressed = np.array((sub_adata[:, genes].X > 0).sum(axis=0)).squeeze()
        cells_in_cluster = sub_adata.shape[0]
        fraction_expressed = n_expressed / cells_in_cluster

        # filter genes by fraction of cells that express them
        threshold_mask = fraction_expressed >= threshold

        result_dict[cluster] = subdf.loc[threshold_mask, :].names.to_list()
    return result_dict


def run_de_pipeline(adata, cluster_key, mod="rna", top_n=5, dotplot=False, **kwargs):
    (
        all_marker_results,
        filtered_marker_results,
        best_markers,
    ) = run_de_analysis(adata, labeling=cluster_key, data_type=mod, **kwargs)
    top_markers = get_top_n_markers(best_markers, n=top_n)
    if dotplot:
        sc.pl.dotplot(
            adata,
            var_names=top_markers,
            groupby=cluster_key,
            standard_scale="var",
            # color_map="Reds",
            # swap_axes=True,
            var_group_rotation=25,
        )
    return all_marker_results, filtered_marker_results, best_markers


def dotplot_rank_genes_threshold(
    adata,
    filtered_marker_results,
    cluster_key,
    threshold=0.2,
    top_n=5,
    figsize=None,
    matrixplot=False,
    as_list=False,
    rot=45,
    invert_cmap=False,
    colors=None,
    swap_axes=True,
    **kwargs,
):
    if "cmap" in kwargs:
        cmap = kwargs["cmap"]
        if invert_cmap:
            inverted_cmap = plt.cm.get_cmap(cmap)
            inverted_cmap = inverted_cmap.reversed()
            kwargs["cmap"] = inverted_cmap
            cmap = inverted_cmap
    else:
        cmap = "viridis"
    best_markers_filtered = rank_genes_with_threshold(
        adata, filtered_marker_results, cluster_key, threshold=threshold
    )
    top_markers = get_top_n_markers(best_markers_filtered, n=top_n)
    if as_list:
        top_markers = [item for sublist in top_markers.values() for item in sublist]
    if matrixplot:
        fig = sc.pl.matrixplot(
            adata,
            var_names=top_markers,
            groupby=cluster_key,
            standard_scale="var",
            # cmap="Reds",
            var_group_rotation=0,
            swap_axes=swap_axes,
            figsize=figsize,
            return_fig=True,
            **kwargs,
        )
        categories = adata.obs[cluster_key].cat.categories
        ax = fig.get_axes()["mainplot_ax"]
        ax.set_xticklabels(
            labels=ax.get_xticklabels(),
            rotation=rot,
            ha="right" if swap_axes else "center",
        )
        # fig.get_axes()["color_legend_ax"].remove()
        cax = fig.get_axes()["color_legend_ax"]
        x0, y0, width, height = cax.get_position().bounds
        cax.set_position(
            [x0 - 0.025 if swap_axes else x0 - 0.015, y0, width / 3.2, height * 2]
        )
        # arange title to the right
        cax.set_title("")
        cax.set_title("Expression", pad=7.5, fontsize=10, loc="left")
        cmap = mpl.cm.get_cmap(cmap)
        norm = mpl.colors.Normalize(vmin=0, vmax=1)

        cax = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cax,
            orientation="vertical",
        )
        # remove ticks
        cax.ax.tick_params(axis="y", which="both", length=0, labelsize=10)
        cax.ax.set_yticks([0.03, 0.5, 0.97])
        cax.ax.set_yticklabels(["0.0", "0.5", "1.0"])

        cmap = mcolors.ListedColormap(colors[0 : len(categories)])
        bounds = np.arange(len(categories) + 1) - 0.5
        norm = mcolors.BoundaryNorm(bounds, cmap.N)
        # Create a divider for the existing axes
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top" if swap_axes else "left", size="1.5%", pad=0.1)
        # Create the categorical colorbar with horizontal orientation
        cbar = plt.colorbar(
            cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cax,
            # ax=ax,
            ticks=np.arange(len(categories)),
            orientation="horizontal" if swap_axes else "vertical",
        )
        if swap_axes:
            cbar.ax.xaxis.set_ticks_position("top")
            cbar.ax.set_xticklabels(categories, rotation=45, ha="left")
            cbar.ax.tick_params(axis="x", which="both", length=0)
            ax.set_xticklabels([])
        else:
            cbar.ax.yaxis.set_ticks_position("left")
            cbar.ax.set_yticklabels(categories[::-1], rotation=0, ha="right")
            cbar.ax.tick_params(axis="y", which="both", length=0)
            ax.set_yticklabels([])

    else:
        sc.pl.dotplot(
            adata,
            var_names=top_markers,
            groupby=cluster_key,
            standard_scale="var",
            # color_map="Reds",
            # swap_axes=True,
            var_group_rotation=25,
            figsize=figsize,
            **kwargs,
        )
    plt.tight_layout()
    return best_markers_filtered


def get_colorbars(fig):
    def check_kids(obj, bars):
        for child in obj.get_children():
            if isinstance(getattr(child, "colorbar", None), Colorbar):
                bars.append(child.colorbar)
            check_kids(child, bars)
        return bars

    return check_kids(fig, [])


def plot_qc(
    adata,
    cell_type_key="cell_type_yu",
    tissue_key="tissue",
    figsize=(13, 9),
    wspace=0.12,
    hspace=0.1,
    save_path=None,
    log_scale=True,
):
    tmp = adata.X.copy()
    adata.X = adata.layers["counts"].copy()
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)

    adata.var["ribo"] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    adata.X = tmp

    # tissue composition
    tissue_composition_per_patient = adata.obs.groupby(["patient"])[
        tissue_key
    ].value_counts(normalize=True)
    tissue_composition_per_patient = pd.DataFrame(tissue_composition_per_patient)
    tissue_composition_per_patient.columns = ["fraction"]
    tissue_composition_per_patient = tissue_composition_per_patient.reset_index()
    # tissue_composition_per_patient

    # tissue composition
    tissue_composition_per_cluster = adata.obs.groupby([cell_type_key])[
        tissue_key
    ].value_counts(normalize=True)
    tissue_composition_per_cluster = pd.DataFrame(tissue_composition_per_cluster)
    tissue_composition_per_cluster.columns = ["fraction"]
    tissue_composition_per_cluster = tissue_composition_per_cluster.reset_index()

    # combined plot
    obs_list = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
    obs_names = ["Number of\ngenes", "Total counts", "% MT\n", "% ribo\n"]

    with sns.axes_style("ticks"):
        # fig, axs = plt.subplots(5, 2, figsize=(13, 8), sharex=False, sharey=False)
        num_patients = len(adata.obs["patient"].unique())
        num_clusters = len(adata.obs[cell_type_key].unique())
        gs = gridspec.GridSpec(
            6,
            2,
            hspace=hspace,
            wspace=wspace,
            width_ratios=[num_patients / num_clusters, 1],
        )  # Change these ratios as needed
        fig = plt.figure(figsize=figsize)
        axs = np.array(
            [[plt.subplot(gs[i, 0]), plt.subplot(gs[i, 1])] for i in range(6)]
        )
        levels = ["patient", cell_type_key]
        for j, level in enumerate(levels):
            if level == "patient":
                tissue_comp = tissue_composition_per_patient
                xlabel = "Patient"
            else:
                tissue_comp = tissue_composition_per_cluster
                xlabel = "Cell type"
            for k, ax in enumerate(axs[:, j]):
                if k in [0, 1, 2, 3]:
                    obs = obs_list[k]
                    sns.violinplot(
                        data=adata.obs,
                        x=level,
                        y=obs,
                        ax=ax,
                        palette="tab20",
                        inner=None,
                        rasterized=True,
                    )
                    sns.stripplot(
                        data=adata.obs,
                        x=level,
                        y=obs,
                        ax=ax,
                        color="black",
                        size=0.1,
                        jitter=0.1,
                        rasterized=True,
                    )
                    ax.set_ylabel(obs_names[k])
                if k in [4]:
                    # plot tissue composition
                    blood_comp = tissue_comp.loc[tissue_comp[tissue_key] == "B", :]
                    kidney_comp = tissue_comp.loc[tissue_comp[tissue_key] == "K", :]
                    colors = ["#65D697", "#FF7163"]
                    sns.barplot(
                        data=blood_comp,
                        x=level,
                        y="fraction",
                        ax=ax,
                        color=colors[0],
                        width=0.9,
                        rasterized=True,
                    )
                    sns.barplot(
                        data=kidney_comp,
                        x=level,
                        y="fraction",
                        ax=ax,
                        color=colors[1],
                        width=0.9,
                        bottom=blood_comp["fraction"],
                        rasterized=True,
                    )
                    ax.set_yticks([0, 0.5, 1])
                    ax.set_yticklabels(["0", "0.5", "1.0"])
                    ax.set_ylabel("Fraction of\ncells")
                if k in [5]:
                    data = pd.DataFrame(adata.obs.groupby(level).size()).reset_index()
                    sns.barplot(
                        data=data,
                        x=level,
                        y=0,
                        ax=ax,
                        color="grey",
                        width=0.9,
                        rasterized=True,
                    )
                    if log_scale:
                        ax.set_yscale("log")
                    ax.set_ylabel("Number of\ncells")
                ax.set_xlabel("")

                if k < axs.shape[0] - 1:
                    ax.set_xticklabels([])
                    ax.set_xticks([])
                    for tick in ax.get_xticklines():
                        tick.set_visible(False)
                else:
                    # rotate x labels
                    for tick in ax.get_xticklabels():
                        tick.set_rotation(90)
                        tick.set_horizontalalignment("center")
                    ax.set_xlabel(xlabel)
                sns.despine(fig, ax, top=True, right=True, left=False, bottom=False)
                if j == 1:
                    ax.set_ylabel("")
        plt.subplots_adjust(hspace=0.1)
        if save_path is not None:
            plt.savefig(
                save_path,
                bbox_inches="tight",
                dpi=400,
                transparent=True,
            )
        plt.show()
