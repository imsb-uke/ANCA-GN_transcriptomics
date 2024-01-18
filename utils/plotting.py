import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import gridspec
import math
import pprint
import matplotlib.colors as mcolors
from source.utils import find_overlap
import matplotlib as mpl
from tqdm import tqdm
from adjustText import adjust_text


def plot_celltype_markers(
    data,
    marker_genes,
    width=3,
    title=None,
    vmax=4,
    name="markers",
    width_fac=3.0,
    ncols=None,
):
    if ncols is not None:
        width = ncols
    if isinstance(marker_genes, list):
        if len(marker_genes) > width:
            plot_new_line = True
        else:
            plot_new_line = False
        marker_genes = {name: marker_genes}
    height = len(marker_genes) * 3
    for ct, markers in marker_genes.items():
        if len(marker_genes) == 1:
            num_markers = len(markers)
            actual_width = min(width, num_markers)
        else:
            actual_width = width
    width_scale = actual_width * width_fac
    fig, axs = plt.subplots(
        len(marker_genes),
        actual_width,
        figsize=(width_scale, height),
        sharex=True,
        sharey=True,
    )
    # if len(marker_genes) == 1:
    #     axs = [axs]
    # else:
    try:
        axs = axs.ravel()
    except:
        axs = [axs]
    axs_reduced = []
    k = 0
    for ct, markers in marker_genes.items():
        for i in range(actual_width):
            try:
                marker = markers[i]
                g = sc.pl.umap(
                    data,
                    color=marker,
                    ax=axs[k],
                    vmin=0,
                    vmax=vmax,
                    sort_order=False,
                    frameon=True,
                    use_raw=False,
                    cmap="Reds",
                    show=False,
                )
                g.set(xticklabels=[])
                g.set(xlabel=None)
                g.set(yticklabels=[])
                g.set(ylabel=None)
                if i == 0:
                    g.set(ylabel=ct)
                axs_reduced.append(axs[k])
            except Exception as e:
                print(e)
                # pprint.pprint(axs[k].properties())
                # axs[k].remove()
                fig.delaxes(axs[k])

            k += 1
    if title is not None:
        plt.suptitle(title, y=0.9)
    fig.tight_layout()
    plt.show()
    if plot_new_line:
        plot_celltype_markers(
            data,
            marker_genes[name][width::],
            width=width,
            vmax=vmax,
            name=name,
            width_fac=width_fac,
        )


# def plot_expression_per_cluster(data, cluster_key, marker_genes, width=3):
#     height = len(marker_genes) * 3
#     width_scale = width * 3.5
#     fig, axs = plt.subplots(len(marker_genes), width, figsize=(width_scale, height), sharex=True, sharey=True)
#     axs = axs.ravel()
#     k = 0
#     for gene in marker_genes:

#         for i in range(width):
#             try:
#                 marker = markers[i]
#                 g = sc.pl.violin(
#                     data,
#                     keys=marker_genes,
#                     groupby=cluster_key,
#                     layer="log_normalized",
#                     use_raw=False,
#                     show=False,
#                     x_label="Leiden",
#                 )
#                 # g.set(xticklabels=[])
#                 # g.set(xlabel=None)
#                 # g.set(yticklabels=[])
#                 # g.set(ylabel=None)
#                 if i == 0:
#                     g.set(ylabel=ct)
#             except:
#                 # delete ax
#                 fig.delaxes(axs[k])
#             k += 1


def plot_specific_cluster(adata, cluster_key, cluster_name):
    adata = adata.copy()
    # adata.obs.loc[adata.obs[cluster_key] == cluster_name, cluster_key] = 1
    # adata.obs.loc[adata.obs[cluster_key] != cluster_name, cluster_key] = np.nan
    # adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
    adata = adata[adata.obs[cluster_key] == cluster_name]
    sc.pl.umap(adata, color=cluster_key, legend_loc="on data")
    plt.show()


def plot_marker_dotplot(adata, markers, title="RNA", groupby="cell_type"):
    fig = sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby=groupby,
        standard_scale="var",
        var_group_rotation=45,
        swap_axes=True,
        # return_fig=True,
        show=False,
    )
    ax = fig["mainplot_ax"]
    # rotate x labels
    labels = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_title(title)
    plt.tight_layout()


def plot_leiden_results(adata, rep_name="X_totalvi", colors=None):
    if colors is None:
        colors = [
            key for key in adata.obs.keys() if "leiden" in key and rep_name in key
        ]
    fig = sc.pl.umap(
        adata,
        color=colors,
        legend_loc="on data",
        return_fig=True,
    )
    plt.tight_layout()
    plt.show()


def visualize_umap(
    adata,
    color,
    title=None,
    legend_loc="on data",
    save=None,
    adjust_legend=False,
    adjust_cbar=False,
    legend_kwargs={"ncols": 2, "bbox_to_anchor": (1.175, 0.75), "fontsize": 5.5},
    return_fig=False,
    ax=None,
    figsize=(6, 5),
    **kwargs,
):
    if not isinstance(color, list):
        color = [color]
    # show celltype labels
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    plot = sc.pl.umap(
        adata,
        color=color,
        legend_loc=legend_loc,
        outline_width=(0.2, 0.05),
        # legend_fontweight="normal",
        ax=ax,
        show=False,
        **kwargs,
    )
    if adjust_legend:
        legend = ax.get_legend()
        handles = legend.legendHandles
        labels = [t.get_text() for t in legend.get_texts()]
        # remove legend
        legend.remove()
        # add new legend
        ax.legend(
            handles,
            labels,
            loc="upper right",
            frameon=False,
            **legend_kwargs,
        )
    if adjust_cbar:
        cax = fig.get_axes()[1]
        x0, y0, width, height = cax.get_position().bounds
        height_scaler = 0.5
        cax.set_position(
            [x0, y0 + 0.5 * height_scaler * height, width, height * height_scaler]
        )

    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    if title is not None:
        ax.set_title(title, fontsize=10)
    if return_fig and ax is not None:
        return fig
    else:
        plt.show()
        return None


def marker_dotplot(adata, var_names, title="RNA", groupby="cell_type", suptitle=None):
    fig = sc.pl.dotplot(
        adata,
        var_names=var_names,
        groupby=groupby,
        # standard_scale="var",
        var_group_rotation=45,
        swap_axes=True,
        # return_fig=True,
        vmax=3,
        show=False,
    )
    ax = fig["mainplot_ax"]
    ax.set_title(title)
    # rotate x labels
    labels = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    if suptitle is not None:
        plt.suptitle(suptitle)
    plt.tight_layout()
    plt.show()


class MarkerPlotter:
    def __init__(
        self,
        rna_data,
        protein_data=None,
        marker_db: dict = None,
        cluster_key: str = None,
    ) -> None:
        self.rna_data = rna_data
        self.protein_data = protein_data
        # expected as dict with one entry per cell type
        # later allow to specify positive and negative markers
        self.marker_db = marker_db
        self.cluster_key = cluster_key

    def adapt_markers(self, adata, markers: list):
        # adapt markers based on the data
        adapted_markers = [marker for marker in markers if marker in adata.var_names]
        removed_markers = list(set(markers) - set(adapted_markers))
        # if len(removed_markers) > 0:
        #     print(f"Removed {len(removed_markers)} markers: {removed_markers}")
        return adapted_markers

    def visualize_markers(
        self,
        markers=None,
        dtype="rna",
        cluster_key=None,
        cell_type=None,
        use_default_plot=False,
        marker_plot=True,
        dotplot=True,
        title_with_dtype=True,
        vmax="p99.9",
        **kwargs,
    ):
        if cell_type is None and markers is None:
            raise ValueError("cell_type or markers must be provided")
        if cell_type is not None:
            markers = self.marker_db[cell_type]

        if dtype == "RNA" or dtype == "rna":
            dtype = "RNA"
            adata = self.rna_data
        elif dtype == "Protein" or dtype == "protein":
            dtype = "Protein"
            if self.protein_data is None:
                raise ValueError("protein data not provided")
            adata = self.protein_data
        else:
            raise ValueError("dtype must be rna or protein")

        adapted_markers = self.adapt_markers(adata, markers)
        print(f"The following markers were retained: {adapted_markers}")
        if len(adapted_markers) == 0:
            return 0

        cluster_key = cluster_key or self.cluster_key

        if dotplot:
            if cluster_key:
                marker_dotplot(adata, adapted_markers, title=dtype, groupby=cluster_key)
            else:
                print("No cluster key provided, cannot plot dotplot!")

        if marker_plot:
            if title_with_dtype:
                title = [f"{marker} ({dtype.lower()})" for marker in adapted_markers]
            else:
                title = adapted_markers
            if cell_type is not None:
                title = [f"{cell_type}: {subtitle}" for subtitle in title]
            if use_default_plot:
                fig = sc.pl.umap(
                    adata,
                    color=adapted_markers,
                    title=title,
                    show=False,
                    vmin=0,
                    vmax=vmax,
                    sort_order=False,
                    frameon=True,
                    use_raw=False,
                    return_fig=True,
                    color_map="Reds",
                    **kwargs,
                )
                # plt.suptitle(dtype)
                # plt.tight_layout()
                return fig
            else:
                if cell_type is not None:
                    name = cell_type
                else:
                    name = ""
                plot_celltype_markers(
                    adata,
                    marker_genes=adapted_markers,
                    title=dtype,
                    name=name,
                    **kwargs,
                )

    def calculate_marker_scores(self, adata, markers):
        # calculate marker scores
        pass

    def plot_marker_scores(
        self, markers=None, dtype="rna", cluster_key="cell_type", **kwargs
    ):
        if cell_types is None:
            cell_types = list(self.marker_db.keys())
        elif isinstance(cell_types, str):
            cell_types = [cell_types]

        # marker score dotplot and feature plot

    def plot_celltype_marker_scores(self, cell_types=None, **kwargs):
        assert self.marker_db is not None
        if cell_types is None:
            cell_types = list(self.marker_db.keys())
        elif isinstance(cell_types, str):
            cell_types = [cell_types]


def plot_cell_type_distribution(adata, cluster_key):
    ax = adata.obs[cluster_key].value_counts().plot.bar()
    # rotate labels
    labels = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    # reduce fontsize
    ax.tick_params(axis="both", which="major", labelsize=7)
    ax.set_ylabel("Number of cells", fontsize=7)
    plt.show()


def plot_cell_type_proportions(
    adata,
    celltype_col="cell_type",
    groupby="tissue",
    figsize=(1.5, 10),
    ylabel="Cell type fractions",
    xlabel="Tissue",
    rot=None,
):
    all_data_tissue = []
    for tissue in adata.obs[groupby].unique():
        mask = adata.obs[groupby] == tissue
        tissue_props = adata.obs.loc[mask, celltype_col].value_counts() / mask.sum(
            axis=0
        )
        local_df = pd.DataFrame()
        # local_df["cluster"] = tissue_props.index.astype(int).tolist()
        local_df["cluster"] = tissue_props.index.tolist()
        local_df["props"] = tissue_props.tolist()
        local_df[groupby] = tissue
        all_data_tissue.append(local_df)
    all_data_tissue = pd.concat(all_data_tissue)
    all_data_tissue = all_data_tissue.sort_values(by="cluster")
    all_data_tissue_wide = all_data_tissue.pivot(
        index="cluster", columns=groupby, values="props"
    )
    fig, ax = plt.subplots(figsize=figsize)
    # generate list of random colors
    colors = []
    num_clusters = len(all_data_tissue_wide.index)
    cm = plt.get_cmap("plasma")
    colors = [cm(i // 1 * 1.0 / num_clusters) for i in range(num_clusters)]
    colors = None
    ax = all_data_tissue_wide.T.plot.bar(
        stacked=True,
        ax=ax,
        rot=90,
        width=0.8,
        color=colors,
        edgecolor="black",
        linewidth=1.0,
        fontsize=11,
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    # rotate labels
    if rot is not None:
        labels = ax.set_xticklabels(ax.get_xticklabels(), rotation=rot, ha="right")

    # change legend heading
    # change legend position
    legend = ax.legend(
        bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0, fontsize=11
    )
    legend.set_title("Cell type")


def plot_tissue_props(adata, cluster_key="cell_type", rot=None, pad=0):
    all_data_tissue = []
    for cluster in adata.obs[cluster_key].unique():
        local_df = pd.DataFrame()
        mask = adata.obs[cluster_key] == cluster
        cluster_props = adata.obs.loc[mask, "tissue"].value_counts() / mask.sum(axis=0)
        local_df["tissue"] = cluster_props.index.tolist()
        local_df["props"] = cluster_props.tolist()
        local_df["cluster"] = cluster
        all_data_tissue.append(local_df)
    all_data_tissue = pd.concat(all_data_tissue)
    # all_data_tissue = all_data_tissue.sort_values(by="cluster", key=lambda x: x.astype(int))
    all_data_tissue_wide = all_data_tissue.pivot(
        index="tissue", columns="cluster", values="props"
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    ax = all_data_tissue_wide.T.plot.bar(
        stacked=True,
        ax=ax,
        rot=90,
        color=["tab:red", "tab:blue"],
        width=0.9,
        edgecolor="black",
        linewidth=1.0,
    )
    ax.set_ylabel("Cell proportions")
    ax.set_xlabel("Cell type")
    if rot is not None:
        labels = ax.set_xticklabels(ax.get_xticklabels(), rotation=rot, ha="right")
    # ax.tick_params(axis='x', pad=pad)
    # change legend heading
    legend = ax.get_legend()
    legend.set_title("Tissue")
    return fig


def plot_violin_scores(
    adata, markers, score_name, groupby="cell_type_joint", ylabel=None
):
    sc.tl.score_genes(
        adata,
        gene_list=markers,
        score_name=score_name,
        ctrl_size=50,
        use_raw=False,
    )
    fig, ax = plt.subplots(figsize=(8, 4))
    if ylabel is None:
        ylabel = score_name
    sc.pl.violin(
        adata, keys=[score_name], groupby=groupby, ylabel=[ylabel], rotation=90, ax=ax
    )


def adjust_bar_spacing(ax, spacings):
    texts = ax.get_yticklabels()
    assert len(ax.patches) == len(spacings)
    positions = []
    widths = []
    for k, (patch, spacing) in enumerate(zip(ax.patches, spacings)):
        # print(patch.get_x(), patch.get_height(), spacing)
        if k == 0:
            pass
        else:
            patch.set_y(positions[-1] + widths[-1] + spacing)

        positions.append(patch.get_y())
        widths.append(patch.get_height())
    ticks = [pos + width / 2 for pos, width in zip(positions, widths)]
    # for text in texts:
    #     text.set_y(text.get_position()[1] + 0.5)
    return positions, ticks


def separate_barplot(
    kidney_cells_new, labels_new, colors_new, widths, blood_cells_new, spacings
):
    with sns.axes_style("ticks"):
        fig, axs = plt.subplots(
            1, 2, figsize=(4, 8), sharey=True, gridspec_kw={"wspace": 0.1}
        )
        # use seaborn barplot and make vertical bar plot

        sns.barplot(
            x=kidney_cells_new,
            y=labels_new,
            ax=axs[0],
            palette=colors_new,
            order=labels_new,
            edgecolor="black",
            width=widths,
            # dodge=0.0,
            linewidth=1,
        )
        sns.barplot(
            x=blood_cells_new,
            y=labels_new,
            ax=axs[1],
            palette=colors_new,
            order=labels_new,
            edgecolor="black",
            width=widths,
            linewidth=1,
        )

        positions, yticks = adjust_bar_spacing(axs[0], spacings)
        positions, yticks = adjust_bar_spacing(axs[1], spacings)
        axs[0].set_yticks(yticks)
        y_lim = [yticks[0] - widths[0] / 2 - 0.1, yticks[-1] + widths[-1] / 2 + 0.1]
        axs[0].set_ylim(y_lim)
        axs[0].set_yticklabels(labels_new)
        # print(axs[0].get_yticklabels())

        x_min = 0
        x_max = 42
        steps = 10
        axs[0].set_xlim([0, x_max])
        axs[1].set_xlim([0, x_max])
        # add x ticks every 10%
        axs[0].set_xticks(np.arange(0, x_max + 1, steps))
        axs[1].set_xticks(np.arange(0, x_max + 1, steps))
        # add x tick labels with % sign
        fontsize = 8
        axs[0].set_xticklabels(
            [str(x) + "%" for x in np.arange(0, x_max + 1, steps)], fontsize=fontsize
        )
        axs[1].set_xticklabels(
            [str(x) + "%" for x in np.arange(0, x_max + 1, steps)], fontsize=fontsize
        )

        axs[0].set_title("Kidney")
        axs[1].set_title("Blood")
        axs[0].grid(axis="both", linestyle="--", linewidth=0.5)
        axs[1].grid(axis="both", linestyle="--", linewidth=0.5)
        # insert xlabel acroos both axes
        fig.text(
            0.56,
            0.1,
            "Relative frequencies",
            ha="center",
            fontsize=10,
            # fontweight="bold",
        )
        yticklabels = axs[0].get_yticklabels()
        # use mathtext
        for ylabel in yticklabels:
            text = ylabel.get_text()
            ylabel.set_text(r"${}$".format(text.replace("+", "^+")))
        sns.despine(top=True, right=True)
        # plt.subplots_adjust(wspace=-3, hspace=0)
        plt.tight_layout()
        plt.show()


def donut_plot(
    adata, cell_type_col, color_map, replacements={}, label_order=None, center_text=None
):
    value_counts = adata.obs[cell_type_col].value_counts()
    value_props = value_counts / value_counts.sum()
    subdata_df = pd.DataFrame(value_props).reset_index()  #
    subdata_df.columns = ["cell_type", "proportion"]
    if label_order is not None:
        subdata_df.set_index("cell_type", inplace=True)
        subdata_df = subdata_df.reindex(label_order)
        subdata_df.reset_index(inplace=True)
    subdata_df["cell_type_new"] = subdata_df["cell_type"].replace(replacements)
    labels = subdata_df["cell_type_new"].tolist()
    # Example data
    fig, ax = plt.subplots(figsize=(5, 4.0), subplot_kw=dict(aspect="equal"))

    # Create the pie chart. Note the `wedgeprops` argument.
    (
        wedges,
        texts,
        # autotexts,
    ) = ax.pie(
        subdata_df["proportion"].tolist(),
        labels=labels,
        # autopct="%1.1f%%",
        startangle=90,
        wedgeprops={"linewidth": 5, "edgecolor": "white", "width": 0.45},
        colors=[color_map[x] for x in labels],
        labeldistance=1.05,  # Increase this value to adjust label positions
        # labeldistance=0.8,
        # pctdistance=0.75,
    )
    kwargs = dict(size=16, fontweight="bold", va="center")
    if center_text is not None:
        ax.text(0, 0, center_text, ha="center", **kwargs)
    # adjust text fontsize
    for text in texts:
        text.set_fontsize(14)

    # for autotext in autotexts:
    #     autotext.set_fontsize(10)

    # Improve appearance: equal aspect ratio ensures pie is drawn as a circle
    ax.axis("equal")

    # Optionally, you can make font size bigger for better clarity
    # for text, autotext in zip(texts, autotexts):
    #     text.set(size=12)
    #     autotext.set(size=12)
    plt.tight_layout()
    # plt.savefig(
    #     "figures/main/exploratory_cd4_donut.pdf",
    #     dpi=300,
    #     bbox_inches="tight",
    # )
    plt.show()


def plot_stacked_bar(
    adata,
    cell_type_col,
    color_map,
    replacements=None,
    label_order=None,
    ylabel=None,
    xticklabels=None,
    figsize=(2, 14),
    width=0.5,
    save_path=None,
):
    if not isinstance(adata, list):
        adatas = [adata]
    else:
        adatas = adata
    if not isinstance(cell_type_col, list):
        cell_type_cols = [cell_type_col]
    else:
        cell_type_cols = cell_type_col
    if not isinstance(color_map, list):
        color_maps = [color_map] * len(adatas)
    else:
        color_maps = color_map
    if not isinstance(replacements, list):
        replacements_list = [replacements] * len(adatas)
    else:
        replacements_list = replacements
    if not isinstance(label_order, list):
        label_order_list = [label_order] * len(adatas)
    else:
        label_order_list = label_order
    if not len(label_order_list) == len(adatas):
        label_order_list = [label_order_list]

    assert (
        len(adatas)
        == len(cell_type_cols)
        == len(color_maps)
        == len(replacements_list)
        == len(label_order_list)
    ), "Please make sure all lists have the same length"

    data_list = []
    labels_list = []

    for k, adata in enumerate(adatas):
        replacements = replacements_list[k]
        label_order = label_order_list[k]

        value_counts = adata.obs[cell_type_cols[k]].value_counts()
        value_props = value_counts / value_counts.sum()
        subdata_df = pd.DataFrame(value_props).reset_index()  #
        subdata_df.columns = ["cell_type", "proportion"]
        subdata_df["cell_type_new"] = subdata_df["cell_type"]
        if replacements is not None:
            subdata_df["cell_type_new"] = subdata_df["cell_type"].replace(replacements)

        if label_order is not None:
            subdata_df.set_index("cell_type_new", inplace=True)
            subdata_df = subdata_df.reindex(label_order)
            subdata_df.reset_index(inplace=True)
        labels = subdata_df["cell_type_new"].tolist()
        proportions = subdata_df["proportion"].tolist()
        # Create a DataFrame from the given proportions and labels for easier plotting
        data = pd.DataFrame({"Labels": labels, "Proportions": proportions})
        # Normalize the proportions so that they sum up to 1 (100%)
        data["Proportions"] /= data["Proportions"].sum()
        data_list.append(data)
        labels_list.append(labels)
    # make sure everything has the same length
    max_len = max([len(x) for x in labels_list])
    # pad with empty strings and zero
    for k, labels in enumerate(labels_list):
        labels_list[k] = labels + [""] * (max_len - len(labels))
        tmp = pd.DataFrame(
            {
                "Labels": [""] * (max_len - len(labels)),
                "Proportions": [0] * (max_len - len(labels)),
            }
        )
        data_list[k] = pd.concat([data_list[k], tmp], axis=0)
    handles = []
    labels_text = []
    num_bars = len(data_list)
    x = np.arange(num_bars)
    with sns.axes_style("ticks"):
        # Set up the matplotlib figure
        fig, ax = plt.subplots(figsize=figsize)

        # Initialize a bottom value to stack the bars.
        bottom = np.array([0 for _ in x], dtype=float)
        # colors = [color_map[x] for x in labels]

        for i in range(max_len):
            labels = [labels_list[k][i] for k in range(num_bars)]
            proportions = [
                data[data["Labels"] == labels_list[k][i]]["Proportions"].iloc[0]
                for k, data in enumerate(data_list)
            ]
            colors = [color_maps[k][labels_list[k][i]] for k in range(num_bars)]
            ax.bar(
                x,
                proportions,
                bottom=bottom,
                color=colors,
                width=width,
            )
            handles.append(
                [
                    plt.Rectangle((0, 0), 1, 1, fc=color, edgecolor="none")
                    for color in colors
                ]
            )
            labels_text.append(labels)

            for j, label in enumerate(labels):
                if label == "":
                    continue
                y_position = bottom[j] + (proportions[j] / 2)
                x_position = x[j]
                # ax.text(
                #     x_position,
                #     y_position,
                #     label,
                #     # ha="center",
                #     ha="center",
                #     va="center",
                #     fontsize=11,
                #     # fontweight="bold",
                #     rotation=90,
                #     color="black",
                #     # color="grey",
                # )
            bottom += np.array(proportions)
        dummy_handle = plt.Rectangle(
            (0, 0), 0, 0, fc="none", edgecolor="none", linewidth=0
        )

        handles = np.array(handles, dtype=object).transpose().flatten().tolist()[::-1]
        labels_text = (
            np.array(labels_text, dtype=object).transpose().flatten().tolist()[::-1]
        )
        split_index = 0
        handles.insert(split_index, dummy_handle)
        labels_text.insert(split_index, "CD4+ EM/RM")
        split_index = 7
        handles.insert(split_index, dummy_handle)
        labels_text.insert(split_index, "\nCD8+ EM/RM")
        anchors = [(1.0, 1.0), (0.85, 0.7)]

        legend = ax.legend(
            handles, labels_text, loc=2, bbox_to_anchor=anchors[0], frameon=False
        )

        # make split text bold
        legend.get_texts()[0].set_fontweight("bold")
        legend.get_texts()[split_index].set_fontweight("bold")

        # for i in range(num_bars):
        #     ax.legend(handles[:, i], labels_text[:, i], loc=2, bbox_to_anchor=anchors[i], frameon=False)
        # ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(0, 1.01)
        ax.set_xticks(x, xticklabels, rotation=0, ha="center")

        # Adjust and Show the plot
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        # despine
        sns.despine(top=True, right=True)

        # plt.subplots_adjust(left=1.0)
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
        plt.show()


def plot_umap(
    adata,
    color,
    title,
    figsize=(6, 5),
    size=30,
    adjust_cbar=False,
    save_path=None,
    axes_fraction=0.2,
    label_dist=0.03,
    legend=True,
    ncol=1,
    plot_clusters_separate=False,
    alpha=0.6,
    **kwargs,
):
    fig, ax = plt.subplots(figsize=figsize)

    if plot_clusters_separate:
        color_values = adata.obs[color].cat.categories.tolist()
        for k in range(len(color_values)):
            kwargs_copy = kwargs.copy()
            if "palette" in kwargs:
                kwargs_copy["palette"] = [kwargs["palette"][k]]
            plot = sc.pl.umap(
                adata[adata.obs[color] == color_values[k]],
                color=color,
                # legend_loc="right",
                size=size,
                title=title,
                # return_fig=True,
                show=False,
                ax=ax,
                alpha=1 if k < 5 else alpha,
                **kwargs_copy,
            )
    else:
        plot = sc.pl.umap(
            adata,
            color=color,
            # legend_loc="right",
            size=size,
            title=title,
            # return_fig=True,
            show=False,
            ax=ax,
            **kwargs,
        )
    # ax = fig.get_axes()[0]
    # Despine right and top
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Adjust the left spine
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin
    new_ymax = ymin + axes_fraction * yrange
    ax.spines["left"].set_bounds(ymin, new_ymax)

    # Adjust the y label position
    ax.yaxis.set_label_coords(-label_dist, axes_fraction / 2)

    # Adjust the bottom spine
    xmin, xmax = ax.get_xlim()
    xrange = xmax - xmin
    new_xmax = xmin + axes_fraction * xrange
    ax.spines["bottom"].set_bounds(xmin, new_xmax)

    # Set spine thickness
    spine_thickness = 2  # you can adjust this value as you like

    ax.spines["left"].set_linewidth(spine_thickness)
    ax.spines["bottom"].set_linewidth(spine_thickness)

    # Adjust the x label position
    ax.xaxis.set_label_coords(axes_fraction / 2, -label_dist)

    if adjust_cbar:
        cax = fig.get_axes()[1]
        x0, y0, width, height = cax.get_position().bounds
        height_scaler = 0.5
        width_scaler = 1.5  # Change this to your desired scaling factor

        # Calculate the new width
        new_width = width * width_scaler

        # Adjust the x0 (left edge) position to accommodate the change in width while keeping the colorbar in the same relative position to the main plot
        # x0 = x0 - (new_width - width) / 2

        cax.set_position(
            [x0, y0 + 0.5 * height_scaler * height, new_width, height * height_scaler]
        )
    # Remove the existing legend
    ax.legend_.remove()
    if legend:
        # Create a new legend with one column
        # Create a new legend without a border and position it outside the plot
        legend = ax.legend(
            ncol=ncol, frameon=False, bbox_to_anchor=(1, 0.5), loc="center left"
        )
        # change alpha values of legend
        for k, lh in enumerate(legend.legendHandles):
            if k > 4:
                lh.set_alpha(alpha)

    # Modify the legend text
    # for text in legend.get_texts():
    #     text.set_text(r"${}$".format(text.get_text().replace("+", "^+")))

    # Adjust plot parameters to make space for the legend
    plt.subplots_adjust(right=0.8)
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight", transparent=True)
    plt.show()


def plot_composition(adata, groupby, key, normalize=True):
    # Assuming 'adata.obs.groupby("leiden_1.4")["patient"].value_counts()' gives a Series
    data = adata.obs.groupby(groupby)[key].value_counts().reset_index(name="counts")

    # Pivot this DataFrame to make 'leiden_1.4' as index, 'patient' as columns and 'counts' as values
    data_pivot = data.pivot(index=groupby, columns=key, values="counts").fillna(0)

    # Normalize the counts by dividing by the sum of counts for each leiden cluster (i.e., row-wise sum)
    if normalize:
        normalized_data_pivot = data_pivot.div(data_pivot.sum(axis=1), axis=0)
    else:
        normalized_data_pivot = data_pivot

    # Plot a stacked bar plot with normalized data
    normalized_data_pivot.plot(kind="bar", stacked=True, figsize=(10, 7))

    # Add labels and title
    plt.xlabel(groupby)
    plt.ylabel(f"Proportion of {key}")
    # plt.title("Normalized Distribution of Patients in each Leiden Cluster")

    # Rotate x-axis labels if they are too long
    plt.xticks(rotation=90)

    # Show legend
    plt.legend(title=key, bbox_to_anchor=(1.05, 1), loc="upper left")

    # Show the plot
    plt.tight_layout()
    plt.show()


def plot_single_stacked_bar(
    adata,
    labels,
    color_map,
    celltype_col,
    ylabel=None,
    xlabel=None,
    save_path=None,
    add_labels=True,
    figsize=(0.6, 10),
    despine_left=False,
    alpha=0.6,
    ylabel_fontsize=17,
):
    value_counts = adata.obs[celltype_col].value_counts()
    value_props = value_counts / value_counts.sum()
    subdata_df = pd.DataFrame(value_props).reset_index()  #
    subdata_df.columns = ["cell_type", "proportion"]
    subdata_df["cell_type"] = subdata_df["cell_type"].replace(
        {
            "CD4+ EM/RM": "uncertain",
        }
    )
    subdata_df.set_index("cell_type", inplace=True)
    subdata_df = subdata_df.loc[labels, :].reset_index()
    subdata_df["cell_type_new"] = subdata_df["cell_type"].replace(
        {
            # "TH1": "TRM1",
            # "TH17": "TRM17",
            # "CD4+ EM": "other CD4+ EM",
        }
    )

    labels = labels[::-1]
    proportions = subdata_df["proportion"].tolist()[
        ::-1
    ]  # You can directly use the proportions list if available

    # Create a DataFrame from the given proportions and labels for easier plotting
    data = pd.DataFrame({"Labels": labels, "Proportions": proportions})

    # Normalize the proportions so that they sum up to 1 (100%)
    data["Proportions"] /= data["Proportions"].sum()
    with sns.axes_style("ticks"):
        # Set up the matplotlib figure
        fig, ax = plt.subplots(figsize=figsize)

        # Initialize a bottom value to stack the bars.
        bottom = 0
        colors = [color_map[x] for x in labels]

        for i, label in enumerate(labels):
            proportion = data[data["Labels"] == label]["Proportions"].iloc[0]
            ax.bar(
                [0],
                proportion,
                bottom=bottom,
                color=colors[i],
                width=0.7,
                rasterized=True,
                alpha=1 if i > len(labels) - 5 else alpha,
            )
            y_position = bottom + (proportion / 2)
            x_position = 0
            x_position = 0.53
            if add_labels:
                ax.text(
                    x_position,
                    y_position,
                    label,
                    # ha="center",
                    ha="left",
                    va="center",
                    fontsize=14,
                    # fontweight="bold",
                    # rotation=90,
                    color="black",
                    # color="grey",
                )
            bottom += proportion
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(0, 1.0)
        ax.set_xticklabels([])
        sns.despine(top=True, right=True, left=despine_left)

        # Adjust and Show the plot
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontdict={"fontsize": ylabel_fontsize})
        if xlabel is not None:
            ax.set_xlabel(xlabel, fontdict={"fontsize": 14})

        if despine_left:
            ax.tick_params(axis="y", colors="white")
            # make ylabel white
            ax.yaxis.label.set_color("white")
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
        plt.show()


def plot_tissue_compostion(
    adata,
    color_map,
    celltype_key="cell_type",
    tissue_key="tissue",
    label_order=None,
    use_rel_values=True,
    save_path=None,
    num_white=5,
    figsize=(5, 4),
):
    kidney_cells = (
        adata[adata.obs[tissue_key] == "K", :].obs[celltype_key].value_counts()
    )
    blood_cells = (
        adata[adata.obs[tissue_key] == "B", :].obs[celltype_key].value_counts()
    )
    labels = kidney_cells.index.tolist()
    assert (
        set(labels)
        == set(blood_cells.index.tolist())
        == set(kidney_cells.index.tolist())
    )
    kidney_cells = kidney_cells.loc[labels].tolist()
    blood_cells = blood_cells.loc[labels].tolist()

    if use_rel_values:
        kidney_cells = np.array(kidney_cells) * 100 / np.sum(kidney_cells)
        blood_cells = np.array(blood_cells) * 100 / np.sum(blood_cells)

    data_df = pd.DataFrame(columns=["cell_type", "kidney_props", "blood_props"])
    data_df[celltype_key] = labels
    data_df["kidney_props"] = kidney_cells
    data_df["blood_props"] = blood_cells

    data_df = (
        data_df.groupby(celltype_key)
        .sum()
        .reset_index()
        .sort_values(by="kidney_props", ascending=False)
    )
    data_df["sort_helper"] = data_df[celltype_key].apply(
        lambda x: 1 if x == "others" else 0
    )
    data_df = data_df.sort_values(
        by=["sort_helper", "kidney_props"], ascending=[True, False]
    )
    if label_order is not None:
        label_order = label_order + [
            label for label in labels if label not in label_order
        ]
        data_df.set_index(celltype_key, inplace=True)
        data_df = data_df.reindex(label_order)
        data_df.reset_index(inplace=True)

    labels_agg = data_df[celltype_key].tolist()
    kidney_cells_agg = data_df["kidney_props"].tolist()
    blood_cells_agg = data_df["blood_props"].tolist()
    colors_agg = [color_map[x] for x in labels_agg]

    x_min = 0
    x_max = 42
    steps = 5
    diff = 10
    linewidth = 0.01
    default_color = "white"

    with sns.axes_style("ticks"):
        gs = gridspec.GridSpec(
            # 1, 2, width_ratios=[x_max / (x_max - diff), 1], wspace=-0.0
            1,
            2,
            width_ratios=[1, x_max / (x_max - diff)],
            wspace=0.0,
        )  # Change these ratios as needed
        fig = plt.figure(figsize=figsize)
        # fig = plt.figure(figsize=(6, 6))
        axs = [plt.subplot(gs[0]), plt.subplot(gs[1])]
        sns.despine(ax=axs[0], top=True, right=True, bottom=False, left=True)
        sns.despine(ax=axs[1], top=True, right=True, left=True)

        # axs[0].axvline(
        #     0, color="black", linestyle="--", linewidth=2
        # )  # adjust color, linestyle, and linewidth as desired
        # #
        widths = [0.5 for label in labels_agg]
        sns.barplot(
            # x=kidney_cells_agg,
            x=blood_cells_agg,
            y=labels_agg,
            ax=axs[0],
            palette=colors_agg,
            order=labels_agg,
            edgecolor=colors_agg,
            # edgecolor="black",
            width=widths,
            # dodge=0.0,
            linewidth=linewidth,
            rasterized=True,
        )
        sns.barplot(
            # x=blood_cells_agg,
            x=kidney_cells_agg,
            y=labels_agg,
            ax=axs[1],
            palette=colors_agg,
            order=labels_agg,
            edgecolor=colors_agg,
            width=widths,
            linewidth=linewidth,
            rasterized=True,
        )
        # Place y-labels inside bars for axs[0]
        bar_heights = [rect.get_height() for rect in axs[0].patches]
        x_positions = [
            0 - (rect_1.get_width() - rect_2.get_width()) / 2
            for rect_1, rect_2 in zip(axs[0].patches, axs[1].patches)
        ]
        num_bars = len(labels_agg)
        num_white = num_bars if num_white == -1 else num_white
        colors_text = num_white * [default_color] + (num_bars - num_white) * ["black"]
        # colors_text = [default_color for label in labels_agg]
        # colors_text[-1] = "black"
        for index, label in enumerate(labels_agg):
            # Adjust -0.05 and 1.02 values to better fit the labels inside the bars
            plt.text(
                x_positions[index],
                # 0,
                index + bar_heights[index] * 0.09,
                label,
                ha="center",
                va="center",
                color=colors_text[index],
                fontsize=10,
            )

        axs[1].set_yticks([])
        axs[1].set_yticklabels([])

        axs[0].set_yticks([])
        axs[0].set_yticklabels([])

        axs[1].set_xlim([0, x_max])
        axs[0].set_xlim([0, x_max - diff])
        # add x ticks every 10%
        range_1 = np.arange(0, x_max + 1, steps)
        range_0 = np.arange(0, x_max + 1 - diff, steps)
        axs[0].set_xticks(range_0)
        axs[1].set_xticks(range_1)
        # add x tick labels with % sign
        fontsize = 9
        x_ticklabels_0 = [str(x) + "" for x in range_0]
        x_ticklabels_0[0] = ""
        axs[0].set_xticklabels(x_ticklabels_0, fontsize=fontsize)
        axs[1].set_xticklabels([str(x) + "" for x in range_1], fontsize=fontsize)
        axs[0].invert_xaxis()
        axs[0].grid(axis="x", linestyle="--", linewidth=0.5)
        axs[1].grid(axis="x", linestyle="--", linewidth=0.5)

        # Try to find and adjust the gridline at x=0
        gridlines = axs[0].get_xgridlines()
        for line in gridlines:
            if line.get_xdata()[0] == 0:  # if the x-position of the gridline is 0
                line.set_visible(False)
                break

        gridlines = axs[1].get_xgridlines()
        for line in gridlines:
            if line.get_xdata()[0] == 0:  # if the x-position of the gridline is 0
                # print(line.get_xdata())
                line.set_color("black")
                line.set_linewidth(1)
                # line.set_xdata((-10,))

        # joint x label
        fig.text(
            0.5,
            0.0,
            "Relative frequencies of blood and kidney cells [%]",
            ha="center",
            va="center",
            rotation="horizontal",
            fontsize=10,
        )
        axs[0].set_title("Blood", fontsize=10)
        axs[1].set_title("Kidney", fontsize=10)

        plt.tight_layout()
        if save_path is not None:
            plt.savefig(
                save_path,
                bbox_inches="tight",
                dpi=300,
                transparent=True,
            )
        plt.show()


def plot_scores(
    adata, score, cmap, title=None, plt_center=0.7, figsize=(5.5, 5), save_path=None
):
    # # score = "cytokine_score"
    # score = "type3_score"
    # plt_center = 0.06

    # score = "type2_score"
    # plt_center = 0.06

    # # score = "type1_score"
    # # center = 0.06

    # score = "cytokine_score"
    # plt_center = 0.7

    # Create a new figure with specified size
    fig, ax = plt.subplots(figsize=figsize)  # You can adjust the size (8, 6) as needed
    size = 0.5
    limit = 1.0
    limit_up = 1.0
    limit_low = -1.0
    colorbar = True
    cbar_kwargs = {"label": "Color Value"}

    cmaps = ["RdBu_r", "RdBu_r"]
    cmaps = [cmap, cmap]

    norm = mcolors.TwoSlopeNorm(vmin=limit_low, vcenter=0.0, vmax=limit_up)
    with sns.axes_style("white"):
        alphas = [1.0, 1.0]
        for i in range(2):
            if i == 0:
                slice = adata[adata.obs[score] <= plt_center, :]
                c = slice.obs[score]
                cmap = cmaps[0]
                alpha = alphas[0]
            else:
                slice = adata[adata.obs[score] > plt_center, :]
                c = slice.obs[score]
                cmap = cmaps[1]
                alpha = alphas[1]

            x = slice.obsm["X_umap"][:, 0]
            y = slice.obsm["X_umap"][:, 1]
            scatter = ax.scatter(
                x,
                y,
                c=c,
                s=size,
                cmap=cmap,
                # vmin=limit_low,
                # vmax=limit_up,
                norm=norm,
                alpha=alpha,
                edgecolor=None,
                facecolor=None,
                rasterized=True,
            )
        # ax.set_title(score.split("_")[0].capitalize(), fontsize=12, pad=-7)
        # Add colorbar for reference
        if colorbar:
            # Create a new axes for the colorbar with [left, bottom, width, height]
            cbar_ax = fig.add_axes(
                [0.86, 0.35, 0.03, 0.3]
            )  # You might need to adjust these values

            # Add colorbar in the new axes
            cbar = fig.colorbar(scatter, cax=cbar_ax)
            # Update the colorbar limits
            # cbar_ax.set_ylim(-limit, limit)

            # cbar.set_label(cbar_kwargs["label"], fontsize=12)
        # Despine right and top
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])
        # Shorten left and bottom spines
        fraction = 1.0  # fraction of original length you want to keep
        label_dist = 0.03

        # Adjust the left spine
        ymin, ymax = ax.get_ylim()
        yrange = ymax - ymin
        new_ymax = ymin + fraction * yrange
        ax.spines["left"].set_bounds(ymin, new_ymax)

        # Adjust the y label position
        ax.yaxis.set_label_coords(-label_dist, fraction / 2)

        # Adjust the bottom spine
        xmin, xmax = ax.get_xlim()
        xrange = xmax - xmin
        new_xmax = xmin + fraction * xrange
        ax.spines["bottom"].set_bounds(xmin, new_xmax)

        # Set spine thickness
        spine_thickness = 2  # you can adjust this value as you like

        ax.spines["left"].set_linewidth(spine_thickness)
        ax.spines["bottom"].set_linewidth(spine_thickness)

        # Adjust the x label position
        ax.xaxis.set_label_coords(fraction / 2, -label_dist)
        ticks = ax.set_xticklabels([])
        ticks = ax.set_yticklabels([])
        plt.subplots_adjust(right=0.8)
        if title is not None:
            ax.set_title(title, pad=-6)
        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
        plt.show()


def dotplot_markers(
    adata,
    celltype_col,
    celltype_order=None,
    marker_dict=None,
    figsize=(11, 5),
    scaled=False,
    vmax=1,
    mod="rna",
    save_path=None,
    title=None,
    pad=0,
    rotation=45,
    return_fig=False,
):
    if scaled:
        layer = "scaled"
        cmap = "RdBu_r"
        standard_scale = None
        vmin = -vmax
    else:
        layer = "log1p" if mod == "rna" else "clr"
        standard_scale = "var"
        cmap = "Reds"
        vmin = None
        vmax = None

    if marker_dict is None:
        marker_dict = {
            "general": [
                # "CD3D",
                "CD3E",
                "CD4",
                "CD8A",
                # "CD8B",
            ],
            "Trm": [
                "CD69",
                "CXCR6",
                # "ITGAE",
                "RGS1",
            ],
            "eff.": [
                "CCR6",
                "CXCR3",
                "IFNG",
                "TNF",
            ],
            # "cytokines": [
            #     "IFNG",
            #     "TNF",
            #     # "IL2",
            #     # "IL4",
            #     # "IL17A",
            #     # "IL17F",
            #     # "IL21",
            #     # "IL22",
            # ],
            "naive/Tcm": [
                "LEF1",
                # "TCF7",
                # "LTB",
                "CCR7",
                "SELL",
                # "KLF2",
            ],
            # "naive": [
            #     "LEF1",
            #     "TCF7",
            #     "LTB",
            # ],
            # "Tcm": [
            #     "CCR7",
            #     "SELL",
            #     "KLF2",
            #     "S1PR1",
            # ],
            "Treg": [
                "FOXP3",
                "IL2RA",
                # "CTLA4",
                # "IKZF2",
                # "TIGIT",
            ],
            "gdT": ["TRDV2", "TRGV9"],
            "MAIT": [
                "TRAV1-2",
                # "KLRB1",
            ],
            "cytotoxic": [
                # "GZMK",
                "GZMA",
                # "GZMB",
                # "GZMH",
                # "GNLY",
                "PRF1",
            ],
            "NK/NKT": [
                "KLRB1",
                "KLRK1",
                "NCAM1",
                "FCGR3A",
                # "B3GAT1",
                # "CD1D",
            ],
            "prolif.": [
                "STMN1",
                "MKI67",
                # "TOP2A",
            ],
            # "stress": [
            #     "LYZ",
            #     "CD14",
            #     "CD33",
            #     "ITGAM",
            #     "FCGR3A",
            #     "FCGR3B",
            #     "CD68",
            #     "CD64",
            #     "S100A8",
            #     "S100A9",
            #     "ICAM1",
            # ],
        }
    adata.obs[celltype_col] = pd.Categorical(
        adata.obs[celltype_col], categories=celltype_order, ordered=True
    )

    # fig, ax = plt.subplots(figsize=figsize)
    fig = sc.pl.dotplot(
        adata,
        # figsize=figsize,
        var_names=marker_dict,
        layer=layer,
        groupby=celltype_col,
        dendrogram=False,
        # ax=ax,
        color_map=cmap,
        standard_scale=standard_scale,
        colorbar_title="Scaled expression",
        vmin=vmin,
        vmax=vmax,
        # return_fig=True,
        var_group_rotation=rotation,
        show=False,
        rasterized=True,
    )
    gene_group_ax = fig["gene_group_ax"]
    for text in gene_group_ax.texts:
        text.set_ha("left")
        x, y = text.get_position()
        text.set_position((x - 0.25, y))
    if title is not None:
        gene_group_ax.set_title(title, pad=pad)
    # change each xticks labels text

    ax = fig["mainplot_ax"]
    ax.set_ylabel("Cell types")
    for text in ax.get_xticklabels():
        # change to italic
        text.set_fontstyle("italic")
    # new_labels = []
    # for text in ax.get_yticklabels():
    #     old_str = text.get_text()
    #     parts = old_str.split("+")
    #     # Encapsulate each part in \mathrm{}
    #     parts_in_mathrm = ["\mathrm{" + part + "}" for part in parts]
    #     new_str = "+".join(parts_in_mathrm)
    #     new_str = r"${}$".format(new_str)
    #     new_str = new_str.replace("+", "^+")
    #     new_str = new_str.replace(" ", "\ ")
    #     text.set_text(new_str)
    #     new_labels.append(text)
    # new_labels = add_superscripts(ax.get_yticklabels())
    # ax.set_yticklabels(new_labels)
    plt.draw()
    # gene_group_ax.set_title("")
    # if title is not None:
    #     plt.title(title)
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.show()
    if return_fig:
        return fig
    else:
        return None


def add_superscripts(labels):
    new_labels = []
    for label in labels:
        # check if label is text object
        try:
            label = label.get_text()
        except:
            pass
        if not "+" in label:
            new_labels.append(label)
            continue
        parts = label.split("+")
        # Encapsulate each part in \mathrm{}
        parts_in_mathrm = ["\mathrm{" + part + "}" for part in parts]
        new_label = "+".join(parts_in_mathrm)
        new_label = r"${}$".format(new_label)
        new_label = new_label.replace("+", "^+")
        new_label = new_label.replace(" ", "\ ")
        new_labels.append(new_label)
    return new_labels


def plot_cluster_composition_per_sample(
    adata,
    colors,
    celltype_col,
    sample_col="patient",
    order=None,
    save_path=None,
    figsize=(8, 5),
    title="Cell type composition per patient",
    title_pad=0,
):
    comp_per_patient = adata.obs.groupby([sample_col])[celltype_col].value_counts(
        normalize=True
    )
    comp_per_patient = pd.DataFrame(comp_per_patient)
    comp_per_patient.columns = ["fraction"]
    comp_per_patient = comp_per_patient.reset_index()
    # pivot table
    comp_per_patient = comp_per_patient.pivot(
        index=sample_col, columns=celltype_col, values="fraction"
    )

    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(comp_per_patient.shape[0])
    if order is None:
        order = adata.obs[celltype_col].unique()
    for k, cell_type in enumerate(order):
        sns.barplot(
            x=comp_per_patient.index,
            y=comp_per_patient[cell_type],
            ax=ax,
            color=colors[k],
            width=0.8,
            rasterized=True,
            bottom=bottom,
        )
        bottom += comp_per_patient[cell_type]
    # rotate xlabels
    labels = ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # construct custom legend with colors for each cell type label
    handles = []
    for k, cell_type in enumerate(order):
        handles.append(
            mpl.lines.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=cell_type,
                markerfacecolor=colors[k],
                markersize=10,
            )
        )
    # add legend without frame
    ax.legend(
        handles=handles, bbox_to_anchor=(0.98, 1.035), loc="upper left", frameon=False
    )
    ax.set_ylabel("Fraction")
    ax.set_xlabel("Patient", labelpad=10)
    ax.set_ylim([0, 1])
    ax.set_title(title, pad=title_pad)
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.show()


def default_scatter(adata, cluster_col, figsize=(5, 5.5), size=4):
    sc.set_figure_params(dpi=100)
    sns.set(style="white")
    fig = plt.figure(figsize=figsize)

    annots = list(adata.obs[cluster_col].cat.categories)
    colors = dict(zip(annots, adata.uns[cluster_col + "_colors"]))

    for ct in tqdm(annots):
        sub = adata[adata.obs[cluster_col] == ct]
        x, y, color = (
            np.array(sub.obsm["X_umap"][:, 0]),
            np.array(sub.obsm["X_umap"][:, 1]),
            colors[ct],
        )
        plt.scatter(x, y, s=size, color=[color] * x.shape[0], label=ct)
    plt.axis("off")
    plt.legend(bbox_to_anchor=(1, 1.1), fontsize=10)
    return fig


def compute_centers(adata, cluster_col):
    df = pd.DataFrame(
        {
            "x": np.array(adata.obsm["X_umap"][:, 0]),
            "y": np.array(adata.obsm["X_umap"][:, 1]),
            "label": np.array(adata.obs[cluster_col]),
        }
    )

    centroids = df.groupby("label").mean()
    left_most = centroids.sort_values(by="x").index[0]
    right_most = centroids.sort_values(by="x").index[-1]
    top_most = centroids.sort_values(by="y").index[-1]
    bottom_most = centroids.sort_values(by="y").index[0]

    considered_centroids = list(set([left_most, right_most, bottom_most, top_most]))

    x, y = centroids.sort_values(by="x").iloc[0]

    left_to_right_bottom = []
    for ct in centroids.sort_values(by="x").index:
        if (
            ct not in considered_centroids
            and centroids.loc[ct][0] > x
            and centroids.loc[ct][1] < y
        ):
            left_to_right_bottom.append(ct)
    considered_centroids += left_to_right_bottom

    left_to_top_left = []
    for ct in centroids.sort_values(by="x").index:
        if (
            ct not in considered_centroids
            and centroids.loc[ct][0] > x
            and centroids.loc[ct][1] > y
        ):
            left_to_top_left.append(ct)
    considered_centroids += left_to_top_left

    x, y = centroids.sort_values(by="y").iloc[0]

    bottom_to_right_right = []
    for ct in centroids.sort_values(by="x").index:
        if (
            ct not in considered_centroids
            and centroids.loc[ct][0] > x
            and centroids.loc[ct][1] > y
        ):
            bottom_to_right_right.append(ct)
    considered_centroids += bottom_to_right_right

    num_located_elements = len(considered_centroids)
    assert num_located_elements == len(
        centroids
    ), f"num_located_elements: {num_located_elements}, len(centroids): {len(centroids)}"

    return (
        centroids,
        left_most,
        right_most,
        top_most,
        bottom_most,
        left_to_right_bottom,
        left_to_top_left,
        bottom_to_right_right,
    )


def custom_umap(
    adata,
    cluster_col,
    figsize=(5, 5.5),
    size=4,
    offset=0.5,
    skip=3,
    auto_adjust=False,
    to_adjust=[],
    poses=[],
    aligns=[],
    dpi=100,
    colors=None,
    use_tex=False,
    str_map={},
    fontsize=14,
    alpha=0.5,
    order=None,
):
    sc.set_figure_params(dpi=dpi)
    sns.set(style="white")
    fig = plt.figure(figsize=figsize)

    if order is not None:
        cat_type = pd.CategoricalDtype(categories=order, ordered=True)
        adata.obs[cluster_col] = adata.obs[cluster_col].astype(cat_type)

    annots = list(adata.obs[cluster_col].cat.categories)
    print(f"Considering {annots} as the order of clusters, {len(annots)} clusters")
    if colors is None:
        colors = dict(zip(annots, adata.uns[cluster_col + "_colors"]))
    else:
        colors = dict(zip(annots, colors))

    (
        centroids,
        left_most,
        right_most,
        top_most,
        bottom_most,
        left_to_right_bottom,
        left_to_top_left,
        bottom_to_right_right,
    ) = compute_centers(adata, cluster_col)

    pos = {}
    for ct in annots:
        pos[ct] = tuple(centroids.loc[ct])

    # annots = (
    #     [left_most]
    #     + left_to_right_bottom
    #     + [bottom_most]
    #     + left_to_top_left
    #     + [top_most]
    #     + bottom_to_right_right
    #     + [right_most]
    # )
    texts = []
    i = j = offset

    sub = adata[adata.obs[cluster_col] == bottom_most]
    x, y, color = (
        np.array(sub.obsm["X_umap"][:, 0]),
        np.array(sub.obsm["X_umap"][:, 1]),
        colors[ct],
    )
    bottom_lim = y.min() + skip

    sub = adata[adata.obs[cluster_col] == top_most]
    x, y, color = (
        np.array(sub.obsm["X_umap"][:, 0]),
        np.array(sub.obsm["X_umap"][:, 1]),
        colors[ct],
    )
    top_lim = y.max() - skip

    if len(to_adjust) > 0:
        assert len(to_adjust) == len(poses) == len(aligns)

        poses = dict(zip(to_adjust, poses))
        aligns = dict(zip(to_adjust, aligns))
    print(len(aligns), len(annots), len(poses))

    for ct in tqdm(annots):
        sub = adata[adata.obs[cluster_col] == ct]
        x, y, color = (
            np.array(sub.obsm["X_umap"][:, 0]),
            np.array(sub.obsm["X_umap"][:, 1]),
            colors[ct],
        )
        plt.scatter(x, y, s=size, label=ct, facecolor=[color] * x.shape[0], alpha=alpha)
        name = ct.replace("+", "^+")
        name = name.replace("like ", "like\ ")
        if len(str_map) > 0:
            for k, v in str_map.items():
                name = name.replace(k, v)
        name_str = r"$\it {}$".format(name)
        with plt.rc_context(
            {
                "text.usetex": use_tex,
                "font.family": "sans-serif",
                "font.sans-serif": "Helvetica",
            }
        ):
            if ct not in to_adjust:
                if ct == left_most:
                    texts.append(
                        plt.text(
                            x.min(),
                            pos[ct][1] + j,
                            name_str,
                            weight="bold",
                            ha="right",
                            va="bottom",
                            color=colors[ct],
                            fontdict={"fontsize": fontsize},
                        )
                    )
                elif ct in left_to_top_left:  # pos[ct][1]
                    texts.append(
                        plt.text(
                            x.min() - j,
                            top_lim - j,
                            name_str,
                            weight="bold",
                            ha="right",
                            va="bottom",
                            color=colors[ct],
                            fontdict={"fontsize": fontsize},
                        )
                    )
                    j = j + 1
                if ct == top_most:
                    texts.append(
                        plt.text(
                            pos[ct][0],
                            pos[ct][1],
                            name_str,
                            weight="bold",
                            ha="right",
                            va="bottom",
                            color=colors[ct],
                            fontdict={"fontsize": fontsize},
                        )
                    )

                elif ct == right_most or ct in bottom_to_right_right:
                    texts.append(
                        plt.text(
                            x.max() + 0.5,
                            pos[ct][1],
                            name_str,
                            weight="bold",
                            ha="left",
                            va="top",
                            color=colors[ct],
                            fontdict={"fontsize": fontsize},
                        )
                    )

                elif ct == bottom_most or ct in left_to_right_bottom:
                    texts.append(
                        plt.text(
                            pos[ct][0] - 0.5,
                            bottom_lim - i,
                            name_str,
                            weight="bold",
                            ha="right",
                            va="bottom",
                            color=colors[ct],
                            fontdict={"fontsize": fontsize},
                        )
                    )
                    i = i + 1
            else:
                texts.append(
                    plt.text(
                        poses[ct][0],
                        poses[ct][1],
                        name_str,
                        weight="bold",
                        ha=aligns[ct][0],
                        va=aligns[ct][1],
                        color=colors[ct],
                        fontdict={"fontsize": fontsize},
                    )
                )

    if auto_adjust:
        adjust_text(None)
    plt.xlim(
        np.array(adata.obsm["X_umap"][:, 0]).min() - 1,
        np.array(adata.obsm["X_umap"][:, 0]).max() + 1,
    )
    plt.axis("off")
    return fig, texts, aligns, annots


def custom_umap_slim(
    adata,
    cluster_col,
    figsize=(5, 5.5),
    size=4,
    auto_adjust=False,
    to_adjust=[],
    poses=[],
    aligns=[],
    dpi=100,
    colors=None,
    use_tex=False,
    str_map={},
    fontsize=10,
    alpha=0.2,
    invert_order=False,
    order=None,
):
    sc.set_figure_params(dpi=dpi)
    sns.set(style="white")
    fig = plt.figure(figsize=figsize)

    if order is not None:
        cat_type = pd.CategoricalDtype(categories=order, ordered=True)
        adata.obs[cluster_col] = adata.obs[cluster_col].astype(cat_type)

    annots = list(adata.obs[cluster_col].cat.categories)
    if invert_order:
        annots = annots[::-1]
        colors = colors[::-1]
    print(f"Considering {annots} as the order of clusters")
    if colors is None:
        colors = dict(zip(annots, adata.uns[cluster_col + "_colors"]))
    else:
        colors = dict(zip(annots, colors))

    texts = []

    if len(to_adjust) > 0:
        assert len(to_adjust) == len(poses) == len(aligns)

        poses = dict(zip(to_adjust, poses))
        aligns = dict(zip(to_adjust, aligns))

    for ct in tqdm(annots):
        print(f"Plotting {ct}")
        sub = adata[adata.obs[cluster_col] == ct]
        x, y, color = (
            np.array(sub.obsm["X_umap"][:, 0]),
            np.array(sub.obsm["X_umap"][:, 1]),
            colors[ct],
        )
        plt.scatter(
            x, y, s=size, color=[color] * x.shape[0], label=ct, alpha=alpha, marker="o"
        )
        name = ct.replace("+", "^+")
        name = name.replace("other", "other\ ")
        for k, v in str_map.items():
            name = name.replace(k, v)
        name_str = r"$\it {}$".format(name)
        with plt.rc_context(
            {
                "text.usetex": use_tex,
                "font.family": "sans-serif",
                "font.sans-serif": "Helvetica",
            }
        ):
            texts.append(
                plt.text(
                    poses[ct][0],
                    poses[ct][1],
                    name_str,
                    weight="bold",
                    ha=aligns[ct][0],
                    va=aligns[ct][1],
                    color=colors[ct],
                    fontdict={"fontsize": fontsize},
                )
            )

    if auto_adjust:
        adjust_text(None)
    plt.xlim(
        np.array(adata.obsm["X_umap"][:, 0]).min() - 1,
        np.array(adata.obsm["X_umap"][:, 0]).max() + 1,
    )
    plt.axis("off")
    return fig, texts
