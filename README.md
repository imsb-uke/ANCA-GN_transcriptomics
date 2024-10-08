## Immune profiling-based targeting of pathogenic T cells with ustekinumab in ANCA-associated glomerulonephritis

The manuscript is available at [Nature Communications](https://www.nature.com/articles/s41467-024-52525-w).

### Instructions for single-cell analysis

Our single-cell analysis consists of two parts one for the exploratory cohort and another for the treatment (Ustekinumab) cohort. Both analysis follow the same general workflow of processing, clustering, annotation, and lastly more specific downstream analysis.

**System requirements and installation**

We ran the code in an Ubuntu 20.04 environment. The corresponding docker image can be obtained by running the following command:

```bash
docker pull imsbuke/dsnb:20211025_jh1.4.2
```
You can then run the docker image with the following command and replace ```/PATH_TO_DATA``` with the path to the data on your local machine:

```bash
docker run -it -v /PATH_TO_DATA:/data --user root --name ustekinumab-test imsbuke/dsnb:20211025_jh1.4.2 bash
```

To setup the single-cell analysis Python environment, please run the following commands:

```bash
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
conda env create -f envs/sc-env.yml
conda activate sc-env
```

To reproduce the figures you also need to install the Arial font. You can do this by running the following command:

```bash
apt-get update
apt-get upgrade -y
echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | sudo debconf-set-selections
echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula seen true" | sudo debconf-set-selections
apt-get install ttf-mscorefonts-installer -y
fc-cache -fv
```

For the R-parts of our analysis you need to install the following package:

```R
setRepositories(ind=1:3)
install.packages("Signac")
```

**Data preparation**

You can find the Links to all raw and aligned single-cell data in Supplementary Table 10 of our manuscript. Please place the aligned data in separate folders in the ```data/single-cell/exploratory/alignment``` or  ```data/single-cell/ustekinumab/alignment```directory.

For ease of use we provide the processed and annotated single-cell data for the both cohorts as MuData objects. This data is deposited at GEO under accession [GSE253633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253633) alongside the new single-cell data of this study. Please place the files in their respective folders in the ```data``` directory.


**Analysis workflow**

The analysis workflow for the single-cell data (CITEseq and scRNAseq) is detailed in the folder ```notebooks\single-cell```. We further split the code for each cohort, namely the exploratory and the ustekinumab treatment cohort. The corresponding code is available in the folders ```notebooks\single-cell\exploratory_cohort``` and ```notebooks\single-cell\ustekinumab_cohort```, respectively. Please use the notebooks in ```notebooks/single-cell/exploratory_cohort/00-preprocessing``` or ```notebooks/single-cell/ustekinumab_cohort/00-preprocessing``` for sample wise quality control and demultiplexing. Then you can proceed with the rest of the annotation workflow and downstream analysis.


Both folders also contain the code to reproduce the figures in the manuscript from the processed and annotated data objects, e.g., see ```notebooks\single-cell\exploratory_cohort\plot-main-figures.ipynb```.

### Instructions for Visium analysis

**Dataset**
The aligned data and processed and annotated data objects are available at [GSE250138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250138). The raw data is available through SRA with bioproject [PRJNA1052093](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1052093).

For reproducibility, the docker image is available at Dockerhub and can be pulled using,
```bash
docker pull imsbuke/dsnb:20211025_jh1.4.2
```
The requirements regarding package versions is available in visium_requirements.txt and can be installed in docker container using,
```bash
pip install -r visium_requirements.txt
```

**Data preparation**
1. Place the cellranger- and spaceranger (v2.0.1)-aligned data in a single folder where each subfolder is a Visium slide_ID (_e.g._ V1_D) and contains the alignment results for the slide_ID. Update the path to the aligned datasets in the first notebook.
2. Place the TIF files corresponding to each spatial sample in the folder ```tif_processed```.

The folders ```annotations_visium*``` contain expert-annotations of the Visium slides.

**Analysis workflow**

Guided code for required to reproduce the results is available in the folder ```notebooks```. The second step of the processing is too big to push to GitHub, this can be downloaded from [02_cluster.ipynb](https://drive.google.com/file/d/11mMGel0VzCgbqmvUIG5L2zHP9qBJGoqu/view?usp=sharing)

**If you encounter any problem, please open an issue.**
