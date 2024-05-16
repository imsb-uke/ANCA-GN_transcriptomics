## Code for the analysis of samples used in the manuscript "Immune profiling-based targeting of pathogenic T cells with ustekinumab in ANCA-associated glomerulonephritis"

### Instructions for single-cell analysis

**System requirements and installation**

We ran the code in an Ubuntu 20.04 environment. To setup the single-cell analysis Python environment, please run the following commands:

```bash
conda env create -f envs/sc-env.yml
conda activate sc-env
```

**Data preparation**

Place the cellranger-aligned data in the folder ```data```. For the correct version please refer to the methods section of the manuscript.

**Analysis workflow**

The analysis workflow for the single-cell data (CITEseq and scRNAseq) is detailed in the folder ```notebooks\single-cell```. We further split the code for each cohort, namely the exploratory and the ustekinumab treatment cohort. The corresponding code is available in the folders ```notebooks\single-cell\exploratory_cohort``` and ```notebooks\single-cell\ustekinumab_cohort```, respectively.

### Instructions for Visium analysis

**Data preparation**
1. Place the cellranger- and spaceranger-aligned data in folder ```data```. For the versions used in the manuscript, please refer to the methods section of the manuscript.
2. Place the TIF files corresponding to each spatial sample in the folder ```tif_processed```

The folders ```annotations_visium*``` contain expert-annotations of the Visium samples

**Analysis workflow**

Guided code for required to reproduce the results is available in the folder ```notebooks```. The second step of the processing is too big to push to GitHub, this can be downloaded from [02_cluster.ipynb](https://drive.google.com/file/d/11mMGel0VzCgbqmvUIG5L2zHP9qBJGoqu/view?usp=sharing)

**If you encounter any problem, please open an issue.**
