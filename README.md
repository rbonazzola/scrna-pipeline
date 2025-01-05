# Single-cell RNA-Seq analysis pipeline
Pipeline for the analysis of single cell RNA-seq data (based on the R pipeline from https://github.com/Bastien-mva/pipeline_seurat_monocle )

## Data
Data for _Arabidopsis thaliana_ inoculated with _Pst_ can be retrieved from: https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/CKVGIN.

## Environment

If you have Conda, you can use it to install the required libraries:

```
# For the Python version
conda install python=3.11
conda install bioconda::scanpy

# For the R version
conda install conda-forge::r-base
conda install conda-forge::r-tidyverse
conda install bioconda::r-seurat
conda install bioconda::r-monocle3

```

Alternatively, you can use a pre-made [Docker image](https://hub.docker.com/r/rbonazzola/scrna) which you can download from DockerHub.

```bash
docker pull rbonazzola/scrna
```
