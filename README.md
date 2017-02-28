# Heterogeneity of human colonic stromal cells in health and inflammatory bowel disease revealed by single cell transcriptome analysis

This repository contains rmarkdown reports and helper functions required to reproduce the analysis included in the manuscript. It can be accessed anonymously.

## Gene count tables

Gene count tables are required to compile the reports. These can be obtained from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95459) via the reviewer link. A gene count table for each of the three subseries is required. The file names are:

GSE95435\_P150057\_full\_gene\_count\_table.txt.gz
GSE95446\_P150409\_full\_gene\_count\_table.txt.gz
GSE95450\_P160099\_full\_gene\_count\_table.txt.gz

These files should be decompressed, and the resulting text files moved to the data subdirectory of the cloned repository prior to compiling the reports.

## R packages

R package versions and their sources are listed below:

 Package | Version | Source 
---------|---------|--------
scater | 1.2.0 | bioconductor
limma | 3.30.6 | bioconductor
scran | 1.2.0 | bioconductor
ggplot2 | 2.2.1 | CRAN
ggbeeswarm | 0.5.2 | CRAN
plyr | 1.8.4 | CRAN
RUVSeq | 1.8.0 | bioconductor
WGCNA | 1.51 | bioconductor
lazyeval | 0.2.0 | CRAN
NMF | 0.22 | CRAN
dendextend | 1.4.0 | CRAN
gplots | 3.0.1 | CRAN
viridis | 0.3.4 | CRAN
reshape2 | 1.4.2 | CRAN
Rtsne | 0.11 | CRAN
matrixStats | 0.51.0 | CRAN
clusterProfiler | 3.2.6 | bioconductor
org.Hs.eg.db | 3.4.0 | bioconductor
xlsx | 0.5.7 | CRAN

CRAN packages are available at <https://cran.r-project.org> and can be installed using the `install.packages("package name")` command.

Bioconductor packages are available at <http://bioconductor.org>. Installation instructions can be found at <http://bioconductor.org/install/>.

R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS Sierra 10.12.3