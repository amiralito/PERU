# Phylogenomics of PERU and FLS2 LRR-RK immune receptors in Solanaceae plant family

Supporting scripts, data, and methods for the "Functional diversification of a wild potato PRR immune receptor at its center of origin".


Resources:

Software                            | Source
------------------------------------| ------------------------------------
*HMMER v3.3.2*                      | (https://github.com/EddyRivasLab/hmmer)
*MAFFT v7.520*                      | (https://github.com/GSLBiotech/mafft)
*FastTree v2.1.11*                  | (http://www.microbesonline.org/fasttree/)
*Dendroscope v3.8.8*                | (https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)

R packages:

```R
install.packages("tidyverse")
install.packages("readr")
install.packages("pheatmap")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```
