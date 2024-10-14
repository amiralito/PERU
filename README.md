# Phylogenetics of PERU and FLS2 LRR-RK immune receptors in the Solanaceae plant family

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8079625.svg)](https://doi.org/10.5281/zenodo.8079625)
[![DOI](https://img.shields.io/badge/science-doi.org/10.1126/science.adg5261-C82426.svg)](https://doi.org/10.1126/science.adg5261)

Supporting scripts, data, and methods for the "Functional diversification of a wild potato PRR immune receptor at its center of origin".

Yerisf C. Torres Ascurra, Lisha Zhang, AmirAli Toghani, Chenlei Hua, Nandeesh Jalahalli Rangegowda, Andres Posbeyikian, Hsuan Pai, Xiao Lin, Pieter J. Wolters h, Doret Wouters, Reinhoud de Blok, Niels Steigenga, Maxence J. M. Paillart, Richard G. F. Visser, Sophien Kamoun, Thorsten Nürnberger, and Vivianne G. A. A. Vleeshouwers


Resources:

Software                            | Source
------------------------------------| ------------------------------------
*HMMER v3.3.2*                      | (https://github.com/EddyRivasLab/hmmer)
*MAFFT v7.520*                      | (https://github.com/GSLBiotech/mafft)
*FastTree v2.1.11*                  | (http://www.microbesonline.org/fasttree/)
*Dendroscope v3.8.8*                | (https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)
*SignalP v6.0*                      | (https://services.healthtech.dtu.dk/services/SignalP-6.0/)
*R*                                 | (https://cran.r-project.org/)

R packages:

*R version 4.2.0 (2022-04-22)*
```R
install.packages("tidyverse")
install.packages("readr")
install.packages("pheatmap")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```

* [Part one - Database preparation](/01_Database_preparation.md)
* [Part two - Search for PK and LRR containing proteins](/02_Search_for_PK_and_LRR_containing_proteins.md)
* [Part three - Phylogenetic tree construction and analyses](/03_Phylogenetic_tree_construction_and_analyses.md)
* [Part four - Data Visualization](/04_Data_Visualization.md)
