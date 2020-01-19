# Specter: A novel tool for clustering large-scale single cell RNA-seq and multi-modal data

Overview
--------

This is a MATLAB Package of Specter. Specter is a novel computational method for clustering large-scale single cell RNA-seq data. In addition, Specter can combine the data from different measurements such as RNA measurements and the antibody-derived tags (collected on the same set of cells). Specter runs in linear time with respect to number of cells, thus it is very suitable for analyzing very big single cell RNA-seq data. On a data set comprising 2 million cells from mouse embryos, Specter requires only 26 minutes to compute the clusters. 

#### Specter enhances cell type identification

Undoubtly, Seurat (or graph-based Louvain algorithm) is a state-of-the-art method for de novo identification of cell type. We showed that Specter outperforms Seurat in term of accurary (ARI) and efficency on a large number of real scRNA-seq datasets. Moreover, Specter can highlight rare cell types in which Seurat might not find.

#### Using Specter with multi-modal data

On multi-modal dataset of 8,617 cord blood mononuclear cells (CBMCs), produced by CITE-seq (Stoeckius et al. (2017)). The authors measure the single cell transcriptomes alongside the expression of 11 surface proteins, whose levels are quantified with DNA-barcoded antibodie. When we combine transtriptomic counts (mRNA) and antibody-derived tags (ADT), Specter reveals CD4 and CD8 T cells, which are quite similar transcriptomically. 

![](img/multimodal.png)


Systems Requirements
--------------------

Specter is independent of operating systems because it is written in Matlab. Basic requirement for running Specter includes MATLAB and the Statistics and Machine Learning Toolbox. 

This Package has been tested using MATLAB 2018a on Linux. 


Usage
-----

Unzip the package. Change the current directory in Matlab to the folder containing the scripts.

This directory includes the following main scripts:
1) Specter_demo.m -- an example run of Specter on a specific dataset
2) preprocessing.m -- do preprocessing of the input data (if applicable) 
3) constructingNetwork.m -- construct a gene-gene co-expression network
4) estimatingscEnergy.m -- estimate the single cell energy (scEnergy) for each cell
5) ECA.m -- prinpipal component analysis of energy matrix
6) clusteringCells.m -- perform unsupervised clustering of single cell data
---------------------------
7) cluster_visualization.m -- visualize cells on two-dimensional space
8) lineage_visualization.m -- display cell lineage hierarchy with transition probability



Please refer to Specter_demo.m for instructions on how to use this code.
Input Data are gene expression data matrix (columns are genes (PCs) and rows are cells). 

If you have any problem or question using the package please contact do@genzentrum.lmu.de

