# HYSC
a clustering algorithm for single-cell data


HYbrid-clustering for Single-cell Categorization (HYSC)

---------------------------------------------------------------------------------------------------------
Three steps to run HYSC:
1) Prepare a .xlsx data file containing expression data. The first column contains the gene symbles and the first row contains cell IDs. The dropouts should be zero-inflated.
2) Prepare the 'HYSC Input.txt' file containing all input parameters, including the file path of data.
3) Put the data file and the 'HYSC Input.txt' files into the same folder of HYSC.m.
4) Run HYSC by typing HYSC in MATLAB commond window.

In this demo, we included two testing datasets stored in the files 'Pollen.xlsx' and 'Treulein.xlsx'. 
The corresponding examples of HYSC parameter files are stored in 'HYSC Input for Pollen.txt' and 'HYSC Input for Treulien.txt'.

---------------------------------------------------------------------------------------------------------

Input parameters:
dataPath                 - The directory of .xlsx file containing expression data. The first column contains the gene symbles and the first row contains cell IDs. The dropouts should be zero-inflated.
ClusterNum               - the given number of clusters. If ClusterNum = 0, the algorithm automatically estimates the number of clusters
Nmax                     - the maximum number of cell clusters. By default, Nmax = 30
minClusterSize_number    - the minimum number of samples (cells) in each cluster. By default, minClusterSize_number = 10
dimClustering            - the number of principle components adopted by k-means clustering. By default, dimClustering = 30
r2cutoff                 - OPLSDA parameters, r2 cutoff value to identify discriminatory variables. By default 0.5
pcutoff                  - OPLSDA parameters, p cutoff value to identify discriminatory variables. By default 0.05 after bonferroni correction
perp                     - perplexity parameter of tSNE. By default, perp = 30
maxHYSCLayer             - the maximum number of HYSC layers. By default, maxHYSCLayer = 5
cores                    - the number of cores used in parallel computation. By default, cores = 8
tSNEScores               - the index of tSNE scores adopted for scores plot. By default, the first two components are illustrated.
geneFiltering_cellCounts - the mimimum number of cells in which a gene must express. By default >0.
geneFiltering_var        - the mimimum variance of a valid gene. By default >0.

All of the input parameter must be given in a .txt file, where the 'dataPath' is compulsory, and the other parameters are optional.
The name of a parameter should be separated from its value by ': '.

---------------------------------------------------------------------------------------------------------

Output:
The output is stored in an excel file with 3 sheets, 
Sheet 1 - the cell categorization indexs
Sheet 2 - the list of variable genes
Sheet 3 - the tSNE scores, cells vs tSNE components

Author Xin Zou, Jie hao, SJTU, China

Copyright Xin Zou, Jie Hao

2018
