HYbrid-clustering for Single-cell Categorization (HYSC)

---------------------------------------------------------------------------------------------------------
To run standalone HYSC without MATLAB:

Download and install the MATLAB Runtime installer: Windows 64-bit version of the MATLAB Runtime for R2017a,  
from the MathWorks Web site by navigating to

   http://www.mathworks.com/products/compiler/mcr/index.html


---------------------------------------------------------------------------------------------------------
Steps to run demo HYSC:
1)  Unzip the file 'Su Fetal liver.rar', open the 'Su Fetal liver.txt' in Excel and save the file into 'Su Fetal liver.xlsx'.
2)  Type in 'HYSC_gui' in the MATLAB commond window or double click the 'HYSC.exe' icon.
3)  Select the data file, 'Su Fetal liver.xlsx', and you MUST click 'Advanced settings' button first to set up parameters.
4)  Check 'Log Transform', 'Gene filtering' and 'Cell Anno', set 'Max HYSC Layers' to '1' and click 'Run'.
5)  Type in 'ViSC_gui' in the MATLAB commond window or double click the 'ViSC.exe' icon.
6)  Select the 'HYSC Output.xlsx' file, and type in the gene symbol of interest. The gene symbol must be contained in the sheet 3 of 'HYSC Output.xlsx'.
7)  Click 'Run'.

---------------------------------------------------------------------------------------------------------
Data file preparation

The expression data should be stored in a .xlsx file with the first column being gene symbols and the first row cell IDs. 
The cell IDs could be either indexs of cells or some a priori information about cell categories, such as time slots. In the latter case, 'Cell Anno' should be checked.
The dropouts should be zero-inflated.

---------------------------------------------------------------------------------------------------------
Input parameters:

Log Transform            - check to perform log transform.
Gene Filtering           - check to perform gene filtering using the parameters specified by 'Filtering No. cells' and 'Filtering Var' in 'AdvancedSetting'.
Cell Anno                - check to incorporate additional information in result plot, e.g., time information of cells, 
                           By default, []. In this case, the annotation of cells will be cluster indexs.
Max HYSC Layers          - the maximum number of HYSC layers.

Advanced Setting:
Cores                    - the number of cores used in parallel computation. By default, cores = 8.
Min Cluster Size         - the minimum number of samples (cells) in each cluster. By default, 10.
Dim Clustering           - the number of principle components adopted by k-means clustering. By default, 5.
r2-cutoff                - OPLSDA parameters, r2 cutoff value to identify discriminatory variables. By default, 0.5.
p-cutoff                 - OPLSDA parameters, p cutoff value to identify discriminatory variables. By default, 0.05 after bonferroni correction.
Perplexity               - perplexity parameter of tSNE. By default, 30.
Filtering No. cells      - the mimimum number of cells in which a gene must express. By default, >2.
Filtering Var            - the mimimum variance of a gene. By default, >0.5.

---------------------------------------------------------------------------------------------------------
Output:
1) 'HYSC Output.xlsx':
Sheet 1 - the cell categorization indexs
Sheet 2 - tSNE scores for scatter plot
Sheet 3 - gene expression data

2) 'HYSC gene markers.xlsx': each sheet contains the gene markers of a cell cluster.

3) A tSNE scatter plot with the cells represented by the indexs of clusters they are assigned to.

4) A gene expression heatmap. The rows are selected gene markers of individual cell clusters and the cell cluster indexs are provided at the bottom of the plot.
The dots above the heatmap indicate cell annotations, and the colors representing different annotation indexs are provided in the colorbar on the right.

---------------------------------------------------------------------------------------------------------
Tips:

Before each run, please click 'Cancel' button to clear memnory.

---------------------------------------------------------------------------------------------------------


Author Jie hao, Xin Zou, SJTU, China

Copyright Jie Hao, Xin Zou

2018
