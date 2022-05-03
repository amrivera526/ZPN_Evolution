# ZPN_Evolution
Files associated with 2022 submitted publication on ZPN Evolution, that included phylogenetics and machine learning analysis.

We outline here the files and code used for both the phylogenetic and machine learning analysis.
Much of the code includes hardcoded file paths and names, that can be changed accordingly.

# Alignments Used
We are including an alignment used for both the phylogenetic inference and machine learning analyses.
Data was processed prior to phylogenetics to remove sequences greater than 90% identical.
For machine learning analyses we only removed completely identical sequences (100% identity cutoff).
We provide both versions of the alignment for repdroducibility. 

Our phylogenetic alignment is : AlgoAln_ZPN_PSI90c_MAFFTProMals_202112140138.fas

Our machine learning alignment is: AlgoAln_ZPN_PSI100C_MAFFTProMals_202203242030_Tent.fas

Our training and testing datasets used: (AlgoAln_ZPN_PSI100c_MAFFTProMals_Train_202203242127.fasta, AlgoAln_ZPN_PSI100c_MAFFTProMals_Test_202203242127.fasta)


# Code
For efficiency reasons our machine learning regression searches were submitted to a cluster.
Then we performed all further data analysis, model regularization, and visualization locally.
We provide one script for each as well as combined version in case it is more convenient for the reader to do this in one go.

The Machine Learning Classification script that was submitted to the cluster: MLClassify_ZPN_ModelSearch.py

The parsimonious model can be selected, analyzed and visualized through the following script: MLClassify_ZPN_AnalyzeRes.py

Combined version of scripts: MLClassify_ZPNPSI_AnalysisRes_CleanCombined_202204292011.py


# Trees
The maximum likelihood tree obtained is included in the newick format. 
It also has nodal support values in TBE (transfer bootstrap expectation):
BestTree_LGG_20211214TBE_Newick.tre


CONTACT:
albertomarcosrivera@gmail.com
