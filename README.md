# Phylogenetic-Dirichlet-multinomial-model-for-microbiome-data
These R codes are for the paper "Phylogenetic Dirichlet-multinomial model for microbiome data" by Yunfan Tang, Ma Li and Dan L. Nicolae. 

The AG raw data (OTU table, metadata and phylogenetic tree) in May 16 2016 batch can be downloaded from

ftp://ftp.microbio.me/AmericanGut/ag-May-16-2016/03-otus/notrim/gg-13_8-97-percent/97_otus.tree
ftp://ftp.microbio.me/AmericanGut/ag-May-16-2016/03-otus/notrim/gg-13_8-97-percent/otu_table.biom
ftp://ftp.microbio.me/AmericanGut/ag-May-16-2016/04-meta/ag-cleaned.txt

Program descriptions:

TailBound.R: Implements the union bound algorithm (Section 4.2)

AmericanGut.R: Two-sample testing (Section 5.1).

LRT.R: Likelihood ratio test for DM vs PhyloDM (Section 5.2). 

ROC.R: Produces ROC curve and power curve by simulation (Section 5.3).

Read_hdf5biom.R: A publically available script to read sparse biom files.

SingleNodeOptimize.R: Source file for LRT.R.

TwoSampleNode.R: Source file for AmericanGut.R, LRT.R and ROC.R.





