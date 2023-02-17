# haploblocks
Genomic prediction with haplotype blocks in wheat

c001.R: Information on the different versions
c002.R: Filter the marker data, creation of splits into training and validation set, random selection of SNP markers for version with reduced marker numbers 
(requires package SelectionTools to be installed from http://population-genetics.uni-giessen.de/~software/)
c003.R: Creation of haplotype blocks with r^2, fixed marker number, fixed chromosome window, creation of files with single SNPs 
(requires package SelectionTools to be installed from http://population-genetics.uni-giessen.de/~software/)
c004.R: Creation of haplotype blocks with HaploBlocker
(requires package HapoloBlocker from Pook et al. 2019)
c005.R: Cross-validation with RR-BLUP and RMLA 
(requires package SelectionTools to be installed from http://population-genetics.uni-giessen.de/~software/)
c006.R: Format the files for GVCHAP 
c007.R: Cross-validation with GVCHAP
(requires software GVCHAP from Prakapenka et al. 2020)
c008.R: Combine results

Note that for the code to function correctly, two additional (empty) folders called "data" and "output" need to be created in the working directory.
