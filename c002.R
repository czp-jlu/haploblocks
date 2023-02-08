############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

############################################################
# Filter the marker data                                   #
# Creation of splits into training and validation set      #
# Random selection of SNP markers for version with reduced #
# marker numbers                                           #
############################################################

# Create a filtered marker data set
###################################

library("SelectionTools")
st.input.dir <- "input"
st.output.dir <- "output"

st.read.marker.data ("haploselekt-marker.mpo", format="m")
st.read.map("haploselekt-map.txt", # Read in linkage map, here in Mbp
            m.stretch=1, # Multiplier for map positions
            format="mcp", # Marker Chromosome Position / "cpms"
            skip=1) # Skip the first line

# Preprocessing of data
st.restrict.marker.data(NoAll.MAX=2) # Maximum number of alleles
st.restrict.marker.data(MaMis.MAX=0.1) # Max missing at a marker
st.restrict.marker.data(ExHet.MIN=0.05) # Minimum gene diversity
st.restrict.marker.data(InMis.MAX=0.1) # Max missing per individual (e.g. deletion of water controls)

st.write.marker.data(mfilename="haploselekt-marker_preselected")
file.copy(from="output/haploselekt-marker_preselected.mpo",
          to="input/haploselekt-marker_preselected.mpo",
          overwrite=T)

# Assignment of lines to training and validation set
####################################################

# The same partitioning of the data set into training and validation set is used
# for each version. The partitioning is stored in a separate file so that
# it is easy to add additional versions without having to run the whole thing
# again to ensure comparability. Additionally, it can be tracked which individuals
# were included in which set in each of the runs.

nn <- 361     # 361 lines in the filtered data set
NRUNS <- 1000 # 1000 cv runs
pp <- 72      # number of individuals in the training/validation set (289/72)

rr <- data.frame(matrix(NA, nrow=361, ncol=1000)) # data frame for the results

set.seed(121)
for (RUN in 1:NRUNS){
  no.p.set <- sample(1:nn, size=pp)
  no.t.set <- c(1:nn)[!(c(1:nn) %in% no.p.set)]
  
  rr[no.p.set, RUN] <- "P" # prediction/validation set
  rr[no.t.set, RUN] <- "T" # training set
}

write.table(rr, "data/assignment_sets.txt", quote=F, row.names=F)


# Selection of SNPs for the versions with reduced numbers of SNP markers
# and creation of marker data files
########################################################################

versions <- read.table("data/version_info.txt", header=T, sep="\t")
versions.reduced <- versions[versions$data=="reduced set of SNPs",]
nn <- unique(versions.reduced$parameters1)

marker <- read.table("input/haploselekt-marker_preselected.mpo")
nm <- nrow(marker)

set.seed(41)
for (NN in nn) {
  
  rr <- data.frame(matrix(FALSE, nrow=nm, ncol=1000))# data frame for the results
  
  for (RR in 1:1000){
        ss <- sample(1:nm, size=NN, replace=F)
    rr[ss, RR] <- TRUE
  }
  
  write.table(rr, file=paste0("data/assignment_snps_", NN, ".txt"), quote=F, row.names=F)

}

