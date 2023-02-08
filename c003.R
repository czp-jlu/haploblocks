############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

#######################################
# Creation of haplotype blocks with   #
# r^2, fixed marker number, fixed     #
# chromosome window                   #  
# Creation of files with single SNPs  #
#######################################

# The whole base population is used for the creation of the
# haplotype blocks.

library ("SelectionTools")

dir.create("blockinfo")

st.input.dir  <- "input"
st.output.dir <- "output"

versions <- read.table("data/version_info.txt", header=T, sep="\t")

versions$nblocks <- NA # empty column for information on how many blocks there are

# Data set for the versions with all the SNP markers
####################################################

st.read.marker.data ("haploselekt-marker_preselected.mpo", format="m")
st.read.map("haploselekt-map.txt", # Read in linkage map, here in Mbp
            m.stretch=1, # Multiplier for map positions
            format="mcp", # Marker Chromosome Position / "cpms"
            skip=1) # Skip the first line

version <- versions[grep("all SNPs", versions$data),]
version <- substr(start=1, stop=3, version[1,1])

fname1 <- paste0("hs-marker-mpo-", version) # first part of file name for marker data
fname2 <- paste0("hs-marker-", version)     # first part of file name for map

st.write.marker.data(mfilename=fname1)
st.write.map(mfilename=fname2)

fname2 <- paste0("hs-marker-", version)     # first part of file name for map
fname3 <- paste0("hs-marker-mpo-", version)     # first part of file name for map

# Copy files to input directory:
file.copy(from = paste0("output/", fname1, ".mpo"), 
          to =   paste0("input/", fname3, ".mpo"), 
          overwrite =T)
file.copy(from = paste0("output/", fname1, ".mpo"), 
          to =   paste0("input/", fname2, ".npo"), # npo is not the correct file extension
          # but it helps to keep the following loops for the cross validation consistent.
          overwrite =T)
file.copy(from = paste0("output/", fname1, ".mmp"), 
          to =   paste0("input/", fname2, ".mmp"),
          overwrite =T)
file.copy(from = paste0("output/", fname1, ".mmp"), 
          to =   paste0("input/", fname3, ".mmp"),
          overwrite =T)

# Data set for the versions with a reduced set of SNP markers
#############################################################

version1 <- version # version with the full set of SNPS from above

version <- versions[grep("reduced set of SNPs", versions$data),]

version2 <- unique(substr(version$version, start=1, stop=3))

fname1 <- paste0("hs-marker-mpo-", version1) # first part of file name for marker data

for (VERSION in version2) {
  
    fname2 <- paste0("hs-marker-", VERSION)     # first part of file name for map
    fname3 <- paste0("hs-marker-mpo-", VERSION)     # first part of file name for map

# Copy files to input directory:
file.copy(from = paste0("output/", fname1, ".mpo"), 
          to =   paste0("input/", fname3, ".mpo"), 
          overwrite =T)
file.copy(from = paste0("output/", fname1, ".mpo"), 
          to =   paste0("input/", fname2, ".npo"), # npo is not the correct file extension
          # but it helps to keep the following loops for the cross validation consistent.
          overwrite =T)
file.copy(from = paste0("output/", fname1, ".mmp"), 
          to =   paste0("input/", fname2, ".mmp"),
          overwrite =T)
file.copy(from = paste0("output/", fname1, ".mmp"), 
          to =   paste0("input/", fname3, ".mmp"),
          overwrite =T)

}

# Data sets for the versions with haplotype blocks with r2
##########################################################

version <- versions[grep("hb with r2", versions$data),]
version1 <- unique(substr(version$version, start=1, stop=3))
version2 <- version$version[seq(1,nrow(version), 3)]

for (VERSION in 1:length(version1)) {
  
 # versions and parameters
  vv1 <- version1[VERSION]
  vv2 <- version2[VERSION]
  par1 <- as.numeric(version$parameters1[version$version==vv2])
  par2 <- version$parameters2[version$version==vv2]
  par3 <- as.numeric(version$parameters3[version$version==vv2])
  
  # file names
  fname1 <- paste0("hs-marker-mpo-", vv1) # first part of file name for marker data
  fname2 <- paste0("hs-marker-", vv1)     # first part of file name for map
  fname3 <- paste0("blockinfo/hs-block-", vv1, ".txt") # file for info on haplotype blocks
  fname4 <- paste0("blockinfo/hs-block-r-", vv1, ".txt") # file for info on haplotype block alleles
    
  st.copy.marker.data(target.data.set="p1",
                      source.data.set="default")
  
  # LD:
  ld <- st.calc.ld ( ld.measure="r2",
                     data.set="p1" )
  
  # haplotype blocks:
  h <- st.def.hblocks ( ld.threshold = par1, # Minimum LD 
                        ld.criterion = par2, # flanking/average
                        tolerance = par3,    # How many markers are allowed that do not fit the criterium
                        data.set="p1" ) 

  # write info on blocks
  write.table(h, file=fname3, quote=F)
  
  # haplotype alleles
  r <- st.recode.hil ( data.set="p1" )
  write.table(r, file=fname4, quote=F)
  
  print(vv1)
  print(table(r$AlleleNr)) # summary of how many alleles there are per haplotype block
  
  versions$nblocks[versions$version==vv2] <- length(unique(r$Block))
  
  # marker data in matrix format:
  st.write.marker.data(format="m", # Matrix
                       mfilename = fname1,
                       data.set="p1")
  
  # marker data in NTSys format:
  st.write.marker.data(format="n", # NTSys
                       nfilename = fname2,
                       data.set="p1")
  
  # Both files are written into the output folder.
  
  # Format NTSys creates a file with one row for each marker
  # and one column for each individual. "marker" in this case
  # means the combination marker/allele, so each allele of the 
  # haplotype block gets its own row.
  
  # So for three alleles:
  
  #    marker
  # i1 1
  # i2 3
  # i3 2
  
  # Individual 1 has allele 1, individual 2 has allele 3 and 
  # individual 3 has allele 2.
  
  # The re-parametrized matrix looks like this:
  
  #    m.1 m.2 m.3
  # i1 1   0   0
  # i2 0   0   1
  # i3 0   1   0
  
  # In NTSys format, the rows and columns change.
  
  # The map needs to be adjusted because it only contains
  # one position per marker and not per marker/allele 
  # combination.
  # Moreover, the names of the genotypes were changed.
  
  map           <- read.table(paste0("output/", fname2, ".nmp"))
  marker.allele <- read.table(paste0("output/", fname2, ".npo"))
  
  marker.allele[1:5,1:5]
  
  # SelectionTools expects 1/2 encoding, not 0/1.
  # 0 is interpreted as missing.
  marker.allele <- marker.allele + 1
  
  # Change names of the genotypes back:
  colnames(marker.allele) <- substr(colnames(marker.allele), 
                                    start=1, stop=11)
  colnames(marker.allele) <- gsub("x", "-", colnames(marker.allele))
  marker.allele[1:5,1:5]
  
  write.table(marker.allele, file=paste0("output/", fname2, ".npo"),
              quote=F)
  
  # Change map:
  marker.allele <- rownames(marker.allele)
  
  head(map)
  head(marker.allele)
  
  markers <- substr(marker.allele, start=1, stop=7)
  head(markers)
  
  new.map <- data.frame()
  
  for (mm in markers) {
    new.map <- rbind(new.map, map[map$V3==mm,])
  }
  
  new.map$V3 <- marker.allele
  
  head(new.map)
  
  write.table(new.map, file=paste0("output/", fname2, ".nmp"),
              quote=F, row.names=F, col.names=F)
  
  # Copy files to input directory:
  file.copy(from = paste0("output/", fname1, ".mpo"), 
            to =   paste0("input/", fname1, ".mpo"),
            overwrite =T)
  file.copy(from = paste0("output/", fname1, ".mmp"), 
            to =   paste0("input/", fname1, ".mmp"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".npo"), 
            to =   paste0("input/", fname2, ".npo"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".nmp"), 
            to =   paste0("input/", fname2, ".nmp"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".ngp"), 
            to =   paste0("input/", fname2, ".ngp"),
            overwrite =T)
  
} # for


# Data sets for the versions with haplotype blocks with fixed numbers of markers
################################################################################

version <- versions[grep("hb with fixed marker number", versions$data),]
version1 <- unique(substr(version$version, start=1, stop=3))
version2 <- version$version[seq(1,nrow(version), 3)]

dir.create("blockinfo")

for (VERSION in 1:length(version1)) {
  
  # versions and parameters
  vv1 <- version1[VERSION]
  vv2 <- version2[VERSION]
  par1 <- as.numeric(version$parameters1[version$version==vv2])
  
  # file names
  fname1 <- paste0("hs-marker-mpo-", vv1) # first part of file name for marker data
  fname2 <- paste0("hs-marker-", vv1)     # first part of file name for map
  fname3 <- paste0("blockinfo/hs-block-", vv1, ".txt") # file for info on haplotype blocks
  fname4 <- paste0("blockinfo/hs-block-r-", vv1, ".txt") # file for info on haplotype block alleles
  
  st.copy.marker.data(target.data.set="p1",
                      source.data.set="default")

  # haplotype blocks:
  h <- st.def.hblocks ( hap = par1,     # number of markers per block
                        hap.unit = 1,   # blocks based on fixed marker number per block
                        data.set="p1" ) 

  # write info on blocks
  write.table(h, file=fname3, quote=F)
  
  # haplotype alleles
  r <- st.recode.hil ( data.set="p1" )
  write.table(r, file=fname4, quote=F)
  
  print(vv1)
  print(table(r$AlleleNr)) # summary of how many alleles there are per haplotype block
  
  versions$nblocks[versions$version==vv2] <- length(unique(r$Block))
  
  # marker data in matrix format:
  st.write.marker.data(format="m", # Matrix
                       mfilename = fname1,
                       data.set="p1")
  
  # marker data in NTSys format:
  st.write.marker.data(format="n", # NTSys
                       nfilename = fname2,
                       data.set="p1")
  
  map           <- read.table(paste0("output/", fname2, ".nmp"))
  marker.allele <- read.table(paste0("output/", fname2, ".npo"))
  
  marker.allele[1:5,1:5]
  
  # SelectionTools expects 1/2 encoding, not 0/1.
  # 0 is interpreted as missing.
  marker.allele <- marker.allele + 1
  
  # Change names of the genotypes back:
  colnames(marker.allele) <- substr(colnames(marker.allele), 
                                    start=1, stop=11)
  colnames(marker.allele) <- gsub("x", "-", colnames(marker.allele))
  marker.allele[1:5,1:5]
  
  write.table(marker.allele, file=paste0("output/", fname2, ".npo"),
              quote=F)
  
  # Change map:
  marker.allele <- rownames(marker.allele)
  
  head(map)
  head(marker.allele)
  
  markers <- substr(marker.allele, start=1, stop=7)
  head(markers)
  
  new.map <- data.frame()
  
  for (mm in markers) {
    new.map <- rbind(new.map, map[map$V3==mm,])
  }
  
  new.map$V3 <- marker.allele
  
  head(new.map)
  
  write.table(new.map, file=paste0("output/", fname2, ".nmp"),
              quote=F, row.names=F, col.names=F)
  
  # Copy files to input directory:
  file.copy(from = paste0("output/", fname1, ".mpo"), 
            to =   paste0("input/", fname1, ".mpo"),
            overwrite =T)
  file.copy(from = paste0("output/", fname1, ".mmp"), 
            to =   paste0("input/", fname1, ".mmp"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".npo"), 
            to =   paste0("input/", fname2, ".npo"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".nmp"), 
            to =   paste0("input/", fname2, ".nmp"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".ngp"), 
            to =   paste0("input/", fname2, ".ngp"),
            overwrite =T)
  
} # for


# Data sets for the versions with haplotype blocks with fixed chromosome windows
################################################################################

# When the map was created, it was assumed that 1 Mbp equals 1 cM.

version <- versions[grep("hb with fixed chromosome window", versions$data),]
version1 <- unique(substr(version$version, start=1, stop=3))
version2 <- version$version[seq(1,nrow(version), 3)]

for (VERSION in 1:length(version1)) {
  
  # versions and parameters
  vv1 <- version1[VERSION]
  vv2 <- version2[VERSION]
  par1 <- as.numeric(version$parameters1[version$version==vv2])
  
  # file names
  fname1 <- paste0("hs-marker-mpo-", vv1) # first part of file name for marker data
  fname2 <- paste0("hs-marker-", vv1)     # first part of file name for map
  fname3 <- paste0("blockinfo/hs-block-", vv1, ".txt") # file for info on haplotype blocks
  fname4 <- paste0("blockinfo/hs-block-r-", vv1, ".txt") # file for info on haplotype block alleles
  
  st.copy.marker.data(target.data.set="p1",
                      source.data.set="default")
  
  # haplotype blocks:
  h <- st.def.hblocks ( hap = par1,     # genetic distance in cM
                        hap.unit = 2,   # blocks based on fixed chromosome windows
                        data.set="p1" ) 

  # write info on blocks
  write.table(h, file=fname3, quote=F)
  
  # haplotype alleles
  r <- st.recode.hil ( data.set="p1" )
  write.table(r, file=fname4, quote=F)
  
  print(vv1)
  print(table(r$AlleleNr)) # summary of how many alleles there are per haplotype block
  
  versions$nblocks[versions$version==vv2] <- length(unique(r$Block))
  
  # marker data in matrix format:
  st.write.marker.data(format="m", # Matrix
                       mfilename = fname1,
                       data.set="p1")
  
  # marker data in NTSys format:
  st.write.marker.data(format="n", # NTSys
                       nfilename = fname2,
                       data.set="p1")
  
  map           <- read.table(paste0("output/", fname2, ".nmp"))
  marker.allele <- read.table(paste0("output/", fname2, ".npo"))
  
  marker.allele[1:5,1:5]
  
  # SelectionTools expects 1/2 encoding, not 0/1.
  # 0 is interpreted as missing.
  marker.allele <- marker.allele + 1
  
  # Change names of the genotypes back:
  colnames(marker.allele) <- substr(colnames(marker.allele), 
                                    start=1, stop=11)
  colnames(marker.allele) <- gsub("x", "-", colnames(marker.allele))
  marker.allele[1:5,1:5]
  
  write.table(marker.allele, file=paste0("output/", fname2, ".npo"),
              quote=F)
  
  # Change map:
  marker.allele <- rownames(marker.allele)
  
  head(map)
  head(marker.allele)
  
  markers <- substr(marker.allele, start=1, stop=7)
  head(markers)
  
  new.map <- data.frame()
  
  for (mm in markers) {
    new.map <- rbind(new.map, map[map$V3==mm,])
  }
  
  new.map$V3 <- marker.allele
  
  head(new.map)
  
  write.table(new.map, file=paste0("output/", fname2, ".nmp"),
              quote=F, row.names=F, col.names=F)
  
  # Copy files to input directory:
  file.copy(from = paste0("output/", fname1, ".mpo"), 
            to =   paste0("input/", fname1, ".mpo"),
            overwrite =T)
  file.copy(from = paste0("output/", fname1, ".mmp"), 
            to =   paste0("input/", fname1, ".mmp"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".npo"), 
            to =   paste0("input/", fname2, ".npo"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".nmp"), 
            to =   paste0("input/", fname2, ".nmp"),
            overwrite =T)
  file.copy(from = paste0("output/", fname2, ".ngp"), 
            to =   paste0("input/", fname2, ".ngp"),
            overwrite =T)
  
} # for

# Write the information on the marker numbers
#############################################

write.table(versions, "data/versions_nblocks_r2.txt", quote=F, sep="\t", row.names=F)
