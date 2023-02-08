############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

#######################################
# Format the files for GVCHAP         #
#######################################

library ("SelectionTools")

st.input.dir  <- "input"
st.output.dir <- "output"

###################################
# Different versions:

versions <- read.table("data/version_info.txt", header=T, sep="\t")
version <- unique(substr(versions$version, start=1, stop=3))

###################################

##############################
# Preparation of marker data #
# p. 11 of the user manual   #
# for GVCHAP                 #
##############################

# GVCHAP always includes the SNP markers in the model, even for
# the models with haplotype blocks.

# Read in marker data and map and reduce to the same number of markers:
mp <- read.table("input/haploselekt-map.txt", stringsAsFactors = T, header=T)
m <- read.table("input/haploselekt-marker_preselected.mpo", stringsAsFactors = T)
mp <- mp[mp$name %in% rownames(m),]
m <- m[rownames(m) %in% mp$name,]

# reorder m according to the map positions:
mp <- mp[order(mp$pos),]
mp <- mp[order(mp$chrom),]
m <- m[rownames(m)[match(mp$name, rownames(m))] , ]
m[1:5,1:5]
m["TA00075951CI04", 1:5]
head(mp)

# GVCHAP expects the following encoding (p. 11 of the user manual):
# 0 for homozygous A1/A1
# 1 for heterozygous A1/A2
# 2 for homozygous A2/A2

m <- t(m)
m <- cbind(rownames(m),m)
m[1:5,"TA00075951CI04"]
m[1:5,1:5]

m <- data.frame(m, stringsAsFactors = T)
str(m)

# Now recode the marker data:

for (ii in 2:ncol(m)){
  
  # If there are NO missing markers, there are two factor levels: one for 
  # A1A1 and one for A2A2 (it is arbitrary which allele is A1 and which
  # is A2 - A1A2 is not present because we have homozygous lines)
  if(length(levels(m[,ii]) == 2)) { 
    m[,ii] <- as.numeric(m[,ii])
    m[m[,ii] == 1,ii] <- 0 # allele 1 (allele 2 is already encoded as 2)
  }
  
  # If there ARE missing markers, there are three factor levels: 
  # level 1 for missing, level 2 for A1A1 and level 3 for A2A2
  if(length(levels(m[,ii]) == 3)) { 
    m[,ii] <- as.numeric(m[,ii])
    m[m[,ii] == 1,ii] <- 4 # missing markers must be encoded with
    # a different number than 0, 1, or 2
    m[m[,ii] == 2,ii] <- 0 # allele 1
    m[m[,ii] == 3,ii] <- 2 # allele 2
  } 
}

# Note that the code above works because the factor levels
# are ordered alphanumerically:
example.levels <- c("2/2", "1/1", "-1/-1", "4/4", "3/3")
sort(example.levels)
# -1/-1 is always at the first position.

m[1:5,1:5]

# Do the marker names now match the names in the map?
sum(colnames(m)[2:ncol(m)] == mp$name) # 16667? yes :)
# The first column name must be "ID":
colnames(m) <- c("ID", colnames(m)[2:ncol(m)])
m[1:5,1:5]
write.table(m, "gvchap/hs-marker-gvchap.dta", sep=" ", quote=F, row.names=F)

# There needs to be one file per chromosome.
# Wheat has 21 chromosomes (7 from each of its three genomes).
unique(mp$chrom)

for (CHR in unique(mp$chrom)) {
  markers <- as.character(mp$name[mp$chrom == CHR])
  marker.subset <- m[,c("ID", markers)]
  print(length(markers))
  
  fname <- paste0("gvchap/hs-marker-gvchap-chrom-", CHR, ".txt")
  write.table(marker.subset, fname, sep="\t", quote=F, row.names=F)
}

#####################################
# Preparation of the haplotype data #
# p. 11 of the user manual          #
# Parameter files, part 1           #
# p. 6 of the user manual           #
#####################################

dir.create("gvchap/parameter")

version <- version[c(2:45, 47:59)]
for (VERSION in version) {
  
  fname1 <- paste0("input/hs-marker-mpo-", VERSION, ".mpo") # file name for marker data
  fname2 <- paste0("input/hs-marker-mpo-", VERSION, ".mmp")     # file name for map
  
  # Read in marker data and map:
  mp <- read.table(fname2, stringsAsFactors = T)
  colnames(mp) <- c("chrom", "pos", "name", "V4")
  m <- read.table(fname1, stringsAsFactors = T)
  mp <- mp[mp$name %in% rownames(m),]
  m <- m[rownames(m) %in% mp$name,]
  
  # transpose:
  m2 <- t(m)
  
  m2[m2=="-1/-1"] <- "0/0" # replace -1/-1 by 0/0

  # There must be two columns per haplotype block, one for each allele.
  # Columns 1 and 2 belong to haplotype block 1, columns 3 and 4 belong
  # to haplotype block 2, etc.
  
  m3 <- NULL # initiate variable for the results

  for (ii in 1:ncol(m2)) {
    a1 <- substr(m2[,ii], start=1, stop=1)
    a2 <- substr(m2[,ii], start=3, stop=3)
    m3 <- cbind(m3, a1, a2)
  }
  m3[1:5, 1:10]
  m3[3,1:200]
  
  m3[m3=="0"] <- -9999 # replace missing allele by -9999
  m3[3,1:200]
  
  # The blocks should be the column names.
  m4 <- data.frame(m3)
  m4[1:5,1:5]
  cn <- colnames(m2)
  cn <- rep(cn, each=2)
  cn <- paste(cn, rep(c(1,2), length.out=length(cn)), sep=".")
  colnames(m4) <- cn
  m4[1:5,1:5]
  
  # There must be a first column with the IDs of the individuals.
  ID <- rownames(m2)
  
  m4 <- cbind(ID, m4)
  m4[1:5,1:5]
  
  # The object h contains the information on which chromosome the
  # blocks are.
  # There must be one file per chromosome.
  
  # information for the parameter file
  code <- readLines("gvchap/parameter.txt")
  ll <- 47:67 # lines in parameter.txt that contain the information on the haplotype block files
  
  for (CHR in unique(mp$chrom)) {
    blocks <- as.character(mp$name[mp$chrom == CHR])
    nb <- length(blocks) # number of blocks on the chromosome
    blocks <- rep(blocks, each=2)
    blocks <- paste(blocks, rep(c(1,2), length.out=length(blocks)), sep=".")
    block.subset <- m4[,c("ID",blocks)]

    fname3 <- paste0("gvchap/", VERSION, "-hs-hblocks-gvchap-chrom-", CHR, ".txt")
    write.table(block.subset, fname3, sep="\t", quote=F, row.names=F)
    
    changed.line <- paste0("geno_hap ", nb, " .\\", "hs-hblocks-gvchap-chrom-", CHR, ".txt")
    code[ll[CHR]] <- changed.line
  }
  
  # create file for the parameter file
  fname4 <- paste0("gvchap/parameter/parameter_", VERSION, ".txt")
  file.create(fname4, overwrite=T)
  
  # write code to new file
  writeLines(code, con = fname4) 
  
  # There is now one parameter file per version with the
  # necessary information on the haplotype blocks.

}

############################
# Phenotypic data          #
# p. 12 of the user manual #
############################

trt2 <- c("yld", "protein", "starch", "hlw", "length", "br", "fus1", "yr", "md", "septoria", "proteinyld")
for (jj in 1:11){
  trt <- trt2[jj]
  fname1 <- paste0("input/hs_", trt, ".txt")
  fname2 <- paste0("gvchap/pheno_hs_", trt, "_gvchap.txt")
  dta <- read.table(fname1, header=T, stringsAsFactors = T, sep=";")
  colnames(dta) <- c("ID", "yy")
  dta <- dta[dta$ID %in% m4$ID,]
  dta$factor <- rep(1, times=nrow(dta))
  write.table(dta, fname2, sep="\t", quote=F, row.names=F)
}

#####################################
# Parameter files, part 2           #
# p. 6 of the user manual           #
#####################################

# The parameter files, as they are so far, are incomplete.
# There is one template per version but there needs to be
# one per version-trait combination and they names of the 
# output files need to be changed as well.

ll <- c(36:41, 45) # lines that need to be changed

version <- c("v01", version)

# There is no parameter file for version 1 yet.
file.copy(from="gvchap/parameter.txt", to="gvchap/parameter/parameter_v01.txt",
          overwrite=T)

for (VERSION in version[1:length(version)]) {
  
  # Read in the file for the version
  fname5 <- paste0("gvchap/parameter/parameter_", VERSION, ".txt")
  code <- readLines(fname5)
  
  # create output directory
  dir.out <- paste0("gvchap/out_", VERSION)
  dir.create(dir.out)
  
  dir.out <- paste0("gvchap\\out_", VERSION)
  
  for (TRT in trt2) {
    
    pheno.file <- paste0("gvchap\\pheno_hs_", TRT, "_gvchap.txt")
    
    l36 <- paste0("out_greml ", ".\\", dir.out, "\\", "greml_", TRT, ".txt")
    l37 <- paste0("out_gblup ", ".\\", dir.out, "\\", "gblup_", TRT, ".txt")
    l38 <- paste0("output_fixed_effect ", ".\\", dir.out, "\\", "fixed_effect_", TRT, ".txt")
    l39 <- paste0("output_hap_heritabilities ", ".\\", dir.out, "\\", "hap_heritabilites_", TRT, ".snpe")
    l40 <- paste0("output_mrk_effect ", ".\\", dir.out, "\\", "mrk_effect_", TRT, ".snpe")
    l41 <- paste0("out_rel ", ".\\", dir.out, "\\", "relationship_", TRT, ".txt")
    l45 <- paste0("phenotype ", ".\\", pheno.file)
  
    code[36] <- l36
    code[37] <- l37
    code[38] <- l38
    code[39] <- l39
    code[40] <- l40
    code[41] <- l41
    code[45] <- l45
    
    # create file for the parameter file
    fname6 <- paste0("gvchap/parameter/parameter_", VERSION, "_", TRT, ".txt")
    file.create(fname6, overwrite=T)
    
    # write code to new file
    writeLines(code, con = fname6) 
    
  }
  
}

# For version 1 (with all SNPs but without haplotype blocks):

for (TRT in trt2) {

fname7 <- paste0("gvchap/parameter/parameter_v01_", TRT, ".txt")
code <- readLines(fname7)
    
    l02 <- "var_snp_a 0.1"
    l05 <- "var_hap_a 0"

    code[2] <- l02
    code[5] <- l05
    
    # create file for the parameter file
    file.create(fname7, overwrite=T)
    
    # write code to new file
    writeLines(code, con = fname7) 
    
  }



