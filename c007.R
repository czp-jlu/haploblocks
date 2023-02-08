############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

################################
# Cross-validation with GVCHAP #
################################

library("sqldf")

# The executable GVCHAP program needs to be in
# the folder gvchap.

###################################

# Number of cross-validation runs
NRUNS <- 1 

# Output folder
dir.create("gvchap/out")

# Data base
dbfile <- "sqldb/hs_cv_gvchap.sqldb"  # database
rfile <- "results"                    # file for the results in the database

# Traits
traits <- c("yld", "protein", "starch", "hlw", "length", "br", "fus1", "yr", "md", "septoria", "proteinyld")

# Assignments of lines to the training and validation set
sets <- read.table("data/assignment_sets.txt", header=T)

rr <- data.frame()        # prepare data frame for the results
nn <- 361                 # number of genotypes in the preselected data set

###################################
# Different versions:

versions.full <- read.table("data/version_info.txt", header=T, sep="\t")
versions <- unique(versions.full$version2)
versions <- versions[c(1:45,47:59)]

EST <- "c" # c is GVCHAP

###################################
# Prepare the parameter files:

for (VERSION in versions) {
  
  for (TRT in traits) {

    # Read in the file for the version
    fname.par1 <- paste0("gvchap/parameter/parameter_", VERSION, "_", TRT, ".txt")
    fname.par2 <- paste0("gvchap/parameter_", VERSION, "_", TRT, ".txt")
    
    code <- readLines(fname.par1)
    
    dir.out <- paste0("gvchap/out")
    pheno.file <- paste0("gvchap/pheno_hs_", TRT, "_gvchap_run.txt")
    
    l36 <- paste0("out_greml ", dir.out, "/greml.txt")
    l37 <- paste0("out_gblup ", dir.out, "/gblup.txt")
    l38 <- paste0("output_fixed_effect ", dir.out, "/fixed_effect.txt")
    l39 <- paste0("output_hap_heritabilities ", dir.out, "/hap_heritabilites.snpe")
    l40 <- paste0("output_mrk_effect ", dir.out, "/mrk_effect.snpe")
    l41 <- paste0("out_rel ", dir.out, "/relationship.txt")
    l45 <- paste0("phenotype ", pheno.file)
    
    code[36] <- l36
    code[37] <- l37
    code[38] <- l38
    code[39] <- l39
    code[40] <- l40
    code[41] <- l41
    code[45] <- l45

    for (LINE in 1:21) {
        ll <- c(11:31)[LINE] # SNP data files
        sspl <- strsplit(code[ll], split=" ")
        n.snps <- as.numeric(sspl[[1]][2])
        if(n.snps == 1016) {sspl[[1]][2] <- 1061 }
        new.line <- paste(sspl[[1]][1], sspl[[1]][2], paste0("gvchap/hs-marker-gvchap-chrom-", LINE, ".txt"), sep=" ")
        code[ll] <- new.line
        
        if (VERSION == "v01") {
            code <- code[1:46]

        } else
            {ll <- c(47:67)[LINE] # haplotype data files
             sspl <- strsplit(code[ll], split=" ")
             new.line <- paste(sspl[[1]][1], sspl[[1]][2], paste0("gvchap/", VERSION, "-hs-hblocks-gvchap-chrom-", LINE, ".txt"), sep=" ")
             code[ll] <- new.line}
    }

    
    # create file for the parameter file
    file.create(fname.par2, overwrite=T)
    
    # write code to new file
    writeLines(code, con = fname.par2) 
  
  }
}

###################################
# ACTUAL CROSS VALIDATION
###################################

# not for single SNPs
versions <- versions[-1]

for (RUN in NRUNS) {
  
  # prediction and training set
  no.p.set <- c(1:nn)[sets[,RUN]=="P"]
  no.t.set <- c(1:nn)[sets[,RUN]=="T"]
  
  for (VERSION in versions) {
    
    for (TRT in traits){
      
      # Read in phenotypic data
      #########################
      
      fname1 <- paste0("gvchap/pheno_hs_", TRT, "_gvchap.txt")
      fname2 <- paste0("gvchap/pheno_hs_", TRT, "_gvchap_run.txt")
      
      # Prepare phenotype file
      ########################
      
      dta <- read.table(fname1, header=T, sep="\t")
      dta2 <- dta # needed for the calculation of the correlation later
      dta$yy[no.p.set] <- -9999 # individuals in the training set are marked as "missing"
      
      write.table(dta, fname2, sep="\t", quote=F, row.names=F)
      
      # run GVCHAP
      ############
      
      parameter.file <- paste0("gvchap/parameter_", VERSION, "_", TRT, ".txt")
      cmd <- paste0("gvchap/gvchap ", parameter.file)
      system(cmd)         # run GVCHAP
      
      # read in predictions
      #####################
      
      gblup <- read.table(paste0("gvchap/out/gblup.txt"), header=T) # predictions
      mu <- read.table(paste0("gvchap/out/fixed_effect.txt"), header=T, skip=4)[1,2] # grand mean
      gv <- mu + gblup$GBLUP_G # GBLUP_G is the genetic effect due to SNPS and/or haploblocks
      cc <- cor(dta2$yy[no.p.set], gv[no.p.set]) # correlation observed/predicted

      rr.cv <- data.frame(version1 = VERSION,
                          version2 = paste0(VERSION, EST),
                          cor = cc,
                          trait = TRT,
                          estimation = EST,
                          run = RUN)
      
      # Save simulation results in data base
      conn <- dbConnect(RSQLite::SQLite(), dbfile)
      dbWriteTable(conn, rfile, rr.cv, append=TRUE)
      dbDisconnect(conn)
        
      } # for TRT
    } # for VERSION
} # for RUN
