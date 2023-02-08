############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################
##########################################
# Cross-validation with RR-BLUP and RMLA #
##########################################

library ("SelectionTools")
library("sqldf")

###################################

# Number of cross-validation runs
NRUNS <- 1 

# Input and output folders
st.input.dir  <- "input"
st.output.dir <- "output"

# Data base
dir.create("sqldb")
dbfile <-"sqldb/hs_cv.sqldb" # database
rfile <- "results"           # file for the results in the database

# Traits
traits <- c("yld", "protein", "starch", "hlw", "length", "br", "fus1", "yr", "md", "septoria", "proteinyld")

# Assignments of lines to the training and validation set
sets <- read.table("data/assignment_sets.txt", header=T)

rr <- data.frame()        # prepare data frame for the results
nn <- 361                 # number of genotypes in the preselected data set
gs.set.num.threads(4) 
st.set.info.level (0)     # suppress messages in SelectionTools

###################################
# Different versions:

versions.full <- read.table("data/version_info.txt", header=T, sep="\t")
versions <- unique(versions.full$version2)

versions.reduced <- versions.full[versions.full$data=="reduced set of SNPs",]
versions.reduced <- unique(versions.reduced$version2)

estimation <- c("a", "b") # a: RR-BLUP, b: RMLA

###################################
# ACTUAL CROSS VALIDATION
###################################

for (RUN in NRUNS) {
  
  # prediction and training set
  no.p.set <- c(1:nn)[sets[,RUN]=="P"]
  no.t.set <- c(1:nn)[sets[,RUN]=="T"]
  
  for (VERSION in versions) {
    
    print(paste0(VERSION, " RUN ", RUN))
    
    fname1 <- paste0(st.input.dir, "/hs-marker-", VERSION) # marker data
    fname2 <- paste0("hs-marker-", VERSION) # marker data
    
    marker <- read.table(paste0(fname1, ".npo"))
    
    if (VERSION %in% versions.reduced){
      
      nm <- as.numeric(unique(versions.full$parameters1[versions.full$version2==VERSION]))
      fname.snps <- paste0("data/assignment_snps_", nm, ".txt")
      selection <- read.table(fname.snps, header=T)
      selection <- selection[,RUN]
      
      marker <- marker[selection,]
      
    } # if VERSION in versions.reduced
    
    p.set <- marker[,no.p.set] # markers only for individuals in the prediction set
    t.set <- marker[,no.t.set] # markers only for individuals in the training set
    
    write.table(p.set, paste0(fname1, "-p.mpo"), quote=F)
    write.table(t.set, paste0(fname1, "-t.mpo"), quote=F)
    
    st.read.marker.data (paste0(fname2, "-p.mpo"), # load marker data for prediction set
                         format="m",
                         data.set = "pp")
    st.read.marker.data (paste0(fname2, "-t.mpo"), # load marker data for training set
                         format="m",
                         data.set = "tt")
    
    # Preprocessing of data
    st.restrict.marker.data(NoAll.MAX=2, data.set="tt") # Maximum number of alleles
    st.restrict.marker.data(MaMis.MAX=0.2, data.set="tt") # Max missing at a marker
    st.restrict.marker.data(ExHet.MIN=0.05, data.set="tt") # Minimum gene diversity      
    
    for (TRT in traits){
      
      # Read in phenotypic data
      fname3 <- paste0("hs_", TRT, ".txt")
      st.read.performance.data(fname3, data.set="tt") # phenotypic data for TS
      st.read.performance.data(fname3, data.set="pp") # phenotypic data for PS
      
      for (EST in estimation){
        
        if (EST == "a"){
          gs.esteff.rr ( method="BLUP", data.set="tt" )
        } else {
          gs.esteff.rr ( method="RMLA", data.set="tt" )
        } # is/else EST
          
        eff <- gs.return.effects(data.set="tt") # read out the marker effects
      
        phe <- gs.predict.genotypes ( training.set   = "tt", # prediction of phenotypes in the PS
                                    prediction.set = "pp")
      
        cc <- cor(phe$y, phe$yhat) # correlation between observed and predicted phenotypes in PS
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
        
        } # EST
      } # for TRT
    } # for VERSION
} # for RUN
