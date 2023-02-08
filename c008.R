############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

###################################
# Data management for the results #
###################################

library("sqldf")

lf <- list.files("sqldb")

rfile <- "results" # file for the results in the database
results <- NULL

for (NFILE in lf) {

dbfile <- paste0("sqldb/", NFILE) # database

conn <- dbConnect(RSQLite::SQLite(), dbfile)
  
  command <- paste0("SELECT * FROM ", rfile)
  Q <- dbGetQuery(conn, command)
  results <- rbind(results, Q)


dbDisconnect(conn)
}

dim(results)
head(results)