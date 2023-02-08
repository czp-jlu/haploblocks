############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

###############################
# Information on the versions #
###############################

n.v <- 68 # number of versions

versions <- data.frame ( version = paste0(rep(sprintf("v%02i", 1:n.v), each=3),
                                          rep(c("a", "b", "c"))),
                         version2 = paste0(rep(sprintf("v%02i", 1:n.v), each=3)),
                         data = rep(c("all SNPs", 
                                      "hb with r2", 
                                      "hb with fixed marker number",
                                      "hb with fixed chromosome window", 
                                      "hb with HaploBlocker",
                                      "reduced set of SNPs", 
                                      "hb with HaploBlocker"), 
                                    times=c(3,108,15,9,3,3*13,3*9)),
                         parameters1 = c(rep(c("-",
                                               as.character(seq(0.1,0.9, 0.1)), 
                                               c("5", "10", "20", "50", "100"), 
                                               c("5", "10", "20"), 
                                               "adaptive=T, coverage=0.9, windowsize=5, overlap=F",
                                               as.character(c(500, seq(1000, 12000, 1000)))), 
                                             times=c(3,rep(3*2*2, 9),rep(3, times=5),rep(3, times=3),3,rep(3, times=13))),
                                         rep(c("adaptive=T, coverage=0.9, windowsize=5, overlap.rm=T",
                                               "adaptive=F, coverage=0.9, windowsize=5, overlap.rm=T",
                                               "adaptive=F, coverage=0.9, windowsize=5, overlap.rm=F",
                                               "adaptive=F, coverage=0.95, windowsize=5, overlap.rm=T",
                                               "adaptive=F, coverage=0.95, windowsize=5, overlap.rm=F",
                                               "adaptive=F, coverage=0.9, windowsize=12, overlap.rm=T",
                                               "adaptive=F, coverage=0.9, windowsize=12, overlap.rm=F",
                                               "adaptive=F, coverage=0.95, windowsize=12, overlap.rm=T",
                                               "adaptive=F, coverage=0.95, windowsize=12, overlap.rm=F"), each=3)),
                         parameters2 = c(rep("-", times=3),
                                         rep(c("flanking", "average"), each=6, length.out=2*6*9),
                                         rep("-", times=n.v*3-(3+2*6*9))),
                         parameters3 = c(rep("-", times=3),
                                         rep(c("0", "1"), each=3, length.out=2*6*9),
                                         rep("-", times=n.v*3-(3+2*6*9))),                   
                         estimation = rep(c("BLUP", "RMLA", "GVCHAP"), length.out=3*n.v) )

write.table(versions, "data/version_info.txt", quote=F, row.names=F, sep="\t")

