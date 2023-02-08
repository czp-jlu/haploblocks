############################################
# HaploSelekt: Cross-validation study      #
# Genomic prediction with haplotype blocks #
############################################

#######################################
# Haplotype blocks with HaploBlocker  #
#######################################

library("HaploBlocker")
dir.create("haploblocker")
dir.create("haploblocker/chromosomes")
dir.create("haploblocker/size")

marker <- read.table("input/haploselekt-marker_preselected.mpo",
                     header=T)
dim(marker)
marker[1:5,1:5]

map <- read.table("input/haploselekt-map.txt", header=T)
head(map)
map <- map[order(map$pos),]
map <- map[order(map$chrom),]
head(map)
map <- map[map$name %in% rownames(marker),]

genotypes <- colnames(marker)

#####################################################
# cut up the marker data into one file per chromosome
#####################################################

cut.chrom <- function(nchr, marker, map){
    marker.names <- map$name[map$chrom == nchr]
    markers <- marker[marker.names,]
    write.table(markers,
                paste0("haploblocker/chromosomes/chrom_", nchr, ".txt"))
    print(nrow(markers))
}

for (NCHR in 1:21) {
    cut.chrom(nchr=NCHR, marker=marker, map=map)
}

######################################################
# define versions for HaploBlocker
######################################################

versions <- data.frame(version = paste0("v", c(46, 60:68)),
                       adaptive = c(TRUE, TRUE, rep(FALSE,8)),
                       coverage = as.numeric(c(0.9,0.9, rep(c(0.9, 0.95), each=2, length.out=8))),
                       window = c(5, 5, rep(c(5, 12), each=4, length.out=8)),
                       overlap = c(FALSE, TRUE, rep(c(TRUE, FALSE), times=4)))
                       
block <- function(nchr, version, haploblocks) {
    marker <- read.table(paste0("haploblocker/chromosomes/chrom_", nchr, ".txt"), header=T)
    blocklist <- block_calculation(marker,
                                   adaptive_mode = versions$adaptive[versions$version==version],
                                   target_coverage = versions$coverage[versions$version==version],
                                   window_size = versions$window[versions$version==version],
                                   overlap_remove = versions$overlap[versions$version==version])
    bm <- block_matrix_construction(blocklist)
    haploblocks <- rbind(haploblocks, bm)
    return(haploblocks)
}

size_block <- function(nchr, version, blocksize) {
    marker <- read.table(paste0("haploblocker/chromosomes/chrom_", nchr, ".txt"), header=T)
    blocklist <- block_calculation(marker,
                                   adaptive_mode = versions$adaptive[versions$version==version],
                                   target_coverage = versions$coverage[versions$version==version],
                                   window_size = versions$window[versions$version==version],
                                   overlap_remove = versions$overlap[versions$version==version])
	bl <- blocklist_startend(blocklist, type="snp")
	blocksize <- rbind(blocksize, bl)
	return(blocksize)
}

for (VERSION in versions$version) {

    haploblocks <- NULL
	blocksize <- NULL

    for (NCHR in 1:21) {

        haploblocks <- block(nchr=NCHR, version=VERSION, haploblocks=haploblocks)
		blocksize <- size_block(nchr=NCHR, version=VERSION, blocksize=blocksize)        

        }

    haploblocks <- haploblocks + 1
    nr <- nrow(haploblocks)
    rownames(haploblocks) <- sprintf("block00%i", 1:nr)
    colnames(haploblocks) <- genotypes
    write.table(haploblocks, file=paste0("input/hs-marker-", VERSION, ".npo"))
	write.table(blocksize, file=paste0("haploblocker/size/blocksize-", VERSION, ".dta"))
        
}



