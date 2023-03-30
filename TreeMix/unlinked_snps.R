## Randomly select unlinked SNPs with SNPRelate
library("SNPRelate")

args <- commandArgs(T)

table <- args[1] # position file
pos <- read.table(table, header = T, row.names = NULL)
colnames(pos) <- c("scaffold", "position")

dist <- args[2] #minimum distance between SNPs
set.seed(12345)
pruned_sites <- snpgdsApartSelection(pos$scaffold, pos$position, min.dist = 10000)

unlinked_snps <- pos[pruned_sites, ]

write.table(unlinked_snps, args[3], sep = "\t", row.names = F)