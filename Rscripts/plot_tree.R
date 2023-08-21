## Plotting FastME Tree ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/01_exploratory/DistTree")

library(ape)
library(phangorn)

trees <- read.tree("petrel_tree", skip=2)
tree <- read.tree("petrel_tree")[[1]]

#view node numbers
plot(tree)
nodelabels()

#root tree
tree_root <- root(tree, node = 183, resolve.root = T)

#plot trees
pdf(file="fastme_tree.pdf", width = 12, height = 10)
plot.phylo(tree_root, type="phylogram", cex=0.4)
invisible(dev.off())

pdf(file="fastme_tree_radial.pdf", width = 12, height = 12)
plotBS(tree_root, trees, type="radial", cex=0.4, bs.adj = c(1.2,1.5), digits = 0, p = 70)
invisible(dev.off())

pdf(file="fastme_tree_bstp.pdf", width = 12, height = 12)
tree_bstp <- extract.clade(tree_root, 171)
plot.phylo(tree_bstp, type="phylogram", cex=0.5, use.edge.length = F)
invisible(dev.off())