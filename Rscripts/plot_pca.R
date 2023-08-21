## Plotting PCAngsd ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/01_exploratory/PCAngsd")
library(ggplot2)
library(ggtree)
library(dplyr)
library(ggrepel)

#Read in the covariance matrix
covmat <- as.matrix(read.table("PCAngsd.cov"))

# for each archipelago from pcadapt
# covmat <- as.matrix(read.table("./pcadapt/pair7.cov"))
# ind <- as.vector(read.table("./pcadapt/pair7.txt"))

#Read in specimen ID
ind <- as.vector(read.table("bampath.txt"))
ID <- gsub(".*BSTP_\\D{2,4}_(\\D*\\d+).*", "\\1", 
           gsub(".*[TOSP|AISP]_(\\d+\\D*?)_.*", "\\1", ind$V1))
site <- gsub(".*BSTP_(\\D+{2,4})_.*", "\\1", 
             (gsub(".*(TOSP).*", "\\1", 
                   gsub(".*(AISP).*", "\\1", ind$V1))))
ind <- as.data.frame(cbind(ID, site))

# ind <- ind[ind$site %in% c("AzV", "AzC", "DesC", "DesH", "SelC", "SelH", "Ber"),]
# ind <- ind[ind$site %in% c("AzH", "AzV", "AzC", "DesC", "DesH", "SelC", "SelH", "Ber", "AsH", "AsC", "SHC", "GalH", "GalC", "HI"),]

#Calculate eigenvalues
e_covmat <- eigen(covmat)

#PC Variance
var <- e_covmat$values/sum(e_covmat$values)

#Convert to data frame
PCs <- as.data.frame(e_covmat$vectors)
PCs <- cbind(ind, PCs)

PCs <- PCs %>% mutate(outgroup = case_when(site=="AISP"|site=="TOSP" ~ "Guad", 
               site=="CVH"|site=="CVC" ~ "CV", .default = "BSTP"))
# PCs <- PCs %>% mutate(outgroup = case_when(site=="AsH"|site=="AsC"|site=="SHC" ~ "South Atlantic", site=="GalH"|site=="GalC"|site=="HI" ~ "Pacific", .default = "BSTP"))

#Plot PCA

#all populations
ggplot(PCs, aes(V1, V2, colour = site)) + 
        geom_point(size = 2, aes(shape = outgroup)) +
        xlab(paste0("PC1 (", round(var[[1]]*100, 1), "% Variance)")) +
        ylab(paste0("PC2 (", round(var[[2]]*100, 1), "% Variance)")) +
        theme_bw() + theme(legend.position = "none")

#for individual archipelagos
ggplot(PCs, aes(V1, V2, colour = site)) + 
        geom_point(size = 2) +
        xlab(paste0("PC1 (", round(var[[1]]*100, 1), "% Variance)")) +
        ylab(paste0("PC2 (", round(var[[2]]*100, 1), "% Variance)")) +
        scale_color_manual(values = c("cyan3","brown1")) +
        theme_bw()

#North Atlantic PCA (separate run)
ggplot(PCs, aes(V1, V2, colour = site, label = ID)) + 
        geom_point(size = 2) + 
        geom_text_repel() +
        xlab(paste0("PC1 (", round(var[[1]]*100, 1), "% Variance)")) +
        ylab(paste0("PC2 (", round(var[[2]]*100, 1), "% Variance)")) +
        theme_bw()

#----

#Read in NJ tree
njtree <- ape::read.tree("PCAngsd.tree")

#Relabel tips
tips <- as.data.frame(cbind(rownames(ind), ind$site))
nj_named <- phylotools::sub.taxa.label(njtree, tips)

#Plot tree
ggtree(nj_named, layout = "circular", branch.length = "none") + geom_tiplab()
