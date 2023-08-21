### FST Plots ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/02_outliers/FST")

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Read fst output from ANGSD
fst <- read.table("./FST_1.2.fst.sw.txt", header = T)
fst <- read.table("./FST_3.4.fst.sw.txt", header = T)
fst <- read.table("./FST_5.6.fst.sw.txt", header = T)
fst <- read.table("./FST_7.8.fst.sw.txt", header = T)
fst <- read.table("./FST_9.10.fst.sw.txt", header = T)
fst <- read.table("./FST_11.12.fst.sw.txt", header = T)
fst <- read.table("./FST_13.14.fst.sw.txt", header = T)

# Format data table
fst.data <- fst %>%
        mutate(ID = NA, 
               CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", chr)), 
               BP = midPos) %>%
        arrange(CHR) %>%
        mutate(ID = paste0(rep("Win", length(rownames(fst))), rownames(fst))) %>%
        select(ID, CHR, BP, Fst) #%>% Azores
        # select(ID, CHR, BP, Fst)

## ID Outliers ----

# Calculate FST outlier threshold (mean FST + 4*SD of FST)
fst.tholds <- function(df, statistic){
        df.tmp <- as.vector(df %>% select(all_of(statistic)))
        thold = mean(pull(df.tmp)) + 4*sd(pull(df.tmp)) #mean FST + 4*SD(FST)
        thold2 = quantile(pull(df.tmp), 0.999) #99.9th percentile
        tholds <- as.vector(cbind(thold, thold2))
        return(tholds)
}

thold <- fst.tholds(fst.data, "Fst")

# Get list of outliers (based on thold2)
outliers <- fst %>% 
        filter(Fst>thold[[2]]) %>% 
        mutate(start = gsub(".+\\((\\d+).+\\(.+", "\\1", region),
               end = gsub(".+\\((\\d+),(\\d+).+\\(.+", "\\2", region),
               sort = as.numeric(gsub(".*_(\\d+)$", "\\1", chr))) %>% 
        arrange(sort) %>%
        select(chr, start, end)

write.table(outliers, file = "out_Guad.bed", row.names = F, col.names = F, quote = F, sep = "\t")

##----

## Manhattan Plots (modified from Burchardlab_Tutorials) ====

# Define Function ====
gg.manhattan <- function(df, statistic, threshold, hlight, col, ylims){
        
        # format df
        df.tmp <- df %>% 
                
                # Compute chromosome size
                group_by(CHR) %>% 
                summarise(chr_len=max(BP)) %>% 
                
                # Calculate cumulative position of each chromosome
                mutate(tot=cumsum(chr_len)-chr_len) %>%
                select(-chr_len) %>%
                
                # Add this info to the initial dataset
                left_join(df, ., by=c("CHR"="CHR")) %>%
                
                # Add a cumulative position of each SNP
                arrange(CHR, BP) %>%
                mutate(BPcum = BP+tot) %>%
                
                # Add highlight and annotation information
                mutate(is_highlight=ifelse(ID %in% hlight, "yes", "no")) %>%
                mutate(is_annotate=ifelse(df %>% select(statistic) > threshold, "yes", "no"))
        
        ggplot(df.tmp, aes_string(x= "BPcum", y= statistic)) +
                # Show all points
                geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
                scale_color_manual(values = rep(col, 22 )) +
                
                # remove space between plot area and x-axis
                scale_y_continuous(expand = c(0, 0), limits = ylims) +
                
                # add plot and axis titles
                labs(x = "Position across Scaffolds",
                     y = "FST") +
                
                # add genome-wide sig and sugg lines
                geom_hline(yintercept = threshold) +
                
                # Add highlighted points
                geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
                
                # Add label using ggrepel to avoid overlapping
                # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(ID), alpha=0.7), size=5, force=1.3) +
                
                # Custom the theme:
                theme_bw(base_size = 18) +
                theme( 
                        plot.title = element_text(hjust = 0.5),
                        legend.position="none",
                        panel.border = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank()
                )
}

# Run Function ====

#get window IDs for the outliers
com_win <- fst %>%
        mutate(ID = NA, 
               CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", chr))) %>%
        arrange(CHR) %>%
        mutate(ID = paste0(rep("Win", length(rownames(fst))), rownames(fst))) %>%
        filter(Fst > thold[[2]]) %>%
        select(ID, CHR, Fst)
hlight <- com_win$ID
hlight <- NULL

# plot
pdf("FST_Azores.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()

pdf("FST_Desertas.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()

pdf("FST_Selvagem.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()

pdf("FST_CV.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()

pdf("FST_SA.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()

pdf("FST_Gal.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()

pdf("FST_Guad.pdf", 24, 6)
gg.manhattan(fst.data, "Fst", thold[[2]], hlight, 
             col = rep(c("black", "grey"), 309), ylims = c(0, 1))
dev.off()
