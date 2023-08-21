### Outlier window plots with fdM, π, and Tajima's D ###
# Libraries ====
library(RcppCNPy) # Numpy library for R
library(dplyr)
library(ggplot2)
library(qvalue)
library(cowplot)

# 1. Setup FST Data ====
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/02_outliers/FST")
fst <- read.table(file.choose(), header = T)
fst <- fst %>% 
        mutate(ID = NA, CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", chr)), BP = midPos) %>%
        arrange(CHR) %>%
        mutate(ID = paste0(rep("Win", length(rownames(fst))), rownames(fst)))
fst.data <- fst %>%
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% #compute scaffold size
        mutate(tot=cumsum(chr_len)-chr_len) %>% #calculate cumulative position of each scaffold
        select(-chr_len) %>%
        left_join(fst, ., by=c("CHR"="CHR")) %>%  #add this info to the initial dataset
        arrange(CHR, BP) %>%
        mutate(BPcum = BP+tot, out = Fst) #add a cumulative position of each window

# 2. Setup PCAdapt Data ====
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/02_outliers/")
pcadapt_format <- function(pair_num) {
        # Load in npy file
        Zscores <<- npyLoad(paste0("./pcadapt/", pair_num, ".pcadapt.zscores.npy"))
        
        #Load in sites file
        sites <<- read.table(paste0("./pcadapt/", pair_num, "_pruned.txt"), 
                             header = T, colClasses = c("marker"="factor"))
        names(sites) <<- c("scaffold","marker")
        
        #Load in alt sites file
        alt_sites <<- read.table(paste0("./pcadapt/", pair_num, ".sites"))
}
pcadapt_format("pair1") #Azores
pcadapt_format("pair2") #Desertas
pcadapt_format("pair3") #Selvagem
source("./pcadapt/pcadapt_script_PC.R")
PCp=read.table("PC.pcadapt.pval.txt")
qval <- qvalue(PCp)$qvalues
sites$V1=alt_sites$V1
sig=sites%>%filter(V1==1)
pcadapt <- data.frame(marker=sig$marker, scaffold=sig$scaffold, q=qval$V1)
pcadapt <- pcadapt %>% mutate(CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", scaffold))) %>% arrange(CHR)
pcadapt.data <- pcadapt %>% 
        group_by(CHR) %>% 
        summarise(chr_len=max(marker)) %>% #compute scaffold size
        mutate(tot=cumsum(chr_len)-chr_len) %>% #calculate cumulative position of each scaffold
        select(-chr_len) %>%
        left_join(pcadapt, ., by=c("CHR"="CHR")) %>%  #add this info to the initial dataset
        mutate(BPcum = marker+tot, BP = marker) %>% #add a cumulative position of each SNP
        filter(q > 0) %>% #remove q-values of 0s
        mutate(logq = -log10(q), out = (logq - min(logq))/(max(logq) - min(logq)))
rm(alt_sites, d2, PCp, qval, sig, sites, zscores, Zscores, args, K)

# 3. Setup fdM Data ====
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/03_introgression")
fd <- read.table(file.choose(), header = T)
fd <- fd %>%
        mutate(ID = NA, 
               CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", chr)), 
               BP = windowStart + ((windowEnd - windowStart)/2)) %>%
        arrange(CHR) %>%
        mutate(ID = paste0(rep("Win", length(rownames(fd))), rownames(fd)))
fd.data <- fd %>% 
        group_by(CHR) %>% 
        summarise(chr_len=max(windowEnd)) %>% #compute chromosome size
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%  #calculate cumulative position of each chromosome
        left_join(fd, ., by=c("CHR"="CHR")) %>% #add this info to the initial dataset
        arrange(CHR, windowEnd) %>%
        mutate(BPcum = windowEnd+tot) #add a cumulative position of each window

# 4. Setup π and Tajima's D Data ====
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/04_overlaps")
cols <- c("Region", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud",
          "fayh", "zeng", "nSites")

#Prepare diversity statistics for the hot season population
nd_hot <- read.table(file.choose(), header = F, col.names = cols)
nd_hot <- nd_hot %>% transmute(ID = rownames(df), 
                               CHR = as.numeric(gsub(".*_(\\d+)", "\\1", Chr)), 
                               BP = WinCenter, tP, Tajima)
nd_hot.data <- nd_hot %>%
        group_by(CHR) %>%
        summarise(chr_len = max(BP)) %>% # compute scaffold size
        mutate(tot = cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>% # calculate cumulative position of each chromosome
        left_join(nd_hot, ., by = c("CHR" = "CHR")) %>% # add values to initial dataset
        arrange(CHR,BP) %>%
        mutate(BPcum = BP + tot) # add cumulative position of each window

#Prepare diversity statistics for the cool season population
nd_cool <- read.table(file.choose(), header = F, col.names = cols)
nd_cool <- nd_cool %>% transmute(ID = rownames(df), 
               CHR = as.numeric(gsub(".*_(\\d+)", "\\1", Chr)), 
               BP = WinCenter, tP, Tajima)
nd_cool.data <- nd_cool %>%
        group_by(CHR) %>%
        summarise(chr_len = max(BP)) %>% # compute scaffold size
        mutate(tot = cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>% # calculate cumulative position of each chromosome
        left_join(nd_cool, ., by = c("CHR" = "CHR")) %>% # add values to initial dataset
        arrange(CHR,BP) %>%
        mutate(BPcum = BP + tot) # add cumulative position of each window

# 5. Plot ====
outlier_plots <- function(outlier_data, fdM_data, scaffold_no, xlims, ylims, ylab){
        
        # Select scaffold of interest
        out.tmp <- outlier_data %>% filter(CHR == scaffold_no)
        fd.tmp <- fdM_data %>% filter(CHR == scaffold_no)
        
        # Outlier/fdM plot
        ggplot(out.tmp, aes(x=BP, y=out)) +
                
                # Add outliers
                geom_line(linewidth = 2, colour = "black") +
                
                # Add fdM values
                geom_line(data = fd.tmp, aes(x=BP, y=f_dM), alpha = 0.8, 
                          linewidth = 2, colour = "orange") +
                
                # remove space between plot area and x-axis
                scale_y_continuous(expand = c(0, 0), limits = ylims) +
                
                # set x-axis limits
                scale_x_continuous(limits = xlims) +
                
                # add plot and axis titles
                labs(x = paste0("Position on Scaffold ", scaffold_no), y = ylab) +
                
                # Custom theme:
                theme_bw(base_size = 18) +
                theme( 
                        plot.title = element_text(hjust = 0.5),
                        legend.position="none",
                        panel.border = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank()
                )
        
        
}
nd_plots <- function(nd_hot.data, nd_cool.data, statistic, scaffold_no, xlims, ylims, ylab){
        
        # Select scaffold of interest
        nd.tmp1 <- nd_hot.data %>% filter(CHR == scaffold_no)
        nd.tmp2 <- nd_cool.data %>% filter(CHR == scaffold_no)
        
        # Plot
        ggplot(nd.tmp1, aes(x=BP, y=.data[[statistic]])) +
                
                # Add values for hot season
                geom_line(linewidth = 2, colour = "orangered") +
                
                # Add values for cool season
                geom_line(data = nd.tmp2, aes(x=BP, y=.data[[statistic]]),
                          linewidth = 2, colour = "turquoise") +
                
                # remove space between plot area and x-axis
                scale_y_continuous(expand = c(0, 0), limits = ylims) +
                
                # set x-axis limits
                scale_x_continuous(limits = xlims) +
                
                # add plot and axis titles
                labs(x = paste0("Position on Scaffold ", scaffold_no), y = ylab) +
                
                # Custom theme:
                theme_bw(base_size = 18) +
                theme( 
                        plot.title = element_text(hjust = 0.5),
                        legend.position="none",
                        panel.border = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank()
                )
}

# Azores ----
# Azores Scaffold 0
p1 <- outlier_plots(pcadapt.data, fd.data, 0, c(14e6, 14.5e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 0, c(14e6, 14.5e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 0, c(14e6, 14.5e6), c(-2.5,2), "Tajima's D")
pdf("Az_Des_scaffold0.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Azores Scaffold 4
p1<-outlier_plots(pcadapt.data, fd.data, 4, c(12.75e6, 13.25e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2<-nd_plots(nd_hot.data, nd_cool.data, "tP", 4, c(12.75e6, 13.25e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3<-nd_plots(nd_hot.data, nd_cool.data, "Tajima", 4, c(12.75e6, 13.25e6), c(-2.5,2), "Tajima's D")
pdf("Az_Des_scaffold4.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Azores Scaffold 6
p1<-outlier_plots(pcadapt.data, fd.data, 6, c(3e6, 3.5e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2<-nd_plots(nd_hot.data, nd_cool.data, "tP", 6, c(3e6, 3.5e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3<-nd_plots(nd_hot.data, nd_cool.data, "Tajima", 6, c(3e6, 3.5e6), c(-2.5,2), "Tajima's D")
pdf("Az_Des_scaffold6.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Azores Scaffold 12
p1<-outlier_plots(pcadapt.data, fd.data, 12, c(5.5e6, 6e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2<-nd_plots(nd_hot.data, nd_cool.data, "tP", 12, c(5.5e6, 6e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3<-nd_plots(nd_hot.data, nd_cool.data, "Tajima", 12, c(5.5e6, 6e6), c(-2.5,2), "Tajima's D")
pdf("Az_Des_scaffold12.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Azores Scaffold 37
p1<-outlier_plots(pcadapt.data, fd.data, 37, c(8.25e6, 8.75e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2<-nd_plots(nd_hot.data, nd_cool.data, "tP", 37, c(8.25e6, 8.75e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3<-nd_plots(nd_hot.data, nd_cool.data, "Tajima", 37, c(8.25e6, 8.75e6), c(-2.5,2), "Tajima's D")
pdf("Az_Des_scaffold37.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Desertas ----
# Desertas Scaffold 9
p1 <- outlier_plots(pcadapt.data, fd.data, 9, c(14.75e6, 15.25e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 9, c(14.75e6, 15.25e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 9, c(14.75e6, 15.25e6), c(-2.5,2), "Tajima's D")
pdf("Des_Sel_scaffold9.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Desertas Scaffold 15
p1 <- outlier_plots(fst.data, fd.data, 15, c(2.5e6, 3e6), c(-0.125,0.5), bquote(F[ST]/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 15, c(2.5e6, 3e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 15, c(2.5e6, 3e6), c(-2.5,2), "Tajima's D")
pdf("Des_Sel_scaffold15.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Desertas Scaffold 18
p1 <- outlier_plots(fst.data, fd.data, 18, c(1e6, 1.5e6), c(-0.125,0.5), bquote(F[ST]/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 18, c(1e6, 1.5e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 18, c(1e6, 1.5e6), c(-2.5,2), "Tajima's D")
pdf("Des_Sel_scaffold18.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Selvagem ----
# Selvagem Scaffold 4
p1 <- outlier_plots(pcadapt.data, fd.data, 4, c(12e6, 12.5e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 4, c(12e6, 12.5e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 4, c(12e6, 12.5e6), c(-2.5,2), "Tajima's D")
pdf("Sel_Des_scaffold4.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Selvagem Scaffold 18
p1 <- outlier_plots(pcadapt.data, fd.data, 18, c(12.25e6, 12.75e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 18, c(12.25e6, 12.75e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 18, c(12.25e6, 12.75e6), c(-2.5,2), "Tajima's D")
pdf("Sel_Des_scaffold18.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()

# Selvagem Scaffold 37
p1 <- outlier_plots(pcadapt.data, fd.data, 37, c(7.75e6, 8.25e6), c(-0.125,0.5), bquote(PCAdapt/f[dM])) + 
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- nd_plots(nd_hot.data, nd_cool.data, "tP", 37, c(7.75e6, 8.25e6), c(0, 50), "Nucleotide Diversity") +
        theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p3 <- nd_plots(nd_hot.data, nd_cool.data, "Tajima", 37, c(7.75e6, 8.25e6), c(-2.5,2), "Tajima's D")
pdf("Sel_Des_scaffold37.pdf", 12, 12)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
dev.off()