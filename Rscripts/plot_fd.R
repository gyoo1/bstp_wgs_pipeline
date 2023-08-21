### Introgression Plots ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/03_introgression")

# Libraries ====
library(dplyr)
library(ggplot2)

# Read data
fd <- read.table(file.choose(), header = T)

# Format data table
fd.data <- fd %>%
        mutate(ID = NA, 
               CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", chr)), 
               BP = windowStart + ((windowEnd - windowStart)/2)) %>%
        arrange(CHR) %>%
        mutate(ID = paste0(rep("Win", length(rownames(fd))), rownames(fd)))

# Calculate introgression outlier threshold
fd.tholds <- function(df, statistic){
        df.tmp <- as.vector(df %>% select(all_of(statistic)))
        thold = quantile(pull(df.tmp), 0.999) #99.9th percentile
        return(thold)
}

thold <- fd.tholds(fd.data, "f_dM")

outliers <- which(fd.data$f_dM > thold)
length(outliers)
Outliers <- fd.data %>% filter(f_dM > thold) %>% 
        select(ID, chr, windowStart, windowEnd, f_dM, BP)
out_bed <- Outliers %>% select(chr, windowStart, windowEnd)
write.table(Outliers, "fdM_outlier.txt", quote = F, sep = "\t")
write.table(out_bed, "fdM.bed", quote = F, sep = "\t", col.names = F, row.names = F)

# Plot fdM
gg.manhattan <- function(df, statistic, threshold, scaffold, col, ylims){
        
        # format df
        df.tmp <- df %>% 
                
                # Compute chromosome size
                group_by(CHR) %>% 
                summarise(chr_len=max(windowEnd)) %>% 
                
                # Calculate cumulative position of each chromosome
                mutate(tot=cumsum(chr_len)-chr_len) %>%
                select(-chr_len) %>%
                
                # Add this info to the initial dataset
                left_join(df, ., by=c("CHR"="CHR")) %>%
                
                # Add a cumulative position of each SNP
                arrange(CHR, windowEnd) %>%
                mutate(BPcum = windowEnd+tot)
        
        ggplot(df.tmp, aes_string(x= "BPcum", y= statistic)) +
                # Show all points
                geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
                scale_color_manual(values = rep(col, 22)) +
                
                # remove space between plot area and x-axis
                scale_y_continuous(expand = c(0, 0), limits = c(-1,1)) +
                
                # add plot and axis titles
                labs(x = "Position across Scaffolds",
                     y = "fdM") +
                
                # add genome-wide sig and sugg lines
                geom_hline(yintercept = threshold) +
                
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

pdf("fdM_DesH_AzH_vNAC.pdf", 24, 6)
gg.manhattan(fd.data, "f_dM", thold, col = rep(c("black", "grey"), 309))
dev.off()

pdf("fdM_DesH_SelH_vNAC.pdf", 24, 6)
gg.manhattan(fd.data, "f_dM", thold, col = rep(c("black", "grey"), 309))
dev.off()

pdf("fdM_DesH_AzH_vSelH.pdf", 24, 6)
gg.manhattan(fd.data, "f_dM", thold, col = rep(c("black", "grey"), 309))
dev.off()

# Plot shared windows
shared <- read.table("fdM_com_DesH_AzH.bed")
fd.data <- fd.data %>% mutate(win = paste0(chr,windowStart,windowEnd))
shared <- shared %>% mutate(win = paste0(V1, V2, V3))
fd.data <- fd.data %>% mutate(com = win %in% shared$win)
df.tmp <- fd.data %>% 
        
        # Compute chromosome size
        group_by(CHR) %>% 
        summarise(chr_len=max(windowEnd)) %>% 
        
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(fd.data, ., by=c("CHR"="CHR")) %>%
        
        # Add a cumulative position of each SNP
        arrange(CHR, windowEnd) %>%
        mutate(BPcum = windowEnd+tot)

p1 <- ggplot(df.tmp, aes(x=BPcum, y=f_dM)) +
        # Show all points
        geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
        scale_color_manual(values = rep(c("black", "grey"), 500)) +
        
        # remove space between plot area and x-axis
        scale_y_continuous(expand = c(0, 0), limits = c(-1,1)) +
        
        # add plot and axis titles
        labs(x = "Position across Scaffolds", y = "fdM") +
        
        # add genome-wide sig and sugg lines
        geom_hline(yintercept = thold) +
        
        geom_point(data=subset(df.tmp, com=="TRUE"), color="orange", size=2) +
        
        # Custom theme:
        theme_bw(base_size = 18) +
        theme( 
                plot.title = element_text(hjust = 0.5),
                legend.position="none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
        )

pdf("fdM_DesH_AzH.pdf", 24, 6)
p1
dev.off()



