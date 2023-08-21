## Plotting fastPCA selection scan results ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/02_outliers/")

# Load Libraries
library(RcppCNPy)
library(dplyr)
library(ggplot2)

# Load Data
pcadapt_format <- function(pair_num) {
        # Load in npy file
        npy_dat <<- npyLoad(paste0("./pcadapt/fastpca_", pair_num, ".selection.npy"))
        
        #Load in sites file
        sites <<- read.table(paste0("./pcadapt/", pair_num, "_pruned.txt"), 
                             header = T, colClasses = c("marker"="factor"))
        names(sites) <<- c("scaffold","marker")
        
        #Load in alt sites file
        alt_sites <<- read.table(paste0("./pcadapt/fastpca_", pair_num, ".sites"))
}

pcadapt_format("pair1")
pcadapt_format("pair2")
pcadapt_format("pair3")
pcadapt_format("pair4")
pcadapt_format("pair5")
pcadapt_format("pair6")
pcadapt_format("pair7")

# Convert test statistics to p-values
pval=1-pchisq(npy_dat,1)
pval=as.data.frame(pval)

#Aligning sites files
sites$V1=alt_sites$V1

#Remove elements of .sites values where V1=0
sig=sites%>%filter(V1==1)

#Create data frame combining sites and pvals
Data=data.frame(marker=sig$marker, scaffold=sig$scaffold, pval=pval$V1)

#Arrange by scaffold number
Data <- Data %>% mutate(CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", scaffold))) %>% arrange(CHR)

#Number of significant hits, and filter non-significant results from data set
threshold <- 0.05/(nrow(Data))
outliers <- which(Data$pval<threshold)
length(outliers)
Outliers=Data%>%filter(pval<threshold)

### Manhattan Plot ----

df.tmp <- Data %>% 
        # Compute chromosome size
        group_by(CHR) %>% 
        summarise(chr_len=max(marker)) %>% 
        
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        dplyr::select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(Data, ., by=c("CHR"="CHR")) %>%
        
        # Add a cumulative position of each SNP
        mutate(BPcum = marker+tot)
        
# # Select every 1000th row
# df.tmp <- df.tmp %>% slice(which(row_number() %% 100 == 1))

plot <- ggplot(df.tmp, aes(x=BPcum, y=-log10(pval))) +
        
                geom_point(aes(color=as.factor(scaffold)), alpha=0.8, size=1) +
                scale_color_manual(values = rep(c("black", "grey"), 750)) +
                
                # scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(df.tmp$pval)))) +
                labs(x = "Position across Scaffolds", y = "-log10(p-value)") +
                
                geom_hline(yintercept=-log10(threshold)) +
                
                theme_bw(base_size = 18) +
                theme( 
                        plot.title = element_text(hjust = 0.5),
                        legend.position="none",
                        panel.border = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank()
                )

pdf("fastpca_Azores.pdf", 24, 6) #pair1
plot
dev.off()

pdf("fastpca_Desertas.pdf", 24, 6) #pair2
plot
dev.off()

pdf("fastpca_Selvagem.pdf", 24, 6) #pair3
plot
dev.off()

pdf("fastpca_SA.pdf", 24, 6) #pair4
plot
dev.off()

pdf("fastpca_Galapagos.pdf", 24, 6) #pair5
plot
dev.off()

pdf("fastpca_Guadalupe.pdf", 24, 6) #pair6
plot
dev.off()

pdf("fastpca_CapeVerde.pdf", 24, 6) #pair7
plot
dev.off()