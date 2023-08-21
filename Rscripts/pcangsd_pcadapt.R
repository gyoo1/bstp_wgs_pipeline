## Plotting pcadapt selection scan results ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/02_outliers/")

# Load Libraries
library(RcppCNPy) # Numpy library for R
library(dplyr)
library(ggplot2)
library(qvalue)

### Load Files ----

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

pcadapt_format("pair1")
pcadapt_format("pair2")
pcadapt_format("pair3")
pcadapt_format("pair4")
pcadapt_format("pair5")
pcadapt_format("pair6")
pcadapt_format("pair7")

### Q-values ----

#convert test statistics to p-values (following pcadapt_script.R)
source("./pcadapt/pcadapt_script_PC.R")
PCp=read.table("PC.pcadapt.pval.txt")

#adjusting p values to q values following Caplins, S. (2021). Marine Genomics 2021. 
#https://baylab.github.io/MarineGenomics/week-7--fst-and-outlier-analysis.html#finding-outliers-using-pcadapt
qval <- qvalue(PCp)$qvalues

#Aligning sites files
sites$V1=alt_sites$V1

#Remove elements of .sites values where V1=0
sig=sites%>%filter(V1==1)

#Create data frame combining sites and pvals
Data=data.frame(marker=sig$marker, scaffold=sig$scaffold, q=qval$V1)

#Arrange by scaffold number
Data <- Data %>% mutate(CHR = as.numeric(gsub(".*_(\\d+)$", "\\1", scaffold))) %>% arrange(CHR)

#Number of significant hits, and filter non-significant results from data set
threshold <- 0.05/(2*1157551241)
outliers <- which(Data$q<threshold)
length(outliers)
Outliers <- Data%>%filter(q<threshold)
out_bed <- Outliers %>% transmute(scaffold, start = pmax(marker-5000, 0), end = marker+5000)

write.table(Outliers, "pcadapt_pair1_out.txt", quote = F, sep = "\t")
write.table(Outliers, "pcadapt_pair2_out.txt", quote = F, sep = "\t")
write.table(Outliers, "pcadapt_pair3_out.txt", quote = F, sep = "\t")
write.table(Outliers, "pcadapt_pair4_out.txt", quote = F, sep = "\t")
write.table(Outliers, "pcadapt_pair5_out.txt", quote = F, sep = "\t")
write.table(Outliers, "pcadapt_pair6_out.txt", quote = F, sep = "\t")
write.table(Outliers, "pcadapt_pair7_out.txt", quote = F, sep = "\t")

write.table(out_bed, "pcadapt_pair1.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(out_bed, "pcadapt_pair2.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(out_bed, "pcadapt_pair3.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(out_bed, "pcadapt_pair4.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(out_bed, "pcadapt_pair5.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(out_bed, "pcadapt_pair6.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(out_bed, "pcadapt_pair7.bed", quote = F, sep = "\t", col.names = F, row.names = F)

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
        mutate(BPcum = marker+tot) %>%

        # Remove q-values of 0s
        filter(q > 0)

# Select every 1000th row
# df.tmp <- df.tmp %>% slice(which(row_number() %% 100 == 1))

plot <- ggplot(df.tmp, aes(x=BPcum, y=-log10(q))) +
        
                geom_point(aes(color=as.factor(scaffold)), alpha=0.8, size=1) +
                scale_color_manual(values = rep(c("black", "grey"), 743)) +
                
                scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(df.tmp$q)))) +
                labs(x = "Position across Scaffolds", y = "-log10(q-value)") +
                
                geom_hline(yintercept=-log10(threshold)) +
                
                theme_bw(base_size = 18) +
                theme( 
                        plot.title = element_text(hjust = 0.5),
                        legend.position="none",
                        panel.border = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank()
                )

pdf("pcadapt_Azores.pdf", 24, 6) #pair1
plot
dev.off()

pdf("pcadapt_Desertas.pdf", 24, 6) #pair2
plot
dev.off()

pdf("pcadapt_Selvagem.pdf", 24, 6) #pair3
plot
dev.off()

pdf("pcadapt_SA.pdf", 24, 6) #pair4
plot
dev.off()

pdf("pcadapt_Galapagos.pdf", 24, 6) #pair5
plot
dev.off()

pdf("pcadapt_Guadalupe.pdf", 24, 6) #pair6
plot
dev.off()

pdf("pcadapt_CapeVerde.pdf", 24, 6) #pair7
plot
dev.off()
