# Load Libraries
library(RcppCNPy) # Numpy library for R
library(dplyr)
library(ggplot2)
library(qvalue)

### Load Files ----

# Load in npy file
Zscores=npyLoad("./PCAngsd/pcadapt/PCAngsd.pcadapt.zscores.npy")

#Load in args file
args=read.csv("./PCAngsd/pcadapt/PCAngsd.args")

#Load in cov file
cov=read.table("./PCAngsd/pcadapt/PCAngsd.cov")

#Load in sites file
sites<-read.table("./PCAngsd/snp.sites", header = T, colClasses = c("marker"="factor"))
names(sites)<-c("scaffold","marker")

#Load in alt sites file
alt_sites=read.table("./PCAngsd/pcadapt/PCAngsd.sites")

###----

### Create QQplot to QC Test Statistics ----

#Run qqchi function
qqchi<-function(x,...){
        lambda<-round(median(x)/qchisq(0.5,1),2)
        qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
        legend("topleft",paste("lambda=",lambda))
}
qqchi(Zscores)

###----

#convert test statistics to p-values (following pcadapt_script.R)
source("./pcadapt_script_PC.R")
PCp=read.table("PC.pcadapt.pval.txt")

#adjusting p values to q values following Caplins, S. (2021). Marine Genomics 2021. 
#https://baylab.github.io/MarineGenomics/week-7--fst-and-outlier-analysis.html#finding-outliers-using-pcadapt
qval <- qvalue(PCp)$qvalues
outliers <- which(qval<0.1)
length(outliers)

#Aligning sites files
sites$V1=alt_sites$V1

#Remove elements of .sites values where V1=0
sig=sites%>%filter(V1==1)

#Create data frame combining sites and pvals
Data=data.frame(marker=sig$marker, scaffold=sig$scaffold, q=qval$V1)

#Number of significant hits, and filter non-significant results from data set
outliers <- which(Data$q<4.89*10^{-9})
length(outliers)
Outliers=Data%>%filter(q<4.89*10^{-9})

### Manhattan Plot ----
df.tmp <- Data %>% 
        
        # Compute chromosome size
        group_by(scaffold) %>% 
        summarise(chr_len=max(marker)) %>% 
        
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(Data, ., by=c("scaffold"="scaffold")) %>%
        
        # Add a cumulative position of each SNP
        arrange(scaffold, marker) %>%
        mutate(BPcum = marker+tot)

ggplot(aes(x=BPcum, y=-log10(q), colour=scaffold), data=df.tmp) +
        geom_point() + theme_bw() + geom_hline(yintercept=-log10(4.89*10^(-9))) +
        xlab("Marker Position") + ylab("-log10(q-value)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) +
        theme(plot.title = element_text(hjust = 0), legend.position = "none")
