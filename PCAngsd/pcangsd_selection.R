# Load Libraries
library(RcppCNPy)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

### Load Files ----

#Load in npy file
npy_dat<- npyLoad("./PCAngsd/PCAngsd.selection.npy")

#Load in .sites file
sites<-read.table("./PCAngsd/snp.sites", header = T, colClasses = c("marker"="factor"))
names(sites)<-c("scaffold","marker")

#Load in alt 0/1 sites file
alt_sites=read.table("./PCAngsd/PCAngsd.sites")

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

#Convert test statistics to p-values
pval=1-pchisq(npy_dat,1)
pval=as.data.frame(pval)

#Align sites files
sites$V1=alt_sites$V1

#Remove elements of .sites values where V1=0
sig=sites%>%filter(V1==1)

#Create data frame combining sites and pvals
Data=data.frame(marker=sig$marker, scaffold=sig$scaffold, pval1=pval$V1)

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

ggplot(aes(x=BPcum, y=-log10(pval1), colour=scaffold), data=df.tmp) +
        geom_point() + theme_bw() + geom_hline(yintercept=-log10(4.89*10^(-9))) +
        xlab("Marker Position") + ylab("-log10(p-value)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) +
        theme(plot.title = element_text(hjust = 0), legend.position = "none")
