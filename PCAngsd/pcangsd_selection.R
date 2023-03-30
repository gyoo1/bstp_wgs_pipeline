# Load Libraries
library(RcppCNPy) # Numpy library for R
library(dplyr)
library(ggplot2)
library(qvalue)

### Load Files ----

# Load in npy file
Zscores=npyLoad("./PCAngsd/PCAngsd.selection.npy")

#Load in args file
args=read.csv("./PCAngsd/PCAngsd.args")

#Load in cov file
cov=read.table("./PCAngsd/PCAngsd.cov")

#Load in sites file
sites<-read.table("./PCAngsd/snp.sites", header = T, colClasses = c("marker"="factor"))
names(sites)<-c("scaffold","marker")

#Load in alt sites file
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
outliers <- which(Data$q<6.5*10^{-5})
length(outliers)
Outliers=Data%>%filter(q<6.5*10^{-5})

### Manhattan Plot ----
ggplot(aes(x=marker, y=-p), data=Data)+
        geom_point()+
        xlab("Marker")+ylab("-log10(q-value)")+ theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))+
        theme(plot.title = element_text(hjust = 0), legend.position = "none")
