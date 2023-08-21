## Estimate Heterozygosity from SFS
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/01_exploratory")

# get list of files
fnames <- list.files(path = "./Heterozygosity", pattern = "SFS_\\d+")
fpaths <- paste0("./Heterozygosity","/",fnames)

# calculate heterozygosity
het <- c()

for(i in 1:length(fnames)){
        tmp <- scan(fpaths[i])
        het <- c(het, tmp[2]/sum(tmp))
}

pop <- c("AzH","AzC","DesH","DesC","SelH","SelC","CVH",
                "CVC","SA_H","SA_C","GalH","GalC","TOSP","AISP")

Het <- as.data.frame(cbind(pop, round(het, 5)))

write.csv(Het, "obs_het.csv")

# plot heterozygosity
library(ggplot2)

colours <- rep(c("orangered","lightblue"), 7)

ggplot(Het, aes(x=pop, y=het, fill=pop)) + 
        geom_bar(stat="identity") +
        coord_flip() +
        scale_fill_manual(values=colours) +
        labs(x="Population", y="Heterozygosity") +
        theme(legend.position = "none",
              panel.background = element_blank())
