## Depth per Population
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/03_depth")
library(dplyr)

depths <- read.table(file="depths.txt")
depths$V1 <- gsub(".*(BSTP_.+)(_.+_).*", "\\1",
                  gsub(".*AISP.*", "AISP",
                       gsub(".*TOSP.*", "TOSP", depths$V1))) #pull out sample locations
depths$V1[109] <- "BSTP_SHC"

# Calculate mean and median depth for each population
depth_by_pop <- depths %>% group_by(V1) %>% 
        summarise(Mean_Depth = mean(V3), Median_Depth = median(V3)) %>% 
        transmute(Sampling_Site = V1, 
                  Mean_Depth = round(Mean_Depth, 2), 
                  Median_Depth = round(Median_Depth,2)) 

write.table(depth_by_pop, "depth_by_pop.txt", sep = "\t", row.names = F)

# Plot a histogram for depth across sites in a single individual
hist_depth <- read.table("hist_site_depth.txt")
hist_depth_clean <- hist_depth %>% 
        transmute(bin = V3, count = V5) %>% arrange(bin)
plot(subset(hist_depth_clean, bin < 60))
