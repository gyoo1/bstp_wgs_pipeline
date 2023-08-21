### Gene Ontology ###
setwd("/Users/gyoo1/Desktop/Queens MSc/Project/Data/Results_Main/04_Overlaps")

library(dplyr)

# Load full annotation file
annot <- read.table("LESP_maker2uni.blastp")
annot <- annot %>% mutate(ID=gsub("sp\\|(.*)\\|.*", "\\1", annot$V2), 
                          sp=gsub("sp\\|.*\\|.*_(.*)", "\\1", annot$V2))

# Extract gene lists for outliers and the whole dataset
genes2GO <- function(genes_of_interest, background, species){
        
        # Select species your genes come from
        back <- background %>% filter(sp==species)
        
        # Load outlier annotations file (.txt)
        outliers <- read.table(genes_of_interest, fill = T)
        
        # Clean up outlier annotations
        goi <- outliers %>% 
                filter(V6=="mRNA") %>% 
                mutate(ID=gsub(".*;Name=(.*mRNA-\\d+);.*", "\\1", V12))
        
        # Match genes to the background dataset
        annot_goi <- back[back$V1 %in% goi$ID,]
        annot_goi_all <- background[background$V1 %in% goi$ID,] %>% 
                mutate(scaffold = gsub(".*-(scaffold_\\d+).*", "\\1", V1)) %>% 
                mutate(sc = as.numeric(gsub(".*(\\d+?)", "\\1", scaffold))) %>%
                arrange(sc) %>%
                select(scaffold, ID, sp, gene = V1)
        
        # Save output
        write.table(annot_goi_all, file = "goi_all.txt", quote = F, sep = "\t", row.names = F)
        cat(paste0(annot_goi$ID, collapse = "\n"), file = "goi.txt")
        cat(paste0(back$ID, collapse = "\n"), file = "back.txt")
}

#Outlier_fdM overlaps ----
genes2GO("./shared_outliers/introgression/annot_fst_fdM_DesH-SelH.txt", annot, "HUMAN")
genes2GO("./shared_outliers/introgression/annot_pcadapt_fdM_DesH-SelH.txt", annot, "HUMAN")
genes2GO("./shared_outliers/introgression/annot_pcadapt_fdM_SelH-DesH.txt", annot, "HUMAN")
genes2GO("./shared_outliers/introgression/annot_pcadapt_common_fdM_AzH-DesH.txt", annot, "HUMAN")

#Outliers ----

#FST Only
genes2GO("shared_outliers/outliers_just_fst/annot_Az.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_fst/annot_Des.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_fst/annot_Sel.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_fst/annot_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_fst/annot_SA.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_fst/annot_Gal.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_fst/annot_Guad.txt", annot, "HUMAN")

#PCAdapt Only
genes2GO("shared_outliers/outliers_just_pcadapt/annot_Az_pcadapt.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_pcadapt/annot_Des_pcadapt.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_pcadapt/annot_Sel_pcadapt.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_pcadapt/annot_CV_pcadapt.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_pcadapt/annot_SA_pcadapt.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_pcadapt/annot_Gal_pcadapt.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_just_pcadapt/annot_Guad_pcadapt.txt", annot, "HUMAN")

#Both ----
genes2GO("shared_outliers/outliers_both/annot_com_Az.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_both/annot_com_Des.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_both/annot_com_Sel.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_both/annot_com_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_both/annot_com_SA.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_both/annot_com_Gal.txt", annot, "HUMAN")
genes2GO("shared_outliers/outliers_both/annot_com_Guad.txt", annot, "HUMAN")

# Parallel outliers ----

#Shared outliers

# With CV
genes2GO("shared_outliers/genetic_parallelism/annot_fst_Des_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_fst_Sel_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_fst_Gal_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_pcadapt_Az_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_pcadapt_Des_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_pcadapt_Sel_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_pcadapt_SA_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_pcadapt_Gal_CV.txt", annot, "HUMAN")
genes2GO("shared_outliers/genetic_parallelism/annot_pcadapt_Guad_CV.txt", annot, "HUMAN")

# Read files back in ----
goi <- read.table("goi.txt", sep = "\t")
goi_all <- read.table("goi_all.txt", sep = "\t")
back <- read.table("back.txt", sep = "\t")
