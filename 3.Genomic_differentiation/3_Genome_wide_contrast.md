# This script runs ANGSD’s score test-based association scan (-doAsso 1) to test for genotype–phenotype associations (binary traits) using genotype likelihoods per lake.
## Step 1: Genome-Wide Contrast (GWC) for Morph Contrasts within Arctic Charr Lakes Using -doAsso 1 in ANGSD
```{r}
#!/bin/bash -l

#SBATCH -A naiss2023-5-221
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 09-00:00:00
#SBATCH --array=1-13
#SBATCH -J asso_no_cov
#SBATCH -e asso_no_cov_%A_%a.err
#SBATCH -o asso_no_cov_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@imbim.uu.se

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load PCAngsd/1.11

# Go to the directory where BAM files are located
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2
cd $BASEDIR

# BAM lists per lake and morph contrast (each list corresponds to one GWAS group)
BAM_LIST_PATH=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/GWAS_list

# SNP sites to be used in the analysis
SITES=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/SITES/sites_sorted.txt

# Reference genome
REFGENOME=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta
REF_INDEXED=${REFGENOME}.fai

# Output directory
OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GWAS


# SELECT FILE FOR ANALYSIS 


# Select BAM list from SLURM array ID (1 file per job)
BAM_LIST=$(ls $BAM_LIST_PATH/assos_*.list | sed -n ${SLURM_ARRAY_TASK_ID}p)
BAM_TARGET=${BAM_LIST/.list/}
OUTPUT=$(basename $BAM_TARGET)

# Get number of individuals in this comparison
x=$(wc -l < $BAM_LIST)

echo "Running association analysis for: $OUTPUT"
echo " - Individuals: $x"
echo " - BAM list: $BAM_LIST"
echo " - Binary phenotype file: $BAM_TARGET.ybin"


# RUN ANGSD GWAS


# Use ANGSD to perform score statistic-based GWAS
# Output: association statistics including p-values per SNP

angsd \
  -yBin $BAM_TARGET.ybin \          # Binary phenotype vector (0/1)
  -doAsso 1 \                       # Perform association using score test
  -out $OUT_DIR/$OUTPUT \          # Output prefix
  -ref $REFGENOME -fai $REF_INDEXED \ # Reference genome
  -doMajorMinor 4 \                # Estimate major/minor alleles using GLs
  -doMaf 1 \                       # Calculate minor allele frequency
  -bam $BAM_LIST -nInd $x \        # Input BAMs and number of individuals
  -uniqueOnly 1 \                  # Only uniquely mapped reads
  -remove_bads 1 \                 # Remove bad reads
  -only_proper_pairs 1 \           # Only properly paired reads
  -minMapQ 30 -minQ 20 \           # Minimum mapping and base quality
  -gl 2 \                          # Use GATK genotype likelihood model
  -Pvalue 1 \                      # Output p-values
  -sites $SITES \                 # Limit analysis to predefined sites
  -nThreads 8                      # Number of threads to use
  
```

## Step 2: R script that visualizes genome-wide contrasts (GWCs) across all morph comparisons (Manhattan plots)
```{r}
#Set up
Load package
require(ggplot2)
require(egg)
require(tidyverse)
require(ggpubr)
require(dplyr)
require(grid)


#Read p-values
dir= '~/Desktop/Comp_UU/REF_SalAlp_UK/GWAS/Data'

#File names
pval_list <- c("Sir_D_vs_N.pvalue",  "Van_D_vs_N.pvalue", "Myv_PCA_SB_PL.pvalue",
                "Thi_DB_vs_LB.pvalue","Thi_DB_vs_Pi.pvalue","Thi_DB_vs_PL.pvalue",
                "Thi_Pi_vs_LB.pvalue", "Thi_PL_vs_LB.pvalue", "Thi_PL_vs_Pi.pvalue")


#Load data 
files <- list.files(dir)

# Find the files in pval_list that exist in the directory
files_to_read <- intersect(files, pval_list)

# Read the desired files
gwas_list <- list()

for(file in files_to_read) {
  file_content <- read.table(file.path(dir, file), header = T)
  name = gsub(".pvalue", "", file)
  # Process the file content as needed
  gwas_list[[name]] <- file_content
}


#Load package
require(ggplot2)
require(egg)
require(tidyverse)
require(ggpubr)
require(dplyr)
require(grid)


orderContigs = c("1","2","3","4","5","6",
                 "7", 
                 "8","9","10","11","12" ,"13","14","15", "16",
                 "17" ,"18","19","20","21",
                 "22", 
                 "23","24","25", "26","27","28",
                 "29","30",
                 "31","32","33","34","35","36","37","38", "39","40","NA")


#Order a data list
order_list <- c("Sir_D_vs_N",  "Van_D_vs_N", "Myv_PCA_SB_PL",
                "Thi_DB_vs_LB","Thi_DB_vs_Pi","Thi_DB_vs_PL",
                "Thi_Pi_vs_LB", "Thi_PL_vs_LB", "Thi_PL_vs_Pi")

gwas_list_order <- gwas_list[order_list]

#Add new names to the list for manhattans
names(gwas_list_order) <- 
                          c("a) Sirdalsvatnet, DB vs. LP",  
                           "b) Vangsvatnet, DB vs. LP",   
                           "c) Mývatn, SB vs. LG",
                           "d) Thingvallavatn, SB vs. LB", 
                           "e) Thingvallavatn, SB vs. Pi", 
                           "f) Thingvallavatn, SB vs. PL", 
                           "g) Thingvallavatn, Pi vs. LB", 
                           "h) Thingvallavatn, PL vs. LB", 
                           "i) Thingvallavatn, PL vs. Pi")


library(ggplot2)
library(dplyr)
library(ggpubr)

# Create a list to store individual plots
ManhattanPlots <- list()


# Run the for loop for the first 8 files and the 9th file separately
for (i in seq_along(gwas_list_order)) {
  file <- names(gwas_list_order)[i]
  gwasResults <- gwas_list_order[[file]]
  
  colnames(gwasResults)[2] <- 'BP'
  colnames(gwasResults)[7] <- 'P'
  colnames(gwasResults)[8] <- 'CHR'
  
  gwasResults$CHR <- factor(gwasResults$CHR, levels = orderContigs)
  
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(as.character(CHR), as.numeric(BP)) %>%
    mutate(BPcum = as.numeric(BP)+ as.numeric(tot) )
  
  # Prepare the X axis
  don$CHR <- factor(don$CHR, levels = orderContigs)
  axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # Threshold bonferroni allowing 10% error
  bonf = -log10(1e-3/nrow(don)) #alpha = 1e3 or 1x10^-3
  max_y = max(-log10(don$P)) + 6
  min_y = bonf-3
  
  # Plot
  manh <- ggplot(don, aes(x = BPcum, y = -log10(P) ) ) +
    # Show all points
    geom_point(aes(color=factor(CHR)), alpha=0.8, size=1) +
    scale_color_manual(values = rep(c("grey", "black"), 9000 )) +
    # Custom X axis
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(3, 25)) +
    # Custom theme
    theme_classic(19) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.subtitle = element_text(size = 20) 
    ) +
    labs(subtitle = paste(file)) +
    geom_hline(yintercept = bonf, linetype="dashed", color = "red", lwd=1.5 )
  
  if (i == 9) {
    manh <- manh +
      theme(
        axis.text.x = element_text(angle = 90)
      ) +
      labs(subtitle = paste(file), x = 'Scaffold')
  }
  
  ManhattanPlots[[paste("Contrast", file)]] <- manh
}

# Arrange and display the plots with common axis labels
final_plot <- ggpubr::ggarrange(plotlist = ManhattanPlots, nrow = 9)

# Save the plot
ggsave('~/Desktop/Comp_UU/Manuscript/Sumbition/Supplementary_Figures/Supplementary_Fig.S5.jpeg', final_plot, width = 16, height = 20, dpi = 300)

```

