#Set up the working directory
dir_adm <- '~/Desktop/Comp_UU/REF_SalAlp_UK/Admixture/Data_rerun'


#Load package
require(ggplot2)
require(egg)
require(tidyverse)

#Load files 
adm_Myv <- read.table(file.path('~/Desktop/Comp_UU/REF_SalAlp_UK/Admixture/Subsampled/Myv_subsampled_2.qopt'), header = F)

#Read bam list with lakes file 
bams <- read.csv("~/Desktop/Comp_UU/REF_SalAlp_UK/BAM_list/all_bam_lakes.list")

#Remove those from Myv list
Myv_rm_list <- c('Myv-12', 'MyvK-20')
bams_Myv <- bams[bams$Lake == "Mývatn" , ]

#Merge the information above 
adm_Myv_list <- cbind(adm_Myv, bams_Myv)
adm_Myv_list <-adm_Myv_list[!adm_Myv_list$SampleID %in% Myv_rm_list, ]

#Make a matrix
adm_Myv_pivot <- 
  adm_Myv_list %>% 
  pivot_longer(cols = c("V1", "V2"),
               names_to =  "AdmComponents", 
               values_to = "AdmProportions")



#Plot admixture with the pcangsd best fit K
library(ggthemes)

Myv <- 
  ggplot(adm_Myv_pivot, 
               aes(x = as.factor( SampleID), y =  AdmProportions, fill = factor(AdmComponents))) +
  geom_col(aes(color = AdmComponents), size = 0.1)+
  facet_grid(~Morph_short, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    
    text=element_text(size=14),
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf"), guide = "none") +
  scale_color_manual(values = c("#ef8a62", "#67a9cf"), guide = 'none')+
  labs(subtitle = "K=2", title = "b")+
  theme(plot.title = element_text(face = "bold", family = "Times")
  )



#Pca plot
dir_pca <- '~/Desktop/Comp_UU/REF_SalAlp_UK/Admixture/Data_rerun'

library(psych)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(egg)
require(ggforce)


#Add pcangsd matrix
pca_Myv_MyvK <-  as.matrix(
  read.table(file.path(dir_pca, paste("Myv.cov")), header = F))
Myv_list <- bams_Myv

#perform the pca using the eigen function.
data <- pca_Myv_MyvK
bamList <- Myv_list

#Run this uniquly
eigen.data <- eigen(data)


#We can then extract the eigenvectors from the pca object and format them into a dataframe for plotting, e.g. using ggplot().
eigenvectors <- as.data.frame(eigen.data$vectors)
eigenvalues <-  eigen.data$values

eigenvectors$Sample <- bamList$SampleID
eigenvectors$Morph_short <- bamList$Morph_short
eigenvectors$Lake <-  bamList$Lake

eigenvectors <- eigenvectors[!eigenvectors$Sample %in% Myv_rm_list, ]

#Get vars
pca.eigenval.sum = sum(eigen.data$values)
varPC1 <- (eigen.data$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (eigen.data$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2


#Set uo color pannel
lake_colors <- c("Mývatn" = "red",
                 "Thingvallavatn" = 'darkorchid','Sirdalsvatnet' = "darkgreen",
                 'Vangsvatnet' = "steelblue")


#All
pcaMyv <- 
  ggplot(data = eigenvectors, 
                 aes(x = V1, y = V2)) +
  geom_point(alpha = 0.8, size = 2.5, shape = 21, color = 'black', aes(fill = Morph_short)) +
  xlab(paste0("PC1: ", round(varPC1,1),"% variance")) +
  ylab(paste0("PC2: ", round(varPC2,1),"% variance")) +
  theme_bw(14) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("SB" = "red", "LG" =  "#67a9cf")) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  ylim(-0.95,0.95)+
  xlim(-0.3, 0.3)+
  labs(fill='Morph', title = "a", subtitle = "Mývatn" )+
  theme(plot.title = element_text(face = "bold", family = "Times")
        
      )


#Arrange all plots together with pca plots 
pc_admix <- ggarrange(pcaMyv,Myv, nrow = 1)
ggsave('~/Desktop/Comp_UU/REF_SalAlp_UK/Admixture/Admix_PC_Myvatn_regroup.pdf',pc_admix, width = 9, height = 3, dpi = 300)

#How many individuals in the pca and admixture
length(eigenvectors$Morph_short[eigenvectors$Morph_short == "LG"])
length(eigenvectors$Morph_short[eigenvectors$Morph_short == "SB"])
length(eigenvectors$Morph_short)
describeBy(adm_Myv_list, "Morph_short")

