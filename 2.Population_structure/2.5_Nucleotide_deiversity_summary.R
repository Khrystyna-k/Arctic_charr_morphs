#Libs
library(tidyverse)
library(patchwork)
library(zoo)
library(data.table)


#Directory
dir <- '~/Desktop/Comp_UU/REF_SalAlp_UK/Nucleotide_diversity/pestPG_data'

#Load data 
files <- list.files(dir)

# Filter files that start with 5kb
desired_files <- grep("theta.5kb.thetasWindow.gz.pestPG$", files, value = TRUE)
col_names <- c('Chr',	'WinCenter',	'tW',	'tP',	'tF',	'tH',	'tL',	'Tajima', 'fuf',
               'fud',	'fayh',	'zeng',	'nSites')


#Look at the list
desired_files

# Read the desired files
thetas_list <- list()

for(file in desired_files) {
  file_content <- read.table(file.path(dir, file), header = F)[,-1]
  colnames(file_content) <- col_names
  # Process the file content as needed
  thetas_list[[file]] <- file_content
}

#Modify names of the data list 
names(thetas_list) <- gsub(".theta.5kb.thetasWindow.gz.pestPG", '', names(thetas_list) )

# List of prefixes to search for
prefixes <- c('maf_MyvDB', 'maf_MyvPi', 'maf_SirDB', 'maf_SirPL', 
              'maf_Thi_LB_SB', 'maf_Thi_Pi_PL', 'maf_ThiLB', 
              'maf_ThiPi', 'maf_ThiPL', 'maf_ThiSB', 'maf_Van_DB',
              'maf_Van_PL', 'Myv', 'Sir', 'Van', 'Thi'
)  # Add more prefixes as needed

# Initialize a new list to store combined data frames
new_data_list <- list()

# Define a function to bind rows of data frames starting with a similar name
for(prefix in prefixes) {
  # Filter data frames starting with the specified prefix
  matching_data_frames <- grep(paste0("^", prefix), names(thetas_list), value = TRUE)
  
  # If there are matching data frames, bind them row-wise
  if (length(matching_data_frames) > 0) {
    combined_data <- do.call(rbind, thetas_list[matching_data_frames])
    new_data_list[[prefix]] <- combined_data
  }
}


# Estimate Pi and Watterson’s theta  -------------------------------------------
l.theta <- list()
l.watt <- list()
l.sum <- list()
l.dfs <- list()

for (pop in seq_along(new_data_list)){
  df.theta <- new_data_list[[pop]]
  pop.name <- names(new_data_list[pop])
  l.theta[[pop.name]] <- mean((as.numeric(df.theta$tP) / as.numeric(df.theta$nSites)), na.rm = T)
  l.watt[[pop.name]] <- mean((as.numeric(df.theta$tW) /as.numeric(df.theta$nSites)), na.rm = T)
  l.sum[[pop.name]] <- sum(as.numeric(df.theta$nSites), na.rm = T)
  df.theta$pop <- pop.name
  l.dfs[[pop.name]] <- df.theta
}

length(l.theta) #as the number of populations
length(l.watt)  #as the number of populations

df.x2 <- data.frame(pop = names(l.theta),
                    pairwise.nuc = round(unlist(l.theta), 4),  
                    wattersons.theta = round(unlist(l.watt), 4),
                    number.sites = round(unlist(l.sum), 4))

mean.theta.nsites <- df.x2
df.theta.nsites <- bind_rows(l.dfs)
names(df.theta.nsites) <- c("Chr", "WinCenter","tW", "tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites", "pop")


#Cleaner summary per lake and per morph
#Filter windows > 10 sites for higher accuracy
df.theta.nsites$nSites <- as.numeric(df.theta.nsites$nSites)
df.theta.nsites$tP <- as.numeric(df.theta.nsites$tP)
df.theta.nsites$tW <- as.numeric(df.theta.nsites$tW)

clean.summary.thetas <- 
  df.theta.nsites %>% 
  filter(nSites > 10) %>%
  group_by(pop) %>%
  summarise(mean.pi = round(mean(tP / nSites),4),
            sd.pi = round(sd(tP / nSites),4),
            mean.wat = round(mean(tW / nSites),4),
            sd.wat = round(sd(tW / nSites),4)) %>%
  mutate(out.pi = paste0(mean.pi,"±",sd.pi),
         out.wat = paste0(mean.wat,"±",sd.wat),)

write.csv(clean.summary.thetas, "~/Desktop/Comp_UU/REF_SalAlp_UK/Nucleotide_diversity/Output/Clean_theta_folded_perpop.csv")

