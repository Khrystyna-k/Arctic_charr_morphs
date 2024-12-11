#Read the file
cov <- read.table('~/Desktop/Comp_UU/REF_SalAlp_UK/Coverage/Coverage_summary.txt' )
cov$Sample_ID <- gsub('_marked_dups_1_2.bam_bamqc', '', cov$V1)
cov$Cov_stats <- as.numeric(gsub('X', '', cov$V5))

#Add extra column
cov$Ref <- "Arctic_charr_ref_fSalAlp1"

#Summary stats
summary(cov)

  
