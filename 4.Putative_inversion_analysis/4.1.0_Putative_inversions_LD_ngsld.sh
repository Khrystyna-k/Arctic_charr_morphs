#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 10
#SBATCH -t 02-00:00:00
#SBATCH --array=1-13
#SBATCH -J LD_inv_B_P
#SBATCH -e LD_inv_B_P_%A_%a.err
#SBATCH -o LD_inv_B_P_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@slu.se


#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load ngsLD/1.1.1

#Directorys
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO_INV
cd $BASEDIR

#Main beagle file
BEAGLE_MAIN=Thin_CharrLowPass_GATKMethod_MinMAF0.05_all_chr.beagle.gz

#Inversion region for 
# Filter lines within the specified range in chromosome and include header
#Chr 4
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "4" && $2 >= 74200000 && $2 <= 77100000'  # Filter data by range
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr4_invLB_SB_ld.beagle.gz

#Chr 5
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "5" && $2 >= 21300000 && $2 <= 23700000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr5_invLB_SB_ld.beagle.gz

#Chr 9
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "9" && $2 >= 60300000 && $2 <= 63100000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr9_invLB_SB_ld.beagle.gz


#Chr 17
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "17" && $2 >= 31400000 && $2 <= 34200000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr17_invLB_SB_ld.beagle.gz

#Contrast for B vs P
#Chr1
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "1" && $2 >= 15300000 && $2 <= 19800000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr1_1_invB_P_ld.beagle.gz

#Chr1
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "1" && $2 >= 18500000 && $2 <= 23200000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr1_2_invB_P_ld.beagle.gz

#Chr3
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "3" && $2 >= 32500000 && $2 <= 36800000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr3_1_invB_P_ld.beagle.gz

#Chr3
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "3" && $2 >= 36300000 && $2 <= 41600000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr3_2_invB_P_ld.beagle.gz

#Chr8
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "8" && $2 >= 28000000 && $2 <= 30900000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr8_invB_P_ld.beagle.gz

#Chr9
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "9" && $2 >= 37400000 && $2 <= 41800000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr9_invB_P_ld.beagle.gz


#Chr14
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "14" && $2 >= 6200000 && $2 <= 7000000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr14_invB_P_ld.beagle.gz



#Chr40
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "40" && $2 >= 15200000 && $2 <= 18000000' 
} | gzip > Thin_CharrLowPass_GATKMethod_MinMAF0.05_Chr40_invB_P_ld.beagle.gz




#Make a list of all *inv*_ld.beagle.gz
ls *inv*_ld.beagle.gz > beagle_ld.list

#Beagle file
BEAGLE=$(cat beagle_ld.list | sed -n ${SLURM_ARRAY_TASK_ID}p)

#Names
NAME=${BEAGLE/.beagle.gz/}

#A beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) can be inputted into ngsLD after the header row and the first three columns (i.e. positions, major allele, minor allele) are removed.

# Preparing the geno file with subsampling for 1 SNP in every 50 snps
#zcat $BEAGLE | awk 'NR % 10 == 0' | gzip > ${NAME}_subsampled.beagle.gz

# Preparing pos file by splitting chromosome and position
zcat $BEAGLE | \
    tail -n +2 |
    awk '{split($1, a, "_"); print a[1] "\t" a[2]}' | \
    gzip > ${NAME}.pos.gz
    

# Counting number of sites
N_SITES=$(zcat ${NAME}.pos.gz | wc -l)

# Print the result information
echo "number of sites $N_SITES"

#Run LD
ngsLD \
--geno $BEAGLE \
--pos ${NAME}.pos.gz \
--probs \
--n_ind 111 \
--n_sites $N_SITES \
--max_kb_dist 0 \
--n_threads 8 \
--out ${NAME}.ld 


#Visualize LD blocks in R
