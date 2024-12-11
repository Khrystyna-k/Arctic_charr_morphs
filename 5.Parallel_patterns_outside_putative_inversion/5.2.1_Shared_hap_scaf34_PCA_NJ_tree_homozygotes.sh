#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 10
#SBATCH -t 0-02:00:00
#SBATCH --array=1-2
#SBATCH -J nj_all
#SBATCH -e nj_all_%A_%a.err
#SBATCH -o nj_all_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@imbim.uu.se


#Load modules
module load bioinfo-tools
module load PCAngsd/1.11

# Khrystyna Kurta, April 2023
######################################################################################

#BEAGLE file dir
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO_SHARED

#STEP 2: Go to the bam files directory
cd $BASEDIR

#Make a list of outside and within shared region
ls Myv_Thin_all_homo.list_shared_34chr_outside_Maf0.05.beagle.gz Myv_Thin_all_homo.list_shared_34chr_18.1_18.5Mb.beagle.gz > list.txt

#Beagle file list
BEAGLE=$(cat list.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
NAME=${BEAGLE/.beagle.gz/}

#Use corresponding sample list
SAMPLE_NAME=${BEAGLE/_shared_34chr_18.1_18.5Mb.beagle.gz/}
SAMPLE_NAME=${SAMPLE_NAME/_shared_34chr_outside_Maf0.05.beagle.gz/}

SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/SHARED_list/Shared_34chr/$SAMPLE_NAME

#Set up directory
OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/PCA_NJ_tree


#Run PCAngsd
pcangsd -b $BEAGLE --tree --tree_samples $SAMPLE_LIST -o $OUT_DIR/$NAME


