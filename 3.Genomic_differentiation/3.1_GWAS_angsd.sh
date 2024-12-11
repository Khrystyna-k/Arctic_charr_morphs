#!/bin/bash -l

#SBATCH -A naiss2023-5-221
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 09-00:00:00
#SBATCH --array=1-11
#SBATCH -J asso_no_cov
#SBATCH -e asso_no_cov_%A_%a.err
#SBATCH -o asso_no_cov_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@imbim.uu.se

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load PCAngsd/1.11


######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION USING ANGSD v.1.11
# Will output the following files per chromosome:
# 1. VCF (with PL field included)
# 2. Beagle (genotype likelihood format)
# 3. MAFs (allele frequencies)

# Khrystyna Kurta, March 2023
######################################################################################
#Go to the bam files directory

#Path to the directory where you have the bam-files
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2
cd $BASEDIR
#Specify which bam list to use
BAM_LIST_PATH=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/GWAS_list

#Sites
SITES=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/SITES/sites_sorted.txt


#Specify which bam to use
BAM_LIST=$(ls $BAM_LIST_PATH/assos_*.list | sed -n ${SLURM_ARRAY_TASK_ID}p)
BAM_TARGET=${BAM_LIST/.list/}
OUTPUT=$(basename $BAM_TARGET)

#Specify the number of individuals
x=`cat $BAM_LIST | wc -l`


#STEP 1: Define paths to Refference genome
# Update this with the name of the ref fasta file
REFGENOME=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta
REF_INDEXED=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta.fai


#Out directory
OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GWAS


#REGION=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_VK-3477/mapped_reads/contig_list_aa

#STEP 2: The association can then be performed on the genotype probabilities using the score statistics

angsd -yBin $BAM_TARGET.ybin -doAsso 1 -out $OUT_DIR/$OUTPUT \
-ref $REFGENOME -fai $REF_INDEXED \
-doMajorMinor 4 -doMaf 1 -bam $BAM_LIST -nInd $x \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-gl 2 \
-Pvalue 1 \ #export P values
-sites $SITES \
-nThreads 8