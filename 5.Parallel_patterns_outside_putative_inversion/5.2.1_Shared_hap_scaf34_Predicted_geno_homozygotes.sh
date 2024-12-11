#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00-07:00:00
#SBATCH --array=1
#SBATCH -J GL_shar
#SBATCH -e GL_shar_%A_%a.err
#SBATCH -o GL_shar_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@imbim.uu.se

#Load modules
module load bioinfo-tools
module load samtools/1.12
module load bamtools/2.5.1
module load ANGSD/0.933


######################################################################################
# CONDUCT GENOTYPE LIKELIHOOD ESTIMATION USING ANGSD v.0.930
# Will output the following files per chromosome:
# 1. VCF (with PL field included)
# 2. Beagle (genotype likelihood format)
# 3. MAFs (allele frequencies)

# A. Khrystyna Kurta, November 2022
######################################################################################

#Index sites 
SITES=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO_SHARED/sites_homo_shared_filter

#angsd sites index $SITES

#Path to the directory where you have the bam-files
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2
cd $BASEDIR

#STEP 1: Define paths to Refference genome
REFGENOME=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta
REF_INDEXED=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta.fai

CHUNK_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/CHUNK_LIST


##STEP 2: Determine chromosome/ or Get all the contig (or scaffold) names from the reference genome fasta file
CHUNK_NAMES=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/CHUNK_LIST

#CHUNK_NAMES_target=$(ls $CHUNK_NAMES/shared_34chr*.region)
#CHUNK_NAMES_target_name=${CHUNK_NAMES_target/.region/}

#STEP 3: Create bam file list
#Text file containing sample bam paths
BAM_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/SHARED_list/Shared_34chr

BAM=$(ls $BAM_LIST/Myv_Thin_all_homo.list | sed -n ${SLURM_ARRAY_TASK_ID}p)
BAM_LIST_NAME=$(basename $BAM)

#ls /proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2/*.bam > all_bam_path.list
#BAM_LIST=all_bam_path.list

OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO_SHARED


#STEP 5: Run ANGSD

for CHUNK in `ls $CHUNK_NAMES/shared_34chr_outside.region`; do

CHUNK_NAMES_target_name=$(basename $CHUNK)
CHUNK_NAMES_target_name=${CHUNK_NAMES_target_name/.region/}

echo "Run angsd for $CHUNK and bam list $BAM"

angsd -b $BAM \
-ref $REFGENOME -fai $REF_INDEXED \
-rf $CHUNK \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-gl 2 -trim 0 -doMajorMinor 4 -domaf 1 -doPost 2 -doGlf 2 -minMaf 0.05 -SNP_pval 1e-6 -docounts 1 -dogeno 2 \
-out $OUT_DIR/${BAM_LIST_NAME}_${CHUNK_NAMES_target_name}_Maf0.05 -P 10 \
-nThreads 10
# -sites $SITES
done

#This argument  -minMaf 0.05 -SNP_pval 1e-6  is not added because the SITE are alredy contain the essential regions

# Explanation of above settings:
# ==============================
# -uniqueOnly = only use uniquely mapped reads (ingnore reads with multiple hits)
# -remove_bads = same as the samtools flags -x which removes read with a flag above 255 (not primary, failure and duplicate reads)
# -only_proper_pairs = include only pairs of reads with both mates (forward and reverse) mapped correctly
# -minMapQ = minimum mapQ quality
# -minQ = minimum base quality score
# -GL = calculate genotype likelihoods (2: using GATK model )
# -doMajorMinor = infer major and minor alleles (1: from GLs)
# -doPost = calculate posterior prob (1: Using frequency as prior)
# -doBcf = output a VCF file (1: yes)
# -doGlf = ouput genotype likelihoods (2: beagle likelihood format)
# -minMaf = minumum minor allele frequency tolerated
# -SNP_pval = significance threshold for determining true polymorphism
#-doGeno 0 1: write major and minor, 2: write the called genotype encoded as -1,0,1,2, -1=not called
# ==============================


