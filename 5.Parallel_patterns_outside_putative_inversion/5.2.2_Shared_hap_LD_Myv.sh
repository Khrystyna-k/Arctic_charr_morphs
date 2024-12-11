#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 10
#SBATCH -t 02-00:00:00
#SBATCH --array=1-2
#SBATCH -J LD_MyvThin
#SBATCH -e LD_MyvThin_%A_%a.err
#SBATCH -o LD_MyvThin_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@slu.se

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load ngsLD/1.1.1

#Directorys
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO
BEAGLE_MAIN=$BASEDIR/Myv_CharrLowPass_GATKMethod_MinMAF0.05_all_chr.beagle.gz

OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO_SHARED
cd $OUT_DIR

#Chr34
{
    zcat $BEAGLE_MAIN | head -n 1  # Extract the header
    zcat $BEAGLE_MAIN | awk -F'[_\t]' '$1 == "34" && $2 >= 18300000 && $2 <= 18500000' 
} | gzip > Myv_CharrLowPass_GATKMethod_MinMAF0.05_shared_34chr_18.3_18.5Mb.beagle.gz

ls Myv_CharrLowPass_GATKMethod_MinMAF0.05_shared_34chr_18.3_18.5Mb.beagle.gz Thin_CharrLowPass_GATKMethod_MinMAF0.05_shared_34chr_18.3_18.5Mb.beagle.gz > shared_b.list

BEAGLE=$(cat shared_b.list | sed -n ${SLURM_ARRAY_TASK_ID}p)

#Names
NAME=${BEAGLE/.beagle.gz/}

#A beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) can be inputted into ngsLD after the header row and the first three columns (i.e. positions, major allele, minor allele) are removed.

# Preparing pos file by splitting chromosome and position
zcat $BEAGLE | \
    tail -n +2 |
    awk '{split($1, a, "_"); print a[1] "\t" a[2]}' | \
    gzip > ${NAME}.pos.gz
    

# Counting number of sites
N_SITES=$(zcat ${NAME}.pos.gz | wc -l)

zcat $BEAGLE | wc -l
zcat ${NAME}.pos.gz | wc -l
echo $N_SITES

# Print the result information
echo "For BEAGLE ${NAME}_subsampled.beagle.gz  positions ${NAME}.pos.gz with number of sites $N_SITES"

#Run LD
ngsLD \
--geno $BEAGLE \
--pos ${NAME}.pos.gz \
--probs \
--n_ind 54 \
--n_sites $N_SITES \
--max_kb_dist 0 \
--n_threads 8 \
--out ${NAME}.ld
#--min_maf 0.1


#Visualize LD blocks in R
