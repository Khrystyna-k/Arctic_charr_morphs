#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 10
#SBATCH -t 02-00:00:00
#SBATCH --array=1-4
#SBATCH -J AdmSubs
#SBATCH -e AdmSubs_%A_%a.err
#SBATCH -o AdmSubs_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@slu.se


#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load PCAngsd/1.11


#Dir
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO
cd $BASEDIR

#OUTDIR 
OUTDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Admixture


#Make a list of files GL per every lake
ls Van_CharrLowPass_GATKMethod_MinMAF0.05_no_unloc_with_header_2.beagle.gz Myv_CharrLowPass_GATKMethod_MinMAF0.05_no_unloc_with_header_2.beagle.gz Sir_CharrLowPass_GATKMethod_MinMAF0.05_no_unloc_with_header_2.beagle.gz Thin_CharrLowPass_GATKMethod_MinMAF0.05_no_unloc.beagle.gz > sub_beagle.list

#Assign ech file to variable
BEAGLE=$(cat subs_beagle.list | sed -n ${SLURM_ARRAY_TASK_ID}p)

#SUBSAMPLE BEAGLE 1 SNP in every 10 SNPs  in the beagle file
zcat $BEAGLE | { head -n 1; awk 'NR > 1 && NR % 10 == 0'; } | gzip > ${BEAGLE}_subsampled
NAME=${BEAGLE/.beagle.gz/}

#Run PCAngsd to select the best supported K automatically
pcangsd -b ${BEAGLE}_subs --admix --threads 10 -o $OUTDIR/${NAME}_subsampled

