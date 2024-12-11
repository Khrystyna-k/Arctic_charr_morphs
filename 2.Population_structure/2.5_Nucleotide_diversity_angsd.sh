#!/bin/bash -l

#SBATCH -A naiss2023-5-222
#SBATCH -p core -n 20
#SBATCH -t 10-00:00:00
#SBATCH --array=1-40
#SBATCH -J nucl_div
#SBATCH -e nucl_div_%A_%a.err
#SBATCH -o nucl_div_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@imbim.uu.se

#Load modules
module load bioinfo-tools
module load ANGSD/0.933
module load samtools

#SET UP CHUNK list

#Specify which bam to use
CHUNK_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/CHUNK_LIST


##STEP 2: Determine chromosome/ or Get all the contig (or scaffold) names from the reference genome fasta file
CHUNK_NAMES=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/CHUNK_LIST/chr_nounloc.list
CHUNK_NAMES_target=$(cat $CHUNK_NAMES | sed -n ${SLURM_ARRAY_TASK_ID}p)
CHUNK_NAMES_target_name=${CHUNK_NAMES_target/.txt/}

#Path to the directory where you have the bam-files
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2

#Text file containing sample bam paths
BAM_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/Nucleotide_diversity_list

#Reference genome
REFGENOME=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta
REF_INDEXED=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta.fai

#Create a directory
OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Nucleotide_diversity

#Step 1: This option creates a fasta.gz file from a sequencing data file (BAM file). For the consensus fasta all letters are capital letters
#angsd -bam $BAM_LIST -dofasta 2 -doCounts 1 -out myFasta #not used in this study

cd $BASEDIR

for POP in `ls $BAM_LIST/pop_*list`; do

OUTPUT=$(basename $POP)

#Step 2: Finding a 'global estimate' of the SFS per each lake
angsd -bam $POP \
-anc $REFGENOME -fai $REF_INDEXED \
-rf $CHUNK_DIR/$CHUNK_NAMES_target \
-doSaf 1 -GL 2 -P 8 -out $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name} \
-doCounts 1 -setMinDepth 15 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 30 -minQ 20 -remove_bads 1 -minInd 10 \
-uniqueOnly 1 -dumpCounts 2

#Step 3: Obtain the folded site frequency spectrum 
realSFS $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.saf.idx -P 8 -fold 1 > $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.sfs

#Step 4: Calculate the thetas for each site
realSFS saf2theta $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.saf.idx  -fold 1 -sfs $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.sfs -outname $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}

#Step 5: Estimate Tajimas D and other statistics do a sliding window analysis by adding -win/-step arguments to the last command
thetaStat do_stat $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.thetas.idx -win 50000 -step 50000 -outnames $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.theta.50kb.thetasWindow.gz
thetaStat do_stat $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.thetas.idx -win 5000 -step 5000 -outnames $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.theta.5kb.thetasWindow.gz
thetaStat do_stat $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.thetas.idx -win 1000 -step 1000 -outnames $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.theta.10kb.thetasWindow.gz

done





#Generate a consensus fasta file - I didn;t use it in this run
#	-doFasta	0
#	1: use a random (non N) base (needs -doCounts 1)
#	2: use the most common (non N) base (needs -doCounts 1)

#Filters
#1) -setMinDepth [int] - Discard site if total sequencing depth (all individuals added together) is below [int]. Requires -doCounts
#2) -setMaxDepth [int] - Discard site if total sequencing depth (all individuals added together) is above [int] -doCounts
#3) -setMinDepthInd [int]
#Discard individual if sequencing depth for an individual is below [int]. This filter is only applied to analysis which are based on counts of alleles i.e. analysis that uses -doCounts
#4)-dumpCounts	0
#	  1: total seqdepth for site	.pos.gz
#	  2: seqdepth persample		.pos.gz,.counts.gz
#	  3: A,C,G,T sum all samples	.pos.gz,.counts.gz
#	  4: A,C,G,T sum every sample	.pos.gz,.counts.gz
#5)-doMajorMinor	0
#	4: Use reference allele as major (requires -ref)

#6)-doMaf	0 (Calculate persite frequencies '.mafs.gz')
#	1: Frequency (fixed major and minor)
#	2: Frequency (fixed major unknown minor)

#Its a 3 step procedure
#1. Estimate an site frequency spectrum. Output is out.sfs file. This is what is being used as the -pest argument in step2.
#2.Calculate per-site thetas. Output is a .thetas.idx/.thetas.gz files. This contains the binary persite estimates of the thetas.
#3.Calculate neutrality tests statistics. Output is a .thetas.idx.pestPG file.