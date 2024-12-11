#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 20
#SBATCH -t 01-00:00:00
#SBATCH --array=1-5
#SBATCH -J sfs
#SBATCH -e sfs_%A_%a.err
#SBATCH -o sfs_%A_%a.out
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
CHUNK_NAMES=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/CHUNK_LIST/chr_ordered_txt.list
CHUNK_NAMES_target=$(cat $CHUNK_NAMES | sed -n ${SLURM_ARRAY_TASK_ID}p)
CHUNK_NAMES_target_name=${CHUNK_NAMES_target/.txt/}

#Path to the directory where you have the bam-files
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2

#Text file containing sample bam paths
BAM_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/FST_list

#Reference genome
REFGENOME=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta
REF_INDEXED=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta.fai

#Create a directory
OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/FST_folded

#Step 1: This option creates a fasta.gz file from a sequencing data file (BAM file). For the consensus fasta all letters are capital letters
#angsd -bam $BAM_LIST -dofasta 2 -doCounts 1 -out myFasta #do I need it?

cd $BASEDIR

for POP in `ls $BAM_LIST/pop_*list`; do

OUTPUT=$(basename $POP)

#Step 2: Finding a 'global estimate' of the SFS per each lake
angsd -bam $POP \
-anc $REFGENOME -fai $REF_INDEXED \
-rf $CHUNK_DIR/$CHUNK_NAMES_target \
-doCounts 1 -setMinDepth 15 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 30 -minQ 20 -remove_bads 1 -minInd 10 \
-uniqueOnly 1 -dumpCounts 2 \
-doSaf 1 -GL 2 -P 8 -out $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}


#Step 3: Obtain the folded site frequency spectrum 
realSFS $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.saf.idx -P 8 -fold 1 > $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.sfs

done

