#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 16
#SBATCH -t 00-12:00:00
#SBATCH --array=1-4
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
CHUNK_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/CHUNK_LIST

CHUNK_NAMES=$CHUNK_DIR/34.txt
CHUNK_NAMES_target=$(basename $CHUNK_NAMES)
CHUNK_NAMES_target_name=${CHUNK_NAMES_target/.txt/}

#Path to the directory where you have the bam-files
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2

#Text file containing sample bam paths
BAM_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/SHARED_list/Shared_34chr

#Reference genome
REFGENOME=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta
REF_INDEXED=/proj/snic2020-2-19/private/arctic_charr/assemblies/fSalAlp1.1/Data_Package_fSalAlp1_assembly_20231024/fSalAlp1.1.hap1.cur.20231016.fasta.fai

#Create a directory
OUT_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Nucleotide_diversity

#Step 1: This option creates a fasta.gz file from a sequencing data file (BAM file). For the consensus fasta all letters are capital letters
#angsd -bam $BAM_LIST -dofasta 2 -doCounts 1 -out myFasta #do I need it?

cd $BASEDIR

POP=$(cat $BAM_LIST/separate_morph_homo | sed -n ${SLURM_ARRAY_TASK_ID}p)
POP_2=$(ls $BAM_LIST/$POP)
OUTPUT=$(basename $POP_2)

#Step 2: Finding a 'global estimate' of the SFS per each lake
angsd -bam $POP_2 \
-anc $REFGENOME -fai $REF_INDEXED \
-rf $CHUNK_DIR/$CHUNK_NAMES_target \
-doSaf 1 -GL 2 -P 8 -out $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name} \
-doCounts 1 -setMinDepth 15 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 30 -minQ 20 -remove_bads 1 -minInd 5 \
-uniqueOnly 1 -dumpCounts 2

#Step 3: Obtain the folded site frequency spectrum 
realSFS $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.saf.idx -P 8 -fold 1 > $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.sfs


#Step 4: Calculate the thetas for each site
realSFS saf2theta $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.saf.idx -fold 1 -sfs $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.sfs -outname $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}


#Step 5: Estimate Tajimas D and other statistics do a sliding window analysis by adding -win/-step arguments to the last command
thetaStat do_stat $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.thetas.idx -win 5000 -step 5000 -outnames $OUT_DIR/${OUTPUT}_${CHUNK_NAMES_target_name}.theta.5kb.thetasWindow.gz

