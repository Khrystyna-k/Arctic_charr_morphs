#!/bin/bash -l

#SBATCH -A naiss2023-5-221
#SBATCH -p core
#SBATCH -n 3
#SBATCH --array=1-285
#SBATCH -t 10:00:00
#SBATCH -J depth
#SBATCH -e depth_%A_%a.err
#SBATCH -o dept_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@imbim.uu.se


#STEP 2: Determine directory 
BAM_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Run_1_2
cd $BAM_DIR



#AVERAGE depth 
output_file=$BAM_DIR/"coverage_summary.txt"

SAMPLE=$(ls *_marked_dups_1_2.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

DIR=$BAM_DIR/${SAMPLE}'_bamqc'

coverage_line=$(grep "mean coverageData =" $DIR/genome_results.txt | awk '{print $1 " " $2 " " $3 " " $4}')


# Append the result to the output file
echo -e "${SAMPLE}_bamqc\t$coverage_line" >> "$output_file"