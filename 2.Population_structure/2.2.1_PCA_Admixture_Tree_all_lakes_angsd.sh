#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 10
#SBATCH -t 02-00:00:00
#SBATCH --array=1-10
#SBATCH -J AdmSubsampledG
#SBATCH -e AdmSubsampledG_%A_%a.err
#SBATCH -o AdmSubsampledG_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@slu.se


#Load modules
module load bioinfo-tools
module load PCAngsd/1.11
module load NGSadmix/32


#Assigns K
K=$SLURM_ARRAY_TASK_ID

#Dir
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO

## Prepare a geno file by subsampling one SNP in every 10 SNPs in the beagle file and removing unlocated scaffolds
zcat $BASEDIR/Charr_RefUK_MAF0.05_all_chr_with_header_2.beagle.gz | sed '/^unloc/d' | gzip > $BASEDIR/Charr_RefUK_MAF0.05_no_unloc_with_header_2.beagle.gz

zcat $BASEDIR/Charr_RefUK_MAF0.05_no_unloc_with_header_2.beagle.gz | { head -n 1; awk 'NR > 1 && NR % 10 == 0'; } | gzip > $BASEDIR/Charr_RefUK_MAF0.05_no_unloc_with_header_2_subsampled.beagle.gz

#OUTDIR 
OUTDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Admixture

#sample list for a nice tree
SAMPLE_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/BEAGLE_list/tree_sample.list

#Run the admix program and used K the number of ancestral popuations defined with the PCA
BEAGLE=$BASEDIR/Charr_RefUK_MAF0.05_no_unloc_with_header_2_subsampled.beagle.gz
BEAGLE_NAME=All_lakes_subsampled_k5

#Allow for K selected by PCAngsd
pcangsd --beagle $BEAGLE --threads 10 --admix --out $OUTDIR/${BEAGLE_NAME}_auto --tree --tree_samples $SAMPLE_LIST

#Set 1-10 K
pcangsd --beagle $BEAGLE --threads 10 --admix --admix_K $K --out $OUTDIR/${BEAGLE_NAME}

# Output file
output_file="summary_PCAangsd_subs_All_lakes.csv"

# Write the header to the output file
echo "file_name,Log-likelihood,Frobenius error" > "$output_file"

# Loop through each .out file
for file in AdmSubsampledG_*.out; do
  # Extract Log-likelihood and Frobenius error using grep and awk
  log_likelihood=$(grep "Log-likelihood" "$file" | awk '{print $2}')
  frobenius_error=$(grep "Frobenius error" "$file" | awk '{print $3}')
  
  # Write the results to the output file
  echo "$file,$log_likelihood,$frobenius_error" >> "$output_file"
done


#Run NGSAdmix with defined K
NGSadmix -likes $BEAGLE -K $K -P 10 -o $OUTDIR/All_lakes_subs_ngxadmix_$K -minMaf 0.05 -minInd 10

#Get errors and log likelihoods
# Output file
output_file="Errors_ngsadmix_subs.csv"

# Write the header to the output file
echo "file_name,Log-likelihood" > "$output_file"

# Loop through each .out file
for file in All_lakes_subs_ngxadmix_$K.log; do
  # Extract Log-likelihood and Frobenius error using grep and awk
  log_likelihood=$(grep "best like" "$file" | awk '{print $2}')
  #frobenius_error=$(grep "Frobenius error" "$file" | awk '{print $3}')
  
  # Write the results to the output file
  echo "$file,$log_likelihood" >> "$output_file"
done
