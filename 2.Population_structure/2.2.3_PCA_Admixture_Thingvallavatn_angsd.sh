#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p node -n 10
#SBATCH -t 1-00:00:00
#SBATCH --array=1-10
#SBATCH -J thin_Adm
#SBATCH -e thin_Adm_%A_%a.err
#SBATCH -o thin_Adm_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@slu.se


#Load modules
module load bioinfo-tools
module load NGSadmix/32
module load PCAngsd/1.11


# Khrystyna Kurta, April 2023
######################################################################################

#BEAGLE file dir
BASEDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/GL_MAF_GENO
cd $BASEDIR

#Beagle file list subsampled in the previous step for 1 SNP in every 10 SNP in the beagle file
BEAGLE=Thin_CharrLowPass_GATKMethod_MinMAF0.05_no_unloc_subsampled.beagle.gz

#Set up directory
OUTDIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Admixture

#Assigns K
K=$SLURM_ARRAY_TASK_ID

#Run NGSadmix
NGSadmix -likes $BEAGLE -K $K -P 10 -o $OUTDIR/Thin_subs_ngxadmix_$K -minMaf 0.05 #-minInd 10

#Get errors and log likelihoods
# Output file
output_file="summary_ngsadmix_Thin.csv"

# Write the header to the output file
echo "file_name,Log-likelihood,Frobenius error" > "$output_file"

# Loop through each .out file
for file in Thin_subs_ngxadmix*.log; do
  # Extract Log-likelihood and Frobenius error using grep and awk
  log_likelihood=$(grep "best like" "$file" | awk '{print $2}')
  #frobenius_error=$(grep "Frobenius error" "$file" | awk '{print $3}')
  
  # Write the results to the output file
  echo "$file,$log_likelihood" >> "$output_file"
done

#Run PCAngsd to select the best fit K
pcangsd -b $BEAGLE --P 10 --admix -o $OUTDIR/Thin_subs_pcangsd_autoK

#Run for every K to score statistics Log-likelihood and Frobenius error
pcangsd -b $BEAGLE --P 10 --admix --admix_K $K -o $OUTDIR/Thin_subs_pcangsd_${K}

# Output file
output_file="summary_pcangsd_Thin.csv"

# Write the header to the output file
echo "file_name,Log-likelihood,Frobenius error" > "$output_file"

# Loop through each .out file
for file in thin_AdmSub_*.out; do
  # Extract Log-likelihood and Frobenius error using grep and awk
  log_likelihood=$(grep "Log-likelihood" "$file" | awk '{print $2}')
  frobenius_error=$(grep "Frobenius error" "$file" | awk '{print $3}')
  
  # Write the results to the output file
  echo "$file,$log_likelihood,$frobenius_error" >> "$output_file"
done
