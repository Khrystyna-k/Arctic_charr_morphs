#!/bin/bash -l

#SBATCH -A naiss2024-5-277
#SBATCH -p core -n 10
#SBATCH -t 03-00:00:00
#SBATCH -J fst_34_M_T_homo
#SBATCH -e fst_34_M_T_homo_%A_%a.err
#SBATCH -o fst_34_M_T_homo_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=khrystyna.kurta@slu.se


#Load modules
module load bioinfo-tools
module load ANGSD/0.933


#Use files from nucleotide diversity directory
NUCL_DIR=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/Nucleotide_diversity

cd $NUCL_DIR

#Two population Fst
#THIS STEP IS FROM NUCLEOTIDE DIVERSITY ALREADT DONE

#STEP 1: first calculate per pop saf for each populatoin -this step was done for nucleotide diversity
#angsd -bam $BAM_LIST -doSaf 1 \
#-ref $REFGENOME -fai $REF_INDEXED \
#-doSaf 1 -GL 1 -P 8 -out ${OUTPUT} \
#-doCounts 1 -setMinDepth 15 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 1 -minQ 20 -remove_bads 1 \
#-uniqueOnly 1 -dumpCounts 2 -doMajorMinor 4 -doMaf 2
##############################################################################################################

#STEP 2:  calculate the 2dsfs prior

# Set the chromosome number using SLURM's array task ID
chr=34


# Assuming FST_LIST is the filename or list containing the paths to the population .list files
POP_LIST=/proj/snic2020-2-19/private/herring/users/khrystyna/Arctic_charr_ref_fSalAlp1/BAM_LISTS/SHARED_list/Shared_34chr/separate_morph_homo

# Convert the list of populations into an array
populations=($(cat $POP_LIST))

# Loop over the populations, avoiding repetition
for ((i = 0; i < ${#populations[@]} - 1; i++)); do
  for ((j = i + 1; j < ${#populations[@]}; j++)); do
      pop_1_name="${populations[$i]}"
      pop_2_name="${populations[$j]}"
      
      
      echo "Comparing pairs: $pop_1_name and $pop_2_name" at chr $chr
      
      # Generate the realSFS .ml file
      realSFS "${pop_1_name}_${chr}.saf.idx" "${pop_2_name}_${chr}.saf.idx" > "${pop_1_name}.${pop_2_name}.${chr}.ml"
      
      # Generate the FST index
      realSFS fst index "${pop_1_name}_${chr}.saf.idx" "${pop_2_name}_${chr}.saf.idx" -sfs "${pop_1_name}.${pop_2_name}.${chr}.ml" -fstout "${pop_1_name}.${pop_2_name}.${chr}"
      
      # Calculate FST stats
      realSFS fst stats "${morph_1_name}.${morph_2_name}.${chr}.fst.idx"
      
      # FST statistics over windows with different sizes and steps
      realSFS fst stats2 "${pop_1_name}.${pop_2_name}.${chr}.fst.idx" -win 20000 -step 10000 > "${pop_1_name}.${pop_2_name}.${chr}_20kb_10kbstep.fst_win"
      realSFS fst stats2 "${pop_1_name}.${pop_2_name}.${chr}.fst.idx" -win 40000 -step 20000 > "${pop_1_name}.${pop_2_name}.${chr}_40kb_20kbstep.fst_win"

  done
done



#STEP 3: prepare the fst for easy window analysis etc
# for morph_1 in $(cat $FST_LIST_1)
# do
#  for morph_2 in $(cat $FST_LIST_2)
#   do
#    if [ "$morph_1" != "$morph_2" ]
#     then
#      echo "Pairs $morph_1 and $morph_2"
#	  realSFS fst index $NUCL_DIR/${morph_1}.saf.idx $NUCL_DIR/${morph_2}.saf.idx -sfs ${morph_1}.${morph_2}.ml -fstout ${morph_1}.${morph_2}
#     fi
#   done
# done


#realSFS fst index pop1.saf.idx pop2.saf.idx -sfs pop1.pop2.ml -fstout out

#STEP 4: get the global estimate
#for fst_files in `ls *fst.idx`; do
#	realSFS fst stats $fst_files
#done


#realSFS fst stats morph_1.morph_2.fst.idx 

#STEP 5: Windowed pairwised Fst
#for fst_files in `ls *fst.idx`; do
#	realSFS fst stats2 $fst_files -win 5000 -step 5000 > ${fst_files/.fst.idx/}.fst_win
#done
	
#realSFS fst stats2 out.fst.idx -win 5000 -step 5000 > fst_win.txt