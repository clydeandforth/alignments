#!/bin/bash

# Specify the file name to write stdout to; %J is replaced by the JOBID
#SBATCH -o Tree_prep.o.%J
#SBATCH -J Biopython
##SBATCH --exclusive
#Specify the file name to write stderr to; %J is replaced by the JOBID
#SBATCH -e Tree_prep.e.%J
##SBATCH -t 0-08:00
##SBATCH -N 1
##SBATCH --ntasks=10
#SBATCH -p dev

res1=$(date +%s.%N)

module purge
module load hpcw
module load clustalw/2.1

#conda activate ashconda
#conda init ETE-toolkit

/scratch/b.bss81c/alignments/multi_complete_script_for_nexus.py 
#/scratch/b.bss81c/alignments/multi_complete_script_with_individual_gene_names.py
/scratch/b.bss81c/alignments/Meld2Nexus/Meld2Nexus -c outfile_padded.nexus outfile_padded_2.nexus outfile_padded_3.nexus outfile_padded_4.nexus outfile_padded_5.nexus outfile_padded_6.nexus -o combined.nex
#sed -i 's/datatype=rna/datatype=dna/' combined.nex
#/scratch/b.bss81c/alignments/other_tree.py
#/scratch/b.bss81c/alignments/ete_temp.py
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
