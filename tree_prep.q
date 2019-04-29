#!/bin/bash

# Specify the file name to write stdout to; %J is replaced by the JOBID
#SBATCH -o Tree_prep.o.%J
#SBATCH -J Biopython
#SBATCH --exclusive
#Specify the file name to write stderr to; %J is replaced by the JOBID
#SBATCH -e Tree_prep.e.%J
##SBATCH -t 2-60:00
##SBATCH -N 1
##SBATCH --ntasks=40
##SBATCH -p dev

res1=$(date +%s.%N)

module purge
module load hpcw
module load clustalw/2.1

/scratch/b.bss81c/alignments/complete_script.py /scratch/b.bss81c/Prokka/concat_file.gbk outfile.fasta 'DNA gyrase subunit B'

res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
