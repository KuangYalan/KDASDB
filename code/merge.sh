#!/bin/bash
#SBATCH -p fat
#SBATCH -J Merge
#SBATCH -N 1
#SBATCH -n 10


#rm Stringtie -r
mkdir Stringtie
cd Stringtie

# find all lncRNA data gtf
#sample1=( `ls *.gtf | xargs -n 1 basename -s .gtf` )
sample1=( `find ../../3-Assembly/ -name '*.gtf' |grep Stringtie` )
for i in "${sample1[@]}"
do
echo $i >> gtf_list.txt
done

sample2=( `find ../../3-Assembly/ -name '*.gtf' |grep transcripts` )
for i in "${sample2[@]}"
do
echo $i >> gtf_list.txt
done

# stringtie merge
stringtie --merge -p 1 -G /public/Database/hg38/Genome/Genes.gtf -o merged.gtf gtf_list.txt
cd ..
