#!/bin/bash
#SBATCH -p cn
#SBATCH -J hisat2
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o hisat2.o
#SBATCH -e hisat2.e

#find -type d | xargs rm -r
genome=hg38

# find all RNA data samples
for study in `ls ~/get_data/Homo`
do
    if [ -d "./$study" ];then
        echo "$study alread exist!"
        continue
    fi
    ReadDIR=/public/home/kuangyl/get_data/Homo/$study
    cd $ReadDIR
    for fastq in `ls ./`
    do
       echo $fastq
       if echo "$fastq" | grep -q -E '\.1.fastq.gz|\.2.fastq.gz'
       then
           end_with_fastq=false
           break
       else
           end_with_fastq=true
           break
       fi
    done 
    echo $end_with_fastq
    if $end_with_fastq
    then
        sample=( `ls *.fastq.gz | xargs -n 1 basename -s .fastq.gz` )
    else
        sample=( `ls *.1.fastq.gz | xargs -n 1 basename -s .1.fastq.gz` )
    fi
    echo $sample
    cd /public/home/kuangyl/RNA-seq/1-Alignment/hisat2
    rm -rf $study
    mkdir $study
    cd $study
    # reads alignment by hisat2
    mkdir logs
    mkdir Bams
    for i in "${sample[@]}"
    do
        if $end_with_fastq
        then
            hisat2 -p 40 -x /public/Database/$genome/Genome/Genome -U $ReadDIR/"$i".fastq.gz -S "$i".sam 1>./logs/"$i".out 2>./logs/"$i".err
        else
            hisat2 -p 40 -x /public/Database/$genome/Genome/Genome -1 $ReadDIR/"$i".1.fastq.gz -2 $ReadDIR/"$i".2.fastq.gz -S "$i".sam 1>./logs/"$i".out 2>./logs/"$i".err
        samtools sort -@40 -o Bams/"$i".bam "$i".sam
        fi
    done

    cd Bams
    ls *.bam | xargs -n 1 -P 40 samtools index
done
