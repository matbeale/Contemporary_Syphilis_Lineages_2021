#!/bin/bash


# Usage: fastq_to_perfect_reads.sh <fasta> <fastq-basename>


if [ $# -ne 2 ]
then
 echo "Usage: shred_fasta_to_perfect_reads.sh <fasta> <fastq-basename>"
 exit 1
fi

module load fastaq/3.17.0-docker3

insertsize=350
readlength=150
coverage=50


myfasta=$1
myfastq=$2

fastaq to_perfect_reads $myfasta $myfastq\.fastqs $insertsize 50 $coverage $readlength

fastaq deinterleave $myfastq\.fastqs $myfastq\_1.fastq $myfastq\_2.fastq

gzip $myfastq\_*.fastq

rm $myfastq\.fastqs
