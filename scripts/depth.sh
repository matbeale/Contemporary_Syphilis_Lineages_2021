#!/bin/bash


# Usage: depth.sh <bamfile>


module load samtools/1.6--h244ad75_4 

outfile=$(echo $1 | perl -pe 's/\.bam//g')

samtools depth -aa $1 > $outfile\.depth


