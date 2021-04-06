#!/bin/bash


# Usage: run_Metagenomic_reads_downsampling_pipeline.sh [options: [-t] [-c]] <sampleid> <xxx_1.fastq.gz> <xxx_2.fastq.gz> <resultsdir>

# Takes raw fastq files and performs read assignment with kraken2, then bins the reads according to defined criteria (e.g. 'Treponema'), followed by downsampling to a specific upper count (note that samples with fewer reads will be labelled as downsampled even if the read counts remain the same). Note that the kraken database is ~36GB (or 8Gb for minikraken), but is not complete - the pipeline will not work well for samples with limited sequence information, highly diverse genomes, or for genomes with novel gene content"

#  Does not need to be bsub-ed!!!


##########################################################
# Specify default options for flags (this was originally written for Treponema)'
readcount=2500000
binning_term='Treponema'
kbinmem=32000
kdbtype='FullKraken'
invert_search='no'
trim='TRUE'
do_kraken='run_kraken'

# Capture command line options and positional arguments
while getopts ":ht:vc:b:dT" opt; do
  case "${opt}" in
    h) echo ""
       echo "Uses kraken2 to classify and bin raw PE fastqs from a species group, then trims and downsamples to a manageable number for assembly/mapping."
       echo "Note: Jobs are submitted to LSF automatically, this does not need to be Bsub-ed"
       echo ""
       echo "Usage: "
       echo "    `basename $0` [options] <sampleid> <xxx_1.fastq.gz> <xxx_2.fastq.gz> <resultsdir>" 
       echo "Options: "
       echo "    -t  -  Specify species to bin reads by [$binning_term]"
       echo "    -v  -  Invert search to remove human reads (no downsampling performed after trimming) [sapiens]"
       echo "    -c  -  Specify max number of reads to retain after binning [$readcount]"
       echo "    -T  -  Binary Switch to turn off read trimming [default:$trim]"
       echo ""
       echo "    -d  -  Binary Switch between databases: FullKraken or MiniKraken8Gb [default:$kdbtype] "
       echo "    -b  -  Specify memory for binning steps [$kbinmem]"
       echo ""
       exit 0
      ;;
    t) binning_term="${OPTARG}" ;;
    c) readcount="${OPTARG}" ;;
    b) kbinmem="${OPTARG}" ;;
    d) kdbtype='MiniKraken8Gb' ;;
    v) invert_search='invert' ;;
    T) trim='FALSE' ;;
    \?)echo "Invalid Option: -$OPTARG" 1>&2
       exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# All 4 positional arguments are mandatory - kill script if any are absent
if [ $# -ne 4 ]
then
 echo "Invalid number of positional arguments "
 echo "Usage: `basename $0`[options] <sampleprefix> <xxx_1.fastq.gz> <xxx_2.fastq.gz> <resultsdir>"
 exit 1
fi

# Print optional flags to screen
#echo "Binning reads classified as $binning_term and downsampling to $readcount reads or less"
#########################################################

#  Does not need to be bsub-ed
sampleid=$1
read1=$2
read2=$3
resultsdir=$4

mkdir -vp $resultsdir

#module load seqtk
module load seqtk/1.3--ha92aebf_0
#module load trimmomatic
module load trimmomatic/0.39--1
#module load kraken2/2.0.8_beta
module load kraken2/2.0.8_beta=pl526h6bb024c_0-c1
#kraken2path=/software/pathogen/external/apps/usr/bin/kraken2
kraken2path=kraken2


if [[ $kdbtype = FullKraken ]]
then
  kraken2db=/lustre/scratch118/infgen/team216/mb29/databases/kraken/kraken2/standard_Kraken2_2019-03-29/standard_Kraken2_2019-03-29.db
  kmem=38000
else
  kraken2db=/lustre/scratch118/infgen/team216/mb29/databases/kraken/kraken2/minikraken2_v2_8GB__2018-10-28
  kmem=10000
fi

# For minikraken database: 
#kraken2db=/lustre/scratch118/infgen/team216/mb29/databases/kraken/kraken2/minikraken2_v2_8GB__2018-10-28
#kmem=10000

# For standard kraken database:
#kraken2db=/lustre/scratch118/infgen/team216/mb29/databases/kraken/kraken2/standard_Kraken2_2019-03-29/standard_Kraken2_2019-03-29.db
#kmem=38000
#########################################################

if [[ $do_kraken = run_kraken ]]
then
  dokraken2=$(bsub -J kraken_$sampleid\_1 -o $resultsdir/k_$sampleid\_1.%J.o -e $resultsdir/k_$sampleid\_1.%J.e -q normal -n 8 -M$kmem -R "span[hosts=1] select[mem>$kmem] rusage[mem=$kmem]" "echo $kraken2db ; $kraken2path --threads 8 --gzip-compressed --paired --use-names --output $resultsdir/$sampleid\.kraken.labels --report $resultsdir/$sampleid\.kraken.report -db $kraken2db $read1 $read2")
  dokraken2_jobID=$(echo $dokraken2 |perl -pe 's/\>.+$//g'|perl -pe 's/^.+\<//g')
  echo "Kraken2 job submitted as jobID $dokraken2_jobID for $sampleid reads 1 & 2 using $kdbtype database and $kmem Mb"
fi


# Combine the kraken binning and downsampling into same job (reduces number of jobs needed) and also include a read trimming step after
# This job will submit with the others, but will remain pending until completion of the kraken runs

# First look at doing inverted search
if [[ $invert_search = invert ]]
then
	  dokbin=$(bsub -w "ended($dokraken2_jobID)" -J kbin_$sampleid -o $resultsdir/kbin_$sampleid.%J.o -e $resultsdir/kbin_$sampleid.%J.e -q normal -n 2 -M$kbinmem -R "span[hosts=1] select[mem>$kbinmem] rusage[mem=$kbinmem]" "
        # Extract reads of interest using kraken outputs and seqtk
        grep -v sapiens $resultsdir/$sampleid\.kraken.labels | awk '{print \$2}' > $resultsdir/$sampleid\.not-sapiens.list ;
        cat $resultsdir/$sampleid\.not-sapiens.list | perl -pe 's/\n/\/1\n/g' > $resultsdir/$sampleid\.not-sapiens.R1.list ;
        cat $resultsdir/$sampleid\.not-sapiens.list | perl -pe 's/\n/\/2\n/g' > $resultsdir/$sampleid\.not-sapiens.R2.list ;
        seqtk subseq $read1 $resultsdir/$sampleid\.not-sapiens.R1.list | gzip -1 > $resultsdir/$sampleid\.not-sapiens.binned_1.fastq.gz ;
        seqtk subseq $read2 $resultsdir/$sampleid\.not-sapiens.R2.list | gzip -1 > $resultsdir/$sampleid\.not-sapiens.binned_2.fastq.gz ;

        # Trim binned reads
        trimmomatic PE -phred33 -threads 2 $resultsdir/$sampleid\.not-sapiens.binned_1.fastq.gz $resultsdir/$sampleid\.not-sapiens.binned_2.fastq.gz $resultsdir/$sampleid\.not-sapiens.paired_1.fastq.gz $resultsdir/$sampleid\.not-sapiens.unpaired_1.fastq.gz $resultsdir/$sampleid\.not-sapiens.paired_2.fastq.gz $resultsdir/$sampleid\.not-sapiens.unpaired_2.fastq.gz ILLUMINACLIP:/nfs/users/nfs_m/mb29/references/AdaptorTrimming/illumina-adaptors.2.fasta:2:30:10 LEADING:3 MINLEN:40 ;
        # Count remaining paired reads, if more than minimum amount downsample
	cp $resultsdir/$sampleid\.not-sapiens.paired_1.fastq.gz $resultsdir/$sampleid\.not-sapiens_1.fastq.gz
	cp $resultsdir/$sampleid\.not-sapiens.paired_2.fastq.gz $resultsdir/$sampleid\.not-sapiens_2.fastq.gz

        # Clean up files
        rm $resultsdir/$sampleid\.not-sapiens.binned_*.fastq.gz
        rm $resultsdir/$sampleid\.not-sapiens.paired_*.fastq.gz
        rm $resultsdir/$sampleid\.list.tmp
        rm $resultsdir/$sampleid\.not-sapiens.R*.list
        rm $resultsdir/$sampleid\.not-sapiens.unpaired_*.fastq.gz ;
        gzip -9 $resultsdir/$sampleid\.not-sapiens.list
        gzip -9 $resultsdir/$sampleid\.kraken.labels
  ")
else

if [[ $trim = TRUE ]]
then
  dokbin=$(bsub -w "ended($dokraken2_jobID)" -J kbin_$sampleid -o $resultsdir/kbin_$sampleid.%J.o -e $resultsdir/kbin_$sampleid.%J.e -q normal -n 2 -M$kbinmem -R "span[hosts=1] select[mem>$kbinmem] rusage[mem=$kbinmem]" "
	# Extract reads of interest using kraken outputs and seqtk
	grep $binning_term $resultsdir/$sampleid\.kraken.labels | awk '{print \$2}' > $resultsdir/$sampleid\.$binning_term\.list ;
	cat $resultsdir/$sampleid\.$binning_term\.list | perl -pe 's/\n/\/1\n/g' > $resultsdir/$sampleid\.$binning_term\.R1.list ;
	cat $resultsdir/$sampleid\.$binning_term\.list | perl -pe 's/\n/\/2\n/g' > $resultsdir/$sampleid\.$binning_term\.R2.list ;
	seqtk subseq $read1 $resultsdir/$sampleid\.$binning_term\.R1.list | gzip -1 > $resultsdir/$sampleid\.$binning_term\.binned_1.fastq.gz ;
	seqtk subseq $read2 $resultsdir/$sampleid\.$binning_term\.R2.list | gzip -1 > $resultsdir/$sampleid\.$binning_term\.binned_2.fastq.gz ;
	
	# Trim binned reads
	trimmomatic PE -phred33 -threads 2 $resultsdir/$sampleid\.$binning_term\.binned_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.binned_2.fastq.gz $resultsdir/$sampleid\.$binning_term\.paired_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.unpaired_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.paired_2.fastq.gz $resultsdir/$sampleid\.$binning_term\.unpaired_2.fastq.gz ILLUMINACLIP:/nfs/users/nfs_m/mb29/references/AdaptorTrimming/illumina-adaptors.2.fasta:2:30:10 LEADING:3 MINLEN:40 ;
	# Count remaining paired reads, if more than minimum amount downsample
	/nfs/users/nfs_m/mb29/scripts/downsample_fastq.sh $resultsdir/$sampleid\.$binning_term $resultsdir/$sampleid\.$binning_term\.paired_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.paired_2.fastq.gz $readcount ;

	# Clean up files	
	rm $resultsdir/$sampleid\.$binning_term\.binned_*.fastq.gz
	rm $resultsdir/$sampleid\.$binning_term\.paired_*.fastq.gz
	rm $resultsdir/$sampleid\.list.tmp
	rm $resultsdir/$sampleid\.$binning_term\.R*.list
	rm $resultsdir/$sampleid\.$binning_term\.unpaired_*.fastq.gz ;
	gzip -9 $resultsdir/$sampleid\.$binning_term\.list
	gzip -9 $resultsdir/$sampleid\.kraken.labels
  ")
else
	echo "Reads will not be trimmed"
   dokbin=$(bsub -w "ended($dokraken2_jobID)" -J kbin_$sampleid -o $resultsdir/kbin_$sampleid.%J.o -e $resultsdir/kbin_$sampleid.%J.e -q normal -n 2 -M$kbinmem -R "span[hosts=1] select[mem>$kbinmem] rusage[mem=$kbinmem]" "
        # Extract reads of interest using kraken outputs and seqtk
        grep $binning_term $resultsdir/$sampleid\.kraken.labels | awk '{print \$2}' > $resultsdir/$sampleid\.$binning_term\.list ;
        cat $resultsdir/$sampleid\.$binning_term\.list | perl -pe 's/\n/\/1\n/g' > $resultsdir/$sampleid\.$binning_term\.R1.list ;
        cat $resultsdir/$sampleid\.$binning_term\.list | perl -pe 's/\n/\/2\n/g' > $resultsdir/$sampleid\.$binning_term\.R2.list ;
        seqtk subseq $read1 $resultsdir/$sampleid\.$binning_term\.R1.list | gzip -1 > $resultsdir/$sampleid\.$binning_term\.binned_1.fastq.gz ;
        seqtk subseq $read2 $resultsdir/$sampleid\.$binning_term\.R2.list | gzip -1 > $resultsdir/$sampleid\.$binning_term\.binned_2.fastq.gz ;
        
        # Count remaining paired reads, if more than minimum amount downsample
        /nfs/users/nfs_m/mb29/scripts/downsample_fastq.sh $resultsdir/$sampleid\.$binning_term $resultsdir/$sampleid\.$binning_term\.binned_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.binned_2.fastq.gz $readcount ;

        # Clean up files        
        rm $resultsdir/$sampleid\.list.tmp
        rm $resultsdir/$sampleid\.$binning_term\.R*.list
        rm $resultsdir/$sampleid\.$binning_term\.unpaired_*.fastq.gz ;
        gzip -9 $resultsdir/$sampleid\.$binning_term\.list
	gzip -9 $resultsdir/$sampleid\.kraken.labels
  ")
fi
fi

dokbin_jobID=$(echo $dokbin |perl -pe 's/\>.+$//g'|perl -pe 's/^.+\<//g')
echo "Kraken read binning and downsampling submitted with $kbinmem Mb as jobID $dokbin_jobID - will start after kraken run finishes"
# Print optional flags to screen
echo "Binning reads classified as $binning_term and downsampling to $readcount reads or less"




