#!/bin/bash


# Usage: count_Ns_in_alignment.sh <myalignment.fasta>

while getopts ":hs" opt; do
  case "${opt}" in
    h) echo ""
       echo "Count number of non-ATGC bases (i.e. N/-) in a multiple sequence alignment"
       echo "Usage: $0 [ options ] <myalignment.fasta>"
#       echo ""
       echo "    -s  -  sort output by proportion of N's"
       echo ""
       exit 0
      ;;
    s) sortoutput='dosort' ;;
    \?)echo "Invalid Option: -$OPTARG" 1>&2
       exit 1
      ;;
  esac
done
shift $((OPTIND -1))


if [ $# -ne 1 ]
then
 echo
 echo "Count number of non-ATGC bases (i.e. N/-) in a multiple sequence alignment"
 echo "Usage: $0 [ options ] <myalignment.fasta>"
 echo "    -s  -  sort output by proportion of N's"
 echo
 exit 1
fi


module load seqtk/1.3--ha92aebf_0


#echo "Sample	N_count	Total_Sites	N_proportion"

#data_out=$(seqtk comp $1 | awk '{x=$3+$4+$5+$6;y=$2;print $1,y-x,y,(y-x)/y}')



if [[ $sortoutput = "dosort" ]]
then
#  data_out=$(sort -k 4 echo $data_out)
  data_out=$(seqtk comp $1 | awk '{x=$3+$4+$5+$6;y=$2;print $1,y-x,y,(y-x)/y}'| sort -k 4 -r )
else
  data_out=$(seqtk comp $1 | awk '{x=$3+$4+$5+$6;y=$2;print $1,y-x,y,(y-x)/y}')
fi

#sort -k 4


echo "Sample    N_count Total_Sites     N_proportion"
echo "$data_out"
