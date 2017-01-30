#!/usr/bin/env bash
DIR=$1
COV=$2
SIM=$3
THR=$4
EXT=$5
OUT=$6


BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

array=()
x=0

if [ $# != 6 ]; then
	echo "***ERROR*** Use: $0 metagenomes_directory coverage similarity threads file_extension outpath"
	exit -1
fi

for elem in $(ls -d $DIR/*.$EXT | awk -F "/" '{print $NF}' | awk -F ".$EXT" '{print $1}')
do
	array[$x]=$elem
	x=`expr $x + 1`
	#echo "X: $elem"
done

for ((i=0 ; i < ${#array[@]} ; i++))
do
	for ((j=i ; j < ${#array[@]} ; j++))
	do
		if [ $i != $j ]; then
			seqX=${array[$i]}
			seqY=${array[$j]}
			#echo "----------${seqX}-${seqY}-----------"
			if [[ ! -f $6/${seqX}-${seqY}.align ]];	then #if the file does not exist

				
				${BINDIR}/IMSAME -query $DIR/${seqX}.$EXT -db $DIR/${seqY}.$EXT -n_threads $THR -coverage $COV -identity $SIM -out $6/${seqX}-${seqY}.align
			fi
			
			

			
			
			if [[ ! -f $6/${seqX}-${seqY}.r.align ]];	then #if the reversed alignment does not exist
				# Compute reverse complement
				${BINDIR}/revComp $DIR/${seqY}.$EXT $DIR/${seqY}.r.${EXT}
				${BINDIR}/IMSAME -query $DIR/${seqX}.$EXT -db $DIR/${seqY}.r.$EXT -n_threads $THR -coverage $COV -identity $SIM -out $6/${seqX}-${seqY}.r.align
			
			fi
			
			
			
		fi
		rm $DIR/${seqY}.r.${EXT}
	done

done
