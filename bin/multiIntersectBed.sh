#!/bin/bash
#PBS -l nodes=1:ppn=4

<<"copyright_info"
multiIntersectBed.sh: intersect genomic coordinates from multiple BED files
Copyright (C) 2015  Sachin Pundhir (pundhir@binf.ku.dk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
copyright_info

#### usage ####
usage() {
    echo
	echo Program: "multiIntersectBed.sh (intersect genomic coordinates from multiple BED files)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: multiIntersectBed.sh -i <files> [OPTIONS]"
    echo "Note: Different from multiIntersectBed (bedtools) since, this only selects the single most observed coordinate among the consecutive overlapping coordinates"
	echo "Options:"
    echo " -i <files>  [input BED files seperated by a comma]"
    echo "[OPTIONS]"
    echo " -j <string> [names to describe each input files seperated by a comma]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:j:h ARG; do
	case "$ARG" in
		i) BEDFILE=$OPTARG;;
        j) NAME=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$BEDFILE" -o "$HELP" ]; then
	usage
fi

###################
#helperfunction
function wait_for_jobs_to_finish {
    for job in `jobs -p`
    do
        echo $job
        wait $job
    done
    echo $1
}
###############

## parse input bam files in an array
IFS=","
BEDFILES=($BEDFILE)
BEDFILES_COUNT=${#BEDFILES[@]}
IFS=" "

## initialize name parameter, if provided
if [ ! -z "$NAME" ]; then
    IFS=","
    NAMES=($NAME)
    NAMES_COUNT=${#NAMES[@]}
    IFS=" "
else
    NAMES_COUNT=0
fi

if [ "$BEDFILES_COUNT" -lt 2 -o ! -z "$NAME" -a "$BEDFILES_COUNT" -ne "$NAMES_COUNT" ]; then
    echo -n "minimum two input bed files are required as input. Also provide name for each input bed file";
    usage
fi

COMMAND_BED=""
COMMAND_NAME=""
for (( i=0; i<$BEDFILES_COUNT; i++ )); do
    TMP_NAME[i]=$RANDOM
    bedtools sort -i ${BEDFILES[$i]} > /tmp/${TMP_NAME[$i]}.txt
    COMMAND_BED="$COMMAND_BED /tmp/${TMP_NAME[$i]}.txt"
    if [ ! -z "$NAME" ]; then
        COMMAND_NAME="$COMMAND_NAME ${NAMES[$i]}"
    fi
done
wait

if [ ! -z "$NAME" ]; then
    bedtools multiinter -i $COMMAND_BED -names $COMMAND_NAME | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_counter!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
else
    bedtools multiinter -i $COMMAND_BED | perl -ane 'if(defined($line)) { if($F[1]==$last_coor) { if($F[3]>$last_counter) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } else { $last_coor=$F[2]; $last_counter=$F[3]; } } elsif($last_counter!=$F[1]) { print "$line"; $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } } elsif(!defined($line)) { $line=$_; $last_coor=$F[2]; $last_counter=$F[3]; } END { print "$line"; }'
fi

## remove temporary files
for (( i=0; i<$BEDFILES_COUNT; i++ )); do
    rm /tmp/${TMP_NAME[$i]}.txt
done
