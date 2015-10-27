#!/bin/bash
#PBS -l nodes=1:ppn=4

<<"copyright_info"
indexBed.sh: create indexing of bed file
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

OUTDIR="."

#### usage ####
usage() {
	echo Program: "indexBed.sh (create index of bed file)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: indexBed.sh -i <file> [OPTIONS]"
	echo " -i <file>     [input map or anno file]"
	echo "[OPTIONS]:"
	echo " -o <dir>      [output directory with index files]"
    echo " -x <string>   [prefix to add to input file name]"
    echo " -s            [seperate index file for each strand]"
    echo " -z            [compress output file]"
    echo " -h            [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:o:x:szh ARG; do
	case "$ARG" in
		i) INFILE=$OPTARG;;
		o) OUTDIR=$OPTARG;;
        x) PREFIX=$OPTARG;;
        s) STRAND=1;;
        z) COMPRESS=1;;
        h) HELP=1;;
	esac
done

if [ ! -f "$INFILE" -o "$HELP" ]; then
	usage
fi

## create output directory
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR;
fi

ID=`echo $INFILE | perl -ane '\$_=~s/\\.[gz|bed]+//gi; \$_=~s/.+\///g; print \$_;'`
if [ -z "$PREFIX" ]; then
    PREFIX=$ID
fi

if [ ! -z "$STRAND" ]; then
    ## split BED file
    if [ `command -v bedSplitOnChrom | wc -l` -eq 1 ]; then
        bedSplitOnChrom $INFILE $OUTDIR -strand
        rename -.bed .N $OUTDIR/chr*
        rename +.bed .P $OUTDIR/chr*
        rename chr $PREFIX.chr $OUTDIR/chr*
    else
        CHR=(`zless $INFILE | cut -f 1 | grep -v "_" | sort | uniq`)

        for chr in "${CHR[@]}";
        do
            zgrep -E "^$chr\s+" $INFILE | grep -w "-" > $OUTDIR/$PREFIX"."$chr".N";
            zgrep -E "^$chr\s+" $INFILE | grep -w "+" > $OUTDIR/$PREFIX"."$chr".P";
        done
    fi
else
    ## split BED file
    if [ `command -v bedSplitOnChrom | wc -l` -eq 1 ]; then
        bedSplitOnChrom $INFILE $OUTDIR -strand
        rename -.bed .N $OUTDIR/chr*
        rename +.bed .P $OUTDIR/chr*
        rename chr $PREFIX.chr $OUTDIR/chr*
    else
        CHR=(`zless $INFILE | cut -f 1 | grep -vE "[\_|\||MT|poly]+" | sort | uniq`)

        for chr in "${CHR[@]}";
        do  
            zgrep -E "^$chr\s+" $INFILE > $OUTDIR/$PREFIX$chr;
        done
    fi
fi

## compress output files
if [ ! -z "$COMPRESS" ]; then
    for file in `ls $OUTDIR/$PREFIX*`; do gzip $file; done
fi
