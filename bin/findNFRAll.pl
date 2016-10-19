#!/usr/bin/env perl

=copyright_info
findNFRAll.pl: determine Nucleosome Free Regions (NFR) using ChIP-seq data for histone marks
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
=cut

use strict;
use warnings;
use Getopt::Long;
use Statistics::Basic qw(:all);
use perlModule;

###############################################################################
## parse input options
use vars qw($readFile $bamFile $outDir $sizeFactor $option $minClusterHeight $minBlockHeight $distance $scale $blockHeight $noMergeOverlapBlocks $minNFRLength $maxNFRLength $extend $genome $fileSuffix $help);
$minClusterHeight=20;
$minBlockHeight=20;
$distance=70;
$scale=0.6;
$blockHeight="abs";
$minNFRLength=20;
$maxNFRLength=1000;
$extend=0;
$genome="mm9";
$fileSuffix="";

GetOptions ("s=s"  => \$readFile,
            "b=s"  => \$bamFile,
            "o=s"  => \$outDir,
            "z=s"  => \$sizeFactor,
            "p=s"  => \$option,
            "c=s"  => \$minClusterHeight,
            "k=s"  => \$minBlockHeight,
            "x=s"  => \$distance,
            "l=s"  => \$scale,
            "g=s"  => \$blockHeight,
            "m"    => \$noMergeOverlapBlocks,
            "n=s"  => \$minNFRLength,
            "v=s"  => \$maxNFRLength,
            "e=s"  => \$extend,
            "y=s"  => \$genome,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$readFile || !$bamFile || !$outDir || !$sizeFactor);

###############################################################################
sub usage {
	print STDERR "\nProgram: findNFRAll.pl (determine Nucleosome Free Regions (NFR) using ChIP-seq data for histone marks)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: findNFRAll.pl -s <file> -b <file> -o <dir> -z <float> [OPTIONS]\n";
	print STDERR " -s <file>         [file containing reads corresponding to regions of interest]\n";
    print STDERR "                   [can be reads in bed format for whole genome]\n";
	print STDERR " -b <file>         [histone ChIP-seq file in BAM format]\n";
	print STDERR " -o <dir>          [directory where output files will be kept]\n";
	print STDERR " -z <float>        [size factor to normalize read expression]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -p <string>       [computation option (default: all)]\n";
    print STDERR "                   a: define blocks and block groups in histone enriched (summit) region\n";
    print STDERR "                   b: define nuclesome free regions\n";
    print STDERR "                   c: determine significant nucleosome free regions\n";
	print STDERR " -c <int>          [mininum number of read in the block group (default: 20)]\n";
	print STDERR " -k <int>          [mininum number of read in the block (default: 20)]\n";
	print STDERR " -x <int>          [maximum distance between the blocks (default: 70)]\n";
	print STDERR " -l <float>        [scale to define blocks (default: 0.6)]\n";
	print STDERR " -g <string>       [relative block height (abs or rel) (default: abs)]\n";
	print STDERR " -m                [do not merge overlapping blocks]\n";
    print STDERR " -n <int>          [minimum length of nucleosome free region (default: 20)]\n";
    print STDERR " -v <int>          [maximum length of nucleosome free region (default: 1000)]\n";
    print STDERR " -e <int>          [extend 3' end of reads by input number of bases (default: 0)]\n";
    print STDERR " -y <string>       [genome (default: mm9)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

my $ID=$bamFile;
$ID=~s/^.*\///g;
$ID=~s/\.gz$//g;
my $start=(); my $end=(); my $coor=(); my @data=();

## create output directory, if does not exist
if ( ! -d $outDir) {
    system("mkdir $outDir");
}

## populate genome file based on input genome
my $GENOME_FILE=`initialize_genome -i $ENV{PAREPATH}/data/annotations/GENOME_FILE -g $genome`;
$GENOME_FILE="$ENV{PAREPATH}/data/annotations/$GENOME_FILE";
if(! -f $GENOME_FILE) {
    print "\ncomputation for $genome is not feasible yet\n";
    print "please add the chromosome size file for $genome at $ENV{PAREPATH}/data/annotations\n";
    print "also update the $ENV{PAREPATH}/data/annotations/GENOME_FILE\n";
    usage();
}

## Step-1: define blocks and block groups
if(!defined($option) || $option=~/[aA]+/) {
    #print("blockbuster.x -minClusterHeight $minClusterHeight -minBlockHeight $minBlockHeight -distance $distance -scale $scale -blockHeight $blockHeight -print 2 $readFile | validateBlockbuster.pl -d $maxNFRLength | perl -ane 'if(\$_=~/^>/) { \$id=\"\$F[1]:\$F[2]-\$F[3]\"; } else { print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$id\\t\$F[4]\\t\$F[5]\\n\"; }' > $outDir/$ID.tmp$fileSuffix\n"); exit;
    system("blockbuster.x -minClusterHeight $minClusterHeight -minBlockHeight $minBlockHeight -distance $distance -scale $scale -blockHeight $blockHeight -print 2 $readFile | validateBlockbuster.pl -d $maxNFRLength | perl -ane 'if(\$_=~/^>/) { \$id=\"\$F[1]:\$F[2]-\$F[3]\"; } else { print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$id\\t\$F[4]\\t\$F[5]\\n\"; }' > $outDir/$ID.tmp$fileSuffix");

    #system("sortBed -i $outDir/$ID.tmp$fileSuffix | bedtools merge -nms -scores sum -i - | perl -ane '\$F[3]=~s/\\,.*//g; \$F[4]=sprintf(\"%0.2f\", \$F[4]); print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[4]\\t+\\n\";' > $outDir/$ID.bg$fileSuffix");
    system("sortBed -i $outDir/$ID.tmp$fileSuffix | bedtools merge -c 4,5 -o distinct,sum -i - -d $minNFRLength | perl -ane '\$F[3]=~s/\\,.*//g; \$F[4]=sprintf(\"%0.2f\", \$F[4]); print \"\$F[0]\\t\$F[1]\\t\$F[2]\\t\$F[3]\\t\$F[4]\\t+\\n\";' > $outDir/$ID.bg$fileSuffix");

    #system("rm $outDir/$ID.tmp$fileSuffix");
}

## Step-2: define nuclesome free regions
if(!defined($option) || $option=~/[bB]+/) {
    if(-e "$outDir/$ID.bg$fileSuffix") {
        @data=openFile("$outDir/$ID.bg$fileSuffix");
        #@data=openFile("$outDir/test");
    }
    else {
        print STDERR "Cannot find $ID.bg$fileSuffix file. Please run the program with option a first\n";
        usage();
    }

    ## create output file for writing
    open(OUTFILE, ">$outDir/$ID.nfr$fileSuffix") || die $!;

    my %bInfo=(); my @t=(); my %NFR=();
    for(my $i=0; $i<scalar(@data)-2; $i++) {
        my @F1=split(/\s+/, $data[$i]);
        my @F2=split(/\s+/, $data[$i+1]);

        if($F1[3]=~/^$F2[3]$/ || ($F2[1]-$F1[2])<$maxNFRLength) {
            #print "DEBUG: $data[$i]\nDEBUG: $data[$i+1]\n";
            ## collect first block group information
            $bInfo{'first'}{'chr'}=$F1[0];
            $bInfo{'first'}{'start'}=$F1[1];
            $bInfo{'first'}{'end'}=$F1[2];
            $bInfo{'first'}{'strand'}=$F1[5];
            $bInfo{'first'}{'length'}=($bInfo{'first'}{'end'}-$bInfo{'first'}{'start'})+1;
            #@t=split(/\,/,$F1[4]);
            #$bInfo{'first'}{'expr'}=0;
            #foreach(@t) {
            #    $bInfo{'first'}{'expr'}+=$_;
            #}
            $coor="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
            $bInfo{'first'}{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
            chomp($bInfo{'first'}{'expr'});

            ## collect second block group information
            $bInfo{'second'}{'chr'}=$F2[0];
            $bInfo{'second'}{'start'}=$F2[1];
            $bInfo{'second'}{'end'}=$F2[2];
            $bInfo{'second'}{'strand'}=$F2[5];
            $bInfo{'second'}{'length'}=($bInfo{'second'}{'end'}-$bInfo{'second'}{'start'})+1;
            #@t=split(/\,/,$F2[4]);
            #$bInfo{'second'}{'expr'}=0;
            #foreach(@t) {
            #    $bInfo{'second'}{'expr'}+=$_;
            #}
            $coor="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
            $bInfo{'second'}{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
            chomp($bInfo{'second'}{'expr'});

            ## define NFR based on first and second block group information
            if($bInfo{'first'}{'strand'}=~/\+/ && $bInfo{'second'}{'strand'}=~/\+/) {
                $NFR{'chr'}=$bInfo{'first'}{'chr'};
                $NFR{'start'}=$bInfo{'first'}{'end'}+1;
                $NFR{'end'}=$bInfo{'second'}{'start'}-1;
                $NFR{'strand'}="+";
                $NFR{'length'}=($NFR{'end'}-$NFR{'start'})+1;
                $coor="$NFR{'chr'}:$NFR{'start'}-$NFR{'end'}";
                $NFR{'startBlock'}="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
                $NFR{'startBlockExpr'}=$bInfo{'first'}{'expr'};
                $NFR{'startBlockLen'}=$bInfo{'first'}{'length'};
                $NFR{'endBlock'}="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
                $NFR{'endBlockExpr'}=$bInfo{'second'}{'expr'};
                $NFR{'endBlockLen'}=$bInfo{'second'}{'length'};
                if($F1[3]=~/^$F2[3]$/) {
                    $NFR{'id'}=$F1[3];
                }
                else {
                    my @t1=split(/[\:\-]+/,$F1[3]);
                    my @t2=split(/[\:\-]+/,$F2[3]);
                    $NFR{'id'}="$t1[0]:$t1[1]-$t2[2]";
                }

                #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
                #print "DEBUG: $bInfo{'first'}{'expr'}\t$bInfo{'second'}{'expr'}\t$NFR{'expr'}\t$NFR{'length'}\n";

                $NFR{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
                chomp($NFR{'expr'});
            }
            elsif($bInfo{'first'}{'strand'}=~/\-/ && $bInfo{'second'}{'strand'}=~/\-/) {
                $NFR{'chr'}=$bInfo{'first'}{'chr'};
                $NFR{'start'}=$bInfo{'second'}{'end'}+1;
                $NFR{'end'}=$bInfo{'first'}{'start'}-1;
                $NFR{'strand'}="-";
                $NFR{'length'}=($NFR{'end'}-$NFR{'start'})+1;
                $coor="$NFR{'chr'}:$NFR{'start'}-$NFR{'end'}";
                $NFR{'startBlock'}="$bInfo{'second'}{'chr'}:$bInfo{'second'}{'start'}-$bInfo{'second'}{'end'}";
                $NFR{'startBlockExpr'}=$bInfo{'second'}{'expr'};
                $NFR{'startBlockLen'}=$bInfo{'second'}{'length'};
                $NFR{'endBlock'}="$bInfo{'first'}{'chr'}:$bInfo{'first'}{'start'}-$bInfo{'first'}{'end'}";
                $NFR{'endBlockExpr'}=$bInfo{'first'}{'expr'};
                $NFR{'endBlockLen'}=$bInfo{'first'}{'length'};
                $NFR{'id'}=$F1[3];

                #print("DEBUG: samtools view -b $bamFile $coor | bedtools bamtobed -i - | perl -ane '\$F[4]=sprintf(\"%0.2f\", 1/$sizeFactor); \$expr+=\$F[4]; END { print \$expr; }'\n");
                #print "DEBUG: $bInfo{'first'}{'expr'}\t$bInfo{'second'}{'expr'}\t$NFR{'expr'}\t$NFR{'length'}\n";

                $NFR{'expr'}=`coor2expr -i $coor -j $bamFile -k $sizeFactor -d -e $extend -g $genome`;
                chomp($NFR{'expr'});
            }
            else {
                print STDERR "WARNING: the first and second block do not have same strand despite having same ID\n";
                print STDERR "--> first block: $data[$i]\n";
                print STDERR "--> second block: $data[$i+1]\n";
            }

            if($NFR{'length'} >= $minNFRLength) {
                ## old scoring and output scheme
                #$NFR{'stddev'}=stddev($bInfo{'first'}{'expr'}, $bInfo{'second'}{'expr'});
                #$NFR{'stddev'}=~s/\,//g;
                #print("DEBUG: (($bInfo{'first'}{'expr'}/$bInfo{'first'}{'length'})+($bInfo{'second'}{'expr'}/$bInfo{'second'}{'length'}))/($NFR{'expr'}/$NFR{'length'}))\n");
                #$NFR{'score'}=((($NFR{'startBlockExpr'}/$NFR{'startBlockLen'})+($NFR{'endBlockExpr'}/$NFR{'endBlockLen'}))/($NFR{'expr'}/$NFR{'length'}));
                #$NFR{'score'}=sprintf("%0.4f", ($NFR{'score'}));
                #$NFR{'score'}=sprintf("%0.4f", log($NFR{'score'}));

                #print OUTFILE "$NFR{'chr'}\t$NFR{'start'}\t$NFR{'end'}\t$NFR{'id'}\t$NFR{'score'}\t$NFR{'strand'}\t$NFR{'startBlock'}\t$NFR{'startBlockExpr'}\t$NFR{'endBlock'}\t$NFR{'endBlockExpr'}\t$NFR{'expr'}\t$NFR{'length'}\n";

                ## new scoring and output scheme
                $NFR{'score'}=((($NFR{'startBlockExpr'}+$NFR{'endBlockExpr'})/($NFR{'startBlockLen'}+$NFR{'endBlockLen'}))-($NFR{'expr'}/$NFR{'length'}));
                $NFR{'score'}=sprintf("%0.4f", ($NFR{'score'}));

                print OUTFILE "$NFR{'chr'}\t$NFR{'start'}\t$NFR{'end'}\t$NFR{'id'}\t$NFR{'score'}\t$NFR{'strand'}\t$NFR{'startBlock'}\t$NFR{'endBlock'}\t$NFR{'startBlockExpr'}\t$NFR{'startBlockLen'}\t$NFR{'endBlockExpr'}\t$NFR{'endBlockLen'}\t$NFR{'expr'}\t$NFR{'length'}\n";
            }
        }
    }
    close OUTFILE;
}

## Step-3: define non-overlapping and significant nuclesome free regions
if(!defined($option) || $option=~/[cC]+/) {
    if(-e "$outDir/$ID.nfr$fileSuffix") {
        @data=`zless $outDir/$ID.nfr$fileSuffix | sort -k 1,1 -k 2n,2 -k 3n,3`;
    }
    else {
        print STDERR "Cannot find $ID.nfr$fileSuffix file. Please run the program with option b first\n";
        usage()
    }

    ## create output file for writing non-overlapping NFR
    open(OUTFILE, ">$outDir/$ID.nfr.uniq$fileSuffix") || die $!;

    my @F1=(); my @F2=();
    for(my $i=0; $i<(scalar(@data)-1); $i++) {
        if(scalar(@F1)==0) {
            chomp($data[$i]);
            @F1=split(/\s+/, $data[$i]);
        }
        @F2=split(/\s+/, $data[$i+1]);

        ## check overlap between current and next NFR
        my $overlap=checkOverlap($F1[1], $F1[2], $F2[1], $F2[2], 4, $F1[0], $F2[0]);

        ## choose the highest scoring one, if overlapped
        if($overlap) {
            if($F1[4]<$F2[4]) {
                @F1=@F2; 
            }
        }
        else {
            ## print the non-overlapping and highest scoring NFR region
            my $NFR=();
            foreach(@F1) { $NFR.="$_\t"; } $NFR=~s/\t$//g;
            print OUTFILE "$NFR\n";

            @F1=();
        }
    }

    close OUTFILE;

    print "$outDir/$ID.nfr.uniq$fileSuffix\n";
}

exit(0);
