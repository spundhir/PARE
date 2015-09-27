#!/usr/bin/perl -w

=copyright_info
commonNFR.pl: determine common Nucleosome Free Regions (NFR) between two replicates
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
use vars qw($nfrFileRep1 $nfrFileRep2 $bamFileRep1 $bamFileRep2 $sizeFactorRep1 $sizeFactorRep2 $outDir $minNFRLength $nfrThreshold $extendRep1 $extendRep2 $fileSuffix $genome $help);

$minNFRLength=20;
$extendRep1=0;
$extendRep2=0;
$genome="mm9";
$fileSuffix="";

GetOptions ("i=s"  => \$nfrFileRep1,
            "j=s"  => \$nfrFileRep2,
            "k=s"  => \$bamFileRep1,
            "l=s"  => \$bamFileRep2,
            "m=s"  => \$sizeFactorRep1,
            "n=s"  => \$sizeFactorRep2,
            "o=s"  => \$outDir,
            "g=s"  => \$minNFRLength,
            "t=s"  => \$nfrThreshold,
            "c=s"  => \$extendRep1,
            "d=s"  => \$extendRep2,
            "y=s"  => \$genome,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$nfrFileRep1 || !$nfrFileRep2 || !$bamFileRep1 || !$bamFileRep2 || !$sizeFactorRep1 || !$sizeFactorRep2 || !$outDir);

###############################################################################
sub usage {
	print STDERR "\nProgram: commonNFR.pl (determine common Nucleosome Free Regions (NFR) between two replicates)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: commonNFR.pl -i <file> -j <file> -k <file> -l <file> -m <float> -n <float> -o <dir> [OPTIONS]\n";
	print STDERR " -i <file>         [input file having NFRs from first replicate]\n";
	print STDERR " -j <file>         [input file having NFRs from second replicate]\n";
	print STDERR " -k <file>         [input BAM file for first replicate]\n";
	print STDERR " -l <file>         [input BAM file for second replicate]\n";
    print STDERR " -m <float>        [size factor to normalize expression of reads from first replicate]\n";
    print STDERR " -n <float>        [size factor to normalize expression of reads from second replicate]\n";
	print STDERR " -o <dir>          [directory where output files will be kept]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -g <int>          [minimum length of nucleosome free region (default: 20)]\n";
    print STDERR " -t <int>          [minimum score to consider a NFR as significant]\n";
    print STDERR " -c <int>          [extend 3' end of reads by input number of bases from first replicate (default: 0)]\n";
    print STDERR " -d <int>          [extend 3' end of reads by input number of bases from second replicate (default: 0)]\n";
    print STDERR " -y <string>       [genome (default: mm9)]\n";
    print STDERR " -f <string>       [a string added at the end of output files. useful when running in parallel]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

## populate genome file based on input genome
my $GENOME_FILE=`initialize_genome -i $ENV{PAREPATH}/data/annotations/GENOME_FILE -g $genome`;
$GENOME_FILE="$ENV{PAREPATH}/data/annotations/$GENOME_FILE";
if(! -f $GENOME_FILE) {
    print "\ncomputation for $genome is not feasible yet\n";
    print "please add the chromosome size file for $genome at $ENV{PAREPATH}/data/annotations\n";
    print "also update the $ENV{PAREPATH}/data/annotations/GENOME_FILE\n";
    usage();
}

my $ID=$nfrFileRep1;
$ID=~s/^.*\///g;
$ID=~s/Rep.*$//g;
$ID=~s/\_$//g;

## create output directory, if does not exist
if ( ! -d $outDir) {
    system("mkdir $outDir");
}

## compute overlapping NFRs
#system("intersectBed -a $nfrFileRep1 -b $nfrFileRep2 > $outDir/$ID.All.nfr");
#system("intersectBed -a $nfrFileRep1 -b $nfrFileRep2 -wo | sort -k 1,1 -k 2n,2 -k 3n,3 -k 17rn,17 | perl -ane 'if((\$F[4]/(\$F[10]+0.0001))>(\$F[16]/(\$F[22]+0.0001))) { print \"\$F[0]\t\$F[1]\t\$F[2]\t\$F[3]\t\$F[4]\t\$F[5]\t\$F[6]\t\$F[7]\t\$F[8]\t\$F[9]\t\$F[10]\t\$F[11]\n\"; } else { print \"\$F[12]\t\$F[13]\t\$F[14]\t\$F[15]\t\$F[16]\t\$F[17]\t\$F[18]\t\$F[19]\t\$F[20]\t\$F[21]\t\$F[22]\t\$F[23]\n\"; }' > $outDir/$ID.All.nfr");

my @data=`intersectBed -a $nfrFileRep1 -b $nfrFileRep2 -wo | sort -k 1,1 -k 2n,2 -k 3n,3`;

## output file for writing overlapping NFR regions
open(OUTFILE, ">$outDir/$ID.All.nfr$fileSuffix") || die $!;
open(SIGFILE, ">$outDir/$ID.All.nfr.sig$fileSuffix") || die $! if(defined($nfrThreshold));

my %NFR=(); my %startBlock=(); my %endBlock=();
foreach my $l(@data) {
    my @F=split(/\s+/,$l);
    @{$startBlock{'rep1Coor'}}=split(/[\:\-]+/,$F[6]);
    @{$endBlock{'rep1Coor'}}=split(/[\:\-]+/,$F[7]);
    @{$startBlock{'rep2Coor'}}=split(/[\:\-]+/,$F[20]);
    @{$endBlock{'rep2Coor'}}=split(/[\:\-]+/,$F[21]);

    $NFR{'overlapCoor'}=returnOverlapCoor($F[1], $F[2], $F[15], $F[16], $F[0], $F[14], "intersect", 0);
    $startBlock{'overlapCoor'}=returnOverlapCoor($startBlock{'rep1Coor'}[1], $startBlock{'rep1Coor'}[2], $startBlock{'rep2Coor'}[1], $startBlock{'rep2Coor'}[2], $startBlock{'rep1Coor'}[0], $startBlock{'rep2Coor'}[0], "union", 1);
    $endBlock{'overlapCoor'}=returnOverlapCoor($endBlock{'rep1Coor'}[1], $endBlock{'rep1Coor'}[2], $endBlock{'rep2Coor'}[1], $endBlock{'rep2Coor'}[2], $endBlock{'rep1Coor'}[0], $endBlock{'rep2Coor'}[0], "union", 2);

    ## check if overlap has been done correctly
    @{$NFR{'overlapCoorSplit'}}=split(/[\:\-]+/,$NFR{'overlapCoor'});
    @{$startBlock{'overlapCoorSplit'}}=split(/[\:\-]+/,$startBlock{'overlapCoor'});
    @{$endBlock{'overlapCoorSplit'}}=split(/[\:\-]+/,$endBlock{'overlapCoor'});

    if($NFR{'overlapCoorSplit'}[1]-$startBlock{'overlapCoorSplit'}[2]>1 || $endBlock{'overlapCoorSplit'}[1]-$NFR{'overlapCoorSplit'}[2]>1) {
        print STDERR "start or end block are not adjacent to NFR\n";
        print STDERR "--> $NFR{'overlapCoorSplit'}[1]\t$startBlock{'overlapCoorSplit'}[2]\t$endBlock{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\n";
        #print "$l\n";
        print STDERR "$F[0]:$F[1]-$F[2]\t$F[6]\t$F[7]\n";
        print STDERR "$F[14]:$F[15]-$F[16]\t$F[20]\t$F[21]\n";
        print STDERR "$NFR{'overlapCoor'}\t$startBlock{'overlapCoor'}\t$endBlock{'overlapCoor'}\n\n";
        exit(-1);
    }

    $NFR{'expr'}=`coor2expr -i $NFR{'overlapCoor'} -j $bamFileRep1,$bamFileRep2 -k $sizeFactorRep1,$sizeFactorRep2 -d -e $extendRep1,$extendRep2 -g $genome`;
    chomp($NFR{'expr'});
    $startBlock{'expr'}=`coor2expr -i $startBlock{'overlapCoor'} -j $bamFileRep1,$bamFileRep2 -k $sizeFactorRep1,$sizeFactorRep2 -d -e $extendRep1,$extendRep2 -g $genome`;
    chomp($startBlock{'expr'});
    $endBlock{'expr'}=`coor2expr -i $endBlock{'overlapCoor'} -j $bamFileRep1,$bamFileRep2 -k $sizeFactorRep1,$sizeFactorRep2 -d -e $extendRep1,$extendRep2 -g $genome`;
    chomp($endBlock{'expr'});

    #$NFR{'stddev'}=stddev($startBlock{'expr'}, $endBlock{'expr'});
    #$NFR{'stddev'}=~s/\,//g;

    $NFR{'length'}=($NFR{'overlapCoorSplit'}[2]-$NFR{'overlapCoorSplit'}[1])+1;
    $startBlock{'length'}=($startBlock{'overlapCoorSplit'}[2]-$startBlock{'overlapCoorSplit'}[1])+1;
    $endBlock{'length'}=($endBlock{'overlapCoorSplit'}[2]-$endBlock{'overlapCoorSplit'}[1])+1;

    ## old scoring scheme
    #$NFR{'score'}=((($startBlock{'expr'}/$startBlock{'length'})+($endBlock{'expr'}/$endBlock{'length'}))/($NFR{'expr'}/$NFR{'length'}));
    #$NFR{'score'}=sprintf("%0.4f", $NFR{'score'});
    #$NFR{'score'}=sprintf("%0.4f", log($NFR{'score'}));

    ## new scoring scheme
    $NFR{'score'}=((($startBlock{'expr'}+$endBlock{'expr'})/($startBlock{'length'}+$endBlock{'length'}))-($NFR{'expr'}/$NFR{'length'}));
    $NFR{'score'}=sprintf("%0.4f", $NFR{'score'});

    if($NFR{'length'} >= $minNFRLength) {
        print OUTFILE "$NFR{'overlapCoorSplit'}[0]\t$NFR{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\t$F[3]\t$NFR{'score'}\t$F[5]\t$startBlock{'overlapCoor'}\t$endBlock{'overlapCoor'}\t$startBlock{'expr'}\t$startBlock{'length'}\t$endBlock{'expr'}\t$endBlock{'length'}\t$NFR{'expr'}\t$NFR{'length'}\n";
        if(defined($nfrThreshold) && $NFR{'score'}>$nfrThreshold) {
            print SIGFILE "$NFR{'overlapCoorSplit'}[0]\t$NFR{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\t$F[3]\t$NFR{'score'}\t$F[5]\t$startBlock{'overlapCoor'}\t$endBlock{'overlapCoor'}\t$startBlock{'expr'}\t$startBlock{'length'}\t$endBlock{'expr'}\t$endBlock{'length'}\t$NFR{'expr'}\t$NFR{'length'}\n";
        }
    }
}

close OUTFILE;
close SIGFILE if(defined($nfrThreshold));

exit(0);
