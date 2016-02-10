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
use Tie::IxHash;

###############################################################################
## parse input options
use vars qw($nfrFiles $bamFiles $sizeFactors $outDir $minNFRLength $nfrThreshold $extends $fileSuffix $genome $help);

$minNFRLength=20;
$genome="mm9";
$fileSuffix="";

GetOptions ("i=s"  => \$nfrFiles,
            "k=s"  => \$bamFiles,
            "m=s"  => \$sizeFactors,
            "o=s"  => \$outDir,
            "g=s"  => \$minNFRLength,
            "t=s"  => \$nfrThreshold,
            "c=s"  => \$extends,
            "y=s"  => \$genome,
            "f=s"  => \$fileSuffix,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$nfrFiles || !$bamFiles || !$sizeFactors || !$outDir);

###############################################################################
sub usage {
	print STDERR "\nProgram: commonNFR.pl (determine common Nucleosome Free Regions (NFR) between two replicates)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: commonNFR.pl -i <file> -k <file> -m <float> -o <dir> [OPTIONS]\n";
	print STDERR " -i <file>         [input file(s) having NFRs]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
	print STDERR " -k <file>         [input BAM file(s)]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
    print STDERR " -m <float>        [size factor(s) to normalize expression of reads]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
	print STDERR " -o <dir>          [directory where output files will be kept]\n";
	print STDERR "[OPTIONS]\n";
    print STDERR " -g <int>          [minimum length of nucleosome free region (default: 20)]\n";
    print STDERR " -t <int>          [minimum score to consider a NFR as significant]\n";
    print STDERR " -c <int>          [value(s) to extend 3' end of reads by input number of bases (default: 0)]\n";
    print STDERR "                   [if multiple, please separate them by a comma]\n";
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

my $ID=$nfrFiles;
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

#my @data=`intersectBed -a $nfrFileRep1 -b $nfrFileRep2 -wo | sort -k 1,1 -k 2n,2 -k 3n,3`;

my @nfrFileArray=split(/\,/,$nfrFiles);
my $SAMPLES=scalar(@nfrFileArray);

## set default value to extends parameter, if not defined
if(!defined($extends)) {
    $extends="";
    foreach(@nfrFileArray) {
        $extends.="0,";
    }
    $extends=~s/\,$//g;
}

$nfrFiles=~s/\,/ /g;
chomp($nfrFiles);
my @data=();
#system("mergeBed.pl -i $nfrFiles | cut -f 1-3 > $outDir/$ID.Common$fileSuffix");
system("multiIntersectBed -i $nfrFiles | perl -ane 'if(\$F[3]=='$SAMPLES') { print \$_; }' | cut -f 1-3 > $outDir/$ID.Common$fileSuffix");
my %info=();
tie %info, 'Tie::IxHash';
my %count=();
foreach my $file(@nfrFileArray) {
    @data=`intersectBed -a $outDir/$ID.Common$fileSuffix -b $file -wo`;
    foreach(@data) {
        #print "$_";
        my @F=split(/\s+/,$_);
        my $key="$F[0]_$F[1]_$F[2]";
        if(!defined($info{$key}{$file})) {
            $info{$key}{$file}{'id'}=$F[6];
            $info{$key}{$file}{'strand'}=$F[8];
            $info{$key}{$file}{'nfr'}="$F[3]:$F[4]-$F[5]";
            $info{$key}{$file}{'startBlock'}=$F[9];
            $info{$key}{$file}{'endBlock'}=$F[10];
            $info{$key}{$file}{'score'}=$F[7];
        } 
        elsif($F[7] > $info{$key}{$file}{'score'}) {
            $info{$key}{$file}{'id'}=$F[6];
            $info{$key}{$file}{'strand'}=$F[8];
            $info{$key}{$file}{'nfr'}="$F[3]:$F[4]-$F[5]";
            $info{$key}{$file}{'startBlock'}=$F[9];
            $info{$key}{$file}{'endBlock'}=$F[10];
            $info{$key}{$file}{'score'}=$F[7];
        }
        $count{$key}++;
    }
}

#print scalar(keys %count);
#exit;

@data=();
foreach my $key(keys(%count)) {
    if($count{$key}>=scalar(@nfrFileArray)) {
        my @nfr=(); my @startBlock=(); my @endBlock=();
        my $id=(); my $strand=();
        foreach my $file(@nfrFileArray) {
            $id=$info{$key}{$file}{'id'};
            $strand=$info{$key}{$file}{'strand'};
            push(@nfr, $info{$key}{$file}{'nfr'});
            push(@startBlock, $info{$key}{$file}{'startBlock'});
            push(@endBlock, $info{$key}{$file}{'endBlock'});
            #print "$info{$key}{$file}{'id'}\t$info{$key}{$file}{'strand'}\t$info{$key}{$file}{'nfr'}\t$info{$key}{$file}{'startBlock'}\t$info{$key}{$file}{'endBlock'}\n";
        }
        #next if($id!~/^chr22\:17348973\-17352265$/);
        #foreach(@nfr) { print "$_\t"; } print "\n";
        my $nfr=returnOverlapCoorArray(\@nfr, "intersect");
        my $startBlock=returnOverlapCoorArray(\@startBlock, "union");
        my $endBlock=returnOverlapCoorArray(\@endBlock, "union");
        #print "$startBlock\t$nfr\t$endBlock\t$id\t$strand\n";
        push(@data, "$startBlock\t$nfr\t$endBlock\t$id\t$strand");
    }
}
## FOR DEBUGGING
#exit;
#print scalar(@data)."\n"; exit;

## output file for writing overlapping NFR regions
open(OUTFILE, ">$outDir/$ID.All.nfr$fileSuffix") || die $!;
open(SIGFILE, ">$outDir/$ID.All.nfr.sig$fileSuffix") || die $! if(defined($nfrThreshold));

my %NFR=(); my %startBlock=(); my %endBlock=();
foreach my $l(@data) {
    chomp($l);
    my ($startBlock, $nfr, $endBlock, $id, $strand)=split(/\s+/,$l);

    ## check if overlap has been done correctly
    @{$startBlock{'overlapCoorSplit'}}=split(/[\:\-]+/,$startBlock);
    @{$NFR{'overlapCoorSplit'}}=split(/[\:\-]+/,$nfr);
    @{$endBlock{'overlapCoorSplit'}}=split(/[\:\-]+/,$endBlock);

    if($NFR{'overlapCoorSplit'}[1]-$startBlock{'overlapCoorSplit'}[2]>1 || $endBlock{'overlapCoorSplit'}[1]-$NFR{'overlapCoorSplit'}[2]>1) {
        print STDERR "start or end block are not adjacent to NFR\n";
        print STDERR "--> $NFR{'overlapCoorSplit'}[1]\t$startBlock{'overlapCoorSplit'}[2]\t$endBlock{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\n";
        print "$l\n";
        print STDERR "$nfr\t$startBlock\t$endBlock\n\n";
        exit(-1);
    }

    $startBlock{'expr'}=`coor2expr -i $startBlock -j $bamFiles -k $sizeFactors -d -e $extends -g $genome`;
    chomp($startBlock{'expr'});
    $NFR{'expr'}=`coor2expr -i $nfr -j $bamFiles -k $sizeFactors -d -e $extends -g $genome`;
    chomp($NFR{'expr'});
    $endBlock{'expr'}=`coor2expr -i $endBlock -j $bamFiles -k $sizeFactors -d -e $extends -g $genome`;
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
        print OUTFILE "$NFR{'overlapCoorSplit'}[0]\t$NFR{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\t$id\t$NFR{'score'}\t$strand\t$startBlock\t$endBlock\t$startBlock{'expr'}\t$startBlock{'length'}\t$endBlock{'expr'}\t$endBlock{'length'}\t$NFR{'expr'}\t$NFR{'length'}\n";
        if(defined($nfrThreshold) && $NFR{'score'}>$nfrThreshold) {
            print SIGFILE "$NFR{'overlapCoorSplit'}[0]\t$NFR{'overlapCoorSplit'}[1]\t$NFR{'overlapCoorSplit'}[2]\t$id\t$NFR{'score'}\t$strand\t$startBlock\t$endBlock\t$startBlock{'expr'}\t$startBlock{'length'}\t$endBlock{'expr'}\t$endBlock{'length'}\t$NFR{'expr'}\t$NFR{'length'}\n";
        }
    }
}

close OUTFILE;
close SIGFILE if(defined($nfrThreshold));
=cut
exit(0);
