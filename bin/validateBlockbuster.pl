#!/usr/bin/perl -w

=copyright_info
validateBlockbuster.pl: validate block groups in blockbuster output file
Copyright (C) 2014  Sachin Pundhir (sachin@rth.dk)

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

##########################################################################################################
#Parse input options
use vars qw($infile $filterThreshold $overThreshold $distThreshold $removeTags $minBGExpr $minBExpr $formatId $maxBGCount $help);

$overThreshold=100;
$distThreshold=1000000000000;
$minBGExpr=10;
$minBExpr=2;
$filterThreshold=0;

GetOptions ("c=s"  => \$infile,
            "o=i"  => \$overThreshold,
            "d=i"  => \$distThreshold,
            "r"    => \$formatId,
            "f=i"  => \$filterThreshold,
            "n=i"  => \$removeTags,
            "k=i"  => \$minBExpr,
            "g=s"  => \$minBGExpr,
            "x=i"  => \$maxBGCount,
            "help" => \$help,
            "h"    => \$help);

usage() if($help);

###########################################################################################################
sub usage {
	print STDERR "\nProgram: validateBlockbuster.pl (validate block groups in blockbuster output file)\n";
	print STDERR "Note: Recommended after blockbuster to sort blocks by start pos. and to make unique block group ids\n";
	print STDERR "Author: University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: sachin\@rth.dk\n";
	print STDERR "Usage: validateBlockbuster.pl -c <file | STDIN> [OPTIONS]\n";
	print STDERR " -c <file>   [block group file in blockbuster format | STDIN]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -o <int>    [minimum percentage of overlap to merge two blocks (default: 100)]\n";
	print STDERR " -d <int>    [minimum distance to split blocks into seperate block groups (default: 1000000000000)]\n";
	print STDERR " -r          [reformat block group Id]\n";
	print STDERR " -f <int>    [minimum relative percentage of tags to include in block groups]\n";
	print STDERR " -n <int>    [remove tags (1: non-unique; 2: dummy)]\n";
	print STDERR " -k <int>    [minimum reads to print a block (default: 2)]\n"; 
	print STDERR " -g <int>    [minimum reads to print a block group (default: 10)]\n";
	print STDERR " -x <int>    [maximum count of block groups to print]\n";
	print STDERR " -h          [help]\n\n";
	exit(-1);
}
###########################################################################################################
## modify overlap threshold, if 0
if($overThreshold==0) { $overThreshold=0.01; }

my $INFILE=();
if(!defined($infile)) {
	$INFILE = *STDIN;
}
elsif($infile=~/\.gz$/) {
	open($INFILE, "gunzip -c $infile |") || die $!;
}
else {
	open($INFILE, $infile) || die $!;
}
my @data=<$INFILE>;
close $INFILE;

my %bgInfo=(); my %bInfo=(); my %tagInfo=();
my $bg=1;

## read block group information
foreach(@data) {
	if(/^>/) {
		exit if(defined($maxBGCount) && $bg>$maxBGCount);
		if(%bInfo && $removeTags) {
			if(removeTags()>0) {
				$bg = validateBG(); $bg++;
			}
		}
		elsif(%bgInfo && $bgInfo{$bg}{'expr'}>=$minBGExpr) {
			$bg = validateBG();
			$bg++;
		}
		#print STDERR $_;
		%bgInfo=(); %bInfo=(); %tagInfo=();
		my @t=split(/\s+/, $_);
		$bgInfo{$bg}{'id'}=$t[0];
		$bgInfo{$bg}{'expr'}=$t[5];
		$bgInfo{$bg}{'tags'}=$t[6];
		$bgInfo{$bg}{'blocks'}=$t[7];
		$bgInfo{$bg}{'name'}=$t[8] if(defined($t[8]));
		$bgInfo{$bg}{'ncrna'}=$t[9] if(defined($t[9]));
		$bgInfo{$bg}{'loci'}=$t[10] if(defined($t[10]));
	}
	else {
		my ($chr, $start, $end, $id, $expr, $strand, $block) = split(/\s+/, $_);
		if(defined($bInfo{$block})) {
			if($start < $bInfo{$block}{'start'}) { $bInfo{$block}{'start'} = $start; }
			if($end > $bInfo{$block}{'end'}) { $bInfo{$block}{'end'} = $end; }
		}
		else { $bInfo{$block}{'start'} = $start; $bInfo{$block}{'end'} = $end; }

		$bInfo{$block}{'bg'} = $bg;
		$bInfo{$block}{'chr'} = $chr;
		$bInfo{$block}{'strand'} = $strand;

		my $key = $id."_".$start."_".$end;
		$tagInfo{$block}{$key}{'chr'} = $chr;
		$tagInfo{$block}{$key}{'start'} = $start;
		$tagInfo{$block}{$key}{'end'} = $end;
		$tagInfo{$block}{$key}{'id'} = $id;
		$tagInfo{$block}{$key}{'expr'} = $expr;
		$tagInfo{$block}{$key}{'strand'} = $strand;
		$tagInfo{$block}{$key}{'block'} = $block;
	}
}

if(%bInfo && $removeTags) {
	if(removeTags()>0) {
		$bg = validateBG();
	}
}
elsif(%bgInfo && $bgInfo{$bg}{'expr'}>=$minBGExpr) {
	$bg = validateBG();
}
exit;

## remove non-unique or dummy tags from the block group
sub removeTags {
	my @bSorted = sort { $bInfo{$a}{'start'} <=> $bInfo{$b}{'start'} } keys(%bInfo);

	my %EXPR=(); $EXPR{'bg'}=0;
	for(my $i=0; $i<(scalar(@bSorted));) {
		$EXPR{'b'}=0;
		foreach my $key(keys(%{$tagInfo{$bSorted[$i]}})) {
			if(($removeTags==1 && $tagInfo{$bSorted[$i]}{$key}{'id'}!~/\_1$/) || ($removeTags==2 && $tagInfo{$bSorted[$i]}{$key}{'id'}=~/dummy/)) {
				delete($tagInfo{$bSorted[$i]}{$key});
			}
			else {
				$EXPR{'b'}+=$tagInfo{$bSorted[$i]}{$key}{'expr'};
			}
		}

		## delete block, if all tags are removed from the block or if expr is less than $minBExpr
		if(keys(%{$tagInfo{$bSorted[$i]}})==0 || $EXPR{'b'}<$minBExpr) {
			delete($bInfo{$bSorted[$i]});
			delete($tagInfo{$bSorted[$i]});
			splice(@bSorted, $i, 1);
		}
		else {
			my $j=0; $EXPR{'bg'}+=$EXPR{'b'};
			foreach my $key(keys(%{$tagInfo{$bSorted[$i]}})) {
				if($j==0) {
					$bInfo{$bSorted[$i]}{'start'}=$tagInfo{$bSorted[$i]}{$key}{'start'};
					$bInfo{$bSorted[$i]}{'end'}=$tagInfo{$bSorted[$i]}{$key}{'end'};
				}
				else {
					if($tagInfo{$bSorted[$i]}{$key}{'start'}<$bInfo{$bSorted[$i]}{'start'}) {
						$bInfo{$bSorted[$i]}{'start'}=$tagInfo{$bSorted[$i]}{$key}{'start'};
					}
					if($tagInfo{$bSorted[$i]}{$key}{'end'}>$bInfo{$bSorted[$i]}{'end'}) {
						$bInfo{$bSorted[$i]}{'end'}=$tagInfo{$bSorted[$i]}{$key}{'end'};
					}
				}
				$j++;
			}
			$i++;
		}
	}

	## return number of blocks left after removing non-unique tags or if block group expression is >=$minBGExpr
	if($EXPR{'bg'}<$minBGExpr) { return 0; }
	else { return scalar(@bSorted); }
}

## validation block group information
sub validateBG {
	my $overlap=();
	my @bSorted = sort { $bInfo{$a}{'start'} <=> $bInfo{$b}{'start'} } keys(%bInfo);

	#foreach(@bSorted) { print "$_\t$bInfo{$_}{'start'}\t$bInfo{$_}{'end'}\t$bInfo{$_}{'bg'}\n"; }

	#if(defined($filterThreshold)) {
	if($filterThreshold>0) {
		## filter low expressed tags from the block group
		for(my $i=0; $i<(scalar(@bSorted));) {
			## determine the expression profile for each block
			for(my $j=$bInfo{$bSorted[$i]}{'start'}; $j<=$bInfo{$bSorted[$i]}{'end'}; $j++) {
				foreach my $key(keys(%{$tagInfo{$bSorted[$i]}})) {
					if($j>=$tagInfo{$bSorted[$i]}{$key}{'start'} && $j<=$tagInfo{$bSorted[$i]}{$key}{'end'}) {
						$bInfo{$bSorted[$i]}{'expProfile'}{$j}+=$tagInfo{$bSorted[$i]}{$key}{'expr'};
					}
				}
			}
			#for(my $j=$bInfo{$bSorted[$i]}{'start'}; $j<=$bInfo{$bSorted[$i]}{'end'}; $j++) { printf("%0.2f,", ($bInfo{$bSorted[$i]}{'expProfile'}{$j}/$bgInfo{$bInfo{$bSorted[$i]}{'bg'}}{'expr'})*100); } print "\n";

			## determine the significant start and end positions for each block
			my $sigStart=(); my $sigEnd=();
			for(my $j=$bInfo{$bSorted[$i]}{'start'}; $j<=$bInfo{$bSorted[$i]}{'end'}; $j++) {
				if(($bInfo{$bSorted[$i]}{'expProfile'}{$j}/$bgInfo{$bInfo{$bSorted[$i]}{'bg'}}{'expr'})*100 >= $filterThreshold) {
					if(!defined($sigStart)) {
						$sigStart=$j;
					}
					else { $sigEnd=$j; }
				}
			}
			#print "$bSorted[$i]\t$sigStart\t$sigEnd\t$bInfo{$bSorted[$i]}{'bg'}\t".(($sigEnd-$sigStart)+1)."\n";
			if(defined($sigStart) && defined($sigEnd) && (($sigEnd-$sigStart)+1) >= 10) {

				## filter low expressed tags
				foreach my $key(keys(%{$tagInfo{$bSorted[$i]}})) {
					my $overlap = checkOverlap($tagInfo{$bSorted[$i]}{$key}{'start'}, $tagInfo{$bSorted[$i]}{$key}{'end'}, $sigStart, $sigEnd);
					## delete tag, if overlapping at less than 90%
					if($overlap < 90) {
						$bInfo{$bSorted[$i]}{'start'}=$sigStart;
						$bInfo{$bSorted[$i]}{'end'}=$sigEnd;
						delete($tagInfo{$bSorted[$i]}{$key});
					}
				}

				## redefine start and end coordinate of block
				foreach my $key(keys(%{$tagInfo{$bSorted[$i]}})) {
					if($tagInfo{$bSorted[$i]}{$key}{'start'}<$bInfo{$bSorted[$i]}{'start'}) {
						$bInfo{$bSorted[$i]}{'start'}=$tagInfo{$bSorted[$i]}{$key}{'start'};
					}
					if($tagInfo{$bSorted[$i]}{$key}{'end'}>$bInfo{$bSorted[$i]}{'end'}) {
						$bInfo{$bSorted[$i]}{'end'}=$tagInfo{$bSorted[$i]}{$key}{'end'};
					}
					
				}

				## delete block, if all tags are removed from the block
				if(keys(%{$tagInfo{$bSorted[$i]}})==0) {
					delete($bInfo{$bSorted[$i]});
					delete($tagInfo{$bSorted[$i]});
					splice(@bSorted, $i, 1);
				}
				else { $i++; }
			}
			else {
				delete($bInfo{$bSorted[$i]});
				delete($tagInfo{$bSorted[$i]});
				splice(@bSorted, $i, 1);
			}
		}

		if(!defined($bSorted[0]) && defined($infile)) { print STDERR "error: no block left for $bgInfo{$bg}{'id'} ($infile) after tag filtering in validateBlockbuster.pl\n"; return $bg; }
		elsif(!defined($bSorted[0])) { print STDERR "error: no block left for $bgInfo{$bg}{'id'} (STDIN) after tag filtering in validateBlockbuster.pl\n"; return $bg; }
	}

	@bSorted = sort { $bInfo{$a}{'start'} <=> $bInfo{$b}{'start'} } keys(%bInfo);
	## initialize start and end coordinates of $bg
	$bgInfo{$bg}{'start'} = $bInfo{$bSorted[0]}{'start'};
	$bgInfo{$bg}{'end'} = $bInfo{$bSorted[0]}{'end'};

	## reinitialize bgInfo parameters before merging and splitting.
	$bgInfo{$bg}{'expr'}=0;
	$bgInfo{$bg}{'tags'}=0;
	$bgInfo{$bg}{'blocks'}=0;

	## merge overlapping blocks and split highly seperated blocks into distinct block groups
	for(my $i=0; $i<(scalar(@bSorted)-1);) {
		for(my $j=$i+1; $j<scalar(@bSorted);) {
			## compute overlap
			$overlap=0;
			for(my $k=$bInfo{$bSorted[$j]}{'start'}; $k<=$bInfo{$bSorted[$j]}{'end'}; $k++) {
				if($k >= $bInfo{$bSorted[$i]}{'start'} && $k <= $bInfo{$bSorted[$i]}{'end'}) {
					$overlap++;
				}
			}
			#print "$bSorted[$i]\t$bInfo{$bSorted[$i]}{'start'}\t$bInfo{$bSorted[$i]}{'end'}\t$bSorted[$j]\t$bInfo{$bSorted[$j]}{'start'}\t$bInfo{$bSorted[$j]}{'end'}\t";
			if((($bInfo{$bSorted[$i]}{'end'}-$bInfo{$bSorted[$i]}{'start'})+1) <  (($bInfo{$bSorted[$j]}{'end'}-$bInfo{$bSorted[$j]}{'start'})+1)) {
				$overlap = sprintf("%0.2f", ($overlap/(($bInfo{$bSorted[$i]}{'end'}-$bInfo{$bSorted[$i]}{'start'})+1))*100);
			}
			else {
				$overlap = sprintf("%0.2f", ($overlap/(($bInfo{$bSorted[$j]}{'end'}-$bInfo{$bSorted[$j]}{'start'})+1))*100);
			}

			## check for overlapping blocks
			if($overlap >= $overThreshold) {
#print "overlapped\n";
				if($bInfo{$bSorted[$j]}{'end'} > $bInfo{$bSorted[$i]}{'end'}) { $bInfo{$bSorted[$i]}{'end'} = $bInfo{$bSorted[$j]}{'end'}; }
				## assign tags of block $j to block $i
				foreach(keys(%{$tagInfo{$bSorted[$j]}})) {
					$tagInfo{$bSorted[$i]}{$_}{'chr'} = delete $tagInfo{$bSorted[$j]}{$_}{'chr'};
					$tagInfo{$bSorted[$i]}{$_}{'start'} = delete $tagInfo{$bSorted[$j]}{$_}{'start'};
					$tagInfo{$bSorted[$i]}{$_}{'end'} = delete $tagInfo{$bSorted[$j]}{$_}{'end'};
					$tagInfo{$bSorted[$i]}{$_}{'id'} = delete $tagInfo{$bSorted[$j]}{$_}{'id'};
					$tagInfo{$bSorted[$i]}{$_}{'expr'} = delete $tagInfo{$bSorted[$j]}{$_}{'expr'};
					$tagInfo{$bSorted[$i]}{$_}{'strand'} = delete $tagInfo{$bSorted[$j]}{$_}{'strand'};
					$tagInfo{$bSorted[$i]}{$_}{'block'} = $bSorted[$i];
				}
				delete($bInfo{$bSorted[$j]}); splice(@bSorted, $j, 1);
			}

			## check for highly seperated blocks
			elsif((($bInfo{$bSorted[$j]}{'start'} - $bInfo{$bSorted[$i]}{'end'})+1) > $distThreshold) {
				#print "splitted\n";
				## update start and end coordinates of $bg
				if(defined($bgInfo{$bg})) {
					if($bInfo{$bSorted[$i]}{'start'} < $bgInfo{$bg}{'start'}) { $bgInfo{$bg}{'start'} = $bInfo{$bSorted[$i]}{'start'}; }
					if($bInfo{$bSorted[$i]}{'end'} > $bgInfo{$bg}{'end'}) { $bgInfo{$bg}{'end'} = $bInfo{$bSorted[$i]}{'end'}; }
				}
				else {
					$bgInfo{$bg}{'start'} = $bInfo{$bSorted[$i]}{'start'};
					$bgInfo{$bg}{'end'} = $bInfo{$bSorted[$i]}{'end'};
				}

				$bg = $bg + 1;
				$i++;
			}
			else {
				#print "fine\n";
				$i++;
			}
			
			# update or initialize coordinates and other information of $bg
			if(defined($bgInfo{$bg})) {
				if($bInfo{$bSorted[$i]}{'start'} < $bgInfo{$bg}{'start'}) { $bgInfo{$bg}{'start'} = $bInfo{$bSorted[$i]}{'start'}; }
				if($bInfo{$bSorted[$i]}{'end'} > $bgInfo{$bg}{'end'}) { $bgInfo{$bg}{'end'} = $bInfo{$bSorted[$i]}{'end'}; }
			}
			else {
				$bgInfo{$bg}{'start'} = $bInfo{$bSorted[$i]}{'start'};
				$bgInfo{$bg}{'end'} = $bInfo{$bSorted[$i]}{'end'};
				$bgInfo{$bg}{'id'}=$bgInfo{$bg-1}{'id'};
				$bgInfo{$bg}{'id'}=~s/\_[0-9]+$/\_$bg/;
				$bgInfo{$bg}{'name'}=$bgInfo{$bg-1}{'name'} if(defined($bgInfo{$bg-1}{'name'}));
				$bgInfo{$bg}{'ncrna'}=$bgInfo{$bg-1}{'ncrna'} if(defined($bgInfo{$bg-1}{'ncrna'}));
				$bgInfo{$bg}{'loci'}=$bgInfo{$bg-1}{'loci'} if(defined($bgInfo{$bg-1}{'loci'}));
			}
			$bInfo{$bSorted[$i]}{'bg'} = $bg;
			$j=$i+1;
		}
	}

	## compute expression, tags and blocks in $bg
	foreach my $block(@bSorted) {
		#print "$block\t$bInfo{$block}{'start'}\t$bInfo{$block}{'end'}\t$bInfo{$block}{'bg'}\t$bgInfo{$bInfo{$block}{'bg'}}{'start'}\t$bgInfo{$bInfo{$block}{'bg'}}{'end'}\t".keys(%{$tagInfo{$block}})."\n";
		foreach(keys(%{$tagInfo{$block}})) { $bgInfo{$bInfo{$block}{'bg'}}{'expr'} += $tagInfo{$block}{$_}{'expr'}; }
		$bgInfo{$bInfo{$block}{'bg'}}{'tags'} += keys(%{$tagInfo{$block}});
		$bgInfo{$bInfo{$block}{'bg'}}{'blocks'}++;
	}

	## print validated block groups
	my %seen=(); my $bCount=(); my $tagExpr=();
	foreach my $block(@bSorted) {
		if(!defined($seen{$bInfo{$block}{'bg'}})) {
			#if(defined($bgInfo{$bInfo{$block}{'bg'}}{'name'}) && !defined($formatId)) {
			if(!defined($formatId)) {
				print "$bgInfo{$bInfo{$block}{'bg'}}{'id'}\t$bInfo{$block}{'chr'}\t$bgInfo{$bInfo{$block}{'bg'}}{'start'}\t$bgInfo{$bInfo{$block}{'bg'}}{'end'}\t$bInfo{$block}{'strand'}\t$bgInfo{$bInfo{$block}{'bg'}}{'expr'}\t$bgInfo{$bInfo{$block}{'bg'}}{'tags'}\t$bgInfo{$bInfo{$block}{'bg'}}{'blocks'}";
				if(defined($bgInfo{$bInfo{$block}{'bg'}}{'name'})) {
					print "\t$bgInfo{$bInfo{$block}{'bg'}}{'name'}\t$bgInfo{$bInfo{$block}{'bg'}}{'ncrna'}\t$bgInfo{$bInfo{$block}{'bg'}}{'loci'}\n";
				}
				else { print "\tn/a\tn/a\tn/a\n"; }
			}
			elsif($formatId) {
				print ">BG_$bInfo{$block}{'bg'}\t$bInfo{$block}{'chr'}\t$bgInfo{$bInfo{$block}{'bg'}}{'start'}\t$bgInfo{$bInfo{$block}{'bg'}}{'end'}\t$bInfo{$block}{'strand'}\t$bgInfo{$bInfo{$block}{'bg'}}{'expr'}\t$bgInfo{$bInfo{$block}{'bg'}}{'tags'}\t$bgInfo{$bInfo{$block}{'bg'}}{'blocks'}";
				if(defined($bgInfo{$bInfo{$block}{'bg'}}{'name'})) {
#if(!defined($bgInfo{$bInfo{$block}{'bg'}}{'ncrna'})) { print STDERR "Undefined ncrna\n"; exit(-1); }
if(!defined($bgInfo{$bInfo{$block}{'bg'}}{'loci'})) { $bgInfo{$bInfo{$block}{'bg'}}{'loci'}="."; }
					print "\t$bgInfo{$bInfo{$block}{'bg'}}{'name'}\t$bgInfo{$bInfo{$block}{'bg'}}{'ncrna'}\t$bgInfo{$bInfo{$block}{'bg'}}{'loci'}\n";
				}
				else { print "\tn/a\tn/a\tn/a\n"; }
			}
			$seen{$bInfo{$block}{'bg'}}=1;
			$bCount=1;
		}
		foreach my $key(keys(%{$tagInfo{$block}})) {
			print "$tagInfo{$block}{$key}{'chr'}\t$tagInfo{$block}{$key}{'start'}\t$tagInfo{$block}{$key}{'end'}\t$tagInfo{$block}{$key}{'id'}\t$tagInfo{$block}{$key}{'expr'}\t$tagInfo{$block}{$key}{'strand'}\t$bCount\n";
		}
		$bCount++;
	}
	return $bg;
}

## subroutine to determine overlap between tag and block coordinates
sub checkOverlap {
	my ($tStart, $tEnd, $bStart, $bEnd) = @_;
	my $overlap=0;
	for(my $k=$tStart; $k<=$tEnd; $k++) {
		if($k>=$bStart && $k<=$bEnd) {
			$overlap++;
		}
	}
	$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	#print "$tStart\t$tEnd\t$bStart\t$bEnd\t$overlap\n";
	return $overlap;
}
