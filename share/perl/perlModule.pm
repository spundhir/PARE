use strict;
use warnings;
package perlModule;

use Carp ();
use Exporter;
use Tie::IxHash;

our $VERSION = 1.0;
our @ISA = qw(Exporter);
our @EXPORT = qw(returnOverlapCoorArray organizeOverlapCoor returnOverlapCoor checkOverlap compReadCount compOverlap readbgInfo compNuclFreq reverse_complement mergeBG parseBED openFile grepBG uniqCoor addRead removeBlock addBlock mapBlock uniqElements reformatbId compbgRelExp printbgInfo checkMapFormat parsePvclust binarySearchCoor checkArrayOverlap commonHashKeys parseDBAoutput checkRevCompExp combinePDF combineMultiPDF coor2annotation);

## return overlapping coordinates between two coordinates
sub returnOverlapCoorArray {
	my ($coor, $mode) = @_;
	return 0 if(!defined($coor));
    my @chr=(); my @start=(); my @end=(); 
    foreach(@{$coor}) {
        #print "$_\n";
        my @F=split(/[\:\-]+/,$_);
        push(@chr, $F[0]);
        push(@start, $F[1]);
        push(@end, $F[2]);
    }

    my @sorted_start=sort { $a <=> $b } @start;
    my @sorted_end=sort { $a <=> $b } @end;
    if($mode=~/intersect/) {
        return "$chr[0]:$sorted_start[scalar(@sorted_start)-1]-$sorted_end[0]";
    }
    else {
        return "$chr[0]:$sorted_start[0]-$sorted_end[scalar(@sorted_start)-1]";
    }
}

## organize overlapping coordinates into three (middle overlapping and flanking coordinates)
sub organizeOverlapCoor {
	my ($tStart, $tEnd, $bStart, $bEnd, $tChr, $bChr) = @_;
	return 0 if(!defined($tStart) || !defined($tEnd) || !defined($bStart) || !defined($bEnd));
	if(defined($tChr) && defined($bChr)) {
		return 0 if($tChr!~/^$bChr$/i);
	}

	my $overlap=0;
	for(my $k=$tStart; $k<=$tEnd; $k++) {
		if($k>=$bStart && $k<=$bEnd) {
			$overlap++;
		}
	}

    ## check if the two coordinates are not completely overlapping 
	if($tEnd-$tStart <= $bEnd-$bStart) {
		$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	}
	else {
		if((($bEnd-$bStart)+1)==0) { print STDERR "Error, start=$bStart and end=$bEnd\n"; exit; }
		$overlap=sprintf("%0.2f", (($overlap/(($bEnd-$bStart)+1))*100));
	}

    ## organize the two coordinates
    my @sorted=sort { $a <=> $b } ($tStart, $tEnd, $bStart, $bEnd);
    my %coor=();
    if($overlap<100) {
        $coor{'left'}=sprintf("%s:%d-%d", $tChr, $sorted[0], $sorted[1]);
        $coor{'mid'}=sprintf("%s:%d-%d", $tChr, ($sorted[1]+1), ($sorted[2]-1));
        $coor{'right'}=sprintf("%s:%d-%d", $tChr, $sorted[2], $sorted[3]);
        if($tStart < $bStart) {
            return($coor{'left'}, $coor{'mid'}, $coor{'right'}, 1);
        }
        else {
            return($coor{'left'}, $coor{'mid'}, $coor{'right'}, 2);
        }
    }
    else {
        return(0, 0, 0, 0);
    }
}

## return overlapping coordinates between two coordinates
sub returnOverlapCoor {
	my ($tStart, $tEnd, $bStart, $bEnd, $tChr, $bChr, $mode, $elseReturn) = @_;
	return 0 if(!defined($tStart) || !defined($tEnd) || !defined($bStart) || !defined($bEnd));
	if(defined($tChr) && defined($bChr)) {
		return 0 if($tChr!~/^$bChr$/i);
	}

	my $overlap=0;
	for(my $k=$tStart; $k<=$tEnd; $k++) {
		if($k>=$bStart && $k<=$bEnd) {
			$overlap++;
		}
	}

    my @sorted=sort { $a <=> $b } ($tStart, $tEnd, $bStart, $bEnd);
    if($overlap>0 && $mode=~/intersect/) {
        return "$tChr:$sorted[1]-$sorted[2]";
    }
    elsif($overlap>0 && $mode=~/union/) {
        return "$tChr:$sorted[0]-$sorted[3]";
    }
    elsif($mode=~/union/ && $elseReturn==1) {
        return "$tChr:$sorted[2]-$sorted[3]";
    }
    elsif($mode=~/union/ && $elseReturn==2) {
        return "$tChr:$sorted[0]-$sorted[1]";
    }
    else {
        print STDERR "Input coordinates do not overlap";
        exit(-1);
    }
}

## check overlap between two coordinates
sub checkOverlap {
	my ($tStart, $tEnd, $bStart, $bEnd, $percentage, $tChr, $bChr, $tStrand, $bStrand) = @_;
	return 0 if(!defined($tStart) || !defined($tEnd) || !defined($bStart) || !defined($bEnd));
	if(defined($tChr) && defined($bChr)) {
		return 0 if($tChr!~/^$bChr$/i);
	}
	if(defined($tStrand) && defined($bStrand)) {
		return 0 if($tStrand!~/\Q$bStrand\E/);
	}

	my $overlap=0;
	for(my $k=$tStart; $k<=$tEnd; $k++) {
		if($k>=$bStart && $k<=$bEnd) {
			$overlap++;
		}
	}

	## take first coordinate as reference
	if(defined($percentage) && $percentage==2) {
		$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	}
	## take second coordinate as reference
	elsif(defined($percentage) && $percentage==3) {
		$overlap=sprintf("%0.2f", (($overlap/(($bEnd-$bStart)+1))*100));
	}
	## take shortest coordinate as reference
	elsif(defined($percentage) && $percentage==1 && $tEnd-$tStart <= $bEnd-$bStart) {
		$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	}
	elsif(defined($percentage) && $percentage==1) {
		if((($bEnd-$bStart)+1)==0) { print STDERR "Error, start=$bStart and end=$bEnd\n"; exit; }
		$overlap=sprintf("%0.2f", (($overlap/(($bEnd-$bStart)+1))*100));
	}
	## take longest coordinate as reference
	elsif(defined($percentage) && $percentage==5 && $tEnd-$tStart >= $bEnd-$bStart) {
		$overlap=sprintf("%0.2f", (($overlap/(($tEnd-$tStart)+1))*100));
	}
	elsif(defined($percentage) && $percentage==5) {
		if((($bEnd-$bStart)+1)==0) { print STDERR "Error, start=$bStart and end=$bEnd\n"; exit; }
		$overlap=sprintf("%0.2f", (($overlap/(($bEnd-$bStart)+1))*100));
	}
	## return 1 if overlap, return 0 otherwise
	elsif(defined($percentage) && $percentage==4 && $overlap>0) {
		$overlap=1;
	}
	#print "$tStart\t$tEnd\t$bStart\t$bEnd\t$overlap\n";
	return $overlap;
}

## compute read count for each coordinate across all samples
## useful for differential expression analysis
sub compReadCount {
	my ($dir, $files, $countFile, $mode, $uniqReadsFilter, $chrFilter, $mapDir) = @_;

	my %inOverlapData=();
	my %outOverlapData=();

	tie %inOverlapData, 'Tie::IxHash';
	tie %outOverlapData, 'Tie::IxHash';

	open(OUTFILE, ">$countFile") || die $!;

	## retrieve coordinates for all blocks of each sample
	print OUTFILE "\t";
	foreach my $f(@{$files}) {
		print OUTFILE "$f\t";
		my %bgInfo=();
		readbgInfo(\%bgInfo, $dir, $f, 1, $uniqReadsFilter, $chrFilter);
		foreach my $bg( sort { $bgInfo{$a}{'start'} <=> $bgInfo{$b}{'start'} } keys(%bgInfo)) {
			## estimate size factor using expression at each block of block group
			if($mode==1) {
				foreach my $b(sort { $bgInfo{$bg}{$a}{'start'} <=> $bgInfo{$bg}{$b}{'start'} } grep(/^[0-9]+$/, keys(%{$bgInfo{$bg}}))) {
					## key is composed of a unique bgId with bId
					my $bgb = $bgInfo{$bg}{'id'}."_".$b;
					$inOverlapData{$f}{$bgb}{'id'} = $bgb;
					$inOverlapData{$f}{$bgb}{'chr'} = $bgInfo{$bg}{'chr'};
					$inOverlapData{$f}{$bgb}{'start'} = $bgInfo{$bg}{$b}{'start'};
					$inOverlapData{$f}{$bgb}{'end'} = $bgInfo{$bg}{$b}{'end'};
					$inOverlapData{$f}{$bgb}{'strand'} = $bgInfo{$bg}{'strand'};
					$inOverlapData{$f}{$bgb}{'reads'} = $bgInfo{$bg}{$b}{'reads'};
				}
			}
			## estimate size factor using expression at each block group
			else {
				## key is composed of a unique bgId
				$inOverlapData{$f}{$bg}{'id'} = $bg;
				$inOverlapData{$f}{$bg}{'chr'} = $bgInfo{$bg}{'chr'};
				$inOverlapData{$f}{$bg}{'start'} = $bgInfo{$bg}{'start'};
				$inOverlapData{$f}{$bg}{'end'} = $bgInfo{$bg}{'end'};
				$inOverlapData{$f}{$bg}{'strand'} = $bgInfo{$bg}{'strand'};
				$inOverlapData{$f}{$bg}{'reads'} = $bgInfo{$bg}{'reads'};
			}
		}
	}
	print OUTFILE "\n";

	## compute overlapping bg or b to given set of coordinates across samples
	compOverlap(\%inOverlapData, \%outOverlapData);

	## print the uniqe coordinates observed in atleast one sample
	my $COUNTER=0;
	foreach my $key(keys(%outOverlapData)) {
		#print "$key\t";
		for(my $i=0; $i<scalar(@{$outOverlapData{$key}{'strand'}}); $i++) {
			print OUTFILE "$outOverlapData{$key}{'chr'}:$outOverlapData{$key}{'start'}-$outOverlapData{$key}{'end'}|$outOverlapData{$key}{'strand'}[$i]\t";
			foreach my $f(@{$files}) {
				my $reads=0; my @data=();
				my $chr = $outOverlapData{$key}{'chr'};
				my $start = $outOverlapData{$key}{'start'};
				my $end = $outOverlapData{$key}{'end'};
				my $strand = $outOverlapData{$key}{'strand'}[$i];

				## retrieve sample id
				my $sampleId=$f;
				$sampleId=~s/\..+//g;
	
				## print read counts from map file, if provided			
				if(defined($mapDir)) {
					if($strand=~/^\+$/) {
						@data=`coor2mapping.pl -c $chr:$start-$end -f $mapDir/$sampleId.*.sorted.$chr.P.gz -s $strand -r segemehl | cut -f 2`;
						if($?!=0) {
							print "Cannot find file, $mapDir/$sampleId.*.sorted.$chr.P.gz\n";
							exit(-1);
						}
					}
					else {
						@data=`coor2mapping.pl -c $chr:$start-$end -f $mapDir/$sampleId.*.sorted.$chr.N.gz -s $strand -r segemehl | cut -f 2`;
						if($?!=0) {
							print "Cannot find file, $mapDir/$sampleId.*.sorted.$chr.N.gz\n";
							exit(-1);
						}
					}
					foreach(@data) { $_=~s/.+\|//g; $reads+=$_; }
					print OUTFILE "$reads\t";
				}
				## print read counts as computed from block groups
				else {
					if(defined($outOverlapData{$key}{$strand}{$f})) {
						print OUTFILE "$outOverlapData{$key}{$strand}{$f}\t";
					}
					else {
						print OUTFILE "0\t";
					}
				}
			}
			print OUTFILE "\n";
		}
		$COUNTER++;
		print STDERR "Done $COUNTER out of ".keys(%outOverlapData)."\n";
	}

	close OUTFILE;
}

## compute block groups (bg) or blocks (b) falling within a coordinate across two or more samples (files)
sub compOverlap {
	my ($inOverlapData, $outOverlapData) = @_;

	## variables
	my %bgInfo=();
	my %exprLoci=();
	my %index=();

	tie %bgInfo, 'Tie::IxHash';

	foreach my $f(keys(%{$inOverlapData})) {
		foreach my $bgb( sort { $inOverlapData->{$f}{$a}{'start'} <=> $inOverlapData->{$f}{$b}{'start'} } keys(%{$inOverlapData->{$f}})) {
			push(@{$bgInfo{$f}{'id'}}, $inOverlapData->{$f}{$bgb}{'id'});
			push(@{$bgInfo{$f}{'chr'}}, $inOverlapData->{$f}{$bgb}{'chr'});
			push(@{$bgInfo{$f}{'start'}}, $inOverlapData->{$f}{$bgb}{'start'});
			push(@{$bgInfo{$f}{'end'}}, $inOverlapData->{$f}{$bgb}{'end'});
			push(@{$bgInfo{$f}{'strand'}}, $inOverlapData->{$f}{$bgb}{'strand'});
			push(@{$bgInfo{$f}{'reads'}}, $inOverlapData->{$f}{$bgb}{'reads'});

			my $key = $inOverlapData->{$f}{$bgb}{'start'}."_".$inOverlapData->{$f}{$bgb}{'end'}."_".$inOverlapData->{$f}{$bgb}{'chr'};
			$exprLoci{$key}{'chr'}=$inOverlapData->{$f}{$bgb}{'chr'};
			$exprLoci{$key}{'start'}=$inOverlapData->{$f}{$bgb}{'start'};
			$exprLoci{$key}{'end'}=$inOverlapData->{$f}{$bgb}{'end'};
			$exprLoci{$key}{'strand'}=$inOverlapData->{$f}{$bgb}{'strand'};
			$exprLoci{$key}{'winLen'}=$inOverlapData->{$f}{$bgb}{'end'}-$inOverlapData->{$f}{$bgb}{'start'};
		}
		
		## initialize start index to 0 for each file
		$index{$f}=0;
	}

	## sort all loci by start and end where a block group is observed in atleast one sample
	#my @exprLociKeys = sort { $exprLoci{$a}{'start'} <=> $exprLoci{$b}{'start'} } keys(%exprLoci);
	#foreach(@exprLociKeys) { print "$_\t$exprLoci{$_}{'start'}\t$exprLoci{$_}{'end'}\t$exprLoci{$_}{'winLen'}\n"; };
	my @exprLociKeysAll = sort { ($exprLoci{$a}{'start'} <=> $exprLoci{$b}{'start'}) || ($exprLoci{$b}{'end'} <=> $exprLoci{$a}{'end'}) } keys(%exprLoci);
	my @exprLociKeys=(); my %lastSeenChr=(); 

	## loop over for each sorted loci to determine unique loci with expression
	for(my $i=0; $i<scalar(@exprLociKeysAll); $i++) {

		## next if loci's start coordinate < last coordinate analyzed for this chromosome
		if(defined($lastSeenChr{$exprLoci{$exprLociKeysAll[$i]}{'chr'}}{'start'}) && $exprLoci{$exprLociKeysAll[$i]}{'start'} <= $lastSeenChr{$exprLoci{$exprLociKeysAll[$i]}{'chr'}}{'start'}+$lastSeenChr{$exprLoci{$exprLociKeysAll[$i]}{'chr'}}{'winLen'}) { next; }

		push(@exprLociKeys, $exprLociKeysAll[$i]);
		#print "$exprLoci{$exprLociKeysAll[$i]}{'chr'}\t$exprLoci{$exprLociKeysAll[$i]}{'start'}\t$exprLoci{$exprLociKeysAll[$i]}{'end'}\n";
		## retrieve block groups overlapping this loci
		## loop over block groups in each sample
		foreach my $f(keys(%bgInfo)) {
			my $start=$index{$f}; my $end=scalar(@{$bgInfo{$f}{'start'}});

			for(my $j=$start; $j<$end; $j++) {
				#print "\t$bgInfo{$f}{'chr'}[$j]\t$bgInfo{$f}{'start'}[$j]\t$bgInfo{$f}{'end'}[$j]\t";
				if($bgInfo{$f}{'start'}[$j] >= $exprLoci{$exprLociKeysAll[$i]}{'start'} && $bgInfo{$f}{'start'}[$j] <= $exprLoci{$exprLociKeysAll[$i]}{'start'}+$exprLoci{$exprLociKeysAll[$i]}{'winLen'} && $bgInfo{$f}{'chr'}[$j]=~/^$exprLoci{$exprLociKeysAll[$i]}{'chr'}$/) {
					if($bgInfo{$f}{'end'}[$j] > $exprLoci{$exprLociKeysAll[$i]}{'end'}) {
					#if(($bgInfo{$f}{'end'}[$j]-$bgInfo{$f}{'start'}) > $exprLoci{$exprLociKeysAll[$i]}{'winLen'}) {
						$exprLoci{$exprLociKeysAll[$i]}{'end'}=$bgInfo{$f}{'end'}[$j];
						$exprLoci{$exprLociKeysAll[$i]}{'winLen'}=$bgInfo{$f}{'end'}[$j]-$exprLoci{$exprLociKeysAll[$i]}{'start'};
					}
				}
				## increment over block groups, if bg start < loci's start coordinate
				if($bgInfo{$f}{'start'}[$j] < $exprLoci{$exprLociKeysAll[$i]}{'start'}) { $index{$f}++; }
				## last if bg start > loci's end coordinate
				if($bgInfo{$f}{'start'}[$j] > $exprLoci{$exprLociKeysAll[$i]}{'start'}+$exprLoci{$exprLociKeysAll[$i]}{'winLen'}) { last; }
			}
		}

		## initialize the start and end coordinate of the last seen loci
		$lastSeenChr{$exprLoci{$exprLociKeysAll[$i]}{'chr'}}{'start'}=$exprLoci{$exprLociKeysAll[$i]}{'start'};
		$lastSeenChr{$exprLoci{$exprLociKeysAll[$i]}{'chr'}}{'winLen'}=$exprLoci{$exprLociKeysAll[$i]}{'winLen'};
		#print "$exprLoci{$exprLociKeysAll[$i]}{'chr'}\t$lastSeenChr{$exprLoci{$exprLociKeysAll[$i]}{'chr'}}{'start'}\n";
	}

	my %bgHits=(); my $found=0; %lastSeenChr=();
	foreach my $f(keys(%bgInfo)) { $index{$f}=0; }

	## loop over for each sorted loci
	for(my $i=0; $i<scalar(@exprLociKeys); $i++) {

		## next if loci's start coordinate < last coordinate analyzed for this chromosome
		if(defined($lastSeenChr{$exprLoci{$exprLociKeys[$i]}{'chr'}}{'start'}) && $exprLoci{$exprLociKeys[$i]}{'start'} <= $lastSeenChr{$exprLoci{$exprLociKeys[$i]}{'chr'}}{'start'}+$lastSeenChr{$exprLoci{$exprLociKeys[$i]}{'chr'}}{'winLen'}) { next; }

		#print "$exprLoci{$exprLociKeys[$i]}{'chr'}\t$exprLoci{$exprLociKeys[$i]}{'start'}\t$exprLoci{$exprLociKeys[$i]}{'end'}\n";
		## retrieve block groups overlapping this loci
		## loop over block groups in each sample
		my @hits=();
		foreach my $f(keys(%bgInfo)) {
			my $start=$index{$f}; my $end=scalar(@{$bgInfo{$f}{'start'}});
			$found=0; my $hits=();

			for(my $j=$start; $j<$end; $j++) {
				#print "\t$bgInfo{$f}{'chr'}[$j]\t$bgInfo{$f}{'start'}[$j]\t$bgInfo{$f}{'end'}[$j]\t";
				if($bgInfo{$f}{'start'}[$j] >= $exprLoci{$exprLociKeys[$i]}{'start'} && $bgInfo{$f}{'start'}[$j] <= $exprLoci{$exprLociKeys[$i]}{'start'}+$exprLoci{$exprLociKeys[$i]}{'winLen'} && $bgInfo{$f}{'chr'}[$j]=~/^$exprLoci{$exprLociKeys[$i]}{'chr'}$/) {
					$found=1;
					#push(@hits, $bgInfo{$f}{'id'}[$j]); last;
					$hits .= "$bgInfo{$f}{'id'}[$j]($bgInfo{$f}{'strand'}[$j]),";

					## compute read count at this coordinate
					if($bgInfo{$f}{'strand'}[$j]=~/^\+$/) {
						$outOverlapData->{$exprLociKeys[$i]}{'+'}{$f}+=sprintf("%0.0f", $bgInfo{$f}{'reads'}[$j]);
					}
					elsif($bgInfo{$f}{'strand'}[$j]=~/^\-$/) {
						$outOverlapData->{$exprLociKeys[$i]}{'-'}{$f}+=sprintf("%0.0f", $bgInfo{$f}{'reads'}[$j]);
					}
				}
				## increment over block groups, if bg start < loci's start coordinate
				if($bgInfo{$f}{'start'}[$j] < $exprLoci{$exprLociKeys[$i]}{'start'}) { $index{$f}++; }
				## last if bg start > loci's end coordinate
				if($bgInfo{$f}{'start'}[$j] > $exprLoci{$exprLociKeys[$i]}{'start'}+$exprLoci{$exprLociKeys[$i]}{'winLen'}) { last; }
			}
			#print "$found\n";
			## put "-" if no block group is observed at the loci for this sample
			if(!$found) { push(@hits, "-");	}
			else { $hits=~s/\,$//; push(@hits, $hits); }
		}

		## results
		$outOverlapData->{$exprLociKeys[$i]}{'chr'}=$exprLoci{$exprLociKeys[$i]}{'chr'};
		$outOverlapData->{$exprLociKeys[$i]}{'start'}=$exprLoci{$exprLociKeys[$i]}{'start'};
		$outOverlapData->{$exprLociKeys[$i]}{'end'}=$exprLoci{$exprLociKeys[$i]}{'end'};

		if(defined($outOverlapData->{$exprLociKeys[$i]}{'+'})) {
			push(@{$outOverlapData->{$exprLociKeys[$i]}{'strand'}}, "+");

		}
		if(defined($outOverlapData->{$exprLociKeys[$i]}{'-'})) {
			push(@{$outOverlapData->{$exprLociKeys[$i]}{'strand'}}, "-");
		}

		foreach(@hits) {
			push(@{$outOverlapData->{$exprLociKeys[$i]}{'hits'}}, $_);
		}

		## initialize the start and end coordinate of the last seen loci
		$lastSeenChr{$exprLoci{$exprLociKeys[$i]}{'chr'}}{'start'}=$exprLoci{$exprLociKeys[$i]}{'start'};
		$lastSeenChr{$exprLoci{$exprLociKeys[$i]}{'chr'}}{'winLen'}=$exprLoci{$exprLociKeys[$i]}{'winLen'};
		#print "$exprLoci{$exprLociKeys[$i]}{'chr'}\t$lastSeenChr{$exprLoci{$exprLociKeys[$i]}{'chr'}}{'start'}\n";
	}
}

## function to read block group information
## merge = 1: merge overlapping, blocks 0: do not merge
sub readbgInfo {
	my($bgInfo, $outDir, $outFile, $merge, $uniqReadsFilter, $chrFilter) = @_;

	if(!defined($merge)) { $merge=0; }
	if(!defined($chrFilter)) { $chrFilter=".+"; }
	if(!defined($uniqReadsFilter)) { $uniqReadsFilter=0; }

	#open(INFILE, "$outDir/$outFile") || die $!;
	my @data = openFile("$outDir/$outFile");

	my $bg=(); my $b=(); my $hit=0; my $uniqReads=0;
	foreach my $l(@data) {
		chomp($l);
		my @t=split(/\s+/, $l);

		## compute total number of uniquely mapped reads in a block group
		$uniqReads = `zless $outDir/$outFile | grep -w '$t[0]' -A $t[6] | grep -v "^>" | cut -f 4,5 | perl -an -F'/\\s+/' -e 'BEGIN { \$sum=0; } if(\$F[0]=~/\\_1\$/) { \$sum+=\$F[1]; } END { print ((\$sum/$t[5])*100); }'` if($l=~/^\>/ && $uniqReadsFilter>0);
		#if($uniqReads>0) { print "$outFile\t$t[0]\t$uniqReads\n"; exit; }

		if($l=~/^\>/ && $t[1]=~/^$chrFilter$/i && $uniqReads>=$uniqReadsFilter) {
			$bg=$t[0];
			$bgInfo->{$bg}{'id'}=$t[0];
			$bgInfo->{$bg}{'chr'}=$t[1];
			$bgInfo->{$bg}{'start'}=$t[2];
			$bgInfo->{$bg}{'end'}=$t[3];
			$bgInfo->{$bg}{'strand'}=$t[4];
			$bgInfo->{$bg}{'reads'}=$t[5];
			$bgInfo->{$bg}{'tags'}=$t[6];
			$bgInfo->{$bg}{'blocks'}=$t[7];
			$bgInfo->{$bg}{'name'}=$t[8];
			$bgInfo->{$bg}{'anno'}=$t[9];
			$bgInfo->{$bg}{'loci'}=$t[10];
			$bgInfo->{$bg}{'origReads'}=$t[5];
			my $sample=$bg; $sample=~s/\_.+//g; $sample=~s/\>//g;
			if($sample!~/^$/ && $sample!~/cluster/ && $sample!~/BG/) { $bgInfo->{$bg}{'sample'}=$sample; }
			else { $bgInfo->{$bg}{'sample'}="NA"; }

			$b=1; $hit=0;
		}
		elsif($l!~/^\s+$/ && $merge && $t[0]=~/^$chrFilter$/i && $uniqReads>=$uniqReadsFilter) {
			if(!defined($bgInfo->{$bg}{$b})) {
				$bgInfo->{$bg}{$b}{'start'}=$t[1];
				$bgInfo->{$bg}{$b}{'end'}=$t[2];
			}
			elsif($t[1]<$bgInfo->{$bg}{$b}{'start'} && $t[2]>=$bgInfo->{$bg}{$b}{'start'}) {
				#if($t[6]!=$b && $t[1]<$bgInfo->{$bg}{$b}{'end'}) { print "$l\n"; $hit=1; }
				$bgInfo->{$bg}{$b}{'start'}=$t[1];
			}
			elsif($t[2]>$bgInfo->{$bg}{$b}{'end'} && $t[1]<=$bgInfo->{$bg}{$b}{'end'}) {
				#if($t[6]!=$b && $t[1]<$bgInfo->{$bg}{$b}{'end'}) { print "$l\n"; $hit=1; }
				$bgInfo->{$bg}{$b}{'end'}=$t[2];
			}
			elsif($t[1]>$bgInfo->{$bg}{$b}{'end'}) {
				#if($hit) { print "\n$l\n"; exit; }
				$b++;
				$bgInfo->{$bg}{$b}{'start'}=$t[1];
				$bgInfo->{$bg}{$b}{'end'}=$t[2];
			}

			$bgInfo->{$bg}{$b}{'reads'}+=$t[4];
			my $lF="$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$b\n";
			push(@{$bgInfo->{$bg}{$b}{'info'}}, $lF);
		}
		elsif($l!~/^\s+$/ && $t[0]=~/^$chrFilter$/i && $uniqReads>=$uniqReadsFilter) {
			if(!defined($bgInfo->{$bg}{$t[6]})) {
				$bgInfo->{$bg}{$t[6]}{'start'}=$t[1];
				$bgInfo->{$bg}{$t[6]}{'end'}=$t[2];
			}
            else {
    			if($t[1]<$bgInfo->{$bg}{$t[6]}{'start'}) {
    				$bgInfo->{$bg}{$t[6]}{'start'}=$t[1];
    			}
    			if($t[2]>$bgInfo->{$bg}{$t[6]}{'end'}) {
    				$bgInfo->{$bg}{$t[6]}{'end'}=$t[2];
	    		}
            }
			$bgInfo->{$bg}{$t[6]}{'reads'}+=$t[4];
			push(@{$bgInfo->{$bg}{$t[6]}{'info'}}, $l);
		}
	}
}

## function to compute nucleotide frequencies
sub compNuclFreq {
	my ($seq, $win, $steps, $entity, $frequency) = @_;

	my @MN=qw(A T G C);

	if($entity=~/mono/) {
		## initialize frequency to 0
		foreach(@MN) { $frequency->{$_}=0; }

		## count mononucleotide frequency in all subsequences for given window size and steps
		for(my $i=0; $i<=(length($seq)-$win); $i+=$steps) {
			my $subSeq=substr($seq, $i, $win);
			for(my $j=0; $j<scalar(@MN); $j++) {
				my @freq=();
				@freq=$subSeq=~m/$MN[$j]/ig;
				$frequency->{$MN[$j]}+=scalar(@freq);
				#print "$subSeq\t".scalar(@freq)."\n";
			}
		}
		#foreach(keys(%{$frequency})) { print "$_\t$frequency->{$_}\n"; }
	}
	elsif($entity=~/di/) {
		## extract all dinucleotides
		my @DN=();
		for (my $i=0; $i<scalar(@MN); $i++) {
			for(my $j=0; $j<scalar(@MN); $j++) {
				push(@DN, "$MN[$i]$MN[$j]");
				## initialize frequency to 0
				$frequency->{"$MN[$i]$MN[$j]"}=0;
			}
		}
		## count dinucleotide frequency in all subsequences for given window size and steps
		for(my $i=0; $i<=(length($seq)-$win); $i+=$steps) {
			my $subSeq=substr($seq, $i, $win);
			for(my $j=0; $j<scalar(@DN); $j++) {
				my @freq=();
				@freq=$subSeq=~m/$DN[$j]/ig;
				$frequency->{$DN[$j]}+=scalar(@freq);
				#print "$subSeq\t".scalar(@freq)."\n";
			}
		}
		#foreach(keys(%{$frequency})) { print "$_\t$frequency->{$_}\n"; }
	}
	elsif($entity=~/tri/) {
		## extract all trinucleotides
		my @TN=();
		for (my $i=0; $i<scalar(@MN); $i++) {
			for(my $j=0; $j<scalar(@MN); $j++) {
				for(my $k=0; $k<scalar(@MN); $k++) {
					push(@TN, "$MN[$i]$MN[$j]$MN[$k]");
					## initialize frequency to 0
					$frequency->{"$MN[$i]$MN[$j]$MN[$k]"}=0;
				}
			}
		}
		## count trinucleotide frequency in all subsequences for given window size and steps
		for(my $i=0; $i<=(length($seq)-$win); $i+=$steps) {
			my $subSeq=substr($seq, $i, $win);
			for(my $j=0; $j<scalar(@TN); $j++) {
				my @freq=();
				@freq=$subSeq=~m/$TN[$j]/ig;
				$frequency->{$TN[$j]}+=scalar(@freq);
				#print "$subSeq\t".scalar(@freq)."\n";
			}
		}
		#foreach(keys(%{$frequency})) { print "$_\t$frequency->{$_}\n"; }
	}
}

sub reverse_complement {
	my $dna = shift;

	my $revcomp=();	
	foreach(split(/\n/, $dna)) {
		if($_=~/^\>/) { $revcomp.="$_\n"; }
		else {
			# reverse the DNA sequence
			$_= reverse($_);

			# complement the reversed DNA sequence
			$_=~ tr/ACGTacgt/TGCAtgca/;
			$revcomp.="$_\n";
		}
	}
	return $revcomp;
}

# function to merge splitted block groups (bg1, bg2)
sub mergeBG {
	my ($bgId, $bgFile, $title) = @_;

	#$bgId="cluster_272232(+),cluster_272233(+),cluster_272234(+)";
	#$file="ERR015544testes.map.clusters.flagged.sig";

	my %bgHeader=(); my @bgReads=(); my @mergedBG=();
	tie %bgHeader, 'Tie::IxHash';
	%bgHeader = ( 'id' => "", 'chr' => "", 'start' => 1000000000000000, 'end' => 0, 'strand' => "",
                'reads' => 0, 'tags' => 0, 'blocks' => 0, 'name' => "n/a", 'anno' => "n/a", 'loci' => ".");

	my $x=0; my $l=(); my $fileFormat=();
	foreach my $id(sort { $a cmp $b } split(/\,/, $bgId)) {
		$fileFormat=$bgFile; $fileFormat=~s/^.+\///g; $fileFormat=~s/\..+//g;
		$id=~s/\(.+//g;
		#print "zgrep -w $id $bgFile | sed 's/^\>//g'";
		my @header = split(/\s+/, `zgrep -w $id $bgFile | sed 's/^\>//g'`);
		if(scalar(@header)==0) { print STDERR "$id not found in $bgFile\n"; exit(1); }
		$bgHeader{'id'} .= $header[0]."_";
		$bgHeader{'chr'} = $header[1];
		$bgHeader{'start'} = $header[2] if($header[2]<$bgHeader{'start'});
		$bgHeader{'end'} = $header[3] if($header[3]>$bgHeader{'end'});
		$bgHeader{'strand'} = $header[4];
		$bgHeader{'reads'} += $header[5];
		$bgHeader{'tags'} += $header[6];
		$bgHeader{'blocks'} += $header[7];
		$bgHeader{'name'} = $header[8] if($header[8]!~/n\/a/);
		$bgHeader{'anno'} = $header[9] if($header[9]!~/n\/a/);
		$bgHeader{'loci'} = $header[10] if($header[10]!~/\./);

		my @reads = `zless $bgFile | grep -w $id -A $header[6] | grep -v "^>"`;
		foreach(@reads) {
			my @t=split(/\s+/, $_);
			$t[6]=$t[6]+$x;
			$l = "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\n";
			push(@bgReads, $l);
		}
		$x = $x + $header[7];
	}
	$bgHeader{'id'}=~s/\_$//g;
	if(defined($title)) { $bgHeader{'id'}=">$fileFormat"."_"."$bgHeader{'id'}"."_"."$title"; }
	else {
		$bgHeader{'id'}=">$fileFormat"."_"."$bgHeader{'id'}";
	}
	$l=(); foreach(keys(%bgHeader)) { $l .= "$bgHeader{$_}\t"; } $l=~s/\t$//g;
	#print "$l\n";
	push(@mergedBG, "$l\n");
	foreach(@bgReads) {
		#print "$_";
		push(@mergedBG, "$_");
	}
	return @mergedBG;
}

# function to retrieve chr, start, end, name, score, strand from BED file
sub parseBED {
	my($l, $BED) = @_;

	my @t=split(/\s+/, $l);

	# differential processing analysis output
	if($t[0]=~/^bg[P|N]+/) {
		($BED->[0], $BED->[1], $BED->[2])=split(/[\:\-]+/,$t[1]);
		$BED->[3]=$t[24];
		$BED->[4]=0;
		if($t[0]=~/^bgP/) { $BED->[5]="+"; }
		else { $BED->[5]="-"; }
	}
	# BED format is "chr start end name score strand"
	elsif(defined($t[5]) && $t[5]=~/[+|-]{1}/) {
		@{$BED}=@t;
	}
	# BED format is "chr start end strand"
	elsif(defined($t[3]) && $t[3]=~/[+|-]{1}/) {
		$BED->[0]=$t[0];
		$BED->[1]=$t[1];
		$BED->[2]=$t[2];
		$BED->[3]="NA";
		$BED->[4]=0;
		$BED->[5]=$t[3];
	}
	## BED format is "chr strat end"
	elsif(scalar(@t)==3 && $t[1]=~/^[0-9]+$/ && $t[2]=~/^[0-9]+$/) {
		$BED->[0]=$t[0];
		$BED->[1]=$t[1];
		$BED->[2]=$t[2];
	}
	else { print STDERR "Error: Incorrect BED format at $l\n"; exit; }
}

# function to open file (including compressed files)
sub openFile {
	my($file)=@_;

	my @data=();
	
	if(defined($file) && $file=~/\.gz$/) {
		open(INFILE, "gunzip -c $file |") || die $!;
		@data=<INFILE>;
		close INFILE;
	}
	elsif(defined($file)) {
		open(INFILE, $file) || die $!;
		@data=<INFILE>;
		close INFILE;
	}
	else {
		my $INFILE = *STDIN;
		@data = <$INFILE>;
	}
	return @data;
}

## function to grep block group(s) from a file
sub grepBG {
	my($id, $file)=@_;

	my @data=openFile($file);

	## grep index corresponding to id
	my @headerIndex = grep { $data[$_]=~/^\>$id\s+/ } 0..$#data;

	## retrieve block group corresponding to the index
	my @bg=();
	foreach my $index(@headerIndex) {
		my @t = split(/\s+/, $data[$index]);
		foreach(@data[$index..($index+$t[6])]) {
			push(@bg, $_);
		}
	}
	return @bg;
}

## function to retrieve coordinates corresponding to non-overlapping set of blocks within block groups at a loci
sub uniqCoor {
	my($C, $bgInfo, $threshold, $print) = @_;

	my %coor=();
	foreach my $bg(keys(%{$bgInfo})) {
		foreach my $b(grep(/^[0-9]+$/, keys(%{$bgInfo->{$bg}}))) {
			my $key="$bgInfo->{$bg}{$b}{'start'}_$bgInfo->{$bg}{$b}{'end'}";
			$coor{$key}{'start'} = $bgInfo->{$bg}{$b}{'start'};
			$coor{$key}{'end'} = $bgInfo->{$bg}{$b}{'end'};
			$coor{$key}{'bLen'}=($bgInfo->{$bg}{$b}{'end'}-$bgInfo->{$bg}{$b}{'start'})+1;
			
		}
	}

	if($print) {
		foreach my $key(sort { $coor{$a}{'start'} <=> $coor{$b}{'start'} } keys(%coor)) { print "$key\t$coor{$key}{'start'}\t$coor{$key}{'end'}\n"; }
	}

	my $i=0; %{$C}=();
	foreach my $key(sort { ($coor{$a}{'start'} <=> $coor{$b}{'start'}) || ($coor{$a}{'bLen'} <=> $coor{$b}{'bLen'}) } keys(%coor)) {
		if(!%{$C}) {
			$C->{$i}{'start'}=$coor{$key}{'start'};
			$C->{$i}{'end'}=$coor{$key}{'end'};
			$C->{$i}{'bLen'}=$coor{$key}{'bLen'};
			$i++;
		}
		else {
			my $maxPerOverlap=0; my $maxPerOverlapIndex=();
			for(my $j=0; $j<$i; $j++) {
				my $overlap = checkOverlap($coor{$key}{'start'}, $coor{$key}{'end'}, $C->{$j}{'start'}, $C->{$j}{'end'}, 1);
				if($overlap > $maxPerOverlap) {
					$maxPerOverlap=$overlap;
					$maxPerOverlapIndex=$j;
				}
			}

			print "\t\t\t$coor{$key}{'start'}\t$coor{$key}{'end'}\t$maxPerOverlap\n" if($print);
			
			if($maxPerOverlap < $threshold) {
				$C->{$i}{'start'}=$coor{$key}{'start'};
				$C->{$i}{'end'}=$coor{$key}{'end'};
				$C->{$i}{'bLen'}=$coor{$key}{'bLen'};
				$i++;
			}
			else {
				if($coor{$key}{'start'}<$C->{$maxPerOverlapIndex}{'start'}) { $C->{$maxPerOverlapIndex}{'start'}=$coor{$key}{'start'}; }
				if($coor{$key}{'end'}>$C->{$maxPerOverlapIndex}{'end'}) { $C->{$maxPerOverlapIndex}{'end'}=$coor{$key}{'end'}; }
			}
		}
	}
	if($print) {
		foreach my $key(sort { $C->{$a}{'start'} <=> $C->{$b}{'start'} } keys(%{$C})) {print "$key\t$C->{$key}{'start'}\t$C->{$key}{'end'}\n";}
	}
}

# function to add a dummy read
sub addRead {
	my ($bgInfo, $bg, $b, $start, $end, $readCount) = @_; 

	my $file = $bgInfo->{$bg}{'id'};
	$file=~s/\_.+//g; $file=~s/\>//g;
	my $readId="$file"."_dummy_1";

	$bgInfo->{$bg}{'reads'}++;
	$bgInfo->{$bg}{'tags'}++;

	## update start and end coordinate of block group
	if($start < $bgInfo->{$bg}{'start'}) { $bgInfo->{$bg}{'start'}=$start; }
	if($end > $bgInfo->{$bg}{'end'}) { $bgInfo->{$bg}{'end'}=$end; }

	if($start<$bgInfo->{$bg}{$b}{'start'}) { $bgInfo->{$bg}{$b}{'start'}=$start; }
	if($end>$bgInfo->{$bg}{$b}{'end'}) { $bgInfo->{$bg}{$b}{'end'}=$end; }   

	$bgInfo->{$bg}{$b}{'reads'}+=$readCount;
	push(@{$bgInfo->{$bg}{$b}{'info'}}, "$bgInfo->{$bg}{'chr'}\t$start\t$end\t$readId\t$readCount\t$bgInfo->{$bg}{'strand'}\t$b");

	$bgInfo->{$bg}{'origReads'}=$bgInfo->{$bg}{'reads'};
}

# function to remove block from a block group
sub removeBlock {
	my ($bgInfo, $bg, $block) = @_;

	#print "$block\n";
	#printbgInfo(\%{$bgInfo}, $bg);
	$bgInfo->{$bg}{'blocks'}--;
	$bgInfo->{$bg}{'reads'}=$bgInfo->{$bg}{'reads'}-$bgInfo->{$bg}{$block}{'reads'};
	$bgInfo->{$bg}{'tags'}=$bgInfo->{$bg}{'tags'}-scalar(@{$bgInfo->{$bg}{$block}{'info'}});
	delete $bgInfo->{$bg}{$block};

	my $i=0;
	foreach my $b(sort { $a <=> $b } grep(/^[0-9]+$/, keys(%{$bgInfo->{$bg}}))) {
		if($b>$block) {
			my @info=();
			foreach(@{$bgInfo->{$bg}{$b}{'info'}}) {
				chomp($_);
				my @t=split(/\s+/, $_);
				$t[6]--;
				push(@info, "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]");
			}
			@{$bgInfo->{$bg}{$b}{'info'}}=@info;
		}

		if($i==0) {
			$bgInfo->{$bg}{'start'}=$bgInfo->{$bg}{$b}{'start'};
			$bgInfo->{$bg}{'end'}=$bgInfo->{$bg}{$b}{'end'};
		}
		else {
			if($bgInfo->{$bg}{$b}{'start'}<$bgInfo->{$bg}{'start'}) {
				$bgInfo->{$bg}{'start'}=$bgInfo->{$bg}{$b}{'start'};
			}
			if($bgInfo->{$bg}{$b}{'end'}>$bgInfo->{$bg}{'end'}) {
				$bgInfo->{$bg}{'end'}=$bgInfo->{$bg}{$b}{'end'};
			}
		}
		$i++;
	}

	if(grep(/^[0-9]+$/, keys(%{$bgInfo->{$bg}}))==0) {
		$bgInfo->{$bg}{'blocks'}=0;
		$bgInfo->{$bg}{'reads'}=0;
		$bgInfo->{$bg}{'tags'}=0;
		$bgInfo->{$bg}{'start'}=100000000000000000000000;
		$bgInfo->{$bg}{'end'}=0;
	}
	#foreach(keys(%{$bgInfo->{$bg}})) { print "$_," if($_=~/^[0-9]+$/); } print "\n";
	#printbgInfo(\%{$bgInfo}, $bg);
}

# function to add block to a block group
sub addBlock {
	my ($bgInfo, $bg, $start, $end, $readCount) = @_;

	$bgInfo->{$bg}{'blocks'}++;
	$bgInfo->{$bg}{'reads'}+=$readCount;
	$bgInfo->{$bg}{'tags'}++;

	if($start < $bgInfo->{$bg}{'start'}) { $bgInfo->{$bg}{'start'}=$start; }
	if($end > $bgInfo->{$bg}{'end'}) { $bgInfo->{$bg}{'end'}=$end; }

	my $readId=$bgInfo->{$bg}{'id'};
	$readId=~s/\>//; $readId=~s/\_.+//g; $readId=$readId."_dummy";

	my $b=$bgInfo->{$bg}{'blocks'};
	$bgInfo->{$bg}{$b}{'start'}=$start;
	$bgInfo->{$bg}{$b}{'end'}=$end;
	$bgInfo->{$bg}{$b}{'reads'}=$readCount;
	$bgInfo->{$bg}{$b}{'info'}->[0]="$bgInfo->{$bg}{'chr'}\t$start\t$end\t$readId\t$readCount\t$bgInfo->{$bg}{'strand'}\t$b";
}

# function to add block to a block group, if read exists in the map file
sub mapBlock {
	my ($bgInfo, $bg, $start, $end, $mapDir, $mapFormat) = @_;

	my $file = $bgInfo->{$bg}{'id'};
	$file=~s/\_.+//g; $file=~s/\>//g;

	my @reads=();
	if($bgInfo->{$bg}{'strand'}=~/^\+$/) {
		@reads=`coor2mapping.pl -c $bgInfo->{$bg}{'chr'}:$start-$end -f $mapDir/$file.*.sorted.$bgInfo->{$bg}{'chr'}.P.gz -s $bgInfo->{$bg}{'strand'} -r $mapFormat`;
	}
	else {
		@reads=`coor2mapping.pl -c $bgInfo->{$bg}{'chr'}:$start-$end -f $mapDir/$file.*.sorted.$bgInfo->{$bg}{'chr'}.N.gz -s $bgInfo->{$bg}{'strand'} -r $mapFormat`;
	}
	#@reads=`coor2mapping.pl -c $bgInfo->{$bg}{'chr'}:$start-$end -f $mapDir/$file.*.sorted.gz -s $bgInfo->{$bg}{'strand'} -r $mapFormat`;

	if(scalar(@reads)>0) { $bgInfo->{$bg}{'blocks'}++; }
	foreach my $l(@reads) {
		my @t=split(/\s+/, $l);
		my($readId, $height) = split(/\|/, $t[1]);
		$readId=$readId."_".$t[14];
		my $readCount=$height/$t[14];
		$bgInfo->{$bg}{'reads'}+=$readCount;
		$bgInfo->{$bg}{'tags'}++;

		if($t[10] < $bgInfo->{$bg}{'start'}) { $bgInfo->{$bg}{'start'}=$t[10]; }
		if($t[11] > $bgInfo->{$bg}{'end'}) { $bgInfo->{$bg}{'end'}=$t[11]; }

		my $b=$bgInfo->{$bg}{'blocks'};
		if(defined($bgInfo->{$bg}{$b}{'start'}) && $t[10]<$bgInfo->{$bg}{$b}{'start'}) {
			$bgInfo->{$bg}{$b}{'start'}=$t[10];
		}
		elsif(!defined($bgInfo->{$bg}{$b}{'start'})) { $bgInfo->{$bg}{$b}{'start'}=$t[10]; }

		if(defined($bgInfo->{$bg}{$b}{'end'}) && $t[11]>$bgInfo->{$bg}{$b}{'end'}) {
			$bgInfo->{$bg}{$b}{'end'}=$t[11];
		}
		elsif(!defined($bgInfo->{$bg}{$b}{'end'})) { $bgInfo->{$bg}{$b}{'end'}=$t[11]; }
		$bgInfo->{$bg}{$b}{'reads'}+=$readCount;
		push(@{$bgInfo->{$bg}{$b}{'info'}}, "$bgInfo->{$bg}{'chr'}\t$t[10]\t$t[11]\t$readId\t$readCount\t$t[9]\t$b");
	}

	$bgInfo->{$bg}{'origReads'}=$bgInfo->{$bg}{'reads'};
	if(scalar(@reads)>0) { return 1; }
	else { return 0; }
}

## function to determine unique elements in an array
sub uniqElements {
	my($arr) = @_;

	my %seen=();
	foreach(@{$arr}) { $seen{$_}=1; }

	return scalar(keys(%seen));
}

# function to reformat block ids in a block group
sub reformatbId {
	my($bgInfo)=@_;

	my %temp=();
	foreach my $bg(keys(%{$bgInfo})) {
		my $i=1;
		foreach my $b(sort { $a <=> $b } grep(/^[0-9]+/, keys(%{$bgInfo->{$bg}}))) {
			$temp{$bg}{$i}{'start'}=$bgInfo->{$bg}{$b}{'start'};
			$temp{$bg}{$i}{'end'}=$bgInfo->{$bg}{$b}{'end'};
			$temp{$bg}{$i}{'reads'}=$bgInfo->{$bg}{$b}{'reads'};
			@{$temp{$bg}{$i}{'info'}}=@{$bgInfo->{$bg}{$b}{'info'}};
			delete $bgInfo->{$bg}{$b};
			$i++;
		}
	}

	foreach my $bg(keys(%{$bgInfo})) {
		foreach $b(keys(%{$temp{$bg}})) {
			$bgInfo->{$bg}{$b}{'start'}=$temp{$bg}{$b}{'start'};
			$bgInfo->{$bg}{$b}{'end'}=$temp{$bg}{$b}{'end'};
			$bgInfo->{$bg}{$b}{'reads'}=$temp{$bg}{$b}{'reads'};
			@{$bgInfo->{$bg}{$b}{'info'}}=@{$temp{$bg}{$b}{'info'}};
		}
	}
}

# function to compute absolute and relative expression in a block group
sub compbgRelExp {
	my($bgInfo, $E, $fArg, $tArg)=@_;

	foreach my $bg(keys(%{$bgInfo})) {
		$E->{$bg}{'reads'}=$bgInfo->{$bg}{'reads'};
		$E->{$bg}{'minRelExp'}=1;
		$E->{$bg}{'minBlockReads'}=$bgInfo->{$bg}{1}{'reads'};
		foreach my $b(keys(%{$bgInfo->{$bg}})) {
			if($b=~/^[0-9]+$/) {
				if(!defined($bgInfo->{$bg}{$b}{'reads'}) || !defined($bgInfo->{$bg}{'reads'})) {
					print STDERR "error: no reads in block ($b) pr block group ($bg)\n";
					exit;
				}
				$E->{$bg}{$b}{'relExp'}=sprintf("%0.2f", $bgInfo->{$bg}{$b}{'reads'}/$bgInfo->{$bg}{'reads'});
				if($E->{$bg}{$b}{'relExp'} < $E->{$bg}{'minRelExp'}) {
					$E->{$bg}{'minRelExp'}=sprintf("%0.2f", $E->{$bg}{$b}{'relExp'});
					$E->{$bg}{'minBlockReads'}=$bgInfo->{$bg}{$b}{'reads'};
				}
			}
		}
	}

	foreach my $bg(keys(%{$E})) {
		${$fArg}.="$E->{$bg}{'reads'},$E->{$bg}{'minBlockReads'},";
		${$tArg}.="$E->{$bg}{'minRelExp'},";
	}
	${$fArg}=~s/\,$//; ${$tArg}=~s/\,$//;
}

# function to print block group
sub printbgInfo {
	my($bgInfo, $outDir, $outFile)=@_;

	my $out2file=0;
	if(defined($outDir) && defined($outFile)) {
		open(OUTFILE, ">$outDir/$outFile") || die $!;
		$out2file=1;
	}
	
	foreach my $bg(keys(%{$bgInfo})) {
		if($out2file==1) {
			print OUTFILE "$bgInfo->{$bg}{'id'}\t$bgInfo->{$bg}{'chr'}\t$bgInfo->{$bg}{'start'}\t$bgInfo->{$bg}{'end'}\t$bgInfo->{$bg}{'strand'}\t$bgInfo->{$bg}{'reads'}\t$bgInfo->{$bg}{'tags'}\t$bgInfo->{$bg}{'blocks'}\t$bgInfo->{$bg}{'name'}\t$bgInfo->{$bg}{'anno'}\t$bgInfo->{$bg}{'loci'}\n";
		}
		else {
			print "$bgInfo->{$bg}{'id'}\t$bgInfo->{$bg}{'chr'}\t$bgInfo->{$bg}{'start'}\t$bgInfo->{$bg}{'end'}\t$bgInfo->{$bg}{'strand'}\t$bgInfo->{$bg}{'reads'}\t$bgInfo->{$bg}{'tags'}\t$bgInfo->{$bg}{'blocks'}\t$bgInfo->{$bg}{'name'}\t$bgInfo->{$bg}{'anno'}\t$bgInfo->{$bg}{'loci'}\n";
		}
		foreach my $b( sort { $a <=> $b } grep(/^[0-9]+$/, keys(%{$bgInfo->{$bg}}))) {
			foreach(@{$bgInfo->{$bg}{$b}{'info'}}) {
				if($out2file==1) { print OUTFILE "$_\n"; }
				else { print "$_\n"; }
			}
		}
	}

	close OUTFILE if($out2file==1);
}

## function to check format of map file
sub checkMapFormat {
	my(@data)=@_;

	foreach my $l(@data) {
		next if($l=~/^\#/);
		my @t=split(/\s+/, $l);
		if($t[0]=~/^Q$/ && $t[9]=~/^[\+\-]+$/ && $t[10]=~/^[0-9]+$/ && $t[11]=~/^[0-9]+$/) {
			return "segemehl";
		}
		elsif($t[1]=~/^[0-9]+$/ && $t[2]=~/^[0-9]+$/ && $t[5]=~/^[\+\-]+$/) {
			return "bed";
		}
		else {
			return "unknown format";
		}
	}
}

## function to parse 'pvclust' cluster file
sub parsePvclust {
	my(@data)=@_;

	my %clusterInfo=(); my $cluster=();
	foreach my $l(@data) {
		last if($l=~/^\$edges/);
		next if($l=~/^$/);
		if($l=~/^\$clusters\[\[[0-9]+\]\]$/) {
			$l=~m/[0-9]+/;
			$cluster=$&;
			tie %{$clusterInfo{$cluster}}, 'Tie::IxHash';
    }
		elsif($l=~/\s*\[[0-9]+\]\s+\"{0,1}/ && defined($cluster)) {
			chomp($l);
			my $description=$l;
			$l=~s/\s*\[[0-9]+\]\s+\"{0,1}//g; $l=~s/[\"\>]+//g; $l=~s/\s+$//g;
			my @t=split(/[\.\:]+/,$l);
			if(scalar(@t)==7) {
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'name'}=$t[1];
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'type'}=$t[2];
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'loci'}=$t[3];
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'coor'}="$t[4]:$t[5]";
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'strand'}=$t[6];
				my @len=split(/\-/,$t[5]);
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'length'}=($len[1]-$len[0]);
			}
			elsif(scalar(@t)==6) {
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'name'}=$t[0];
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'type'}=$t[1];
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'loci'}=$t[2];
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'coor'}="$t[3]:$t[4]";
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'strand'}=$t[5];
				my @len=split(/\-/,$t[4]);
				$clusterInfo{$cluster}{"$t[0].$t[1]"}{'length'}=($len[1]-$len[0]);
			}
			else {
				print STDERR "Error: input pvclust file is not in correct format at $l\n";
				exit(-1);
			}
			$clusterInfo{$cluster}{"$t[0].$t[1]"}{'description'}=$description;
			#print "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\n";
    }
  }
	return %clusterInfo;
}

## determine overlapping coordinates by binary sort (coordinate file should be sorted)
sub binarySearchCoor {
	my($start, $end, $mid, $posmin, $posmax, $mapFormat, $line)=@_;
	my $result=();
	my @t=split(/\s+/, $line);
	if($mapFormat=~/segemehl/) {
		if($t[10]>$end) { ${$posmax}=$mid; $result=0; }
		elsif($t[11]<$start) { ${$posmin}=$mid; $result=0; }
		else { $result=1; }
	}
	elsif($mapFormat=~/bed/) {
		if($t[1]>$end) { ${$posmax}=$mid; $result=0; }
		elsif($t[2]<$start) { ${$posmin}=$mid; $result=0; }
		else { $result=1; }
	}

	if((${$posmax}-${$posmin})==1) { $result=1; }
	return $result;
}

## determine the common elements between two arrays
sub checkArrayOverlap {
	my($fArray, $sArray)=@_;

	my @common_array=();

	## sort the array elements
	my @sort_fArray = sort { $a cmp $b } @{$fArray};
	my @sort_sArray = sort { $a cmp $b } @{$sArray};

	## determine common elements
	foreach my $f(@sort_fArray) {
		foreach my $s(@sort_sArray) {
			push(@common_array, $f) if($f=~/$s/);
		}
	}

	## remove duplicates in the common_array
	my %uniq=();
	foreach(@common_array) {
		$uniq{$_}=1;
	}

	## print common elements
	#foreach(keys(%uniq)) { print "$_\n"; }
	
	## return number of common elements
	return keys(%uniq);
}

## determine the common keys shared between two or more hashes
## usage: commonHashKeys(\@allKeys, scalar(keys(%hash)))
sub commonHashKeys {
	my($keys, $count) = @_;

	## read all keys
	my %seen=();
	foreach(@{$keys}) {
		$seen{$_}++;
	}

	## determine common keys
	my @common_keys=();
	foreach(keys(%seen)) {
		push(@common_keys, $_) if($seen{$_}==$count);
	}

	## return common keys
	return @common_keys;
}

## parse deepBlockAlign output
sub parseDBAoutput {
	my($l) = @_;
	my %bgInfo=();

	$l=~s/^\>//g;
	my @t=split(/\|/, $l);
	$bgInfo{'id'}=$t[0];
	$bgInfo{'name'}=$t[1];
	$bgInfo{'type'}=$t[2];
	$bgInfo{'class'}=$t[3];
	($bgInfo{'chr'}, $bgInfo{'start'}, $bgInfo{'end'}, $bgInfo{'strand'})=split(/[\:\-\(]+/, $t[scalar(@t)-1]);
	$bgInfo{'strand'}=~s/[\(\)]+//g;
	
	return %bgInfo;
}

# function to check for reverse complement expression at a loci
sub checkRevCompExp {
	my ($l) = @_;
	my @lF=(); my $identifier=();
    if($l=~/cluster\_/) { $identifier="cluster"; }
    elsif($l=~/BG\_/) { $identifier="BG"; }

	if($l=~/\(\+\)/ && $l=~/\(\-\)/) {
		my @bgId = split(/\s+/, $l);
		my $bgPos = "$bgId[0]\t"; my $bgNeg = "$bgId[0]\t";
		shift(@bgId);
		foreach(@bgId) {
			if(grep(/$identifier\_[0-9]+\(\+\)/, $_)) {
				while($_=~m/$identifier\_[0-9]+\(\+\)/g) {
					$bgPos .= "$&,";
				}
				$bgPos=~s/\,$/\t/;
			}
			else { $bgPos .= "-\t"; }
			if(grep(/$identifier\_[0-9]+\(\-\)/, $_)) {
				while($_=~m/$identifier\_[0-9]+\(\-\)/g) {
					$bgNeg .= "$&,";
				}
				$bgNeg=~s/\,$/\t/;
			}
			else { $bgNeg .= "-\t"; }
		}
		push(@lF, $bgPos);
		push(@lF, $bgNeg);
		return (1, \@lF);
	}
	else { push(@lF, $l); return (0, \@lF); }
}

## combine two PDF images as single PDF file using latex
sub combinePDF {
	my($pdf1, $pdf2, $pdfOut, $orientation, %info)=@_;
	open(OUTFILE, ">$pdfOut.tex") || die $!;

print OUTFILE <<LATEX;
\\documentclass{standalone}
\\usepackage[hmargin=0.5cm, vmargin=1cm]{geometry}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{tikz}
\\usepackage{subfigure}
\\usepackage{adjustbox}
\\usepackage{multirow}
\\usepackage{rotating}
\\definecolor{dark-red}{RGB}{100,0,0}
\\definecolor{dark-green}{RGB}{0,100,0}
\\definecolor{dark-blue}{RGB}{0,0,100}
\\begin{document}
LATEX

	## orientation (1: one above another; 2: one next to other)
	if(defined($orientation) && $orientation==1) {
		print OUTFILE "\\begin{adjustbox}{width=30 cm,height=6 cm,keepaspectratio}\n";
		print OUTFILE "\\begin{tabular}{cc}\n";
		print OUTFILE "{\\begin{sideways}\\hspace{1.5cm}\\textbf{FIRST SAMPLE}\\end{sideways}} & \\includegraphics[type=pdf, ext=.pdf, read=.pdf, scale=1]{$pdf1} \\\\\n";
		print OUTFILE "{\\begin{sideways}\\hspace{1.5cm}\\textbf{SECOND SAMPLE}\\end{sideways}} & \\includegraphics[type=pdf, ext=.pdf, read=.pdf, scale=1]{$pdf2} \\\\\n";
		print OUTFILE "\\end{tabular}\n";
		print OUTFILE "\\end{adjustbox}\n";
	}
	else {
		$pdf1=~s/\.pdf//g; $pdf2=~s/\.pdf//g;
		print OUTFILE "{\\Large\\begin{tabular}{|l|l|l} \\hline\n";
		print OUTFILE " & Replicate 1 & Replicate 2 \\\\\\hline\n";
		print OUTFILE " & \\includegraphics[type=pdf,ext=.pdf,read=.pdf, scale=1]{$pdf1} & \\includegraphics[type=pdf,ext=.pdf,read=.pdf, scale=1]{$pdf2} \\\\\\hline\n";
		if(defined($info{'rep1'}{'id'}) && $info{'rep2'}{'id'}) {
			$info{'rep1'}{'id'}=~s/\_/\\_/g; $info{'rep2'}{'id'}=~s/\_/\\_/g;
			print OUTFILE "Identifier & $info{'rep1'}{'id'} & $info{'rep2'}{'id'} \\\\\\hline\n";
		}
		if(defined($info{'rep1'}{'score'}) && defined($info{'rep2'}{'score'})) {
			print OUTFILE "Cluster score & $info{'rep1'}{'score'} & $info{'rep2'}{'score'} \\\\\\hline\n";
		}
		if(defined($info{'rep1'}{'clusterInfo'}) && defined($info{'rep2'}{'clusterInfo'})) {
			print OUTFILE "Clusters & $info{'rep1'}{'clusterInfo'} & $info{'rep2'}{'clusterInfo'} \\\\\\hline\n";
		}
		if(defined($info{'overlap'})) {
			print OUTFILE "Overlap & $info{'overlap'} & \\\\\\hline\n";
		}
		print OUTFILE "\\end{tabular}}\n";
	}

print OUTFILE <<LATEX;
\\end{document}
LATEX

	close OUTFILE;
	system("pdflatex $pdfOut.tex > /dev/null");
	system("rm $pdfOut.aux");
	system("rm $pdfOut.log");
	system("rm $pdfOut.tex");
	system("rm $pdf1.pdf");
	system("rm $pdf2.pdf");
}

## combine multiple PDF images as single PDF file using latex
sub combineMultiPDF {
	my($pdfIn, $pdfOut, $title)=@_;
	open(OUTFILE, ">$pdfOut.tex") || die $!;

	my $WIDTH=30; my $HEIGHT=scalar(@{$pdfIn})*3;

print OUTFILE <<LATEX;
\\documentclass{standalone}
%\\usepackage[hmargin=0.5cm, vmargin=1cm]{geometry}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{tikz}
\\usepackage{subfigure}
\\usepackage{adjustbox}
\\usepackage{multirow}
\\usepackage{rotating}
\\definecolor{dark-red}{RGB}{100,0,0}
\\definecolor{dark-green}{RGB}{0,100,0}
\\definecolor{dark-blue}{RGB}{0,0,100}
\\begin{document}
\\begin{adjustbox}{width=$WIDTH cm,height=$HEIGHT cm,keepaspectratio}
LATEX

	print OUTFILE "\\begin{tabular}{c}\n";
foreach(@{$pdfIn}) {
	print OUTFILE "\\includegraphics[scale=1]{$_} \\\\\n";
}
	if(defined($title) && $title!~/^$/) {
		$title=~s/\_/\\_/g;
		print OUTFILE "$title \\\\\n";
	}
	print OUTFILE "\\end{tabular}\n";

print OUTFILE <<LATEX;
\\end{adjustbox}
\\end{document}
LATEX

	close OUTFILE;

	system("pdflatex $pdfOut.tex > /dev/null");
	system("rm $pdfOut.aux");
	system("rm $pdfOut.log");
	system("rm $pdfOut.tex");
}

## retrieve non-coding RNA annotation overlapping to the input coordinate (sorted by start coordinate in BED format)
sub coor2annotation {
	my %coor=();
	my $annotation=();
	my $minPerOverlap=();
	($coor{'coor'}, $coor{'strand'}, $annotation, $minPerOverlap)=@_;

	($coor{'chr'}, $coor{'start'}, $coor{'end'}) = split(/[\:\-]+/, $coor{'coor'});

	if(!defined($coor{'chr'}) || !defined($coor{'start'}) || !defined($coor{'end'})) { print STDERR "Incorrect coordinate\n"; return "None"; }

	if(scalar(@{$annotation})==0) {
		return "None";
	}

	my $posmin=0;
	#my $posmax=$#annotation;
	my $posmax=scalar(@{$annotation})-1;
	my $result=0;
	while($result==0) {
		my $mid = int (($posmin + $posmax) /2);
		#print "$annotation->[$posmin]\n$annotation->[$mid]\n$annotation->[$posmax]\n\n";
		$result = binarySearchCoor($coor{'start'}, $coor{'end'}, $mid, \$posmin, \$posmax, "bed", $annotation->[$mid]);
		$result=1 if($posmin==$posmax);
	}

	$coor{'annoOverlap'}=$minPerOverlap;
	$coor{'annoLength'}=10000000000000000000000000000000000000000000000;

	#foreach my $l(@{$annotation[$posmin..$posmax]}) {
	for(my $i=$posmin; $i<=$posmax; $i++) {
		my @anno=split(/\s+/, $annotation->[$i]);
		my $overlap=(); my $length=();

		if(defined($coor{'strand'})) {
			$overlap=checkOverlap($anno[1], $anno[2], $coor{'start'}, $coor{'end'}, 1, $anno[0], $coor{'chr'}, $anno[5], $coor{'strand'});
			#print "$anno[1]\t$anno[2]\t$coor{'start'}\t$coor{'end'}\t$anno[0]\t$coor{'chr'}\t$anno[5]\t$coor{'strand'}\t$anno[3]\t$overlap\t";
		}
		else {
			$overlap=checkOverlap($anno[1], $anno[2], $coor{'start'}, $coor{'end'}, 1, $anno[0], $coor{'chr'});
		}
		$length=($anno[2]-$anno[1])+1;
		#print "$length\n";

		if($overlap >= $coor{'annoOverlap'} && $length < $coor{'annoLength'}) {
			$coor{'anno'}=$annotation->[$i];
			$coor{'annoOverlap'}=$overlap;
			## determine the length of the overlapping annotation
			$coor{'annoLength'}=($anno[2]-$anno[1])+1;
			#print "$coor{'anno'}\n";
		}
	}

	## print the annotation overlapping to the input coordinate
	if(defined($coor{'anno'})) {
		return "$coor{'anno'}";
	}
	else {
		return "None";
	}
}

1;
