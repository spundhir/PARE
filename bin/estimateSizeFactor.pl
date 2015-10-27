#!/usr/bin/perl -w

=copyright_info
estimateSizeFactor.pl: estimate size factor for normalization across samples
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
use Tie::IxHash;
use Getopt::Long;
use perlModule;

#############################################################################################################
## parse input options
use vars qw($option $bgDir $countFile $mapDir $bamFile $coorFile $replicate $uniqReadsFilter $chrFilter $fileFilter $progDir $extend $genome $programToUse $countMultiMap $help); 
$uniqReadsFilter=0;
$chrFilter=".+";
$fileFilter="flagged";
$replicate="*";
$genome="mm9";
$programToUse=2;
$countMultiMap=0;

GetOptions ("o=s"  => \$option,
            "d=s"  => \$bgDir,
            "r=s"  => \$countFile,
            "m=s"  => \$mapDir,
            "b=s"  => \$bamFile,
            "x=s"  => \$coorFile,
            "t=s"  => \$replicate,
            "u=i"  => \$uniqReadsFilter,
            "c=s"  => \$chrFilter,
            "f=s"  => \$fileFilter,
            "p=s"  => \$progDir,
            "e=s"  => \$extend,
            "g=s"  => \$genome,
            "s=s"  => \$programToUse,
            "l=s"  => \$countMultiMap,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$option || !$countFile);

#############################################################################################################
sub usage {
	print STDERR "\nProgram: estimateSizeFactor.pl (estimate size factor for normalization across samples)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: estimtateSizeFactor.pl -o <option> [OPTIONS]\n";
	print STDERR "          -o <computation option>\n";
	print STDERR "          -h <help>\n";
	print STDERR "Option a: compute read count per block, block group or gene coordinate\n";
	print STDERR "          -d <dir>    [input directory with block groups in blockbuster format]\n";
	print STDERR "          -r <file>   [output file with read count at each block or block groups across all samples]\n";
	print STDERR "          -m <dir>    [directory with indexed map files (optional)]\n";
	print STDERR "                      [if not given: use read count corresponding to block groups]\n";
	print STDERR "                      [if given: use acutal number of mapped reads]\n";
	print STDERR "          -x <file>   [input BED file with coordinates corresponding to which read count to compute (optional)]\n";
    print STDERR "Option b: compute read count per BED file region using HTSeq-count\n";
    print STDERR "          -b <file>   [input bam file]\n";
    print STDERR "                      [if multiple, seperate them by a comma]\n";
	print STDERR "          -x <file>   [input BED file with coordinates corresponding to which read count to compute]\n";
	print STDERR "          -r <file>   [output file with read count at each block or block groups across all samples]\n";
    print STDERR "          -e <int>    [extend 3' end of reads by input number of bases (default: 0; useful for ChIP-seq)]\n";
    print STDERR "                      [if multiple, seperate them by a comma]\n";
    print STDERR "          -g <string> [genome (default: mm9)]\n";
    print STDERR "          -s <int>    [use program 1:htseq-count or 2:featureCounts (default: 2)]\n";
    print STDERR "          -l <int>    [count multiple mapping reads, only applicable when -s is 2 (default: 0)]\n";
	print STDERR "Option c: estimate size factor using read count file computed at step a or b\n";
	print STDERR "          -r <file>   [input file with read count at each block or block groups across all samples]\n";
	print STDERR "          -t <string> [compute size factor for given replicate (Rep1, Rep2 etc)]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR "          -p <dir>    [directory having dependent programs]\n";
	print STDERR "          -u <int>    [compute read count for block groups with at least given number of unique reads]\n";
	print STDERR "          -c <string> [compute read count for block groups of given chromosome only]\n";
	print STDERR "          -f <string> [consider only those block group files that match the given name (default: flagged)]\n\n";
	exit(-1);
}

#############################################################################################################

## populate genome file based on input genome
my $GENOME_FILE=`initialize_genome -i $ENV{PAREPATH}/data/annotations/GENOME_FILE -g $genome`;
$GENOME_FILE="$ENV{PAREPATH}/data/annotations/$GENOME_FILE";
if(! -f $GENOME_FILE) {
    print "\ncomputation for $genome is not feasible yet\n";
    print "please add the chromosome size file for $genome at $ENV{PAREPATH}/data/annotations\n";
    print "also update the $ENV{PAREPATH}/data/annotations/GENOME_FILE\n";
    usage();
}

if($option=~/[aA]+/) {
	usage if(!$bgDir);

    opendir(INDIR, $bgDir) || die $!;
    my @files=(); my %seen=();

    ## read block group information (coordinate, identifier, chr, strand) from each block group file
    foreach my $bgFile(sort { $a cmp $b } readdir(INDIR)) {
    	if($bgFile=~/.*$fileFilter.*/) {
				my $file=$bgFile; $file=~s/\..+$//g;
				if(!defined($seen{$file})) {
					#print "$bgFile\n";
    			push(@files, $bgFile);
					$seen{$file}=1;
				}
    	}
    }
    close INDIR;

    if(defined($coorFile) && defined($mapDir)) {
			## open count file for writing
			open(OUTFILE, ">$countFile") || die $!;

			for(my $i=0; $i<scalar(@files); $i++) {
				$files[$i]=~s/\..+$//g;
				print OUTFILE "\t$files[$i]";
 			} print OUTFILE "\n";

			## read unique chromosomes in the coordinate file
			my @chr=`zless $coorFile | cut -f 1 | grep -v "_" | sort | uniq`;
			my $count_field=`zless $coorFile | perl -an -F'/\\t+/' -e 'print scalar(\@F)."\\n";'| sort | uniq`;
			chomp($count_field); $count_field++;

			foreach my $chr(@chr) {
				my %readCount=(); my %coor=();
				chomp($chr);
				system("zgrep -E \"$chr\\s+.*\\s+\\+\\s+\" $coorFile > $coorFile.$chr.P");
				@{$coor{'P'}}=`zless $coorFile.$chr.P | perl -ane 'print "\$F[0]:\$F[1]-\$F[2]|\$F[5]\n";'`;
				system("zgrep -E \"$chr\\s+.*\\s+\\-\\s+\" $coorFile > $coorFile.$chr.N");
				@{$coor{'N'}}=`zless $coorFile.$chr.N | perl -ane 'print "\$F[0]:\$F[1]-\$F[2]|\$F[5]\n";'`;

				foreach my $file(@files) {
					$file=~s/\..+$//g;
					if(-e "$mapDir/$file.bed.sorted.$chr.P") {
						@{$readCount{'P'}{$file}}=`zless $coorFile.$chr.P | intersectBed -c -s -sorted -a stdin -b $mapDir/$file.bed.sorted.$chr.P | cut -f $count_field`;
					}
					else { foreach(@{$coor{'P'}}) { push(@{$readCount{'P'}{$file}}, 0); } }
					if(-e "$mapDir/$file.bed.sorted.$chr.N") {
						@{$readCount{'N'}{$file}}=`zless $coorFile.$chr.N | intersectBed -c -s -sorted -a stdin -b $mapDir/$file.bed.sorted.$chr.N | cut -f $count_field`;
					}
					else { foreach(@{$coor{'N'}}) { push(@{$readCount{'N'}{$file}}, 0); } }
				}

				## print results
				for(my $i=0; $i<scalar(@{$coor{'P'}}); $i++) {
					chomp($coor{'P'}[$i]);
					print OUTFILE "$coor{'P'}[$i]\t";
					foreach my $file(@files) {
						chomp($readCount{'P'}{$file}[$i]);
						print OUTFILE "$readCount{'P'}{$file}[$i]\t";
					} print OUTFILE "\n";
				}
				for(my $i=0; $i<scalar(@{$coor{'N'}}); $i++) {
					chomp($coor{'N'}[$i]);
					print OUTFILE "$coor{'N'}[$i]\t";
					foreach my $file(@files) {
						chomp($readCount{'N'}{$file}[$i]);
						print OUTFILE "$readCount{'N'}{$file}[$i]\t";
					} print OUTFILE "\n";
				}
				
				system("rm $coorFile.$chr.P");
				system("rm $coorFile.$chr.N");
			}
			close OUTFILE;
    }
    elsif(defined($mapDir)) {
    	## compute read count for each coordinate across all samples, if mapDir is provided
    	compReadCount($bgDir, \@files, $countFile, 1, $uniqReadsFilter, $chrFilter, $mapDir);
    }
    else {
    	## compute read count for each coordinate across all samples, if mapDir is not provided
    	compReadCount($bgDir, \@files, $countFile, 1, $uniqReadsFilter, $chrFilter);
    }
}
elsif($option=~/[bB]+/) {
	usage if(!$bamFile || !$coorFile);

    my @files=split(/\,/,$bamFile);
    if(!defined($extend)) {
        $extend="";
        foreach(@files) { $extend.="0,"; }
        $extend=~s/\,$//g;
    }

    my @extends=split(/\,/,$extend);

    usage() if(scalar(@files) != scalar(@extends));

    my $count=1;
    my $header=();

    unless(-e("$coorFile.gtf")) {
        system("zless $coorFile | perl -ane 'chomp(\$F[5]); if(\$F[0]=~/\\_/) { next; } \$count++; \$F[3]=~s/^.*\\\///g; print \"\$F[0]\\tregion\\tall\\t\$F[1]\\t\$F[2]\\t.\\t\$F[5]\\t.\\tgene_id \\\"\$F[3]\\\"\\n\";' > $coorFile.gtf");
    }
    foreach my $file(@files) {
        if(-e "$file") {
            if($programToUse==1) {
                ## uncomment to debug
                #print("htseq-count -f bam -s no -t all $file <\(zless $coorFile | perl -ane '\$count++; \$F[3]=~s/^.*\\\///g; print \"\$F[0]\\tregion\\tall\\t\$F[1]\\t\$F[2]\\t.\\t\$F[5]\\t.\\tgene_id \\\"\$F[3]\\\"\\n\";'\)\n");
                system("bedtools bamtobed -i $file | bedtools slop -i stdin -g $GENOME_FILE -s -l 0 -r $extends[$count-1] | bedtools bedtobam -i stdin -g $GENOME_FILE | samtools view - | htseq-count -f sam -s no -t all - $coorFile.gtf > $countFile.$count");
                #system("htseq-count -f bam -s no -t all $file $coorFile.gtf > $countFile.$count && rm $coorFile.gtf");
            }
            else {
                if($countMultiMap==0) {
                    system("featureCounts -t all -a $coorFile.gtf -o $countFile.$count -T 5 --readExtension3 $extends[$count-1] $file");
                } else {
                    system("featureCounts -t all -a $coorFile.gtf -o $countFile.$count -T 5 --readExtension3 $extends[$count-1] -M $file");
                }
            }
            $file=~s/^.*\///g;
            $header.="\\t$file";
            $count++;
        }
        else {
            print STDERR "Cannot open $file\n"; exit(-1);
        }
    }
    #system("rm $coorFile.gtf");
    system("echo -e \"$header\" > $countFile");
    if($programToUse==1) {
        system("paste $countFile.* | perl -ane 'next if(\$_=~/^\\_\\_/); \$line=(); \$line.=\"\$F[0]\\t\"; \$index=1; foreach(\@F) { \$line.=\"\$F[\$index]\\t\"; \$index+=2; } \$line=~s/\\t\$//g; print \"\$line\\n\";' >> $countFile");
    }
    else {
        system("rm -f $countFile.*.summary");
        system("paste $countFile.* | perl -ane 'next if(\$_=~/^[\\#|Geneid]+/); \$line=(); \$line.=\"\$F[0]\\t\"; \$index=6; foreach(\@F) { \$line.=\"\$F[\$index]\\t\"; \$index+=7; } \$line=~s/[\\t\\r]+\$//g; print \"\$line\\n\";' >> $countFile\n");
    }
    #system("rm -i $countFile.*");
}
elsif($option=~/[cC]+/) {
	## compute size factor for all samples
	#print "R --no-save --vanilla --slave < ~/software/myScripts/deseq.R --args $countFile 1 $replicate\n"; exit;
	my @result=();
	#print "R --no-save --vanilla --slave < deseq.R --args $countFile 1 \"$replicate\"\n";
	if(defined($progDir)) {
		@result=`R --no-save --vanilla --slave < $progDir/deseq.R --args $countFile 1 \"$replicate\" 2>/dev/null`;
	}
	else {
		@result=`R --no-save --vanilla --slave < ~/software/myScripts/RNAseq/deseq.R --args $countFile 1 \"$replicate\" 2>/dev/null`;
	}

	## retrieve size factor for each sample
	my %sizeFactor=(); my @key=();
	tie %sizeFactor, 'Tie::IxHash';
	foreach my $l(@result) {
		chomp($l);
		if($l=~/^\s*[A-Za-z]+/) {
			$l=~s/^\s+//g;
			@key=split(/\s+/, $l);
			for(my $i=0; $i<scalar(@key); $i++) {
				$key[$i]=~s/^\s*//g;
				$key[$i]=~s/\..+//g;		
			}
		}
		else {
			$l=~s/^\s+//g;
			my @t=split(/\s+/, $l);
			for(my $i=0; $i<scalar(@key); $i++) {
				$t[$i]=~s/\s+//g;
				$sizeFactor{$key[$i]}=sprintf("%0.2f", $t[$i]);
			}
		}
	}

	## print size factor for each sample
	foreach my $key(keys(%sizeFactor)) {
		print "$key\t$sizeFactor{$key}\n";
	}
}

exit(0);
