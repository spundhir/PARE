#!/usr/bin/env perl

=copyright_info
nfr2ucsc.pl: create UCSC track file for input NFR regions
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

###############################################################################
## parse input options
use vars qw($nfrFile $help);

GetOptions ("i=s"  => \$nfrFile,
            "help" => \$help,
            "h"    => \$help);

usage() if($help || !$nfrFile);

###############################################################################
sub usage {
	print STDERR "\nProgram: nfr2ucsc.pl (create UCSC track file for input NFR regions)\n";
	print STDERR "Author: BRIC, University of Copenhagen, Denmark\n";
	print STDERR "Version: 1.0\n";
	print STDERR "Contact: pundhir\@binf.ku.dk\n";
	print STDERR "Usage: nfr2ucsc.pl -i <file> [OPTIONS]\n";
	print STDERR " -i <file>         [input file having NFRs]\n";
	print STDERR "[OPTIONS]\n";
	print STDERR " -h                [help]\n\n";
	exit(-1);
}

###############################################################################

## create output file for writing UCSC tracks of NFR (all)
open(INFILE, $nfrFile) || die $!;
my @data=<INFILE>;

print "track name=\"Predicted NFR ($nfrFile)\" description=\"Predicted NFR ($nfrFile)\" itemRgb=\"On\"\n";

foreach my $l(@data) {
    my @F=split(/\s+/,$l);
    my @startBlock=split(/[\:\-]+/,$F[6]);
    my @endBlock=split(/[\:\-]+/,$F[7]);
    my $expr_up=sprintf("%0.4f", $F[8]/$F[9]);
    my $expr_down=sprintf("%0.4f", $F[10]/$F[11]);
    my $expr_nfr=sprintf("%0.4f", $F[12]/$F[13]);
    print "$startBlock[0]\t$startBlock[1]\t$startBlock[2]\t$F[6]\t$expr_up\t.\t$startBlock[1]\t$startBlock[2]\t0,255,0\n";
    print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$expr_nfr\t$F[5]\t$F[1]\t$F[2]\t255,0,0\n";
    print "$endBlock[0]\t$endBlock[1]\t$endBlock[2]\t$F[7]\t$expr_down\t.\t$endBlock[1]\t$endBlock[2]\t0,255,0\n";
}

close INFILE;
exit(0);
