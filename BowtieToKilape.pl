#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  BowtieToSam.pl
#
#        USAGE:  ./BowtieToSam.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  YOUR NAME (), 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  06/14/12 14:09:22
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;



use Getopt::Long;

my $contig_file;
my $lastread;
GetOptions ('c=s'=>\$contig_file, 'l=s'=>\$lastread);

unless ($contig_file && $lastread) { 
	print STDERR "Usage: BowtieToSam.pl -c <reference file> -l <last read number>\n";
	exit;
}

my %lengths;
my @length_accs;
open F, $contig_file;
my $acc;
while (my $line = <F>) { 
	chomp $line;
	if ($line =~ />/) { 
		$line =~ s/>//;
		$acc = $line;
		$lengths{$acc} = 0;
		push @length_accs, $acc;
	}
	else { $lengths{$acc} += length($line); }
}

foreach my $acc (@length_accs) {
	if ($lengths{$acc}) {
		print "\@SQ\tSN:$acc\tLN:$lengths{$acc}\n";
	}
}

my @trash_array = (0,0,"*", "*");

my @lastarray;
my $lastnum = 0;
while (my $line = <STDIN>) { 
	my $oline;
	chomp $line;
	my @array = split /\s+/, $line;
	if ($array[1] eq "+") { $array[1] = 0; }
	else { $array[1] = 16; }
	@array = ($array[0], $array[1], $array[2], $array[3] + 1);

	$oline = join "\t", @array;
	if ($array[0] =~ /\D/) { 
		print "$oline\n";
	}
	else {
		while ($array[0] > ($lastnum + 1)) {
			$lastnum++;
			$trash_array[0] = $lastnum;
			$oline = join "\t", @trash_array;
			print "$oline\n";
		}
		$lastnum = $array[0];
		$oline = join "\t", @array;
		print "$oline\n";
	}
}
			
$lastnum++;
while ($lastnum <= $lastread) { 
	$trash_array[0] = $lastnum;
	my $oline = join "\t", @trash_array;
	print "$oline\n";
	$lastnum++;
}


		
