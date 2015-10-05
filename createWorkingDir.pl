#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  createWorkingDir.pl
#
#        USAGE:  ./createWorkingDir.pl  
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
#      CREATED:  05/16/11 08:32:41
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Parallel::ForkManager;
use File::Path;
use Getopt::Long;

my $scaffold_file;
my $contig_file;
my @contig_used;
my $DEBUG;
my $fill_gap_only;
my $fix = 0;

GetOptions ('d' => \$DEBUG, 's=s' => \$scaffold_file, 'c=s' => \$contig_file, 'f'=>\$fill_gap_only, 'fix'=>\$fix);
unless ($contig_file) { 
	print STDERR "Usage: createWorkingDir.pl -s <scaffold file> -c <assembly file>\n";
	exit; 
}

my $num_scaffolds = 0;
if ($scaffold_file) { 
	open (S, "<$scaffold_file") or die "Cannot open the file \'$scaffold_file\': $!\n";
	while (<S>) { 
		my @contigs = split /\s+/;	
		foreach my $contig (@contigs) {
			$contig_used[$contig] = 1;
		}
		$num_scaffolds++;	
	}
	close S;
}

unless ($fill_gap_only) { 
	my $last_contig = `tail -n 2 $contig_file | head -n 1`;
	$last_contig =~ s/>//; 
	chomp $last_contig;

	if ($fix) { $num_scaffolds = $last_contig; }
	else { 
		for (my $i = 1; $i <= $last_contig; $i++) { 
			unless ($contig_used[$i]) { $num_scaffolds++; }
		}
	}
}

#my $num_dirs = $ARGV[0];
#unless ($num_dirs) { print STDERR "Usage: createWorkingDir.pl <directory #>\n"; exit; }
$num_scaffolds++;
print "Calculated number of scaffolds as: $num_scaffolds\n" if ($DEBUG);

open F, ">scaffolds.count";
print F "$num_scaffolds\n";
close F;

if (-d "working") { 
	my $i = 1;
	my $dir = "working_1";
	while (-d $dir) { 
		$i++;
		$dir = "working_" . $i;
	}
	rename "working", $dir;
}
