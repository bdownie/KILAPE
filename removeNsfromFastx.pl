#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  removeNsfromFastx.pl
#
#        USAGE:  ./removeNsfromFastx.pl  
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
#      CREATED:  08/18/11 12:38:28
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Getopt::Long;

my $min_length = 30;
my $qual = 0;
my $renumber = 0;
my $is_fastq = 0;
my $quiet = 0;
GetOptions('l=i'=>\$min_length, 'q=i'=>\$qual, '-quiet'=>\$quiet);
my $file = shift @ARGV;
unless ($file) { print STDERR "Usage: removeNsFromFastx.pl [ -l <min read length> -q <quality offset> ] <Fastq/a file>\n\n"; exit; } 
open F, $file;
my $guess = <F>;
my $reps = "NNNNN\n";
my $repq = "BBBBB\n";
if ($guess =~ /^@/) { 
	$is_fastq = 1;
	unless ($quiet) { 
		print STDERR "Expecting fastq file.\n";
	}
	<F>; <F>; 
	$guess = <F>; 
	unless ($qual) { 
		unless ($quiet) { 
			print STDERR "Guessing quality offset\n";
		}
		chomp $guess;
		$qual = 64;
		my @guess_array = split "", $guess;
		foreach my $qualval (@guess_array) { 
			if (ord($qualval) < 64) {
				#print "($qualval)\n";
				$qual = 33;
				$repq = "!!!!!\n";
			}
		}
	}
	unless ($quiet) { 
		print STDERR "Using quality offset of $qual\n";
	}
}
elsif ($guess =~ /^>/) { 
	print STDERR "Expecting fasta file\n";
}
else { 
	print STDERR "Bad file format in $file\n";
	exit;
}
seek(F,0,0);


unless ($quiet) { 
	print STDERR "Min read length: $min_length\n";
}
my $i = 1;

my $line3;
my $line4;
while (my $line = <F>) { 
	my $line2 = <F>; 
	if ($is_fastq) { 
		$line3 = <F>; 
		$line4 = <F>; 
	}
	
	chomp $line2;
	my @line_array = split /N+/, $line2;
	#my @line_array = split /N+/, $line2;
	#print "(@line_array)($line2)\n";
	my $max_len = 0;
	my $best_val = "";
	foreach my $val (@line_array) { 
		my $len = length($val);
		if ($len > $max_len) { $best_val = $val; $max_len = $len; }
	}
	if ($max_len < $min_length) { 
		$line2 = $reps; 
		if ($is_fastq) { 
			$line4 = $repq;
		}
	}
	else { 
		$line2 =~ /(.*)$best_val(.*)/;
		my $strip1 = $1;
		my $strip2 = $2; 
		my $strip1_len = length($strip1);
		my $strip2_len = length($strip2);
		if ($is_fastq) { 
			chomp $line4; 
			$line4 =~ s/.{$strip1_len}(.*).{$strip2_len}/$1/;
			$line4 = $line4 . "\n";
			$line2 = $best_val . "\n"; 
		}
	}

	print "$line$line2";
	if ($is_fastq) { print "$line3$line4"; }
}


		
