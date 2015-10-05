#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  make_histogram.pl
#
#        USAGE:  ./make_histogram.pl  
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
#      CREATED:  06/01/11 08:53:55
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


my $file = $ARGV[0];

open F, $file;
my @results;

while (my $line = <F>) { 
	chomp $line;
	$line =~ s/^>//;
	unless ($results[$line]) { $results[$line] = 0; }
	$results[$line]++;
}

for (my $i = 0; $i < $#results; $i++) { 
	print "$i\t$results[$i]\n";
}
