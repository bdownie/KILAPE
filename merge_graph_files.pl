#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  merge_graph_files.pl
#
#        USAGE:  ./merge_graph_files.pl  
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
#      CREATED:  10/25/11 12:14:29
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my %dists;
my %edges;
my %direct1s;
my %direct2s;
my %buildscaffolds_dist;
my %avg_insert_sizes;

if ($#ARGV < 1) { print STDERR "Usage: merge_graph_files.pl <contigs.fasta> [<graph file> <insert size>] \n"; exit; }

my $contigs = shift @ARGV;
my @contig_sizes;
my $contig_fasta;
my $contig_acc;
#open F, $contigs;
#while (my $line = <F>) {
#	if ($line =~ />(\d+)/) {
#		if ($contig_acc) {
#			$contig_sizes[$contig_acc] = length($contig_fasta);
#		}
#		$contig_acc = $1;
#		$contig_fasta = "";
#	}
#	else {
#		$contig_fasta .= $line;
#	}
#}
#$contig_sizes[$contig_acc] = length($contig_fasta);
#close F;
my $biggest_insert_size = 0;

while (@ARGV) { 
	my $file = shift @ARGV;
	my $insert_size = shift @ARGV;
	if ($insert_size > $biggest_insert_size) { $biggest_insert_size = $insert_size; }
	print STDERR "($file)($insert_size)\n";

	open F, $file;
	while (my $edge = <F>) { 
		my $dist = <F>; 
		my $direct1 = <F>; 
		my $direct2 = <F>; 

		$edge =~ m/edge:\s+(\d+)\s+(\d+)\s+(\d+)/;
		my $edge1 = $1;
		my $edge2 = $2;
		my $edge_count = $3;

#		next if (($contig_sizes[$edge1] + $contig_sizes[$edge2]) < $insert_size);
	
		my $edge_string = $edge1 . " " . $edge2;
		
		$dist =~ m/dist\s+(\d+)/;
		$dist = $1;
		$direct1 =~ m/1direct\s+(\d+)/;
		$direct1 = $1;
		$direct2 =~ m/2direct\s+(\d+)/;
		$direct2 = $1;
	
		my $scaff_dist = $dist;
		#$scaff_dist -= $insert_size;
		$scaff_dist = $insert_size - $scaff_dist;
		$scaff_dist *= $edge_count;
	#	print "($scaff_dist)($dist)($insert_size)($edge_count)($edge_string)\n";

	#	$dist -= $insert_size;
		$dist *= $edge_count;
			
		unless($edges{$edge_string}) { $edges{$edge_string} = 0; }
		unless($dists{$edge_string}) { $dists{$edge_string} = 0; }
		unless($direct1s{$edge_string}) { $direct1s{$edge_string} = 0; }
		unless($direct2s{$edge_string}) { $direct2s{$edge_string} = 0; }
		unless($avg_insert_sizes{$edge_string}) { $avg_insert_sizes{$edge_string} = 0; }
		unless($buildscaffolds_dist{$edge_string}) { $buildscaffolds_dist{$edge_string} = 0; }

		$edges{$edge_string} += $edge_count;
		$dists{$edge_string} += $dist;
		$direct1s{$edge_string} += $direct1;
		$direct2s{$edge_string} += $direct2;
		$buildscaffolds_dist{$edge_string} += $scaff_dist;
		$avg_insert_sizes{$edge_string} /= $insert_size;
	}
	close F;
}

open F, ">BuildScaffolds.graph";
#my @array;
#foreach my $key (keys %edges) { 
#	push @array, $edges{$key};
#}
#@array = sort { $a <=> $b } @array;
#print STDERR "Median count: " . $array[$#array*.9]. "\n";
open G, ">discarded_edges.graph";
foreach my $key (keys %edges) {
	my $print_edge = 0;
	if ($direct1s{$key} > (0.9 * $edges{$key})) { 
		$direct1s{$key} = $edges{$key}; 
		$print_edge++;
	}
	elsif ($direct1s{$key} < (0.1 * $edges{$key})) { 
		$direct1s{$key} = 0;
		$print_edge++;
	}
	if ($direct2s{$key} > (0.9 * $edges{$key})) { 
		$direct2s{$key} = $edges{$key}; 
		$print_edge++;
	}
	elsif ($direct2s{$key} < (0.1 * $edges{$key})) { 
		$direct2s{$key} = 0;
		$print_edge++;
	}
	if ($print_edge >= 2) { 
	#if ((($direct1s{$key} == $edges{$key}) || ($direct1s{$key} == 0)) && 
	#	(($direct2s{$key} == $edges{$key}) || ($direct2s{$key} == 0))) { 
		my $total_dist = int($dists{$key} / $edges{$key});
		my $buildscaffolds_distance = int($buildscaffolds_dist{$key} / $edges{$key});
	#	my $avg_insert_size = int($avg_insert_sizes{$key}/ $edges{$key});
	#	print "($avg_insert_size)\n";
		#$buildscaffolds_distance = $buildscaffolds_dist{$key} / $edges{$key};
	#	print "($buildscaffolds_dist{$key})($edges{$key})\n";
		#if ($total_dist < ($biggest_insert_size * 1.2)) { 
		if ($buildscaffolds_distance > (0.8 * (0-$biggest_insert_size))) { 
			print F "edge: $key $edges{$key}\n";
			print F "dist $buildscaffolds_distance\n";
			print F "1direct $direct1s{$key}\n";
			print F "2direct $direct2s{$key}\n";
		#}

		#if ($total_dist < 0) { $total_dist = 2; }
			print "edge: $key $edges{$key}\n";
			print "dist $total_dist\n";
			print "1direct $direct1s{$key}\n";
			print "2direct $direct2s{$key}\n";
		}
		else { 
			print G "edge: $key $edges{$key}\n";
			print G "dist $buildscaffolds_distance\n";
			print G "1direct $direct1s{$key}\n";
			print G "2direct $direct2s{$key}\n";
		}
	}
	else { 
		my $total_dist = int($dists{$key} / $edges{$key});
		my $buildscaffolds_distance = int($buildscaffolds_dist{$key} / $edges{$key});
		print G "edge: $key $edges{$key}\n";
		print G "dist $buildscaffolds_distance\n";
		print G "1direct $direct1s{$key}\n";
		print G "2direct $direct2s{$key}\n";
	}
}
