#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  FinishScaffolder_2.pl
#
#        USAGE:  ./FinishScaffolder_2.pl  
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
#      CREATED:  02/03/11 17:28:13
#     REVISION:  ---
#===============================================================================

my $PREFIX = "/home/lakatos/bdownie/bin/kilape";
use Parallel::ForkManager;
use strict;
use warnings;
use Getopt::Long;

use File::Copy;

my $velvetkmer = 31;
my $ins_length;
my $ins_length_sd;
my $exp_coverage;
my $num_duplicates = 10;
my $MAX_PROCESSES;
my $DEBUG;
my $firstdir;
my $lastdir;
#$MAX_PROCESSES=5;
#my $min_contig_length = 1000;

my $current_dir = `pwd`;
chomp $current_dir;
#chdir "working";
#my @files = `find working | grep "/contigs.fasta"`;
#$num_scaffolds = $ARGV[0];



#my $lastdir = $ARGV[0];
my $min_size = 0;
#my $ratio = $ARGV[1];

my $kill_lc = 0;
my $bso = 0;
my $jb;
my $min_ratio = 0;
GetOptions ('d' => \$DEBUG, 'il=i' => \$ins_length, 'ilsd=i' => \$ins_length_sd, 'min_size=i' => \$min_size,
			'exp_cov=i' => \$exp_coverage, 'dir_count=i' => \$lastdir, 'fd=i' => \$firstdir, 'min_ratio=f'=>\$min_ratio,
			'kill_lc' => \$kill_lc, 'bso'=>\$bso, 'just_best'=>\$jb);

unless ($MAX_PROCESSES) { $MAX_PROCESSES=25; }
unless ($lastdir) { print STDERR "Usage: -dir_count <# dirs>\n"; exit; }
my @files;
my $i = 1;
if ($firstdir) { $i = $firstdir; }
if ($lastdir < 1) { print STDERR "Dir count must be larger than 0\n"; exit; }
for (; $i <= $lastdir; $i++) { 
	my $string = "working/" . $i %10 . "/". $i %100 . "/". $i %1000 . "/". $i %10000 . "/". $i %100000 . "/". $i;
	#my $string = "working/" . $i %10 . "/". $i %100 . "/". $i %1000 . "/". $i %10000 . "/". $i %100000;
	push @files, $string;
}

#open F, "contigs.found";
#while (my $a = <F>) { 
#	chomp $a; 
#	push @files, $a;
#}
#close F;

if (-e "scaffolds.final2.fasta") { system("rm ./scaffolds.final2.fasta"); }


#$MAX_PROCESSES=15;
my $pm = new Parallel::ForkManager($MAX_PROCESSES);

my $num_threads = 0;
my $max_threads = 8;
my $start_kmer = 31;
my $end_kmer = 71;
my $kmer_step = 10;

my $return_val;
my $children = 0;
my $progress = 0;
my $last_progress_length = 1;
my $numfiles = $#files + 1;
my $num_lines = 0;
$| = 1;
my $lastlen = 1;
if ($firstdir) { 
	$progress = $firstdir; 
	$last_progress_length = length($firstdir);
}
$progress = 0;
$last_progress_length = 1;
print STDERR "\n";

#@files = ();
#$i = 1;
#if ($firstdir) { $i = $firstdir; }
#for (; $i <= $lastdir; $i++) { 
#	my $string = "working/" . $i %10 . "/". $i %100 . "/". $i %1000 . "/". $i %10000 . "/". $i %100000 . "/". $i %1000000;
#	push @files, $string;
#}

#if ($DEBUG) { exit; }
#my $n50 = 0;
#open SCAFF, ">$current_dir/scaffolds.final.fasta";
print STDERR "Consolidating: $progress";
my $also_contigs = 0;
if (-e "$files[0]/celera.ctg.fasta") { 
	open CONTIGS, ">contigs.final.fasta";
	$also_contigs = 1;
}
foreach my $dir (@files) {
	$progress++;
	if (($progress %1000) == 0) { 
		my $backspace = "\b" x $lastlen;
		$lastlen = length($progress);
		print STDERR $backspace;
		print STDERR $progress;
	}
	my $total_length = 0;
	next if (!-e "$dir/contigs.fasta");
	if ($jb) { 
		my $kmer = 0;
		if ((!-z "$dir/kmer") && (-e "$dir/kmer")) {
			open F, "$dir/kmer";
			my $line;
			while ($line = <F>) { 
				if ($line =~ /(\d+)/) { $kmer = $1; }
			}
			close F;
		}
		if ($kmer) { 
			if ((!-z "$dir/velvet_$kmer.contigs.fa") && (-e "$dir/velvet_$kmer.contigs.fa")) {
				open ASSEMBLY, "$dir/velvet_$kmer.contigs.fa";
				while (my $line = <ASSEMBLY>) {
					print $line;
				}
				close ASSEMBLY;
			}
			elsif ((!-z "$dir/velvet_$kmer/contigs.fa") && (-e "$dir/velvet_$kmer/contigs.fa")) {
				open ASSEMBLY, "$dir/velvet_$kmer/contigs.fa";
				while (my $line = <ASSEMBLY>) {
					print $line;
				}
				close ASSEMBLY;
			}
		}
		else { 
			if ((!-z "$dir/contigs.fasta") && (-e "$dir/contigs.fasta")) {
				open ASSEMBLY, "$dir/contigs.fasta";
				while (my $line = <ASSEMBLY>) {
					print $line;
				}
				close ASSEMBLY;
			}
		}
	}
	else { 

		if ((-z "$dir/best.fasta") || (!-e "$dir/best.fasta")) { 
			open F, "$dir/contigs.scaff.fasta"; 
			while (my $line = <F>) {
				print $line;
			}
			
			close F;
		}
		else {
			open F, "$dir/best.fasta";
	
			my $output_sequence = "";
			my $best_length = 0;
			my $tmp_sequence;
			while (my $line = <F>) { 
				if ($kill_lc) { 
					unless ($line =~ /^>/) {
						$line =~ tr/a/N/;
						$line =~ tr/c/N/;
						$line =~ tr/g/N/;
						$line =~ tr/t/N/;
					}
				}
				$output_sequence .= $line ;
				if ($line =~ /^>/) {
					if ($tmp_sequence) { 
						my $tmp_len = length($tmp_sequence);
						if ($tmp_len > $best_length) { 
							$best_length = $tmp_len;
						}
					}
					$tmp_sequence = "";
				}
				else {
					chomp $line;
					$tmp_sequence .= $line;
					#$best_length += length($line);
				}
				#print "$line\n"; 
			}
			close F;

			if ($min_ratio) { 
				open G, "$dir/contigs.scaff.fasta";
				my $contigs_line = "";
				my $contigs_length = 0;
				while (my $line = <G>) { 
					$contigs_line .= $line;
					unless ($line =~ /^>/) {
						chomp $line;
						$contigs_length += length($line);
					}
				}
				#print "($contigs_length)($min_ratio)($best_length)\n";
				if (($contigs_length * $min_ratio) > $best_length) {
					$output_sequence = $contigs_line;
				}
				close G;
			}
			print $output_sequence;
					
			#print "($also_contigs)($dir/celera.ctg.fasta)\n";
				#if (!(-e "$dir/celera.ctg.fasta")) { 
				#	print "$dir\n";
				#}
			if ($also_contigs) { 
				if (-e "$dir/celera.ctg.fasta") { 
					open G, "$dir/celera.ctg.fasta";
				}
				else {
					open G, "$dir/contigs.fasta";
				}
				while (my $line = <G>) { 
					print CONTIGS $line;
				}
				close G;
			}
		}
	}
}

my $backspace = "\b" x $lastlen;
$lastlen = length($progress);
print STDERR $backspace;
print STDERR $progress;
print STDERR "\n";
