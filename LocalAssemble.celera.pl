#!/usr/bin/perl 

# LocalAssembles need to be joined together, but that represents a major project which must come later.

use Cwd;
use File::Path;
use File::Copy;
#Path to executables

my $wgs_path = $ENV{'WGS_PATH'};
if (-e "$wgs_path/runCA") { $wgs_path .= "/"; }
unless (-e "${wgs_path}runCA") { 
	my $tmp1 = `which runCA`;
	my $tmp2 = `which fastqToCA`;
	chomp $tmp1;
	chomp $tmp2;
	$tmp1 =~ s/runCA//;
	$tmp2 =~ s/fastqToCA//;
	if ($tmp1 && ($tmp1 eq $tmp2)) { 
		$wgs_path = $tmp1;
	}
	else { 
		print STDERR "Couldn't find runCA. Make sure the WGS_PATH environmental variable is set or include it in your config file.\n"; 
		exit; 
	} 
}



push @INC, "/home/lakatos/bdownie/lib/perl5"; ## better would be a global "export PERL5LIB=/..." command beforehand
use Parallel::ForkManager;
use strict;
use warnings;
use Getopt::Long;
use File::Spec;

## Get the path where this program is located. All related programs should also be there.
#my $prefix = dirname(File::Spec->rel2abs($0));
my $prefix = "/home/lakatos/bdownie/bin/kilape/";

my $ins_length;
my $ins_length_sd;
my $exp_coverage;
my $num_duplicates;
my $MAX_PROCESSES;
my $lastdir;
my $low_kmer;
my $high_kmer;
my $kmer_step;
my $firstdir;
my $final_assembly;
my $clean;
my $reference = 0;
my $strict = 1;
my $build_scaffolds_only = 0;
my $mp = 0;
my $fa = 0;
my $unmasked_assemble = 0;
my $complete = 0;
my $last = 0;

my $DEBUG;

my $bwa = 0;
my $fix = 0;
my $seq454 = 0;

my $trim = 0;
my $ramdisk;
my $ec = 0;
my $l_value = 0;
my $m_value = 0;
my $noassemble;
my $velvet_correction = 0;
my $keep_low = 0;
my $num_read_libs = 0;
my $cluster_params;
my $contigs = 0;
my $start_dir = `pwd`;
my $clean_working_dir = 0;
chomp $start_dir;

GetOptions ('debug' => \$DEBUG,'dir_count=i' => \$lastdir, 't=i'=>\$MAX_PROCESSES, 
			'fd=i' => \$firstdir, 'ramdisk=s'=>\$ramdisk, 'contigs'=>\$contigs, '454'=>\$seq454, 'clean'=>\$clean_working_dir); 

unless ($lastdir) { 
	print STDERR "Usage: LocalAssemble.celera.pl -dir_count <#> [-k <kmer size> ] [-t <number of threads (25)] [-reference]\n";
	exit;
}
if ($DEBUG) { $MAX_PROCESSES = 1; }

unless ($MAX_PROCESSES) { $MAX_PROCESSES = 25; }
unless ($firstdir) { $firstdir = 1; }

if ($DEBUG) { print STDERR "In Debug mode...\n"; }
unless (-d "working") { print STDERR "No working directory found!\n"; exit; }


my $dirs_per_fork = 50;
my @files;
my $i = $firstdir;


my $total_iterations = 1+ int(($lastdir - $i + 1)/$dirs_per_fork);


$dirs_per_fork = 1;
for (; $i <= $lastdir; $i+=$dirs_per_fork) { 
	my $string;
	for (my $j = 0; $j < $dirs_per_fork; $j++) { 
		my $k = $j + $i;
		next if ($k > $lastdir);
		$string .= "working/" . $k %10 . "/". $k %100 . "/". $k %1000 . "/". $k %10000 . "/". $k %100000 . "/". $k . ":";
	}
	chop $string;
	push @files, $string;
}

if (-e "scaffolds.final.fasta") { unlink "scaffolds.final.fasta"; }

my @change_insert_after_line;
my @insert_sizes;


my $pm;
unless ($DEBUG || ($MAX_PROCESSES == 1)) { 
	$pm = new Parallel::ForkManager($MAX_PROCESSES);
}

my $num_threads = 0;
my $max_threads = 8;

my $return_val;
my $children = 0;
my $progress = 0;
my $numfiles = $#files + 1;
$| = 1;
print STDERR "Progress ($lastdir): 0";
$progress = $firstdir;
my $last_progress = 0;
my $last_progress_length = 1;
my $tmp = 1;
my $start_time = time;
my $orig_max_procs = $MAX_PROCESSES;
# First loop needs to be put into a commonly used file or method.
foreach my $dirset (@files) {
 	$progress+= $dirs_per_fork;
	my $progress_meter = 100 * int($progress/ 100);
	if ($progress_meter > $last_progress) { 
		for (my $i = 0; $i < $last_progress_length; $i++) { 
			print STDERR "\b";
		}
		print STDERR $progress_meter;
 		$last_progress_length = length($progress_meter);
		$last_progress = $progress_meter;
	}
	#if ($cluster_params) { sleep 1; }
 	unless ($DEBUG || ($MAX_PROCESSES == 1)) { my $pid = $pm->start and next; }

	my @directories = split /:/, $dirset;

	my @min_contig_lengths;
	my %contig_lengths;
	my %total_contig_lengths;
	my %num_contigs;
	my %backup_lengths;
	my %exp_covs;
	my %contigs_in_scaffolds;
	my %scaffold_coverage;
	my %reads_in_dir;
	my %max_read_size;

	foreach my $dir (@directories) {
		next if ((-z "$dir/contigs.fasta") || (!-e "$dir/contigs.fasta"));

		my $total_length = 0;
		my $total_length_scaff = 0;
		my $total_length_contig = 0;
		my @final_assembly = ();
		my $longest_contig_length = 0;
		my $smallest_contig_length = 0;
		my $count = 0;
		my @read_files_tmp = <$dir/*>;
		my @read_files = ();

		my @lengths = ();
		my $total_ns = 0;



		if ($DEBUG) { print "starting $dir\n"; }
		next unless (-e "$dir/contigs.fasta" && (!-z "$dir/contigs.fasta")); 
		#open LOG, ">$dir/kilape.log";
		open F, "$dir/contigs.fasta";
		my $contigs;
		my $has_n = 0;
		while (my $line = <F>) { 
			if (($line !~ /^>/) && ($line =~ /N/)) { $has_n = 1; }
			$contigs .= $line;
		}
		close F;
		
		next unless ($has_n);

		#open H, ">$dir/contigs.noNs.fasta";
		#open J, ">$dir/contigs.noNs.qual";
		#open SCAFFQUAL, ">$dir/contigs.scaff.qual";
		#while (my $line = <F>) { 
		#	my $line2 = <F>; 
		#	chomp $line;
		#	chomp $line2;
		#	my @array = split /N+/, $line2;
		#	my $z = 1;
		#	my @sequence = split "", $line2;
		#	foreach my $nuc (@sequence) { 
		#		if ($nuc eq "N") { $nuc = 0; }
		#		else { $nuc = 10; }
		#	}
		#	my $qual = join " ", @sequence;
		#	print SCAFFQUAL "$line\n$qual\n";
		#	
#
#			foreach my $seq (@array) {
#				chomp $seq;
#				if (length($seq) > 20) { 
#					my $tmp = "10 " x length $seq;
#					chop $tmp;
#					print H "$line-$z\n$seq\n";
#					print J "$line-$z\n$tmp\n";
#					$z++;
#				}
#			}
#		}
#		close F;
#		close H;
#		close J;
#		close SCAFFQUAL;


		my $max_insert_size = 0;
		#push @read_files_tmp, "$dir/contigs.noNs.fasta";
		foreach my $elem (@read_files_tmp) { 
			chomp $elem;
			next if ($elem =~ /unmatched/);
			next if ($elem =~ /masked/);
			if (($elem =~ /reads.*.fast./) || ($elem eq "contigs.noNs.fasta")) { 
				my $is = 0;
				my $issd = 0;
				if ($elem =~ /s\.(\d+)\.(.*)\.fast./) { 
					$is = $1; 
					$issd = $2;
				}
				if ($is > 0) { 
					my $mate_name = $elem;
					$mate_name =~ s/fasta/mates/; 
					$mate_name =~ s/fastq/mates/; 
					open X, ">$mate_name";
				}
				open Z, $elem || next;
				$num_read_libs++;
				my $qual_file = $elem;
				if ($elem =~ /fasta/) { 
					$qual_file =~ s/fasta/qual/;
					open Y, ">$qual_file";
					while (my $line = <Z>) {
						my $a = <Z>; 
						chomp $a; 
						$a = "40 " x length($a);
						print Y "$line$a\n";
						if ($is > 0) { 
							$line =~ s/>//; chomp $line;
							print X $line;
							$line = <Z>;
							$a = <Z>; 
							chomp $a; 
							$a = "40 " x length($a);
							print Y "$line$a\n";
							$line =~ s/>//; 
							print X " $line";
						}
					}
				}
				elsif (($elem =~ /fastq/) && ($is > 0)) { 
					while (my $line = <Z>) {
						<Z>; <Z>; <Z>; 
						my $line2 = <Z>; <Z>; <Z>; <Z>;
						$line =~ s/\@//; 
						chomp $line;
						$line2 =~ s/\@//; 
						print X "$line $line2";
					}
				}
				close Z;
				close Y;
				push @read_files, $elem;
			}
		}

		unless (@read_files) {
			open F, "$dir/contigs.fasta";
			open G, ">$dir/best.fasta";
			while (my $line = <F>) { 
				print G $line;
			}
			close F;
			close G;
			next;
		}

		$reads_in_dir{$dir} = join ":", @read_files;
	}
	if ($build_scaffolds_only) { 
		if ($MAX_PROCESSES != 1) {  $pm->finish; }
	}
	else { 
		# Much of this logic can be put into a commonly used subroutine (with other LocalAssembles)
		my $failed_long_assembly = 1;
		my $best_kmer = 0;
		if ($DEBUG) { print "Starting to process directories\n"; }
		foreach my $dir (@directories) {
			if ($DEBUG) { sleep 1; }
			my $best_line = "";
			unless ($reads_in_dir{$dir}) { 
				if (-e "$dir/contigs.fasta" && (!-z "$dir/contigs.fasta")) { 
					open F, "$dir/contigs.fasta";
					open G, ">$dir/best.fasta";
					while (my $line = <F>) { 
						print G $line; 
						$best_line .= $line;
					}
					close F;
					close G;
				}
			}
			else { 
				my @reads = split ":", $reads_in_dir{$dir};
				my $cur_dir = cwd();
				my $working_dir;
				my $tmp_dir;
				my $source_dir = "$cur_dir/$dir";
				system "rm -f $source_dir/*masked* $source_dir/best.fasta";
				if ($ramdisk) { 
					if ($dir =~ /(\d+)$/) {
						$tmp_dir = $1;
					}
					$working_dir = "$ramdisk/$tmp_dir";
					if ($DEBUG) { 
						print "mkdir -p $working_dir\n";
						print "cp $source_dir/*fastq $source_dir/*fasta $source_dir/*qual $source_dir/*mates $working_dir 2> /dev/null\n";
						print "chdir $working_dir\n";
					}
					else { 
						mkpath($working_dir);
						system "cp $source_dir/*fastq $source_dir/*fasta $source_dir/*qual $source_dir/*mates $working_dir 2> /dev/null";
						chdir "$working_dir";
					}
				}
				else { 
					$working_dir = $source_dir;
					if ($DEBUG) { 
						print "chdir $source_dir\n";
						print "rm -rf celera_dir\n";
					}
					else { 
						chdir $source_dir;
						rmtree("celera_dir");
					}
				}
	
				#my $runCA_cmd = "ionice -c 3 $wgs_path/runCA -d $working_dir/celera_dir -p celera_prefix unitigger=bogart cnsMinFrags=1000";;
				my $runCA_cmd = "ionice -c 3 " . $wgs_path . "runCA -d $working_dir/celera_dir -p celera_prefix unitigger=bogart";
				my $frg_exec;
				#push @reads, "contigs.0.0.fasta";
				foreach my $read (@reads) { 
					$read =~ s/$dir//; 
					$read =~ s/\///; 
					if ($read =~ /reads\.(\d+)\.(.*)\.fast./) { 
						my $is = $1;
						my $issd = $2;
						if ($issd =~ /\D/) { 
							$issd = int($is * 0.2);
						}
						#$issd *=2;
						my $lib_name = $read;
						$lib_name =~ s/\.fast.//;
						$lib_name =~ s/reads\.//; 
						my $qual_name;
						if ($read =~ /fastq/) { 
							my $is_illumina = 1;
							open TEST_QUAL, "$working_dir/$read";
							for (my $line_num = 0; ($line_num < 100) && $is_illumina; $line_num++) { 
								<TEST_QUAL>; <TEST_QUAL>; <TEST_QUAL>;
								my $line = <TEST_QUAL>;
								last unless ($line);
								chomp $line;
								my @tmp_array = split "", $line;
								foreach my $qual_val (@tmp_array) { 
									if (ord($qual_val) < 59) { 
										$is_illumina = 0; 
									}
								}
							}
							close TEST_QUAL;
							my $qual_type = "sanger";
							if ($is_illumina) { $qual_type = "illumina"; } 

							$frg_exec .=  "ionice -c 3 " . $wgs_path . "fastqToCA -libraryname $lib_name -type $qual_type";
							if ($is) { $frg_exec .= " -mates $read -insertsize $is $issd"; }
							else { $frg_exec .= " -reads $read"; }
						}
						else {
							$qual_name = $read; 
							$qual_name =~ s/fasta/qual/; 
							$frg_exec .=  "ionice -c 3 " . $wgs_path . "fastaToCA -l $lib_name -s $read -q $qual_name";
							if ($is) { 
								my $mates_name = $read; 
								$mates_name =~ s/fast./mates/; 
								my $issd2 = $issd * 2;
								$frg_exec .= " -m $mates_name -mean $is -stddev $issd2";
							}
						}
						my $frg_out = $read;
						$frg_out =~ s/fast./frg/; 
						$frg_exec .= " > $working_dir/$frg_out 2>> celera.log ;";
						$runCA_cmd .= " $frg_out";
					}
				}
				if ($contigs) { 
					#unless ($DEBUG) { 
					open F, "contigs.noNs.fasta";
					open G, ">contigs.noNs.qual";
					#open F, "contigs.fasta";
					#open G, ">contigs.qual";
					while (my $line = <F>) { 
						my $a = <F>; 
						chomp $a;
						my @tmp_array = split "", $a; 
						foreach my $nuc (@tmp_array) { 
							if ($nuc eq "N") { $b .= "0 "; }
							else { $b .= "20 "; }
						}
						chop $b;
						$b .= "\n"; 
						print G "$line$b\n";
					}
					close F;
					close G;
					#}
	
					$frg_exec .=  "ionice -c 3 " . $wgs_path . "fastaToCA -l contigs -s contigs.noNs.fasta -q contigs.noNs.qual > contigs.frg 2>> celera.log ; ";
					#$frg_exec .=  "ionice -c 3 $wgs_path/fastaToCA -l contigs -s contigs.fasta -q contigs.qual > contigs.frg 2>> celera.log ; ";
		
					$runCA_cmd .= " contigs.frg";
				}

				$runCA_cmd .= "	>> $working_dir/celera.log";
				$runCA_cmd .= " 2>&1"; 
				$runCA_cmd = $frg_exec . " " . $runCA_cmd;
				if ($DEBUG) { print "$runCA_cmd\n"; }
				else { 

					rmtree("$working_dir/celera_dir");
					system $runCA_cmd; 
					if ((-e "celera_dir/9-terminator/celera_prefix.scf.fasta") 
					 && (!-z "celera_dir/9-terminator/celera_prefix.scf.fasta")) { 
						move "celera_dir/9-terminator/celera_prefix.scf.fasta", "$source_dir/best.fasta";
					}
					else { copy "$source_dir/contigs.scaff.fasta", "$source_dir/best.fasta"; }
					if ((-e "celera_dir/9-terminator/celera_prefix.ctg.fasta") 
					 && (!-z "celera_dir/9-terminator/celera_prefix.ctg.fasta")) { 
						move "celera_dir/9-terminator/celera_prefix.ctg.fasta", "$source_dir/celera.ctg.fasta";
					}
					if ($ramdisk && ((!-e "$source_dir/celera.log") || (-z "$source_dir/celera.log"))) { 
						copy "$working_dir/celera.log", "$source_dir/celera.log";
					}
					my $to_del;
					if ($ramdisk) { 
						$to_del = "$working_dir";
					}
					else { 
						$to_del = "$source_dir/celera_dir";
					}
					chdir $source_dir;
					rmtree("$to_del");
					#if ($clean_working_dir) { rmtree("$to_del"); }
				}
			}
	
			if ($DEBUG) { sleep 1; }
		}
	
		#close LOG;
		unless ($MAX_PROCESSES == 1) { 
   		 	$pm->finish; # Terminates the child process
		}
	}
}
unless ($MAX_PROCESSES == 1) { 
	$pm->wait_all_children;
}


for (my $i = 0; $i < $last_progress_length; $i++) { 
	print STDERR "\b";
}
print STDERR $progress;
print STDERR "\n";


chdir "..";	
sub rev_comp { 
	my $seq = shift;
	if ($seq =~ /\d/) { return $seq; }

	$seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;

	return $seq;
}

