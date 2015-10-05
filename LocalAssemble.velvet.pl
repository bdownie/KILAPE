#!/usr/bin/perl 

#Path to executables
#my $path_to_velvet="/path/to/velvet_without_binary";
#my $bwa_path = "/path/to/bwa";
#my $bowtie_path = "/path/to/bowtie";
my $velvet_path = $ENV{'VELVET_PATH'};
if (-e "$velvet_path/velveth") { $wgs_path .= "/"; }
unless (-e "${velvet_path}velveth") { print STDERR "Couldn't find velvet. Make sure the VELVET_PATH environmental variable is set or include it in your config file.\n"; exit; } 

#push @INC, "/path/to/libs"; ## better would be a global "export PERL5LIB=/..." command beforehand
use Parallel::ForkManager;
use strict;
use warnings;
use Getopt::Long;
use File::Path;
use File::Copy;
#use File::Path;
#use File::Basename;
use File::Spec;
my $min_contig_ratio = .8;
my $NX_value = 80;

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
my $cov_cutoff = "auto";

my $current_dir = `pwd`;
chomp $current_dir;
my $DEBUG;

my $bwa = 0;
my $bowtie = 0;
my $bowtie2 = 0;
my $fix = 0;

my $trim = 0;
my $ramdisk;
my $ec = 0;
my $l_value = 0;
my $m_value = 0;
my $noassemble;
my $velvet_correction = 0;
my $keep_low = 0;
my $num_read_libs = 0;
my $velvet_correct = 1;

GetOptions ('debug' => \$DEBUG,'lk=i' => \$low_kmer, 'hk=i' => \$high_kmer, 'ks=i' => \$kmer_step,  
			'exp_cov=i' => \$exp_coverage, 'dir_count=i' => \$lastdir, 't=i'=>\$MAX_PROCESSES, 'fd=i' => \$firstdir, 'final_assembly'=>\$fa,
			'clean'=>\$clean, 'reference'=>\$reference,'strict=i'=>\$strict,  'mp=i'=>\$mp,'c'=>\$complete,
			'bwa'=>\$bwa, 'bowtie'=>\$bowtie, 'bowtie2'=>\$bowtie2,'last'=>\$last, 'ramdisk=s'=>\$ramdisk, 
			'noassemble'=>\$noassemble, 'velv_corr'=>\$velvet_correction,'keep_low'=>\$keep_low, 'cov_cutoff=s'=>\$cov_cutoff );

unless ($lastdir) { 
	print STDERR "Usage: LocalAssemble.velvet.pl -dir_count <#> [-k <kmer size> ] [-t <number of threads (25)] [-reference]\n";
	exit;
}
if ($bowtie2) { 
	$bwa = 0; 
	$bowtie = 0;
}
elsif ($bowtie) { 
	$bwa = 0; 
	$bowtie2 = 0;
}
elsif ($bwa) { 
	$bowtie = 0; 
	$bowtie2 = 0;
}

unless ($low_kmer) { $low_kmer = 31; }
unless ($high_kmer) { $high_kmer = 71; }
unless ($kmer_step) { $kmer_step = 10; }
if ($DEBUG) { $MAX_PROCESSES = 1; }

if ($fix) { $low_kmer = 25; $high_kmer = 37; $kmer_step = 2; }


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

if (-e "scaffolds.final.fasta") { system("rm ./scaffolds.final.fasta"); }

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
my $start_dir = `pwd`;
chomp $start_dir;

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
 	unless ($DEBUG || ($MAX_PROCESSES == 1)) { my $pid = $pm->start and next; }

	my @directories = split /:/, $dirset;

	my @min_contig_lengths;
	my %contig_lengths;
	my %total_contig_lengths;
	my %num_contigs;
	my @tmp_directories = ();
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
		my @read_files = ();
		chdir $dir;
		my @read_files_tmp = <*>;

		open CONTIGS, "contigs.fasta";
		my @lengths = ();
		my $total_ns = 0;
		while (my $line = <CONTIGS>) { 
			push @final_assembly, $line; 
			if ($line =~ />/) { 
				$count++;
			}
			else {
				my $len = length($line);
				$total_ns += $line =~ tr/N//;
				$total_ns += $line =~ tr/n//;
				push @lengths, $len;
				$total_length_contig += $len;
			}
		}
		close CONTIGS;

		next if ($total_ns == 0);


		if ($DEBUG) { print "starting $dir\n"; }
		next unless (-e "contigs.fasta" && (!-z "contigs.fasta")); 
		open LOG, ">kilape.log";
		open H, ">contigs.noNs.fasta";
		my $z = 1;
		foreach my $line (@final_assembly) { 
			next if $line =~ />/;
			chomp $line;
			my @seqs = split /N{100}N+/, $line;
			foreach my $seq (@seqs) {
				print H ">$z\n$seq\n";
			}
		}
		close H;


		if (-e "reads.graph") { 
			open G, "reads.graph";
			open H, ">reads.graph2";
			$mp = 0;
			while (my $a = <G>) {
				chomp $a;
				if ($a =~ /edge:\s\d+\s\d+\s(\d+)/) { 
					my $tmp = $1; 
					if ($mp == 0) { $mp = $tmp; } 
					elsif ($tmp < $mp) { $mp = $tmp; }
				}
				if ($a =~ /\dedge/) { $a =~ s/edge/\nedge/; }
				if ($a =~ /dist -/) { $a = "dist 2"; }
				unless ($a eq "") { print H "$a\n"; }
			}
			close H;
			close G;
			print LOG "Finished making graph2\n";
		}

		 
		open G, "reads.scaff";
		while (my $line = <G>) { 
			chomp $line;
			if ($line eq "") { 
				open H, "contigs.fasta";
				while (my $line2 = <H>) { 
					if ($line2 =~ />(.*)/) {
						$contigs_in_scaffolds{$dir} = $1;

					}
				}
			}
			else { 
				$contigs_in_scaffolds{$dir} = $line;
			}
		}
		close G;
		print LOG "Finished parsing scaffold file: reads.scaff\n";
			
		

		my $max_insert_size = 0;
		foreach my $elem (@read_files_tmp) { 
			chomp $elem;
			next if ($elem =~ /unmatched/);
			next if ($elem =~ /masked/);
			if ($elem =~ /reads\.(\d+).*\.fast/) { 
				my $is = $1;
				if ($is < $total_length_contig) { 
					push @read_files, $elem;  
					if ($max_insert_size < $is) {
						$max_insert_size = $is; 
					}
				}
			}
		}
		$max_read_size{$dir} = $max_insert_size;
		unless (@read_files) {
			open F, "contigs.fasta";
			open G, ">best.fasta";
			while (my $line = <F>) { 
				print G $line;
			}
			close F;
			close G;
			next;
		}
		my $total_read_length = 0;
		$scaffold_coverage{$dir} = 0;
		my @read_files_tmp2 = ();
		foreach my $read_file (@read_files) {
			chomp $read_file;
			push @read_files_tmp2, $read_file;
			my $read_length = 0;
			open Z, $read_file || next;
			$num_read_libs++;
			while (my $line = <Z>) {
				unless ($line =~ /^>/) { 
					chomp $line;
					$total_read_length += length($line);
					$read_length += length($line);
					unless ($read_file =~ /-1/)  { $scaffold_coverage{$dir} += length($line); }
				}
			}
			close Z;
			if ($unmasked_assemble) { 
				if ($DEBUG) { print "$prefix/matchSeedswithContigs -c contigs.fasta -f <placeholder> -s 21 -l $read_length > $read_file.unmatched\n"; }
				push @read_files_tmp2, "$read_file.unmatched";
			}
		}	
		print LOG "Found read files: (@read_files)\n";
		@read_files = @read_files_tmp2;
		$reads_in_dir{$dir} = join ":", @read_files;

		#if (!$fix && !$last && !$reference) { 
		#	if ($count == 1) { 
		#		open G, ">best.fasta";
		#		foreach my $line (@final_assembly) { 
		#			print G $line; 
		#		}
		#		unless ($MAX_PROCESSES == 1) { $pm->finish; }
		#	}
		#	close G;
		#}
		my $exp_cover = int($total_read_length/$total_length_contig) + 1;
		$scaffold_coverage{$dir} = int($scaffold_coverage{$dir}/$total_length_contig);
		$exp_covs{$dir} = $exp_cover;

		@lengths = sort {$a <=> $b} @lengths;
		$longest_contig_length = $lengths[$#lengths] + 1;
		my $lrg_backup_lengths;
		if ($count > 1) { $lrg_backup_lengths = $lengths[0]  + $lengths[1]; }
		else { $lrg_backup_lengths = $lengths[0]; }
		$backup_lengths{$dir} = "$lengths[0]:$lrg_backup_lengths";

		my $mean_contig_length = $total_length_contig / $count;


		$contig_lengths{$dir} = int($mean_contig_length);
		$total_contig_lengths{$dir} = $total_length_contig;
		push @tmp_directories, $dir; 
		chdir $start_dir;
	}

	@directories = @tmp_directories;
	if ($build_scaffolds_only) { 
		if ($MAX_PROCESSES != 1) {  $pm->finish; }
	}
	else { 

	my $failed_long_assembly = 1;
	my $best_kmer = 0;
	if ($DEBUG) { print "Starting to process directories\n"; }
	foreach my $dir (@directories) {
		if ($DEBUG) { sleep 1; }
		open LOG, ">>$dir/kilape.log";
		print LOG "Starting to execute local assemblies\n";
		open MASTERLOG, ">>velvet.log";
		print MASTERLOG "Processing $dir\n";
		close MASTERLOG;
		my @reads = split ":", $reads_in_dir{$dir};
		my $best_line = "";
		if ($#reads < 0) { 
			open F, "$dir/contigs.fasta";
			open G, ">$dir/best.fasta";
			while (my $line = <F>) { 
				print G $line; 
				$best_line .= $line;
			}
			close F;
			close G;
			print LOG "Error: Couldn't find reads\n";
		}
		else { 
			my $working_dir;
			my $source_dir = $start_dir . "/" . $dir;
			my $total_length = $total_contig_lengths{$dir};
			my $kmer;
			my $mcl = $contig_lengths{$dir};
			my $exp_cov = $exp_covs{$dir};
			if ($ramdisk) { 
				$dir =~ /(\d+)$/;
				$working_dir = $ramdisk . "/" . $1;
				if ($DEBUG) { 
					print "rm -rf $working_dir\n";
					print "mkdir -p $working_dir\n";
					print "cp -r $source_dir $working_dir\n";
					print "cd $working_dir\n";
				}
				else { 
					system "rm -rf $working_dir";
					system "mkdir -p $working_dir";
					chdir $source_dir;
					foreach my $read_file (@reads) {
						my $new_read_file = "$working_dir/" . $read_file;
						open A, $read_file;
						open B, ">$new_read_file";
						while (my $line = <A>) {
							print B $line;
						}
						close A;
						close B;
					}
					#foreach my $read_file (@reads) {
					#	my $masked_read_file = $read_file;
					#	$masked_read_file =~ s/fast/masked.fast/;
					#	my $new_masked_read_file = $masked_read_file;
					#	$new_masked_read_file =~ s/$dir/$working_dir/;
					#	open A, $masked_read_file;
					#	open B, ">$new_masked_read_file";
					#	while (my $line = <A>) {
					#		print B $line;
					#	}
					#	close A;
					#	close B;
					#}

					open A, "contigs.scaff.fasta";
					open B, ">$working_dir/contigs.scaff.fasta";
					while (my $line = <A>) {
						print B $line;
					}
					close A;
					close B;
					open A, "contigs.fasta";
					open B, ">$working_dir/contigs.fasta";
					while (my $line = <A>) {
						print B $line;
					}
					close A;
					close B;
					open A, "contigs.noNs.fasta";
					open B, ">$working_dir/contigs.noNs.fasta";
					while (my $line = <A>) {
						print B $line;
					}
					close A;
					close B;
				}
			}
			else { $working_dir = $source_dir; }
			chdir $working_dir;
			my @velveth_cmds = build_velveth_cmd($low_kmer, $high_kmer, $kmer_step, $working_dir, $mcl, $total_length, $exp_cov, @reads);
			my @velvetg_cmds = build_velvetg_cmd($low_kmer, $high_kmer, $kmer_step, $working_dir, $mcl, $total_length, $exp_cov, @reads);
			my $best_assembly_ratio = 1;
			my $velvethash = shift @velveth_cmds;
			my @velveth_cmds_tmp = @velveth_cmds;
			my @velvetg_cmds_tmp = @velvetg_cmds;

			for (my $kmer = $high_kmer; $kmer >= $low_kmer; $kmer-= $kmer_step) {
				mkdir "$working_dir/velvet_$kmer";
				if (-e "$working_dir/velvet_$kmer/Sequences") { unlink "$working_dir/velvet_$kmer/Sequences"; }
				symlink "$working_dir/Sequences", "$working_dir/velvet_$kmer/Sequences";
			}
		#	unless ($DEBUG) { 
		#		open G, ">$working_dir/allreads.fasta";
		#		foreach my $read_file (@reads) { 
		#			$read_file =~ s/$dir/$working_dir/;
		#			open F, "$read_file";
		#			while (my $line = <F>) { 
		#				print G $line;
		#			}
		#			close F;
		#		}
		#		close G;
		#	}
			#open G, ">$working_dir/allreads.masked.fasta";
			#foreach my $read_file (@reads) { 
			#	if ($ramdisk) { $read_file =~ s/$dir/$ramdisk/; }
			#	$read_file =~ s/fast/masked.fast/; 
			#	if ((-e "$read_file") && (!-z "$read_file")) { 
			#		open F, "$read_file";
			#		while (my $line = <F>) { 
			#			print G $line;
			#		}
			#		close F;
			#	}
			#}

			my $best_substitute_kmer = 0;
			my $best_length = 0;
			my $run_again= 1;
			my $smallest_assembly_length = 0;

			if ($ec) { 
				open F, "$working_dir/contigs.scaff.fasta";
				while (my $line = <F>) { 
					$best_line .= $line;
				}
				$best_kmer = 1;
			}
			elsif ($noassemble) { 
				my $best_file;
				open F, "$source_dir/kmer";
				while (my $line = <F>) {
					if ($line =~ /Using kmer (\d+)/) { 
						$best_kmer = $1;
						$best_file = "velvet_$best_kmer.contigs.fa";
					}
					if ($line =~ /contigs.scaff.fasta/) {
						$best_file = "contigs.scaff.fasta";
					}

				}
				close F;
				open G, $best_file;
				while (my $line = <G>) {
					$best_line .= $line;
				}
			}

			else { 
				#my $tmp = `pwd`; print $tmp;
				if ($DEBUG) { print "$velvethash\n"; }
				else { 
					system $velvethash; 
				}
				$best_line = "";

				my $velvetcmd = join " ", @velveth_cmds;
				$velvetcmd .= join " ; ", @velvetg_cmds;
				if ($DEBUG) { print "$velvetcmd\n"; }
				else { 
					system $velvetcmd; 
					my $least_contig_count = -1;
					my $best_sequence ="";
					my @acc = ();
					my @fasta = ();
					my $total_length = 0;
					open SCAFF, "$working_dir/contigs.scaff.fasta";
					while (my $scaff_line = <SCAFF>) { 
						next if ($scaff_line =~ /^>/);
						chomp $scaff_line;
						$total_length += length($scaff_line);
					}
					close SCAFF; 

					$total_length *= $min_contig_ratio;

					my $best_n_length = -1;
					my $best_no_n_length = 0;
					my $node_count = 0;
					for (my $i = $low_kmer; $i <= $high_kmer; $i += $kmer_step) {
						@fasta = ();
						@acc = ();
						open F, "$working_dir/velvet_$i/contigs.fa";
						while (my $line = <F>) { 
							chomp $line;
							if ($line =~ /^>/) { push @acc, $line; }
							else { $fasta[$#acc] .= $line; }
						}
						close F;
						my @sorted_fasta = sort { length $b <=> length $a } @fasta;
						my $tmp_length = 0;
						my $n_length = 0;
						my $no_n_length = 0;
						foreach my $seq (@sorted_fasta) { 
							$tmp_length += length($seq);
							$n_length += $seq =~ tr/N//; 
							$no_n_length += length($seq);
							$node_count++;
							if ($tmp_length > $total_length) {
								if (($best_n_length == -1) || (($best_n_length > $n_length) && ($best_no_n_length < $no_n_length))) {
									$least_contig_count = $node_count;
									$best_n_length = $n_length;
									$best_no_n_length = $no_n_length;
									$best_kmer = $i;
									$best_line = "";
									for (my $j = 0; $j <= $#acc; $j++) {
										$best_line .= $acc[$j] . "\n" . $fasta[$j] . "\n";
									}
								}
								last;
							}
						}
					}
				}
			}
	
			# New velvet requires sam entry...
				if ($best_kmer == 0) { 
					$best_line = "";
					if ($last && -e "$source_dir/contigs.scaff.fasta" && !-z "$source_dir/contigs.scaff.fasta") { 
						open F, "$source_dir/contigs.scaff.fasta";
						open G, ">kmer"; print G "Unable to assemble, using contigs.scaff.fasta\n"; close G; 
					}
					else { 
						open F, "$source_dir/contigs.fasta";
						open G, ">kmer"; print G "Unable to assemble, using contigs.fasta\n"; close G; 
					}
					while (my $line = <F>) { 
						$best_line .= $line;
					}
					close F;
				}
				else { open F, ">kmer"; print F "Using kmer $best_kmer\n"; close F; }
	
			if ($best_line)  { 
				if ($l_value < 0) { 
					open G, ">best.fasta";
					print G $best_line;
					close G;
				}
				else {
					open G, ">$working_dir/best.fasta";
					print G $best_line;
					close G;
					my $validate_cmd;
					my $mapping_reads;
					if ($trim) { 
						$mapping_reads = "allreads.masked.fasta";
					} 
					else { 
						$mapping_reads = "allreads.fasta";
						$mapping_reads = "allreads.masked.fasta";
					}
					my $to_correct;
					if ($ec) { 
						$to_correct = "contigs.scaff.fasta";
					}
					else { 
						$to_correct = "best.fasta";
					}
		#			if ($bwa) { 
		#				$validate_cmd = "$bwa_path/bwa index $working_dir/best.fasta 2> /dev/null; $bwa_path/bwa aln $working_dir/best.fasta $working_dir/allreads.fasta > $working_dir/best.sai 2> /dev/null; $bwa_path/bwa samse  $working_dir/best.fasta $working_dir/best.sai  $working_dir/$mapping_reads 2> /dev/null | $prefix/bin/BreakBadScaffolds -o $working_dir/contigs.scaff.fasta -i $max_read_size{$dir} > best.fasta";
		#				}
		#			elsif ($bowtie) { 
		#				$validate_cmd = "$bowtie_path/bowtie-build --quiet $working_dir/$to_correct $working_dir/$to_correct.bowtie 2> /dev/null; $bowtie_path/bowtie -S --quiet -f  $working_dir/$to_correct.bowtie $working_dir/$mapping_reads | $prefix/bin/BreakBadScaffolds ";
		#				if ($l_value) { $validate_cmd .= "-l $l_value "; }
		#				if ($m_value) { $validate_cmd .= "-m $m_value "; }
		#				$validate_cmd .= " -o $working_dir/$to_correct -i 0 ";
		#				if ($trim) { $validate_cmd .= " -t "; }
		#				$validate_cmd .= " > $working_dir/best.fasta2; cp $working_dir/best.fasta2 best.fasta";
		#			}
		#			if ($validate_cmd) { 
		#				if ($DEBUG) { 
		#					print("$validate_cmd\n");
		#				}
		#				else { 
		#					unless ($working_dir eq $dir) { 
		#						system("cp $working_dir/best.fasta $dir/best.fasta");
		#					}
		#				}
		#			}
					if ($ramdisk) { 
						if ($DEBUG) { 
							print "cp -r $working_dir/velvet* $source_dir\n";
							print "cp -r $working_dir/best.fasta $source_dir\n";
						}
						elsif (!$ec) { 
							for (my $i = $low_kmer; $i <= $high_kmer; $i+=$kmer_step) { 
								unless (-z "$working_dir/velvet_$i/contigs.fa" || (!-e "$working_dir/velvet_$i/contigs.fa")) { 
									open A, "$working_dir/velvet_$i/contigs.fa";
									open B, ">$source_dir/velvet_$i.contigs.fa";
									while (my $line = <A>) {
										print B $line;
									}
									close A;
									close B;
								}
							}
							system "cp -r $working_dir/best.fasta $source_dir\n";
						}
						if ($DEBUG) { 
							print "rm -rf $working_dir";
						}
						else { 
							rmtree("$working_dir");
							#system "rm -rf $working_dir";
						}
					}
				}
			}
				

		}
		if ($DEBUG) { sleep 1; }
	}
	
	close LOG;
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

sub create_markers { 
	my $orig_file = shift;

	my @return_markers;

	open Z, $orig_file;
	my @sequences;
	my $seq;
	while (my $line = <Z>) { 
		chomp $line;
		unless ($line =~ />/) {
			$seq = $line;
		}
	}
	close Z;
	if (length($seq) < 1300) { 
		return ();
	}

	my $marker_loc = 100;
	my $len = length($seq);
	while ($len >= 200) { 
		my $tmp_marker =  substr($seq,100,100);
		$tmp_marker =~ s/N*$//;
		$tmp_marker =~ s/^N*//;
		my @tmp_markers = split /N+/, $tmp_marker;
		my $longest = shift @tmp_markers;
		foreach my $tmp (@tmp_markers) {
			if (length($tmp) > length($longest)) { $longest = $tmp; }
		}
		$tmp_marker = $longest;
		unless ($tmp_marker) { $tmp_marker = ""; }

		if (length($tmp_marker) > 50) {
			push @return_markers, $tmp_marker;
			push @return_markers, $marker_loc;
		}
		if ($len < 700) { $len = 0; }
		else { 
			$seq = substr($seq,500,$len - 500);
			$len = length($seq);
		}
		$marker_loc += 500;
	}
	return @return_markers;
}


sub has_inversion {
	my $file = shift;
	my @markers = @_;
	return -1 unless ($markers[0]);
	

	my @inverted_markers;
	foreach my $marker (@markers) {
		push @inverted_markers, rev_comp($marker);
	}

	open G, $file;
	my @lines = ();
	my $i = 0;
	while (my $line = <G>) { 
		if ($line) { 
			chomp $line;
			if ($line =~ /^>/) { 
				if ($lines[$i]) { $i++; }
				$lines[$i] = $line;
				$i++;
			}
			else { 
				chomp $line;
				$lines[$i] .= $line;
			}
		}
	}
	close G;
		
	$i = 0;
	my $j = 0;
	my $lastpos = 0;
	my $inversion_count = 0;
	my $good_markers = 0;
	foreach my $line (@lines) {  
		if ($line =~ />/) {
			chomp $file;
			chomp $line;
		}
		else { 
			$line = uc($line);
			my $last_marker_loc;
			for (my $index = 0; $index < $#markers; $index+=2) { 
				my $marker = $markers[$index];
				my $marker_loc = $markers[$index+1];
				if (($line =~ /$marker/) && $last_marker_loc) { 
					my $putative_dist = abs($marker_loc - $last_marker_loc);
					my $actual_dist = abs($-[0] - $lastpos);
					my $discrepency = abs($putative_dist - $actual_dist);

					$inverted_markers[$index] = "XXX";
					if ($lastpos) { 
						if ($discrepency > 3000) { 
							$inversion_count++;
						}
						else { $good_markers++; }
					}
					$lastpos = $-[0];
					$i++;
				}
				$last_marker_loc = $marker_loc;
			}
			foreach my $marker (@inverted_markers) {
				next if ($marker =~ /\d/);
				next if ($marker =~ /^X/);
				if ($line =~ /$marker/) { 
					$j++;
					$inversion_count++;
				}
			}
		}
	}
	unless ($inversion_count) { 
		$good_markers = 0 - $good_markers;
		return $good_markers;
	}

	close G;
	return $inversion_count;
}

		
		

	
sub restore_Ns {
	my $orig_file = shift;
	my $file = shift;
	my $full_line;
	my $return_line;
	my @lines;
	my $i = 0;
	my @starts;
	my @ends;
	if ((!-e $file) || (-z $file)) { 
		open F, $orig_file;
		my $line;
		while (my $tmp = <F>) { $line .= $tmp; } 
		return $line;
	}
	open F, $file;
	open G, $orig_file;
	while (my $line = <G>) { 
		$line = <G>;
		$line = uc $line;
		while ($line =~ /(.{20})N+(.{20})/) {
			my $start = $1;
			my $end = $2;
			push @starts, $start;
			push @ends, $end;
			$line =~ s/.*${start}N+//;
		}
	}
	close G;


	while (my $line = <F>) {
		if ($line) { 
			chomp $line;
			if ($line =~ /^>/) { 
				if ($lines[$i]) { $i++; }
				$lines[$i] = $line;
				$i++;
			}
			else { 
				chomp $line;
				$lines[$i] .= $line;
			}
		}
	}
	close F;

	foreach my $line (@lines) { 
		if ($line) { 
			if ($line =~ /^>/) { $line .= "\n"; }
			else { 
			$line = uc $line;
			while (@starts) { 
				my $start = shift @starts;
				my $end = shift @ends;
				my $rev_start = rev_comp($end);
				my $rev_end = rev_comp($start);
				if ($line =~ /${start}([A]+)${end}/) { 
					my $count = length($1);
					my $ns = "N" x $count;
					$line =~ s/${start}[A]+$end/$start$ns$end/;
				}
				elsif ( $line =~ /${start}([T]+)${end}/) { 
					my $count = length($1);
					my $ns = "N" x $count;
					$line =~ s/${start}[T]+$end/$start$ns$end/;
				}
				elsif ($line =~ /${rev_start}([A]+)${rev_end}/) { 
					my $count = length($1);
					my $ns = "N" x $count;
					$line =~ s/${rev_start}[A]+$rev_end/$rev_start$ns$rev_end/;
				}
				elsif ($line =~ /${rev_start}([T]+)${rev_end}/) { 
					my $count = length($1);
					my $ns = "N" x $count;
					$line =~ s/${rev_start}[T]+$rev_end/$rev_start$ns$rev_end/;
				}
			}
				while ($line =~ /([aA]{100}[aA]*)/) {
					my $as = $1;
					my $len = length($as);
					my $ns = "N" x $len;
					$line =~ s/$as/$ns/;
				}
				while ($line =~ /([tT]{100}[tT]*)/) {
					my $as = $1;
					my $len = length($as);
					my $ns = "N" x $len;
					$line =~ s/$as/$ns/;
				}
				$line .= "\n";
			}
			$return_line .= $line;
		}
	}
	return $return_line;
}

sub rev_comp { 
	my $seq = shift;
	if ($seq =~ /\d/) { return $seq; }

	$seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;

	return $seq;
}

sub build_velveth_cmd {
	my ($low_kmer, $high_kmer, $kmer_step, $dir, $min_contig_length, $length, $exp_coverage, @assembly_reads) = @_;

	my @velveth_cmds;
	my $velvet_cmd = "";

	$velvet_cmd .= $velvet_path . "velveth . $high_kmer -noHash "; 
	#open F, "contigs.noNs.fasta";
	#open G, ">$dir/allreads.withcontigs.fasta";
	#while (my $line = <F>) { 
	#	print G $line;
	#}
	#close F;
	#$velvet_cmd .= "-long $dir/allreads.withcontigs.fasta ";

	my $i = 1;
	my $last_was_fastq = 0;
	foreach my $read (@assembly_reads) { 
		chomp $read;
		if ($read =~ /q$/) { 
			unless ($last_was_fastq) {
				$velvet_cmd .= "-fastq ";
				$last_was_fastq = 1;
			}
		}
		elsif ($last_was_fastq) {
			$velvet_cmd .= "-fasta ";
			$last_was_fastq = 0;
		}
	#	open F, $read;
	#	while (my $line = <F>) { 
	#		print G $line;
	#	}
	#	close F;
		my $isize;
		my $isd;
		if ($read =~ /reads.(\d+)\.-1..*fast./) { 
			$isize = $1;
			$isd = -1;
		}
		elsif ($read =~ /reads.(\d+)\.(\d+)..*fast./) { 
			$isize = $1;
			$isd = $2;
		}

		if ($isize < 20) { 
			if ($num_read_libs == 1) { 
				$velvet_cmd .= "-long $read ";
			}
			else { 
				$velvet_cmd .= "-short $read ";
			}
		}
		else {
			if ($num_read_libs == 1) { 
				$velvet_cmd .= "-longPaired $read ";
			}
			else { 
				$velvet_cmd .= "-shortPaired$i $read ";
				$i++;
			}
		}
	}
	$velvet_cmd .=  "> velveth.out 2>&1 ; "; 
	push @velveth_cmds, $velvet_cmd;
	close G;

	for (my $kmer = $high_kmer; $kmer >= $low_kmer; $kmer -= $kmer_step) { 
		$velvet_cmd = "";
		$velvet_cmd = $velvet_path . "velveth velvet_$kmer $kmer -reuse_Sequences >> velveth.out 2>&1 ;";
		push @velveth_cmds, $velvet_cmd;
	}

	return @velveth_cmds;
}

sub build_velvetg_cmd {
	my ($low_kmer, $high_kmer, $kmer_step, $dir, $min_contig_length, $length, $exp_coverage, @assembly_reads) = @_;

	my $velvet_cmd;
	my @velvet_cmds;
	my $mp = 0;
	if ((!-e "reads.graph2") || (-z "reads.graph2")) {
		$mp = 1;
	}
	else { 
		open MP, "reads.graph2";
		while (my $line = <MP>) {
			if ($line =~ /edge:\s\d+\s\d+\s(\d+)/) {
				my $tmp_mp = $1;
				if (!$mp) { $mp = $tmp_mp; }
				elsif ($tmp_mp < $mp) { $mp = $tmp_mp; }
			}
		}
		close MP;
	}


	for (my $kmer = $high_kmer; $kmer >= $low_kmer; $kmer -= $kmer_step) { 
  		$velvet_cmd .= $velvet_path . "velvetg velvet_$kmer ";
		my $i = 1;
		foreach my $read (@assembly_reads) { 
			chomp $read;
			my $isize;
			my $isd;
			if ($read =~ /reads.(\d+)\.-1..*fast./) { 
				$isize = $1;
				$isd = -1;
			}
			elsif ($read =~ /reads.(\d+)\.(\d+)..*fast./) { 
				$isize = $1;
				$isd = $2;
			}
	
	
			next if ($isize < 20);
			if (($isd < 0) || ($isize > $length)) { 
				$i++;
				next; 
			}
			elsif ($isize >= 20) {
				if ($num_read_libs == 1) { 
					$velvet_cmd .= "-ins_length_long $isize -ins_length_long_sd $isd ";
				}
				else { 
					if ($isize > 1000) { $velvet_cmd .= "-shortMatePaired$i yes "; }
					$velvet_cmd .= "-ins_length$i $isize -ins_length${i}_sd $isd ";
				}
				$i++;
			}
		}
	
		my $working_min_contig_length = int($min_contig_length * 0.01);
		if ($fix) { 
			$min_contig_length = 200; 
			$velvet_cmd .= "-scaffolding no ";
		}
		else { 
			$velvet_cmd .= "-long_mult_cutoff 0 -scaffolding yes -conserveLong yes "; 
		}

			
			#$mp = 5;
			#unless ($keep_low) { $velvet_cmd .= "-exp_cov auto "; }
			$velvet_cmd .= "-exp_cov $exp_coverage "; 
			if ($kmer == $high_kmer) { 
				$velvet_cmd .= "-cov_cutoff $cov_cutoff -min_contig_lgth $working_min_contig_length -min_pair_count $mp -clean yes  > velvetg.out 2>&1 "; 
			}
			else { 
				$velvet_cmd .= "-cov_cutoff $cov_cutoff -min_contig_lgth $working_min_contig_length -min_pair_count $mp -clean yes  >> velvetg.out 2>&1 "; 
			}




		push @velvet_cmds,$velvet_cmd;
		$velvet_cmd = "";
	}
		
	return @velvet_cmds;
}

sub get_NX {
	my ($file, $orig_length, $NX) = @_;

	open M, "$file";
	#while (my $tmp = <M>) { print $tmp; }
	my @contig_lengths;
	my $length = 0;
	while (my $line = <M>) {
#		print "($file)($line)\n";
		chomp $line;

		if ($line =~ /^>/) { 
			if ($length) { 
				push @contig_lengths, $length;
				$length = 0;
			}
		}
		else { 
			$length += length($line); }
#			print "($length)($file)\n";
	}
	push @contig_lengths, $length;

	my @sorted_contig_lengths = sort { $b <=> $a } @contig_lengths;

	my $threshold = $orig_length * $NX;
	if ($NX > 1) { $threshold = int($threshold/100); }
	foreach my $len (@sorted_contig_lengths) {
		$threshold -= $len;
		if ($threshold <= 0) { return $len; }
	}
	return 0;
}


