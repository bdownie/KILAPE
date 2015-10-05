#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  wrapper.pl
#
#  DESCRIPTION:  Script to execute a single run of KILAPE. 
#			     
#				 Execute with -h flag for help running.
#				 
#				 
#
#       AUTHOR:  Bryan Downie
#  AFFILIATION:  Leibniz Institute for Age Research, Jena Germany
#      VERSION:  1.2
#      CREATED:  04/07/11 14:55:34
#
#        NOTES: This script still needs to be extensively cleaned up and
#		 	 	commented
#===============================================================================

#use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use File::Spec;
use File::Basename;
use Pod::Usage qw(pod2usage);

my $man = 0;
my $help = 0;

#my $bwa_path = "/path/to/bwa";
#my $bowtie_path = "/path/to/bowtie";
my $bwa_path = $ENV{'BWA_PATH'};
my $bowtie_path = $ENV{'BOWTIE_PATH'};
my $bowtie2_path = $ENV{'BOWTIE2_PATH'};
my $wgs_path = $ENV{'WGS_PATH'};
my $velvet_path = $ENV{'VELVET_PATH'};
unless ($bowtie2_path) { $bowtie2_path = ""; }
unless ($bowtie_path) { $bowtie_path = ""; }
unless ($bwa_path) { $bwa_path = ""; }

my $prefix = dirname(File::Spec->rel2abs($0));

# Initializing variables set by options. Many are obsolete.
my $DEBUG;
my $MAX_TASKS = 1000;
#my $final_assembly;
#my $allow_repeats = 0;
#my $strict = 0;
#my $delete_all = 0;
#my $soap = 0;
my $clean = 0;
my $reference = 0;
my $threshold = 0;
my $mp = 2;
my $absmp = 2;
my $bwa = 0;
#my $LOCAL_ASSEMBLY = 1;
my $sub_threads = 0;
my $fgo = 0; # Fill gap only option
my $bso = 0; # Build scaffolds only option
#my $process_unmasked = 0;
#my $map = 0;
my $bowtie = 0;
my $last = 0;
my $longest = 0;
#my $complete = 0;
my $nomap = 0;
my $noindex = 0;
#my $mismatch = 1; # Default number of mismatches per read mapping
my $fix = 0;
#my $trim = 0;
my $ec = 0;
my $ramdisk = "";
#my $l_value;
#my $m_value;
my $high_threshold = 0;
my $psl = 0;
my $s_value = -5;
my $bowtie2 = 0;
my $celera=0;
my $velvet=0;
my $rd_override = 0;
my $sge = 0;
my $ionice = 1;
my $no_local = 0;
my $cov_cutoff;
my $gapfiller = 0;


# Many of these options are obsolete and will be removed with code cleanup.
GetOptions ('d' => \$DEBUG,'a=s' => \$assembly,   't=i' => \$threads, 'fd=i'=>\$first_dir, 'threshold=i'=>\$threshold,'bwa'=>\$bwa,'reference'=>\$reference, 'm=i'=> \$sub_threads, 'bso'=>\$bso, 'fgo'=>\$fgo, 'bowtie'=>\$bowtie, 'last'=>\$last, 'longest'=>\$longest, 'fix' =>\$fix, 'nomap'=>\$nomap, 'noindex'=>\$noindex, 'ec'=>\$ec, 'ramdisk=s'=>\$ramdisk, 'ht=i'=>\$high_threshold, 'psl'=>\$psl,'s'=>\$s_value, 'bowtie2'=>\$bowtie2, 'celera'=>\$celera, 'velvet'=>\$velvet,'sge'=>\$sge, 'absmp=i'=>\$absmp, 'no_local'=>\$no_local, 'cov_cutoff=s'=>\$cov_cutoff, 'gapfiller'=>\$gapfiller,

# These options are so that multiple libraries can be used simultaneously in one run.
# 
# form -libX <unmasked library> <masked library> <insert size> <insert size standard deviation> < 1: use for scaffold/0: use only for local assembly >
# The last element for -libX can be removed because of the -bso option.
'lib1=s{5}'=>\@mainlib,
'lib2=s{5}'=>\@lib2,
'lib3=s{5}'=>\@lib3,
'lib4=s{5}'=>\@lib4,
'lib5=s{5}'=>\@lib5,
'lib6=s{5}'=>\@lib6,
'lib7=s{5}'=>\@lib7,
'lib8=s{5}'=>\@lib8,
'lib9=s{5}'=>\@lib9,
'lib10=s{5}'=>\@lib10,
'lib11=s{5}'=>\@lib11,
'lib12=s{5}'=>\@lib12,
'lib13=s{5}'=>\@lib13,
'lib14=s{5}'=>\@lib14,
'lib15=s{5}'=>\@lib15,
'lib16=s{5}'=>\@lib16,
'lib17=s{5}'=>\@lib17,
'lib18=s{5}'=>\@lib18,
'lib19=s{5}'=>\@lib19,
'lib20=s{5}'=>\@lib20,
'lib21=s{5}'=>\@lib21,
'lib22=s{5}'=>\@lib22,
'lib23=s{5}'=>\@lib23,
'lib24=s{5}'=>\@lib24,
'lib25=s{5}'=>\@lib25,
'lib26=s{5}'=>\@lib26,
'lib27=s{5}'=>\@lib27,
'lib28=s{5}'=>\@lib28,
'lib29=s{5}'=>\@lib29,
'lib30=s{5}'=>\@lib30,

'mp=i'=>\$mp,

'h|help|?'=>\$help,man=>\$man) or pod2usage(2);

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
## If no arguments were given, then allow STDIN to be used only
## if it's not connected to a terminal (otherwise print usage)
#pod2usage("$0: No files given.")  if ((@ARGV == 0) && (-t STDIN));


if ($help) {
	print STDERR "Usage: ";
	exit;
}

# Set defaults (but many can be moved up to initialization section or eliminated).
unless ($mp) { $mp = 2; }
#if ($bso) { $LOCAL_ASSEMBLY = 0; }
else {$last = 1; }
#unless ($exp_cov) { $exp_cov = 50; }
if ($psl) { $noindex = 1; }

# Let the person running know that they forgot to include a -t option
unless ($threads) {
	print STDERR "No number of threads specified. Continuing in 5 seconds single threaded...\n";
	sleep 5;
	$threads = 1;
}
# Let the person running know that they forgot to include a -m option (which will be obsoleted in
# future KILAPE version anyway).
unless ($sub_threads && !$bso) { 
	print STDERR "No subthreads specified, will process one library at a time for PrepareLocalAssemblies\n";
	$sub_threads = $threads;
}

# If skipping mapping step, also skip the indexing step.
if ($nomap) { $noindex = 1; }

# Put all the library information into @extra_libs. This part can also be cleaned up and made
# more flexible, but it works. Will have to reconsider this data structure if users want to
# use more than 10 libraries simultaneously.
($original_reads, $masked_reads, $insert_size, $insert_size_sd) = @mainlib;
if (@lib2) { push @extra_libs, [ @lib2]; }
if (@lib3) { push @extra_libs, [ @lib3]; }
if (@lib4) { push @extra_libs, [ @lib4]; }
if (@lib5) { push @extra_libs, [ @lib5]; }
if (@lib6) { push @extra_libs, [ @lib6]; }
if (@lib7) { push @extra_libs, [ @lib7]; }
if (@lib8) { push @extra_libs, [ @lib8]; }
if (@lib9) { push @extra_libs, [ @lib9]; }
if (@lib10) { push @extra_libs, [ @lib10]; }
if (@lib11) { push @extra_libs, [ @lib11]; }
if (@lib12) { push @extra_libs, [ @lib12]; }
if (@lib13) { push @extra_libs, [ @lib13]; }
if (@lib14) { push @extra_libs, [ @lib14]; }
if (@lib15) { push @extra_libs, [ @lib15]; }
if (@lib16) { push @extra_libs, [ @lib16]; }
if (@lib17) { push @extra_libs, [ @lib17]; }
if (@lib18) { push @extra_libs, [ @lib18]; }
if (@lib19) { push @extra_libs, [ @lib19]; }
if (@lib20) { push @extra_libs, [ @lib20]; }
if (@lib21) { push @extra_libs, [ @lib21]; }
if (@lib22) { push @extra_libs, [ @lib22]; }
if (@lib23) { push @extra_libs, [ @lib23]; }
if (@lib24) { push @extra_libs, [ @lib24]; }
if (@lib25) { push @extra_libs, [ @lib25]; }
if (@lib26) { push @extra_libs, [ @lib26]; }
if (@lib27) { push @extra_libs, [ @lib27]; }
if (@lib28) { push @extra_libs, [ @lib28]; }
if (@lib29) { push @extra_libs, [ @lib29]; }
if (@lib30) { push @extra_libs, [ @lib30]; }
push @extra_libs, [ @mainlib ];

# When multithreaded (obsolete), set number of threads
$MAX_THREADS = $#extra_libs + 1;


if ($ramdisk eq "override") { $rd_override = 1; $ramdisk = ""; } 
# For local assembly, make sure user knows that ramdisk should be used to reduce disk I/O if
# they provide no ramdisk location. 
unless (($bso) || ($rd_override) || ($ramdisk)) { 
	print STDERR "Use of a ramdisk is strongly recommended! If you're sure you don't want to use one, use the option '-ramdisk override'\n";
	exit;
}


# This part will be reworked when removing last argument of -libX.
$max_insert_size = 0;
foreach my $lib_tmp (@extra_libs) { 
	my @lib = @{ $lib_tmp };
	my $lib_size = $lib[2];
	my $lib_sd = $lib[3];
	my $scaffold = $lib[4];
	if (($lib_sd <= 0) && $scaffold) { 
		print STDERR "Can't scaffold pairs without valid insert size standard deviation!\n";
		exit;
	}
	if ($lib_sd < 0) { $lib_sd = 0; }
	$lib_size += $lib_sd;
	if ($scaffold && ($lib_size > $max_insert_size)) { 
		$max_insert_size = $lib_size;
	}
}
print "Max insert size available: $max_insert_size\n";
	
# Obsolete
#unless ($kmer_size) { $kmer_size = 41; }
#unless ($min_contig_length) { $min_contig_length = 200; }
unless ($first_dir) { $first_dir = 1; }

# Can only use one (bwa, bowtie, or bowtie2). Should change this to alert user if they provide
# two programs.
if ($bwa) { $soap = 0; $bowtie = 0;}
if ($bowtie) { $soap = 0; $bwa = 0;}

# If important files are missing, tell user how to use program.
unless ($assembly && $original_reads && $masked_reads) { 
	printUsage();
	exit; 
}

# Tell user if files are missing.
unless (-e $assembly) { 
	print STDERR "Couldn't find $assembly\n"; 
	exit;
}
unless (-e $original_reads) { 
	print STDERR "Couldn't find $original_reads\n";
	exit;
}
unless (-e $masked_reads) { 
	print STDERR "Couldn't find $masked_reads\n";
	exit;
}
# If important files are missing, tell user how to use program.\
# Need to combine with above.
unless (@mainlib) { 
	printUsage();
	exit; 
}


# Initialize file names to be used by KILAPE. All are based on the original assembly name.
# Need to include validation that the assembly name ends with .faXXX (otherwise, it will
# semi-break).
$joined_fasta = $assembly;
$joined_fasta =~ s/\.fa.{0,3}$/.join/;

$sai = $assembly;
$sai =~ s/fasta/sai/;
$sam = $assembly;
$sam =~ s/fasta/sam/;
$filtered_sam = $assembly;
$filtered_sam =~ s/\.fa.{0,3}$/.filtered.sam/;



$graph = $assembly;
$graph =~ s/fasta/graph/;
$prefiltergraph = $assembly;
$prefiltergraph =~ s/\.fa.{0,3}$/.pre.graph/;
$scaff = $assembly;
$scaff =~ s/\.fa.{0,3}$/.scaff/;

# Initialization should be done in a single line here.
my $max_distance = 0;
$max_distance = $insert_size;


# This should be performed later, right before scaffolder_int system call.
$scaffolder_int = "$prefix/bin/scaffolder_int -s BuildScaffolds.graph -o tmp.scaff > scaff.output";
# When running multiple libraries, nice them to keep them from getting in each other's way.
# Should also make sure ionice is used.
my $niceval = 0;

my $assembly_ends = "$assembly.ends";
unless ($noindex) { 
	# Each @contigs entry should contain accession (integer) and sequence of assembly
	# seperated by a single newline (ie, no newlines found in the sequence).
	my @contigs = join_and_renumber($assembly);
	# Overwrite the assembly file. This should be changed in the future so that all work is done on a copy
	# (e.g. symbolic links)
	open F, ">$assembly";
	# If we're building scaffolds only, then we should only index the ends of the contigs.
	if ($bso) { 
		open G, ">$assembly_ends";
	}
	foreach my $contig (@contigs) {
		my ($acc,$fasta) = split /\n/, $contig;
		chomp $acc;
		chomp $fasta;
		print F "$acc\n$fasta\n"; 
		# Output the ends of the contigs to the assembly_ends file.
		# Take as much sequence as the largest insert size used. This value could 
		# possibly be optimized.
		#
		# Note that .a and .b references in the SAM file will be eliminated later.
		if ($bso) { 
			my $len = length($fasta);
			$max_frag_length = $max_insert_size + 200;
			if ($len > (2 * $max_frag_length)) { 
				my $seqa = substr($fasta,0,$max_frag_length);
				my $seqb = substr($fasta,$len - $max_frag_length, $max_frag_length);
				print G "$acc.a\n$seqa\n$acc.b\n$seqb\n";
			}
			else { print G "$acc\n$fasta\n"; }
		}
	}
	close F; 
	if ($bso) { 
		close G;
	}
}

		



# Use Parallel::ForkManager to perform multiple alignments.
my $pm = new Parallel::ForkManager($MAX_THREADS);

my @graph_files;
my @graph_insert_sizes;

# Reporting to the user. Change this to STDERR
$date = localtime;
print "Starting mapping at: $date\n";

# If we're not building scaffolds only, then index the assembly file.
unless ($bso) { $assembly_ends = $assembly; }

my $index;
# Soap is not currently implemented
#if ($soap) { 
#	$index = "2bwt-builder $assembly_ends >& soap.output";
#}
# BWA has seperate indexing options based on the size of the genome.
if ($bwa) { 
	my $assembly_size = -s $assembly;
	if ($assembly_size > 1000000000) { 
		$index = $bwa_path . "bwa index -a bwtsw $assembly_ends >& bwa.output";
	}
	else { 
		$index = $bwa_path . "bwa index -a is $assembly_ends >& bwa.output";
	}
}
# Bowtie options
elsif ($bowtie) { 
	 $index = $bowtie_path . "bowtie-build --quiet $assembly_ends $assembly_ends.bowtie";
}
# Bowtie2 options
elsif ($bowtie2) { 
	 $index = $bowtie2_path . "bowtie2-build --quiet $assembly_ends $assembly_ends.bowtie";
}
# Tell the user what we're doing (maybe more explicit)
if ($index) { 
	print "$index\n";
}

# Don't actually index the assembly if -noindex (or -nomap)
unless ($noindex) { 
	system $index;
}

open F, ">libraries";
my $merge = "$prefix/merge_graph_files.pl $assembly ";
my $max_lib_size = 0;
my $max_lib_graph;
my $at_least_one_scaffold_library = 0;
my $no_adjust = 0;
# For each library passed via -libX
foreach my $lib_tmp  (@extra_libs) { 
	my $aln = "";
	# Ionice option to reduce HD thrashing
	if ($ionice) { 
		$aln = "ionice -c 3 ";
	}
	@lib = @{ $lib_tmp };
	my $unmasked = $lib[0];
	my $masked = $lib[1];

	my $last_unmasked_acc;
	my $last_masked_acc;
	# Check to see if the library is FASTA or FASTQ
	open TMP_FILE, $masked;
	my $tmp_fastq = <TMP_FILE>;

	if ($tmp_fastq =~ /^\@/) {
		# Lib is FASTQ
		# Get the last masked/unmasked accession.
		$last_unmasked_acc = `tail -n 4 $unmasked | head -n 1`;
		$last_masked_acc = `tail -n 4 $masked | head -n 1`;
		$last_unmasked_acc =~ s/\@//; 
		$last_masked_acc =~ s/\@//; 
	}
	elsif ($tmp_fastq =~ /^>/) {
		# Lib is FASTA
		# Get the last masked/unmasked accession.
		$last_unmasked_acc = `tail -n 2 $unmasked | head -n 1`;
		$last_masked_acc = `tail -n 2 $masked | head -n 1`;
		$last_unmasked_acc =~ s/>//; 
		$last_masked_acc =~ s/>//; 
	}
	close TMP_FILE;
	chomp $last_unmasked_acc;
	chomp $last_masked_acc;
	# Make sure accessions are integers (only).
	if ($last_unmasked_acc =~ /\D/) {
		print STDERR "Accessions must have sequential integers as accessions..\n";
		exit;
	}
	# Make sure that the accessions of both masked and unmasked match each other.
	if ($last_unmasked_acc != $last_masked_acc) {
		print STDERR "Final accession numbers don't match $masked $unmasked!\n";
		exit;
	}
	my $lib_size = $lib[2];
	my $scaffold = $lib[4];
	# Can be obsoleted with -bso
	if ($scaffold) {  $at_least_one_scaffold_library = 1; }
	my $lib_sam = $assembly;
	$lib_sam =~ s/\.fa.{0,3}$/.$lib[2].sam/;

	my $lastread = `tail -n 2 $masked | head -n 1 | sed 's/>//'`;
	chomp $lastread;
	#my $soap_threads = $threads;
	# This was used for future blat implementation, but will probably be obsoleted.
	if ($psl) {
		my $lib_psl = $assembly;
		$lib_psl =~ s/\.fa.{0,3}/.$lib[2].psl/;
	}
	# Soap not supported
	#elsif ($soap) { 
	#	my $lib_soap = $lib_sam;
	#	$lib_soap =~ s/sam/soap/;
	#	$aln = "nice -n $niceval soap -D $assembly_ends.index -a $masked -o $lib_soap -r 0 -p $soap_threads -v $mismatch >> soap.$lib[2].output 2>&1; nice -n $niceval cat $lib_soap | $prefix/soap_merge_with_sam.pl $assembly_ends $lastread > $lib_sam "; 
	#}
	elsif ($bwa) { 
		my $lib_sai = $lib_sam;
		$lib_sai =~ s/sam/sai/;
		$aln = "nice -n $niceval " . $bwa_path . "bwa aln -k 0 -n 0 -o 0 -t $threads $assembly_ends $masked  2>> bwa.output |  " . $bwa_path . "bwa samse $assembly_ends - $masked 2>> bwa.output2 > $lib_sam";
	}
	elsif ($bowtie) { 
		# BowtieToKilape necessary because sam output of bowtie is SLOW. Better to use bowtie2
		# Might have to verify quality score for bowtie + FASTQ (see bowtie2 section below)
		$aln = "nice -n $niceval " . $bowtie_path . "bowtie --best -v 0 -f $assembly_ends.bowtie $masked | $prefix/BowtieToKilape.pl  -c $assembly_ends -l $lastread ";
		$aln .= " > $lib_sam";
	}
	elsif ($bowtie2) { 
		open TMP_FILE, $masked;
		my $score_string = "";
		#my $score_string = "--very-sensitive";
		#if ($scaffold) { $score_string = "--score-min L,0,0"; }
		my $tmp_fastq = <TMP_FILE>;
		if ($tmp_fastq =~ /^\@/) { 
			seek(TMP_FILE, 0, 0);
			# If file is fastq, check the first thousand entries to guess which phred offset to use.
			my $phred_option = "--phred64";
			for (my $i = 0; $i < 1000; $i++) { 
				<TMP_FILE>; <TMP_FILE>; <TMP_FILE>;
				$tmp_fastq = <TMP_FILE>;
				chomp $tmp_fastq;
				@qualities = split "", $tmp_fastq;
				# Assume phred64
				if ($tmp_fastq =~ /\^454/) { 
					$phred_option = "--phred33";
				}
				foreach my $qual (@qualities) {
					# If any of the quality scores are below 66 ("B"), guess phred33.
					if (ord($qual) < 66) { $phred_option = "--phred33"; }
				}
			}
			#$aln = "nice -n $niceval $bowtie_path/bowtie2 --score-min L,0,0 --quiet --reorder -p $threads $phred_option -q $assembly_ends.bowtie $masked > $lib_sam";
			$aln = "nice -n $niceval " . $bowtie2_path . "bowtie2 $score_string --quiet --reorder -p $threads $phred_option -x $assembly_ends.bowtie -q $masked > $lib_sam";
		}
		else {
			$aln = "nice -n $niceval " . $bowtie2_path . "bowtie2 $score_string --quiet --reorder -p $threads -x $assembly_ends.bowtie  -f $masked > $lib_sam";
		}
		close TMP_FILE;
	}


	my $lib_graph = $assembly;
	$lib_graph =~ s/\.fa.{0,3}/.$lib[2].graph/;

	push @graph_files, $lib_graph;
	push @graph_insert_sizes, $lib_size;

	if ($lib_size > $max_lib_size) { 
		$max_lib_size = $lib_size; 
		$max_lib_graph = $lib_graph;
	}

	my $final_sam = $lib_sam;
	$final_sam =~ s/sam/final.sam/;

	if ($scaffold) {
		#  Validate was in a previous iteration and may be implemented in the future, but is relatively slow and
		# of arguable value.
		my $vunmasked = $unmasked;
		my $vlib_sam = $lib_sam;
		$vunmasked =~ s/\.fa.{0,3}/.validate.fasta/;
		$vlib_sam =~ s/sam/validate.sam/;
		my $pre_sam = $lib_sam;
		$pre_sam =~ s/sam/pre.sam/;
		#my $vlib_soap = $lib_soap;
	
		if ($bso) { 
			if ($psl) { $final_sam = $lib_sam; }
			else { 
				# Need to restore the proper locations in the SAM when building scaffolds
				$aln .= "; $prefix/bin/RestoreFullContigs -s $lib_sam -c $assembly | $prefix/bin/countSamSites >  $final_sam";
			}
			# Build scaffold graph. 
			#$aln .= "; nice -n $niceval $prefix/bin/filterSamOnCount -s $final_sam -c sites.count | $prefix/bin/BuildGraphFromSAM -c $assembly -i $lib[2] -d $lib[3] > $lib_graph";  
			$aln .= ";cat $final_sam |  $prefix/bin/BuildGraphFromSAM -c $assembly -i $lib[2] -d $lib[3] > $lib_graph";  
		}
		# Redundant. Either scaffold or local assembly
		elsif ($lib[2]) { 
			$aln .= "; nice -n $niceval $prefix/bin/BuildGraphFromSAM -s $lib_sam -c $assembly -i $lib[2] -d $lib[3] > $lib_graph";  
		}
		else { print STDERR "Libraries to be used for scaffolding need insert size and insert size standard deviation\n"; }
		$aln .= " ";
		print F "$lib[0] $masked $final_sam $lib[2] $lib[3]\n";
		$merge .= "$lib_graph $lib[2] ";
	}
	else { 
		$aln =~ s/$lib_sam/$final_sam/;
		$aln .= "; mv $final_sam $lib_sam";
		#print  "$aln\n"; exit;
		print F "$lib[0] $masked $lib_sam $lib[2] $lib[3]\n";
	}
	push @cmds, $aln;
	$niceval+=10;
	if ($niceval > 20) { $niceval = 20; }
}
# Actually run the alignment commands.
foreach my $cmd (@cmds) { 
	sleep 1;
	my $pid = $pm->start and next;
	print "Executing: $cmd\n";
	unless ($nomap) { 
		system ($cmd);
	}
	$pm->finish;
}

# Wait for alignments to finish.
$pm->wait_all_children;
$date = localtime;
print "Mapping finished at: $date\n";
my @graph_edges;

# Future implementation (break contigs on bad alignments).
if ($ec) { 
	my $breakbad_cmd = "$prefix/bin/BreakBadContigs -f libraries  -c $assembly -s $s_value > scaffolds.final.fasta.tmp";
	print ("Executing: $breakbad_cmd");
	system ($breakbad_cmd);
	open F, "scaffolds.final.fasta.tmp";
	open G, ">scaffolds.final.fasta";
	while (my $line = <F>) {
		my $line2 = <F>;
		if (length($line2) >= 200) { print G "$line$line2"; }
	}
	exit;
}


# To rework as $bso
if ($at_least_one_scaffold_library) { 

	unless ($fix) { 
		$merge_output = $prefiltergraph;
		$merge_output =~ s/pre.graph/merge.graph/;
	
		$scaffolder_int = "$prefix/bin/scaffolder_int -s $graph -o $scaff -f $assembly > scaff.output";
		$merge .= "> $merge_output";
		print "Executing $merge\n";
			# Merge graph files from aligning multiple libraries (weight distance based on pairing)
			system $merge;
			if (-z $merge) { 
				print STDERR "Failed to merge graph files!\n";
				exit;
			}
	close F;
}

# Filter the graph based on $absmp (default to 2). Contig connections need at least 2 pairs to pass filter.
open F, "BuildScaffolds.graph";
my @graph_lines;
my @graph_line_counts;
while (my $line = <F>) { 
	$line .= <F>; 
	$line .= <F>;
	$line .= <F>;
	if ($line =~ /^edge:\s\S+\s\S+\s(\d+)/) {
		my $count = $1; 
		if ($absmp <= $count) { 
			push @graph_line_counts, $count;
			push @graph_lines, $line;
		}
	}
}
close F;

# Obsolete.
if ($high_threshold == -1) {
	my @graph_line_counts_sorted = sort { $a <=> $b } @graph_line_counts;
	$high_threshold = 2* $graph_line_counts_sorted[int($#graph_line_counts_sorted * 0.5)];
}

# Obsolete.
open G, ">BuildScaffolds.filt.graph";
for (my $i = 0; $i < $#graph_lines; $i++) {
	if ($high_threshold) {
		if ($graph_line_counts[$i] <= $high_threshold) {
			print G $graph_lines[$i];
		}
	}
	else {
		print G $graph_lines[$i];
	}
}
close G;

# Will need to look into this. Should be obsolete
if ($no_adjust) { 
	$scaffolder_int = "$prefix/bin/scaffolder_int -s BuildScaffolds.filt.graph -o $scaff -f $assembly > scaff.output";
	print "Executing $scaffolder_int\n";
	system "$scaffolder_int";
	system "cp BuildScaffolds.filt.graph $graph";
}
elsif ($at_least_one_scaffold_library) { 
	print "\n";
	if ($threshold) { 
		$mp = $threshold;
	}
	else { 
		unless ($total_contigs) { 
			$total_contigs = 0;
			open F, $assembly;
			while (my $line = <F>) { 
				if ($line =~ />/) { $total_contigs++;}
			}
			$fewest_scaffold_number = $total_contigs;
		}
		$best_i = 0;
		if ($longest) { 
			open F, $assembly;
			while (my $line = <F>) { 
				my $line2 = <F>; 
				$line =~ s/>//; 
				chomp $line;
				next if ($line =~ /\D/);
				$length_array[$line] = length($line2);
			}
			close F;
		}
		$max_evaluation_size = $mp + 1;
		#print "($mp)\n";
		
		# Choose the mp value which gives the fewest number of scaffolds (see adjust_graph documentation)
		for (my $i = $mp; $i <= $max_evaluation_size; $i++) { 
			my $total_length = 0;
			# adjust_graph removes discrepencies in the scaffolding graph.
			my $adjust_graph = "$prefix/bin/adjust_graph -g BuildScaffolds.filt.graph -p $i -f $assembly > $graph";
			print "$adjust_graph\n";
			system $adjust_graph;
			if (-z $graph) { 
				print STDERR "Empty graph file!\n";
				$i = $max_evaluation_size + 1; 
			}
			else { 
	
				print "$scaffolder_int\n";
				system "$scaffolder_int";
				open F, $scaff;
				my $contigs = $total_contigs;
				my $scaffolds = 0;
				while (my $line = <F>) { 
					@a = split /\s+/, $line;
					if ($longest) { 
						foreach my $contig (@a) { 
							if ($length_array[$contig]) { 
								$total_length += $length_array[$contig];
							}
						}
					}
					else { 
						$scaffolds++;
						$contigs -= $#a + 1;
					}
				}
				close F;
				$scaffolds += $contigs;
				if ($longest) { 
					if ($total_length > $best_total_length) { 
						$best_i = $i;
						$best_total_length = $total_length;
						$max_evaluation_size++;
					}
				#	print "($i)($total_length)($best_total_length)\n";
				}
				else {
				#	print "($i)($fewest_scaffold_number)($scaffolds)\n";
					if (($scaffolds <= $fewest_scaffold_number ) || (($i == 1) && $fewest_scaffold_number )) {
						$best_i = $i;
						$fewest_scaffold_number = $scaffolds;
						$max_evaluation_size++;
					}
				}
				sleep 1;
				print "here\n";
			}
		}
		unless ($best_i) { 
	#	if (-z "BuildScaffolds.filt.graph" || !(-e "BuildScaffolds.filt.graph")) { 
			open F, "start.fasta";
			open G, ">scaffolds.final.fasta";
			while (my $line = <F>) {
				print G $line;
			}
			close F;
			close G;
			return 1;
		}
		$threshold = $best_i;
		if ($mp == -1) { $mp = $threshold; }
	}


	print "Using $threshold\n";

	# This could be reworked by just keeping a temporary file of the "best". But that can be done later.
	my $adjust_graph = "$prefix/bin/adjust_graph -g BuildScaffolds.filt.graph -p $threshold -f $assembly > $graph";
	print STDERR "Executing $adjust_graph\n";
	system "$adjust_graph";

	print STDERR "Executing $scaffolder_int\n";
	system "$scaffolder_int";
}

# Extra curly brace at end? Need to fix whitespace
}


# Obsolete-ish. This now just changes working to working_X. Can be brought into wrapper.pl
$create_working_dir = "$prefix/createWorkingDir.pl -c $assembly"; 

if ($at_least_one_scaffold_library) {
	$create_working_dir .= " -s $scaff"; 
	unless ($last) { $create_working_dir .= " -f"; }
}
print "Executing $create_working_dir\n";
system "$create_working_dir ";


# Need better variable name than finalfinish. Also, remove sub_threads possibility.
if ($sub_threads) { 
	$finalfinish = "$prefix/bin/PrepareLocalAssemblies -c $assembly -t 1 -m $sub_threads";
}
else { 
	$finalfinish = "$prefix/bin/PrepareLocalAssemblies -c $assembly -t $threads";
}

if ($at_least_one_scaffold_library) { $finalfinish .= " -i $scaff -g $graph"; }
if (!$bso) { 
	$finalfinish .= " -r libraries"; 
	if ($fgo ) { $finalfinish .= " -f"; }
}
# If -bso, PrepareLocalAssemblies will look for overlap and put appropriate # of Ns between scaffolds.
else { $finalfinish .= " > scaffolds.final.fasta"; }
if (-e "contigs.orphan.fasta") { system "rm contigs.orphan.fasta"; }

my $time = localtime;
print "[$time] Executing $finalfinish\n";
system "$finalfinish";
# Put if !-bso here, exit
if ($no_local) { exit; }
# wgs-assembler (supported)
if ($celera) { 
	my $local_assemble;
	$count = `cat scaffolds.count`;
	chomp $count;
	$count--;
	unless ($count) { print STDERR "No count file found.\n"; exit; }
	# Each assembly can be distributed across a cluster (not completely supported)
	if ($sge) { 
		my $cur_dir = `pwd`; chomp $cur_dir;
		$local_assemble = "echo \"ionice -c 3 $prefix/LocalAssemble.celera.cluster.pl -cluster $cur_dir -fd $first_dir -dir_count $count";
		if ($ramdisk) { $local_assemble .= " -ramdisk $ramdisk"; }
		$local_assemble .= " \"";
	}
	else { 
		$local_assemble_threads = $threads;
		#$local_assemble = "$prefix/LocalAssemble.celera.pl -dir_count $count -t $local_assemble_threads -fd $first_dir -contigs ";
		$local_assemble = "$prefix/LocalAssemble.celera.pl -dir_count $count -t $local_assemble_threads -fd $first_dir ";
		if ($ramdisk) { $local_assemble .= " -ramdisk $ramdisk"; }
	}
	$time = localtime;
	print STDERR "[$time] Executing $local_assemble\n";
	system "$local_assemble"; 
	# Brings all the successful assemblies together into resulting file.
	my $consolidate = "$prefix/BuildAndConsolidateScaffolds.pl -dir_count $count > scaffolds.final.fasta";
	print STDERR "Executing $consolidate\n";
	system "$consolidate";
}
elsif ($gapfiller) { 
	my $local_assemble;
	$count = `cat scaffolds.count`;
	chomp $count;
	$count--;
	unless ($count) { print STDERR "No count file found.\n"; exit; }
	$local_assemble_threads = $threads;
	$local_assemble = "$prefix/LocalAssemble.gapfiller.pl -dir_count $count -t $local_assemble_threads -fd $first_dir ";
	if ($ramdisk) { $local_assemble .= " -ramdisk $ramdisk"; }
	$time = localtime;
	print STDERR "[$time] Executing $local_assemble\n";
	system "$local_assemble"; 
	# Brings all the successful assemblies together into resulting file.
	my $consolidate = "$prefix/BuildAndConsolidateScaffolds.pl -dir_count $count > scaffolds.final.fasta";
	print STDERR "Executing $consolidate\n";
	system "$consolidate";
}
# Velvet section (working, but has some issues.)
elsif ($velvet || !$bso) {
	$count = `cat scaffolds.count`;
	chomp $count;
	$count--;
	unless ($count) { print STDERR "No count file found.\n"; exit; }

	if ((-e "contigs.orphan.fasta") && (!-z "contigs.orphan.fasta") && !$bso) { 
		my $dist_orphans = "$prefix/bin/DistributeOrphanedContigs -i $scaff -c contigs.orphan.fasta -g orphan.graph  > contigs.unused.fasta";
		print STDERR "Executing: $dist_orphans\n"; 
		unless ($DEBUG) { system $dist_orphans; }
	}
	$local_assemble_threads = $threads;
	$local_assemble = "$prefix/LocalAssemble.velvet.pl -dir_count $count -t $local_assemble_threads -fd $first_dir ";
	if ($cov_cutoff) { $local_assemble .= " -cov_cutoff $cov_cutoff"; }
	if ($ramdisk) { $local_assemble .= " -ramdisk $ramdisk"; }
	if ($bwa) { $local_assemble .= " -bwa"; }
	elsif ($bowtie) { $local_assemble .= " -bowtie"; }
	elsif ($bowtie2) { $local_assemble .= " -bowtie2"; }
	if ($reference) { 
		$local_assemble .= " -reference"; 
	}
	#	if ($no_validate) { $local_assembly .= " -no_validate"; }
	#	if ($l_value) { $local_assembly .= " -l $l_value"; }
	#	if ($m_value) { $local_assembly .= " -l $m_value"; }
	#	if ($ec) { $local_assemble .= " -ec"; }
	#	if ($last) { $local_assemble .= " -last"; }
	#	if ($bso) { $local_assemble .= " -bso"; }
	#	else { $local_assemble .= " -c"; }
	#	if ($clean) { $local_assemble .= " -clean"; }
		#elsif ($bwa) { $local_assemble .= " -bwa"; }
		#elsif ($bowtie) { $local_assemble .= " -bowtie"; }
		#if ($map) { $local_assemble .= " -map"; }
		#if ($final_assembly) { $local_assemble .= " -final_assembly"; }
		#if ($trim) { $local_assemble .= " -trim"; }
	#}
	$consolidate = "$prefix/BuildAndConsolidateScaffolds.pl -dir_count $count > scaffolds.final.fasta";



	print STDERR "Executing $local_assemble\n";
	system "$local_assemble"; 
	print STDERR "Executing $consolidate\n";
	system "$consolidate";
	if ($bso) { 
		system "cat contigs.orphan.fasta >> scaffolds.final.fasta";
	}
	# Obsolete
	elsif ($fgo) { 
		system "cat contigs.orphan.fasta >> scaffolds.final.fasta";
	}
	# Obsolete
	elsif (!$last) { 
		system "cat contigs.orphan.fasta >> scaffolds.final.fasta";
	}
}



# Subroutines

# Remove newlines from fasta sequence and change accessions to integers.
sub join_and_renumber { 
	my $assembly_file = shift;
	open F, $assembly_file;
	my @contigs;
	my $contig;
	while (my $line = <F>) { 
		chomp $line;
		if ($line =~ /^>/) { 
			if ($contig) { 
				$contig .= "\n";
				push @contigs, $contig;
				$contig = "";
			}
		}
		else { 
			$contig .= $line;
		}
	}
	$contig .= "\n";
	push @contigs, $contig;
	close F;

	@contigs = sort { length $b <=> length $a } @contigs;
	my $i = 1;
	for (my $j = 0; $j <= $#contigs; $j++) { 
		if ($contigs[$j]) { 
			$contigs[$j] = ">$i\n" . $contigs[$j];
			$i++;
		}
	}

	return @contigs;
}

# Must be updated
#sub printUsage {
#	print STDERR "Usage: wrapper.pl -a <assembly fasta> lib<#> <unmasked> <masked> <insert size> <insert size stdev> <use library for scaffold> [ -no_local_assemble ] [ -exp_cov <expected coverage> ] [ -t <threads>] \n";
#}

__END__

=head1 NAME

sample - Using GetOpt::Long and Pod::Usage

=head1 SYNOPSIS

sample [options] -a <file> -lib1 [file ...]

 Required:
   -a				assembly file to be scaffolded/gap-filled
 Options:
   -t				simultaneous threads to use (e.g. bowtie/bwa)
   -m				number of sub-threads available to use during local assembly (e.g. by velvet)
   -no_local		Don't perform local assembly/gap-filling
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut

