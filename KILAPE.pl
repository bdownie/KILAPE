#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  KILAPE.pl
#
#        USAGE:  ./KILAPE.pl  
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
#      CREATED:  05/07/13 13:59:44
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use File::Basename;

use Pod::Usage qw(pod2usage);

my $help=0;
my $man=0;

my $kilape_prefix = dirname(File::Spec->rel2abs($0));
my $config_file = "kilape.conf";
my $contig_file = "start.fasta";
my $prefix = "kilape_working";
my $scaffold_only = 0;
my $preprocess_only = 0;
my $BSO = 1;
my $LOCAL_ASSEMBLY = 0;
my $skip_until_after_jellyfish = 0;
my $min_pair_count = 0;
my $abs_min_pair_count = 0;

require "$kilape_prefix/RepARK.pl";

GetOptions('help|h|?'=>\$help, 'man'=>\$man,'c=s'=>\$config_file, 'a=s'=>\$contig_file, '-prepare_lib'=>\$preprocess_only,
		   'scaffold'=>\$scaffold_only, 'o=s'=>\$prefix, 'nojellyfish'=>\$skip_until_after_jellyfish);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

unless (-e $config_file) { 
	print STDERR "No config file given or found\n\tUse -h for help\n"; 
	exit;
}
unless (-e $contig_file) { 
	print STDERR "No configuration file given or found\n\tUse -h for help\n"; 
	exit;
}

open CONFIG, $config_file;
my $cur_lib;
my @libraries;
my %config_hash;
my %sizes_hash;
my @command_list;

# Set default values if not in config file
my $jellyfish_hash_size = 100000000;
my $jellyfish_kmer_size = 17;
my $thread_count = 1;
my $mapper = "BOWTIE2";
my $ramdisk = "override";
my $expected_seq_coverage = "NORMAL";

my @kilape_order;
my @la_order;

while (my $line = <CONFIG>) { 
	chomp $line;
	$line =~ s/\#.*//; 

	next unless ($line);
	if ($line =~ /<LIB=(\S+)>/) {
		$cur_lib = $1;
		push @libraries, $cur_lib;
	}
	else {
		my ($ATTR, $VALUE) = split /=/, $line;
		$VALUE =~ s/\s//g; 
		next unless ($ATTR);
		if ($ATTR eq "SCAFFOLD") { 
			push @kilape_order, $VALUE; 
			push @la_order, 0;
		}
		elsif ($ATTR eq "LA_VELVET") { 
			push @kilape_order, $VALUE; 
			push @la_order, 1; 
		}
		elsif ($ATTR eq "LA_CELERA") { 
			push @kilape_order, $VALUE; 
			push @la_order, 2; 
		}
		elsif ($ATTR eq "LA_CELERA_CLUSTER") { 
			push @kilape_order, $VALUE; 
			push @la_order, 3; 
		}
		elsif ($ATTR eq "LA_GAPFILLER") { 
			push @kilape_order, $VALUE; 
			push @la_order, 4; 
		}
		elsif ($ATTR eq "JELLYFISH_HASH_SIZE") { $jellyfish_hash_size = $VALUE; }
		elsif ($ATTR eq "JELLYFISH_KMER_SIZE") { $jellyfish_kmer_size = $VALUE; }
		elsif ($ATTR eq "MAX_THREADS") { $thread_count = $VALUE; }
		elsif ($ATTR eq "ALIGN") { $mapper = $VALUE; }
		elsif ($ATTR eq "RAMDISK") { $ramdisk = $VALUE; }
		#elsif ($ATTR eq "EXPECTED_SEQUENCE_COVERAGE") { $ramdisk = $VALUE; }
		elsif ($ATTR =~ /_PATH/) { $ENV{$ATTR} = $VALUE . "/"; }
		else { 
			$config_hash{$cur_lib}{$ATTR} = $VALUE;
			if ($ATTR eq "INSERT_SIZE") { 
				if ($sizes_hash{$VALUE}) { 	
					print STDERR "Size: $VALUE for lib: $cur_lib is already in use by $sizes_hash{$VALUE}\n";
					print STDERR "Try increasing the insert size by one and run again.\n";
					exit;
				}
				$sizes_hash{$VALUE} = $cur_lib;
			}
		}
	}
}

#print $ENV{'WGS_PATH'} . "\n"; exit;
my @kilape_libraries;
my @masked_libraries;
mkdir $prefix;
mkdir "$prefix/kilape_libraries";
unless ($scaffold_only) {
	foreach my $lib (@libraries) { 
		my $loc = $config_hash{$lib}{'LOCATION'};
		my $loc2 = $config_hash{$lib}{'LOCATION2'};
		unless ($loc) { 
			$loc = $config_hash{$lib}{'LOCATION1'};
		}
		if (!(-e $loc) || -z $loc) { print STDERR "File $loc doesn't exist.\n\n"; exit; }
		my $orientation = $config_hash{$lib}{"ORIENTATION"};
		push @kilape_libraries, preprocess_reads($loc, $loc2, $orientation, $lib);
	}
	unless ($skip_until_after_jellyfish) { 
		run_jellyfish($jellyfish_hash_size, $jellyfish_kmer_size, $thread_count, @kilape_libraries);
	}
	@masked_libraries = run_kmasker(@kilape_libraries);
}

if ($contig_file) { 
	symlink "../../$contig_file", "$prefix/kilape_libraries/start.fasta";
}
else {
#	run_velvet_on_masked_libraries(@masked_libraries);
	symlink "velvet/contigs.fa", "$prefix/kilape_libraries/start.fasta";
}


my $iteration_count = 1;
my $last_iteration_string = "kilape_libraries";
chdir "$prefix";

my $last_kilape_event = 0;

my $final_assembly_name = "scaffolds.final.fasta";
for (my $i = 0; $i <= $#kilape_order; $i++, $iteration_count++) { 
	my $read_libraries  = $kilape_order[$i];
	my $kilape_event = $la_order[$i];
	mkdir "iteration$iteration_count";
	chdir "iteration$iteration_count";

	if ($iteration_count == 1) { 
		open F, "../$last_iteration_string/start.fasta";
	}
	else { 
		open F, "../$last_iteration_string/$final_assembly_name";
	}

	open G, ">start.fasta";
	while (my $line = <F>) { print G $line; }
	close F;
	close G;


	my @libs = ("0");
	push @libs , split(/,/,$read_libraries);

	my $nomap = 0;
	if ($last_kilape_event > 0) { 
		$nomap = 1; 
	}
	$last_kilape_event = $kilape_event;

	my $wrapper_cmd = build_wrapper_cmd($nomap, $mapper,$kilape_event,$ramdisk, $thread_count, \%config_hash,\@libs);
	print STDERR "Executing $wrapper_cmd\n"; 
	system $wrapper_cmd;
	if ($kilape_event == 1) { 
		rename "scaffolds.final.fasta", "scaffolds.velvet.fasta";
		$final_assembly_name = "scaffolds.velvet.fasta";
	}
	elsif (($kilape_event == 2) || ($kilape_event == 3)) { 
		rename "scaffolds.final.fasta", "scaffolds.celera.fasta";
		$final_assembly_name = "scaffolds.celera.fasta";
	}
	$last_iteration_string = "iteration$iteration_count";
	if (-d "working") {
		print STDERR "Archiving working directory (This might take a while)\n";
		system "ionice -c 3 tar -cf - working | gzip -1 > working.tar.gz";
		system "rm -rf working";
	}
	chdir "..";
}
chdir "..";

symlink "kilape_working/$last_iteration_string/$final_assembly_name", $final_assembly_name;  
		

############################################################################################
#
# Subroutines from here on

sub run_velvet_on_masked_libraries {
	my @libs_to_assemble  = @_;

	my $velveth_cmd = "velveth $prefix/kilape_libraries/velvet $jellyfish_kmer_size";
	foreach my $lib (@libs_to_assemble) { 
		my $velvet_prefix = " -fastq -short ";
		open F, $lib;
		my $line = <F>;
		if ($line =~ /^>/) {
			$velvet_prefix = " -fasta -short ";
		}
		$velveth_cmd .= $velvet_prefix . $lib;
	}
	my $velvetg_cmd = "velvetg $prefix/kilape_libraries/velvet -scaffolding no -exp_cov auto -cov_cutoff auto";

	print STDERR "Executing: $velveth_cmd\n";
	system($velveth_cmd);
	print STDERR "Executing: $velvetg_cmd\n";
	system($velvetg_cmd);
}

sub build_wrapper_cmd {
	my ($do_not_map, $map_prog, $event, $rd, $threads, $conf_hash_ref, $libraries_ref) = @_;
	my %conf_hash = %{ $conf_hash_ref };
	my @libraries = @{ $libraries_ref };

	$map_prog = lc($map_prog);
	my $return_string = "$kilape_prefix/wrapper.pl -a start.fasta -$map_prog -t $threads ";
	if ($do_not_map) { $return_string .= "-nomap "; }
	my $scaffold_libs = 0;
	if ($event == 0) { 
		$return_string .= "-bso "; 
		$scaffold_libs = 1;
	}
	elsif ($event == 1) { $return_string .= "-velvet -ramdisk $rd "; }
	elsif ($event == 2) { $return_string .= "-celera -ramdisk $rd "; }
	elsif ($event == 3) { $return_string .= "-celera -sge -ramdisk $rd "; }
	elsif ($event == 4) { $return_string .= "-gapfiller -ramdisk $rd "; }
	my %sizes_hash;
	for (my $i = 1; $i <= $#libraries; $i++) { 
		my $insert_size = $conf_hash{$libraries[$i]}{'INSERT_SIZE'};

		my $insert_size_sd = $conf_hash{$libraries[$i]}{'INSERT_SIZE_SD'};
		my $mp = $conf_hash{$libraries[$i]}{'MIN_PAIR_COUNT'};
		my $abs_mp = $conf_hash{$libraries[$i]}{'ABS_MIN_PAIR_COUNT'};
		if ($mp) { $return_string .= "-mp $mp "; }
		if ($abs_mp) { $return_string .= "-absmp $abs_mp "; }
		my $unmasked_reads = $libraries[$i];
		$unmasked_reads =~ s/MP/PE/i;
		if (-e "../kilape_libraries/$unmasked_reads.fq") { $unmasked_reads .= ".fq"; }
		elsif (-e "../kilape_libraries/$unmasked_reads.fa") { $unmasked_reads .= ".fa"; }
		else { print STDERR "Couldn't find library $unmasked_reads.fq or .fa\n"; exit; }
		my $masked_reads = $unmasked_reads . ".masked";
		symlink "../kilape_libraries/$unmasked_reads", $unmasked_reads;
		symlink "../kilape_libraries/$masked_reads", $masked_reads;
		$return_string .= "-lib$i $unmasked_reads $masked_reads $insert_size $insert_size_sd $scaffold_libs ";
	}
	$return_string .= "> kilape.out 2>&1";
	return $return_string;
}

sub run_kmasker {
	my @libraries_to_mask = @_;

	my $low;
	my $high;
	my $peak;
	my $coverage;

	($low,$high) = run_calc_threshold();
	#$high *= 2;
	rename "jf_RepARK.kmers","$prefix/kilape_libraries/jf_kilape.kmers";
	rename "jf_RepARK.db","$prefix/kilape_libraries/jf_kilape.db";
	rename "jf_RepARK.histo","$prefix/kilape_libraries/jf_kilape.histo";
	unlink "jf_RepARK.repeat.kmers";
	#print  "mv jf_RepARK.db $prefix/kilape_libraries/jf_kilape.kmers\n";
	#rename "$prefix/jf_RepARK.histo","$prefix/kilape_libraries/jf_kilape.histo";
	
	my $kmasker_cmd = "$kilape_prefix/bin/Kmasker -l $low -u $high -k $prefix/kilape_libraries/jf_kilape.kmers ";
	$kmasker_cmd .= join " ", @libraries_to_mask;
	print STDERR "Executing $kmasker_cmd\n";
	system ($kmasker_cmd);


	my @return_libs;
	foreach my $lib (@libraries_to_mask) {
		if (-e $lib && !-z $lib) { 
			my $purify_reads_cmd = "$kilape_prefix/removeNsfromFastx.pl $lib.k$low-$high -quiet > $lib.masked";
			print STDERR "Purifying: $purify_reads_cmd\n";
			system $purify_reads_cmd;
		#	my $tmp_lib = $lib;
		#	$tmp_lib =~ s/.*\///;
		#	symlink "$tmp_lib.k$low-$high", "$lib.masked";
			push @return_libs, "$lib.masked";
		}
		else { print STDERR "Couldn't find $lib to purify!\n"; }
	}
	return @return_libs;
}
sub preprocess_reads {
	my ($location, $location_rev_read, $orient, $lib_name) = @_;


	my @return_libraries;
	unless ($location && -e $location) { 
		print STDERR "File $location doesn't exist\n";
		exit;
	}
	open LIB, $location;
	my $is_interleaved = 1;
	if ($location_rev_read) { 
		open LIB_REV, $location_rev_read;
		$is_interleaved = 0;
	}
	my $line = <LIB>;
	my $suffix = ".fa";
	my $is_fastq = 0;
	my $acc_header = ">";
	if ($line =~ /^\@/) { 
		$suffix = ".fq"; 
		$is_fastq = 1;
		$acc_header = "@";
	}
	elsif ($line =~ /^>/) { $suffix = ".fa"; }
	else { print STDERR "$location is not FASTA/Q\n\n"; exit; }
	my $new_lib_name = $lib_name;
	$new_lib_name =~ s/MP/PE/i;
	push @return_libraries, "$prefix/kilape_libraries/$new_lib_name$suffix";
	if ($skip_until_after_jellyfish) { 
		return @return_libraries;
	}
	open K_LIB, ">$prefix/kilape_libraries/$new_lib_name$suffix";
	seek(LIB,0,0);
	
	my $index = 1;
	while ($line = <LIB>) {
		my $seq = <LIB>;
		chomp $seq;
		my $qual;
		if ($is_fastq) { 
			<LIB>;
			$qual = <LIB>;
			chomp $qual;
		}
		if ($orient eq "MP") { 
			$seq = rev_comp($seq);
			if ($is_fastq) { 
				$qual = reverse($qual);
			}
		}
		print K_LIB "$acc_header$index\n";
		print K_LIB "$seq\n";
		if ($is_fastq) { print K_LIB "+\n$qual\n"; }
		$index++;
		unless ($is_interleaved) { 
			$line = <LIB_REV>;
			$seq = <LIB_REV>;
			chomp $seq;
			if ($is_fastq) { 
				<LIB_REV>;
				$qual = <LIB_REV>;
				chomp $qual;
			}
			if ($orient eq "MP") { 
				$seq = rev_comp($seq);
				if ($is_fastq) { 
					$qual = reverse($qual);
				}
			}
			print K_LIB "$acc_header$index\n";
			print K_LIB "$seq\n";
			if ($is_fastq) { print K_LIB "+\n$qual\n"; }
			$index++;
		}
	}

	return @return_libraries;
}
		
sub rev_comp { 
	my $seq = shift;
	if ($seq =~ /\d/) { return $seq; }

	$seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;

	return $seq;
}



__END__

=head1 NAME

Run full KILAPE pipeline

=head1 SYNOPSIS

KILAPE.pl [ -c <config file> ] [ -a <contig file> ] [ -scaffold ] [ -nojellyfish ] [ -prepare_lib ]

=head1 ARGUMENTS

=over 4

=item B<-a>

assembly file to be scaffolded/gap-filled (default: start.fasta)

=item B<-c>

config file to be read (default: kilape.conf)

=back 

=head1 OPTIONS

=over 4

=item B<-scaffold>

Perform scaffolding only

=item B<-nojellyfish>

Start kilape without preprocessing libraries or running jellyfish

=item B<-prepare_lib>

Only prepare libraries

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

Please read the README for details

=cut

