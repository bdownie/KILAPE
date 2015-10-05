#!/usr/bin/perl
 
# Description
#
# Autor: Philipp Koch
# modified date: 11.04.2012

# TODO
#
# at the end warnings off
#
######

####
# HowTo: extending this script by one additional entry:
# - add new empty variable at the beginning section : my $var = "";
#
# - add the prompting and checking at your favourite position in the 'create-section'
#
# - add new variable at the beginning of the 'modify-section' to parse the current value out of the existing config file: my $var_old = $Config{'var'};
# - add the prompting and checking at your favourite position in the 'modify-section'
#
# - add one line for printig the key-value pair to the last section: print CONF "var = $var\n";
# 
# - check for bugs
####


use strict;
use warnings;
use File::Basename;
use File::Spec;
use POSIX qw/strftime/;


use FindBin;
use lib "$FindBin::RealBin/Modules";
use KmScTools;


my $base_name = basename($0);
my $usage = <<EOF;
USAGE: $base_name <project name>
  This script helps to create or modify a configuration file for the
  project <project name>.
INPUT:	
	project name   	a custom name for the project (no spaces)
OUTPUT:
	config file   	this file is required for the main program PipelineX.pl
EOF

if ($#ARGV < 0){die $usage;}

my $project_name = $ARGV[0];

#my $fastq = $ARGV[0];
#my $basename = $fastq;
#$basename =~ s/\.fastq//;
my $config_file = "$project_name.conf";


my $decision1 = 0;
my $decision2 = 0;
my $decision3 = 0;
my $decision4 = 0;

my $qual = "";
my $procs = "";
my $k = "";
#!rm!# my $seqsize = "";
my $estgenomesize = "";
my $sensitivity = "";

my $dict = "";
my $jf_path = "";
my $gp_path = "";

my @read_libs;
my $num_of_read_libs = "";

my $assembly_file = "";
my $bwa_path = "";

my $soap_path = "";
my $bwt_path = "";

my $backup = "";

## get the path where this program and all related programs are located
my $program_dir = dirname(File::Spec->rel2abs($0))."/";

## check if there exists a config file
if (-e "$config_file"){
	print "a configuration file was found: $config_file\nwhat do you want to do now? \n";
	while ( !($decision1 =~ m/^[umc]$/) ){
		$decision1 = &promptUser("use it (u), modify it (m) or create a new file and overwrite the old one (c)");
	}
}


###############################################
##
## the case, if we want to use the existing config file
##
###############################################

if ($decision1 eq "u"){
	print "keeping the existing config file\n";
}


###############################################
##
## the case, if we want to modify the existing config file
##
###############################################

if ( $decision1 eq "m" ){

	print "Note:\n";
	print "In the following steps all values in square brackets [] represent the content of the existing config file.\n";
	print "In case of a jellyfish and gnuplot execution file was found, these paths will also be provided in [].\n"; 


	## first lets do a backup of the existing config file
	system ("cp $config_file $config_file.bak");
	print "backup of existing file $config_file to $config_file.bak done\n(NOTE: next call of Configure.pl with 'm' will overwrite this .bak again!)\n";
	
	
	### parse the config file
	my %Config;
	&parse_config_file ($config_file, \%Config);

	### asign all values from the config file to the local variables
#	my $qual_old = $Config{'qual'};
	my $procs_old = $Config{'procs'};
	my $k_old = $Config{'k'};
#!rm!#	my $seqsize_old = $Config{'seqsize'};
	my $estgenomesize_old = $Config{'estgenomesize'};
	my $sensitivity_old = $Config{'sensitivity'};

	my $jf_path_old = $Config{'jf_path'};

	my $dict_old = $Config{'dict'};
	my $gp_path_old = $Config{'gp_path'};
	
	my $assembly_file_old = $Config{'assembly_file'};
#	my $bwa_path_old = $Config{'bwa_path'};
	
	my $soap_path_old = $Config{'soap_path'};
	my $bwt_path_old = $Config{'2bwt_path'};
		
	my @read_libs_old;
	for (my $i= 0; $i <= 20; $i++){
		if ($Config{"read_lib_$i"}){
			$read_libs_old[$i] = $Config{"read_lib_$i"};
		}	
	}
	
	### start asking the user for changes
	
	## ask the user for the jellyfish path
	my $jf_test = "";
	while ( !((-e $jf_path) && ($jf_path =~ m/^(\/.*jellyfish)$/)) ) { 
		$jf_test = `which jellyfish`;
		chomp $jf_test;
		if ($jf_test =~ m/^(\/.*jellyfish)$/){ # found a jellyfish installation
			$jf_path = &promptUser("absolute path to your jellyfish execution file", $1);
		}
		else { # no jellyfish found
			$jf_path = &promptUser("absolute path to your jellyfish execution file", $jf_path_old);
		}
	}
	
	
	## ask if the user want to set the gnuplot path
	while ( !($decision2 =~ m/^[yn]$/) ){
		$decision2 = &promptUser("do you want to set the gnuplot path? (y/n)", "y");
	}	
	if ($decision2 eq "y"){	
		my $gp_test = "";
		while ( !((-e $gp_path) && ($gp_path =~ m/^(\/.*gnuplot)$/)) ) { 
			$gp_test = `which gnuplot`;
			chomp $gp_test;
			if ($gp_test =~ m/^(\/.*gnuplot)$/){ # found a gnuplot installation
				$gp_path = &promptUser("absolute path to your gnuplot execution file", $1);
			}
			else { # no gnuplot found
				$gp_path = &promptUser("absolute path to your gnuplot execution file", $gp_path_old);
			}
		}
	}


	## ask if the user wants to provide a kmer dictionary file
	while ( !($decision3 =~ m/^[yn]$/) ){
		$decision3 = &promptUser("do you want to provide a kmer dictionary file? (y/n)", "n");
	}
	if ($decision3 eq "y"){	
		while ( ! -e $dict) { 
			$dict = &promptUser("kmer dictionary file", $dict_old);
		}
	}
	

	
	print "the existing config file contains information about the following read libraries:\n";

	for (my $i = 0; $i < @read_libs_old; $i++){
		print "\t".$read_libs_old[$i]."\n";
	}
	
	while ( !($decision4 =~ m/^[yn]$/) ){
		$decision4 = &promptUser("do you want to modify them? (y/n)", "n");
	}
	
	if ($decision4 eq "y"){
		print "the existing entries of all read libraries will be deleted.\n";
		print "they have to be set again like in the following style \"<reads> <pe|mp> <IS> <SD> <NT>\".\n";
		print "with\t<reads>\tfastq file including read pairs preprocessed by \"PreProcessReadLibs.pl\"\n";
		print "\t<pe|mp>\tare these pairs 'paired end' (pe) or 'mate pairs' (mp)?\n";
		print "\t<IS>\taverage insert size\n";
		print "\t<SD>\tstandard deviation of the insert size[in bp]\n";
		print "\t<NT>\ttotal number of nucleotides\n";
		
		while (!($num_of_read_libs =~ m/^\d+$/)){
			$num_of_read_libs = &promptUser("how many read libraries will you provide in the following? (max. 20)", "1");
		}
	
		for (my $i = 0; $i < $num_of_read_libs; $i++){
			while ( !($read_libs[$i] =~ m/.+\s+(mp|pe)\s+\d+\s+\d+\s+\d+/) ){
				$read_libs[$i] = &promptUser("please enter the string for read library ".($i+1));
			}
		}
	}
	else {
		@read_libs = @read_libs_old;
	}

#$$$$$$$ SCAFFOLDING STUFF
	## ask for the assembly file needed for scaffolding
	while ( (! -e $assembly_file) && ($assembly_file ne "none")) { 
		$assembly_file = &promptUser("initial assembly file for the scaffolding step", $assembly_file_old);
	}
	
# currently there is no bwa support	
#	## ask the user for the bwa path
#	my $bwa_test = "";
#	while ( !((-e $bwa_path) && ($bwa_path =~ m/^(\/.*bwa)$/)) ) { 
#		$bwa_test = `which bwa`;
#		chomp $bwa_test;
#		if ($bwa_test =~ m/^(\/.*bwa)$/){ # found a bwa installation
#			$bwa_path = &promptUser("absolute path to your bwa execution file", $1);
#		}
#		else { # no bwa found
#			$bwa_path = &promptUser("absolute path to your bwa execution file", $bwa_path_old);
#		}
#	}
	
##	## ask the user for the soap path
	my $soap_test = "";
	while ( !((-e $soap_path) && ($soap_path =~ m/^(\/.*soap)$/)) ) { 
		$soap_test = `which soap`;
		chomp $soap_test;
		if ($soap_test =~ m/^(\/.*soap)$/){ # found a soap installation
			$soap_path = &promptUser("absolute path to your soap execution file", $1);
		}
		else { # no bwa found
			$soap_path = &promptUser("absolute path to your soap execution file", $soap_path_old);
		}
	}
	
##	## ask the user for the 2bwt-builder path
	my $bwt_test = "";
	while ( !((-e $bwt_path) && ($bwt_path =~ m/^(\/.*2bwt-builder)$/)) ) { 
		$bwt_test = `which 2bwt-builder`;
		chomp $bwt_test;
		if ($bwt_test =~ m/^(\/.*2bwt-builder)$/){ # found a 2bwt-builder installation
			$bwt_path = &promptUser("absolute path to your 2bwt-builder execution file", $1);
		}
		else { # no 2bwt-builder found
			$bwt_path = &promptUser("absolute path to your 2bwt-builder execution file", $bwt_path_old);
		}
	}
	
	
#$$$$$$$		

# currently not needed
#	while (!($qual =~ m/^\d+$/)){
#		$qual = &promptUser("quality score for error correction", $qual_old);
#	}

	while (!($procs =~ m/^\d+$/)){
		$procs = &promptUser("number of processes", $procs_old);
	}
	
	while (!($k =~ m/^\d+$/)){
		$k = &promptUser("kmer length", $k_old);	
	}
	
#!rm!#	while (!($seqsize =~ m/^\d+$/)){
#!rm!#		$seqsize = &promptUser("number of bases of the input fasta file (in bp)", $seqsize_old);	
#!rm!#	}
	
	while (!($estgenomesize =~ m/^\d+$/)){
		$estgenomesize = &promptUser("estimated number of bases of the genome size (in bp)", $estgenomesize_old);	
	}
	
	while ( !(($sensitivity =~ m/^[01]\.\d{2}$/) && ($sensitivity < 1.01) && ($sensitivity > 0.00)) ){
		$sensitivity = &promptUser("sensitivity of the threshold computation algorithm (range 0.01..1.00)", $sensitivity_old);	
	}

	
}


###############################################
##
## the case, if we want to create a new config file
##
###############################################

if ( !(-e "$config_file") || $decision1 eq "c"){
	
	## ask the user for the jellyfish path
	my $jf_test = "";
	while ( !((-e $jf_path) && ($jf_path =~ m/^(\/.*jellyfish)$/)) ) { 
		$jf_test = `which jellyfish`;
		chomp $jf_test;
		if ($jf_test =~ m/^(\/.*jellyfish)$/){ # found a jellyfish installation
			$jf_path = &promptUser("absolute path to your jellyfish execution file", $1);
		}
		else { # no jellyfish found
			$jf_path = &promptUser("absolute path to your jellyfish execution file");
		}
	}
	
	
	## ask if the user want to set the gnuplot path
	while ( !($decision2 =~ m/^[yn]$/) ){
		$decision2 = &promptUser("do you want to set the gnuplot path? (y/n)", "y");
	}	
	if ($decision2 eq "y"){	
		my $gp_test = "";
		while ( !((-e $gp_path) && ($gp_path =~ m/^(\/.*gnuplot)$/)) ) { 
			$gp_test = `which gnuplot`;
			chomp $gp_test;
			if ($gp_test =~ m/^(\/.*gnuplot)$/){ # found a gnuplot installation
				$gp_path = &promptUser("absolute path to your gnuplot execution file", $1);
			}
			else { # no gnuplot found
				$gp_path = &promptUser("absolute path to your gnuplot execution file");
			}
		}
	}


	## ask if the user wants to provide a kmer dictionary file
	while ( !($decision3 =~ m/^[yn]$/) ){
		$decision3 = &promptUser("do you want to provide a kmer dictionary file? (y/n)", "n");
	}
	if ($decision3 eq "y"){	
		while ( ! -e $dict) { 
			$dict = &promptUser("kmer dictionary file");
		}
	}


	print "please provide at least one read library like in the following style \"<reads> <pe|mp> <IS> <SD> <NT>\"\n";
	print "with\t<reads>\tfastq file including read pairs preprocessed by \"PreProcessReadLibs.pl\"\n";
	print "\t<pe|mp>\tare these pairs 'paired end' (pe) or 'mate pairs' (mp)?\n";
	print "\t<IS>\taverage insert size\n";
	print "\t<SD>\tstandard deviation of the insert size[in bp]\n";
	print "\t<NT>\ttotal number of nucleotides\n";
	
	while (!($num_of_read_libs =~ m/^\d+$/)){
		$num_of_read_libs = &promptUser("how many read libraries will you provide in the following? (max. 20)", "1");
	}
	
	for (my $i = 0; $i < $num_of_read_libs; $i++){
		while ( !($read_libs[$i] =~ m/.+\s+(mp|pe)\s+\d+\s+\d+/) ){
			$read_libs[$i] = &promptUser("please enter the string for read library ".($i+1));
		}
	}
	
#$$$$$$$ SCAFFOLDING STUFF
	## ask for the assembly file needed for scaffolding
	while ( (! -e $assembly_file) && ($assembly_file ne "none")) { 
		$assembly_file = &promptUser("initial assembly file for the scaffolding step", "none");
	}

# currently there is no bwa support
#	## ask the user for the bwa path
#	my $bwa_test = "";
#	while ( !((-e $bwa_path) && ($bwa_path =~ m/^(\/.*bwa)$/)) ) { 
#		$bwa_test = `which bwa`;
#		chomp $bwa_test;
#		if ($bwa_test =~ m/^(\/.*bwa)$/){ # found a bwa installation
#			$bwa_path = &promptUser("absolute path to your bwa execution file", $1);
#		}
#		else { # no bwa found
#			$bwa_path = &promptUser("absolute path to your bwa execution file");
#		}
#	}
	
##	## ask the user for the SOAP path
	my $soap_test = "";
	while ( !((-e $soap_path) && ($soap_path =~ m/^(\/.*soap)$/)) ) { 
		$soap_test = `which soap`;
		chomp $soap_test;
		if ($soap_test =~ m/^(\/.*soap)$/){ # found a bwa installation
			$soap_path = &promptUser("absolute path to your soap execution file", $1);
		}
		else { # no soap found
			$soap_path = &promptUser("absolute path to your soap execution file");
		}
	}
	
##	## ask the user for the 2bwt-builder path
	my $bwt_test = "";
	while ( !((-e $bwt_path) && ($bwt_path =~ m/^(\/.*2bwt-builder)$/)) ) { 
		$bwt_test = `which 2bwt-builder`;
		chomp $bwt_test;
		if ($bwt_test =~ m/^(\/.*2bwt-builder)$/){ # found a 2bwt-builder installation
			$bwt_path = &promptUser("absolute path to your 2bwt-builder execution file", $1);
		}
		else { # no soap found
			$bwt_path = &promptUser("absolute path to your 2bwt-builder execution file");
		}
	}
#$$$$$$$

# currently not needed	
#	while (!($qual =~ m/^\d+$/)){
#		$qual = &promptUser("quality score for error correction", "10");
#	}

	while (!($procs =~ m/^\d+$/)){
		$procs = &promptUser("number of processes", "1");
	}
	
	while (!($k =~ m/^\d+$/)){
		$k = &promptUser("kmer length", "19");	
	}
	
#!rm!#	while (!($seqsize =~ m/^\d+$/)){
#!rm!#		$seqsize = &promptUser("number of bases of the input fasta file (in bp)");	
#!rm!#	}
	
	while (!($estgenomesize =~ m/^\d+$/)){
		$estgenomesize = &promptUser("estimated number of bases of the genome size (in bp)");	
	}
	
	while ( !(($sensitivity =~ m/^[01]\.\d{2}$/) && ($sensitivity < 1.01) && ($sensitivity > 0.00)) ){
		$sensitivity = &promptUser("sensitivity of the threshold computation algorithm (range 0.01..1.00)", "0.01");	
	}
	
}



####################################
##
## if we got new input from the user, write the config file
##
####################################
if ( $decision1 ne "u" ){
	open (CONF, ">$config_file") or die "can not create file $config_file: $!";
		print CONF "# file created on ".strftime('%a %b %d %H:%M:%S %Y',localtime)."\n";
		print CONF "project_name = $project_name\n";
		$program_dir =~ s/\/$//;
		print CONF "program_dir = $program_dir\n";
		print CONF "jf_path = $jf_path\n";	
		print CONF "gp_path = $gp_path\n";
# curr no supp	print CONF "qual = $qual\n";
		print CONF "procs = $procs\n";
		print CONF "k = $k\n";
#!rm!#		print CONF "seqsize = $seqsize\n";
		print CONF "dict = $dict\n";
		print CONF "estgenomesize = $estgenomesize\n";
		print CONF "sensitivity = $sensitivity\n";
		
		for (my $i = 0; $i < @read_libs; $i++){
			print CONF "read_lib_$i = ".$read_libs[$i]."\n";
		}
		print CONF "assembly_file = $assembly_file\n";
# curr no supp	print CONF "bwa_path = $bwa_path\n";
		print CONF "soap_path = $soap_path\n";
		print CONF "2bwt_path = $bwt_path\n";
		
	close CONF;
}

print "content of your configuration file:\n-----------------------------------\n";
system ("cat $config_file");
print "-----------------------------------\n";



print "end of configuration\n";
print "now you can run the main program by using the following command:\n";
print "Kmask9.pl -c $config_file\n";





## später warnings abschalten



####################################
## sub routines
####################################

sub promptUser {
  my($prompt, $default) = @_;
  my $defaultValue = $default ? "[$default]" : "";
  print "$prompt $defaultValue: ";
  chomp(my $input = <STDIN>);
  return $input ? $input : $default;
}
