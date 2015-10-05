#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  CalculateThreshholds.pl
#
#  DESCRIPTION: Performs a linear fit to the ascending/descending portions of
#				k-mer frequency histogram curve to predict frequencies at which
#				error/repetitive k-mers will occur. Provides lower and upper bounds.
#				Bounds can be provided to Kmasker to mask appropriate k-mers from
#				sequence data.
#
#       AUTHOR:  Bryan Downie
#  AFFILIATION:  Leibniz Institute for Age Research, Jena Germany
#===============================================================================
#
# Modifications PK:
# removed STDERR output -> all information is just sent to STDOUT



#use strict;
#use warnings;

# This is most likely unnecessary
use Statistics::LineFit;
#use List::Util qw(sum);

use Getopt::Long;
use File::Basename;

my $script = basename($0);
my $usage = <<EOT;
Usage:  $script -c size [ <options> ] <histogram file>

	-c : Expected coverage
    -s / --sensitivity <value>: Sensitivity of repeat detection
                                (default is 0.01)	
    -es / --errorsensitivity <value>: Sensitivity of error detection
                                (default is 0.1)	


    Note on sensitivity: 
		* To aggressively include all "good" kmers, 
		  set sensitivity to 0.
		* To exclude as many repetitive kmers as possible, 
		  set sensitivity to 1.

        Author: Bryan Downie (bdownie\@fli-leibniz.de).
 
		 	 			    
EOT
if(!scalar(@ARGV)){printf "$usage";exit(0);}

GetOptions ('sensitivity|s=f' => \$s, 'errorsensitivity|es=f'=>\$es, 'c=i'=>\$c);

# Coverage is necessary for smoothing of data. Smoothing window size determined later.
if (!defined($c)) { print STDERR "Pass expected coverage using -c\n"; exit (0); }
if (defined($s)) { $sensitivity = $s; }
else { $sensitivity = 0.01; }
if (defined($es)) { $errorsensitivity = $es; }
else { $errorsensitivity = 0.1; }

my @raw;
$kmer_count = 0;
for (my $i = 1; $i < 10000; $i++) {
	$x_values[$i] = $i;
	$raw[$i] = 0;
}

open FILE, $ARGV[0] or die "Couldn't open file $ARGV[0]";
while (<FILE>) { 
	@line = split /\s+/;
	$raw[$line[0]] = $line[1];
	$kmer_count += $line[1] * $line[0];
}


my $coverage = $c;
print "Using coverage: $coverage\n";

# Remove the last entry. Jellyfish puts everything > 10000 into one bin which causes a peak
# at the end of the curve.
pop @raw;
pop @x_values;

# Smoothing window is 1 for each 20x coverage above 10.
# e.g. coverage 30 -> smoothing window 1, coverage 70 -> window 3
$smooth_param = int($coverage/10);
if ($smooth_param < 3) { $smooth_param = 1; }
elsif ($smooth_param < 5) { $smooth_param = 3; }
elsif (($smooth_param %2 ) != 1) { $smooth_param += 1; }

# Offset is to exclude the lowest values after taking the derivative.
my $smooth_offset = int($smooth_param/2);

# Valley assigned elsewhere. This is probably obsolete.
my $valley = 0;
my $peak = 0;
my $best_peak_val = 0;
my $is_ascending = 0;
my $last_value = $raw[1];
for ($i = 2; $i < $#raw; $i++) {
	if ($is_ascending) { 
		if ($raw[$i] > $best_peak_val) { 
			$peak = $i;
			$best_peak_val = $raw[$i];
		}
	}
	elsif ($last_value < $raw[$i]) {
		$is_ascending = 1;
	}
	else { $last_value = $raw[$i]; }
}


# Do smoothing;
@smoothed_raw = smooth_array($smooth_param, @raw);
for ($i =  1+$smooth_offset; $i < $#smoothed_raw - 1; $i++) {
	if ($smoothed_raw[$i-1]) { 
		$first_deriv[$i] = ($smoothed_raw[$i] - $smoothed_raw[$i-1])/($x_values[$i] - $x_values[$i-1]);
	}
}

my $last_value = $raw[$peak];
for ($i = $peak; $i > 0; $i--) { 
	if ($smoothed_raw[$i] > $last_value) { 
		$valley = $i + 1;
		$i = 0;
	}
	$last_value = $smoothed_raw[$i];
}
print "Peak: $peak\nValley: $valley\n";


# Guess genome size based on the total number of k-mers divided by the peak frequency.
# If there are excessive error k-mers, The genome size will be over-estimated.  
$genome_size = int($kmer_count/$peak);
print "Predicted genome size: $genome_size\n";

# Smooth again (we need the third derivative.
@first_deriv_smoothed = smooth_array($smooth_param,@first_deriv);

for ($i = 2+$smooth_offset; $i < $#first_deriv_smoothed - 2; $i++) {
	if ($first_deriv_smoothed[$i-1]) { 
		$second_deriv[$i] = ($first_deriv_smoothed[$i] - $first_deriv_smoothed[$i-1])/($x_values[$i] - $x_values[$i-1]);
	}
}

# Final smoothing
@second_deriv_smoothed = smooth_array($smooth_param,@second_deriv);
# TO CHANGE: 0 should be $valley here
for ($i = 0; $i < $peak; $i++) {
	# Left regression start/end is between the valley and peak
	if (($second_deriv_smoothed[$i] * $second_deriv_smoothed[$i-1]) < 0) { 
		$left_regression_end = $x_values[$i+1];
	}
}

# Smooth again so that we know when the rate of curving changes (3rd derivative)
for ($i = 3+$smooth_offset; $i < $#second_deriv_smoothed - 3; $i++) {
	if ($second_deriv_smoothed[$i-1]) { 
		$third_deriv[$i] = ($second_deriv_smoothed[$i] - $second_deriv_smoothed[$i-1])/($x_values[$i] - $x_values[$i-1]);
	}
}

@third_deriv_smoothed = smooth_array($smooth_param,@third_deriv);

# These are backup values in case the next section fails.
# But doesn't fit the logic exactly. Will need to adjust this.
# Should include whether or not the curve is up or down (second derivative)
$left_regression_start = $valley;
$left_regression_end = $peak;
for ($i = $valley + 1; $i < $peak; $i++) { 
	if ($third_deriv_smoothed[$i] == 0) { $third_deriv_smoothed[$i] = $third_deriv_smoothed[$i-1]; }
	elsif (($third_deriv_smoothed[$i] * $third_deriv_smoothed[$i-1]) < 0) { 
		if ($third_deriv_smoothed[$i] < 0) {
			if (!defined($left_regression_start)) {
				# This is where the errors end
				$left_regression_start = $x_values[$i];
				if ($left_regression_start > $peak) {
					$left_regression_start = $valley;
				}
			}
			else { $left_regression_end = $i - 1; }
		}
	}
}
print "Left regression start: $left_regression_start\nLeft regression end: $left_regression_end\n";

# Determine window for right regression fit.
for ($i = $peak; $i < $#third_deriv_smoothed - 3; $i++) { 
	if ($third_deriv_smoothed[$i] == 0) { $third_deriv_smoothed[$i] = $third_deriv_smoothed[$i-1]; }
	elsif (!$right_regression_start) { 
		if (($third_deriv_smoothed[$i] * $third_deriv_smoothed[$i-1]) < 0) { 
			$right_regression_start = $i;
		}
	}
	else {
		if (($third_deriv_smoothed[$i] * $third_deriv_smoothed[$i-1]) < 0) { 
			$right_regression_end = $i - 1;
			if ($right_regression_end == ($right_regression_start + 1)) { $right_regression_end++; }
			last;
		}
	}
}

print "Right regression start: $right_regression_start\nRight regression end: $right_regression_end\n";


$total = 0;
$count = 0;
# Take the average of the slope between each point between left regression start/end
for ($i = $left_regression_start; $i <= $left_regression_end; $i++) { 
	$total += $raw[$i] - $raw[$i-1];
	$count++;
}
$leftslope = $total/$count;

# Determine the intercept of the left regression.
$total = 0;
$count = 0;
for ($i = $left_regression_start; $i <= $left_regression_end; $i++) { 
	$count++;
	$total += int($raw[$i] - ($i * $leftslope));
}
$left_intercept = int($total/$count);



# Take the average of the slope between each point between right regression start/end
$total = 0;
$count = 0;
for ($i = $right_regression_start; $i <= $right_regression_end; $i++) { 
	$count++;
	$total += $raw[$i] - $raw[$i-1];
	$tot = $raw[$i] - $raw[$i-1];
}
$rightslope = $total/$count;

# Determine the intercept of the right regression.
$total = 0;
$count = 0;
if ($right_regression_end < ($right_regression_start + $smooth_param)) { 
	$right_regression_end = $right_regression_start + $smooth_param;
}
for ($i = $right_regression_start; $i <= $right_regression_end; $i++) { 
	$count++;
	$total += int($raw[$i] - ($i * $rightslope));
}
$right_intercept = int($total/$count);

# Predict where errors start after adjusting for error sensitivity
for ($i = $left_regression_start; $i > 0; $i--) { 
	$error_sum = $raw[$i];
	$extrapolated_sum = int(($leftslope * $x_values[$i]) + $leftintercept);
	if ($extrapolated_sum < 0) { $extrapolated_sum = 0; }
	if ($error_sum == 0) { 
		$left_threshold = $x_values[$i];
		print "($left_threshold)($x_values[$i])($i)\n";
		last;
	}
	else { 
		$percent = $extrapolated_sum/$error_sum;
		if ($percent > $errorsensitivity) {
			$left_threshold = $x_values[$i];
			last;
		}
	}
}
	
$error_sum = 0;
$extrapolated_sum = 0;
$percent = 0;

$right_threshold = 1 + int((0- $right_intercept)/$rightslope);
#$x0_intercept = 1 + int((0- $right_intercept)/$rightslope);

# Predict where errors start after adjusting for sensitivity
#if ($sensitivity > 0.9) { $sensitivity = 0.9; }
#for ($i = $x0_intercept ; $i  >= $right_regression_start; $i--) {
#	if ($percent > $sensitivity) {
#		$right_threshold = $i;
#		last;
#	}
#	
#	$actual_val = $raw[$i];
#	$extrapolated_val = int(($rightslope * $x_values[$i]) + $right_intercept);
#	$percent = $extrapolated_val/$actual_val;
#}
unless ($right_threshold) { 
	print "Couldn't calculate right threshold. Try a different sensitivity value\n";
	exit; 
	$right_threshold = $right_regression_start; 
	}
unless ($left_threshold) { 
	$left_threshold = $left_regression_end; 
	print "Couldn't calculate left threshold. Try adjusting the coverage or check your histogram file!\n";
	exit;
	#exit;
}

if ($right_threshold < $right_regression_end) { 
	$right_threshold = $right_regression_end;
}
if ($valley) { 
	$left_threshold = $valley;
}
#$right_threshold *=2;

print "Left Threshold: $left_threshold\n";
print "Right Threshold: $right_threshold\n";


# averages values in a array with its neighbors (window size determined by $number)
sub smooth_array {
	my $number = shift;
	my @array = @_;
	my @return_array;
	
	unless (($number %2) == 1) { print STDERR "Array smoothing requires odd integer.\n"; exit (0); }
	$div = int($number/2);
	for ($i = $div; $i < $#array; $i++) {
		$count = 0;
		$value = 0;
		next if !defined($array[$i]);
		for ($j = $i - $div; $j <= $i + $div; $j++) {
			if ($j < 0) { $j = 0; }
			if ($j <= $#array) { 
				if ($array[$j]) { 
					$value += $array[$j];
					$count++;
				}
			}
		}
		if ($count) { 
			$return_array[$i] = int($value/$count);
		}
	}
	return @return_array;
}

# Add an array together (better as eval)
sub sum {
	my @array = @_;

	my $val = 0;

	foreach $num (@array) {
		if ($num =~ /\D/) { return "error"; }
		$val += $num;
	}

	return $val;
}
