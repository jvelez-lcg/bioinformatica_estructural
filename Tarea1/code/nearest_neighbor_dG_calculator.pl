#!/usr/bin/perl -w
# V1.0
# Nearest Neighbor dG calculator.
# Salazar Diana, Velez Jesus.
# dsalazar@lcg.unam.mx, jvelez@lcg.unam.mx

# Returns:
# 	Print the following data in tabular format and if a field contains more than one value, these will be separated by a comma.
#	seq_name	dGs	positive_index	energy_difference	starts	ends	type_segment	tag
#	

# Examples of use:
# perl nearest_neighbor_dG_calculator.pl -i data/K12_400_50_sites -t 37 -w 15 -s 200 -o 25 -c 3.4,-15.99 -r -400,50 -d -200,50 -e test1
# perl nearest_neighbor_dG_calculator.pl --inputfile data/K12_400_50_sites --temperature 37 --window_length 15 --size_compare 200 --size_overlap 25 --cutoffs 3.4,-15.99 --region_seq -400,50 --delimited_region -200,50 --experiment_tag test1
# Combinations of short or long forms of parameters are also valid.
# Notes:
# 	For now there are no default values so all must be specified.
#	All the necessary validations have not been verified so you should be careful in the parameters that you enter make sense.

use strict;
use List::Util qw(sum);
use Getopt::Mixed;

Getopt::Mixed::init(
 "i=s t=i w=i s=i o=i c=s r=s d=s e=s
 inputfile>i temperature>t window_length>w size_compare>s size_overlap>o cutoffs>c region_seq>r delimited_region>d tag>e");
my ($infile, $T, $windowL, $size_of_seqs_to_compare, $size_overlap, @cutoffs, @region_of_seq, @delim_region_to_check, $experiment_tag);
while( my( $option, $value, $pretty ) = Getopt::Mixed::nextOption())
{
    OPTION: {
      $option eq 'i' and do {
        $infile = $value;
        last OPTION;
      };
      $option eq 't' and do {
        $T = $value;
        last OPTION;
      };
      $option eq 'w' and do {
        $windowL = $value;
        die "Window Size must be greater or equal to 1." if $windowL < 1;
        last OPTION;
      };
      $option eq 's' and do {
        $size_of_seqs_to_compare = $value;
        die "Size of sequence to compare must be greater or equal to 1." if $size_of_seqs_to_compare < 1;
        last OPTION;
      };
      $option eq 'o' and do {
        $size_overlap = $value;
        die "Size of overlap must be greater or equal to 1." if $size_overlap < 1;
        last OPTION;
      };
      $option eq 'c' and do {
        @cutoffs = split(/\,/,$value);
        die "Cutoffs array must be of length 2." if (@cutoffs != 2);
        last OPTION;
      };
      $option eq 'r' and do {
        @region_of_seq = split(/\,/,$value);
        die "Region array must be of length 2." if (@region_of_seq != 2);
        die "Region of sequence: First region must be lower than second one." if ($region_of_seq[0] > $region_of_seq[1]);
        last OPTION;
      };
      $option eq 'd' and do {
        @delim_region_to_check = split(/\,/,$value);
        die "Region array must be of length 2." if (@delim_region_to_check != 2);
        die "Delim region to check: First region must be lower than second one." if ($delim_region_to_check[0] > $delim_region_to_check[1]);
        last OPTION;
      };
	  $option eq 'e' and do {
		$experiment_tag = $value;
		last OPTION;
	  }
    }
};
Getopt::Mixed::cleanup();

print "# Parameters:\n";
print "# Local date and time: " . localtime() . "\n";
print "# Input file: $infile\n";
print "# Temperature = $T°C\tWindow Size = $windowL\n";
print "# Size of sequence to compare = $size_of_seqs_to_compare\tSize of Overlaping to check = $size_overlap\n";
print "# Cut-off of DeltaGs = $cutoffs[0]\tCut-off of E1 = $cutoffs[1]\n";
print "# Sequence region = $region_of_seq[0]:$region_of_seq[1]\n";
print "# Region to check for true positives = $delim_region_to_check[0]:$delim_region_to_check[1]\n";
print join("\t",("seq_name","dGs","positive_index", "energy_difference","starts","ends","type_segment","tag")) . "\n";

my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/k�mol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H',   0, 'S',-1.4 });

my ($dGs,$positives_idxs, $starts, $ends, $promoter_type);
open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>)
{   
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		my ($name,$seq) = ($1,$2); 
	
		my ($indexs_Ref, $deltaGs_Ref) = get_indexs_and_deltaGs_of_windows($seq, $windowL);
		my ($e1s_Ref, $e2s_Ref, $differences_Ref, $possible_regions_Ref) = calculate_e1_and_e2($indexs_Ref,$deltaGs_Ref);
		my @positives_idxs = extract_indexs_of_positive_signals($possible_regions_Ref, $differences_Ref, $e1s_Ref);
		my @overlaps = check_for_segments_that_overlap(@positives_idxs);
		my ($starts_Ref, $ends_Ref) = get_range_of_segments(\@overlaps,\@positives_idxs);
		my @trues_or_falses = check_true_false_segment($starts_Ref, $ends_Ref);

		my @starts_Ref = (na_if_empty_else_get_strings_versions($starts_Ref));
		my @ends_Ref = (na_if_empty_else_get_strings_versions($ends_Ref));

		my @results = na_if_empty_else_get_strings_versions($deltaGs_Ref, \@positives_idxs, $differences_Ref, \@starts_Ref, \@ends_Ref, \@trues_or_falses);
		print(join("\t",$name, @results, $experiment_tag) . "\n");
	}
}
close(SEQ);

# Calculate NN free energy of a DNA duplex , dG(t) = (1000*dH - t*dS) / 1000
# Parameters: 1) DNA sequence string; 2) Celsius temperature.
# Returns: 1) Free energy scalar.
# Notes: Uses global hash %NNparams.
sub duplex_deltaG
{
   	my ($seq,$tCelsius) = @_; 
	
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);
	my @sequence = split(//,uc($seq));
	my $tK = 273.15 + $tCelsius;
	
	sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] }
	
	# add dG for overlapping dinculeotides
	for(my $n=0;$n<$#sequence;$n++) 
	{
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.
				complement($sequence[$n].$sequence[$n+1]);
			
			if(!defined($NNparams{$DNAstep}))
			{
				$DNAstep = reverse($DNAstep);
			}
			
			$dG = ((1000*$NNparams{$DNAstep}{'H'}) -
					($tK*$NNparams{$DNAstep}{'S'})) / 1000 ;
			
			$total_dG += $dG; 
	}
	
	# add correction for helix initiation
	$nt = $sequence[0]; # first pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'})) / 1000; 
	
	$nt = $sequence[$#sequence]; # last pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) }
	$total_dG += ((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'})) / 1000;
					
	# please complete for symmetry correction
	
	return $total_dG;
}

# Calculate NN free energy of a DNA duplex for each window in it.
# Parameters: 1) DNA sequence string; 2) Celsius temperature; 3) Size of windows to evaluate dG.
# Returns: 1) Indexs of Windows; 2) dGs array of windows.
# Notes: Uses global hash %NNparams and scalar $T.
sub get_indexs_and_deltaGs_of_windows
{
	my ($seq, $windowL) = @_;
	my $start = ($windowL -1) / 2;
	my $end = (length $seq) - $start - 1;
	my @indexs = $start..$end;
	my @deltaGs = map {duplex_deltaG(substr($seq, $_ - $start, $windowL), $T)} @indexs; 
	return(\@indexs, \@deltaGs);
}

# Calculates e1 and e2 as well as its differeces.
# Parameters: 1) Reference array to indexs of windows; 2) Refence array of dGs.
# Returns: 1) Array e1; 2) Array e2; 3) Array of differences (e1 -e2); 4) Array of possible regions of promoter.
# Notes: Uses global scalar $size_of_seqs_to_compare and then it is assigned to scalar $size.
sub calculate_e1_and_e2
{
	my ($indexs_Ref, $deltaGs_Ref) = @_;
	my (@e1s, @e2s, @differences, @possible_regions, @seq);
	my ($e1, $e2);
	my $size = $size_of_seqs_to_compare;
	my $deltaGsLength = @$deltaGs_Ref;

	die "Invalid Size." unless ($size >= 1 and $size % 4 == 0);

	foreach my $index (@{$indexs_Ref}){
		@seq = $index..($index + $size);

		last if ($seq[-1] > $deltaGsLength);

		$e1 = sum(@{$deltaGs_Ref}[@seq[0..((1/4 * $size) - 1)]]) / (1/4 * $size);
		$e2 = sum(@{$deltaGs_Ref}[@seq[(1/2 * $size)..($size-1)]]) / (1/2 * $size);

		push @e1s, $e1;
		push @e2s, $e2;
		push @differences, $e1 - $e2;
		push @possible_regions, $seq[1/4 * $size];
	}
	return(\@e1s, \@e2s, \@differences, \@possible_regions);
}

# Extract indexs of  positive signals if dGs > cutoff[0] and e1 > cutoff[2].
# Parameters: 1) Referece array of possible regions of promoter; 2) Reference array of differences between e1 and e2; 3)Reference array of e1.
# Returns: 1) Indices that meet the condition extracted from the array of possible_regions.
# Notes: Uses global array @cutoffs.
sub extract_indexs_of_positive_signals
{
	my ($possible_regions_Ref, $differences_Ref, $e1s_Ref) = @_;
	my (@positives_idxs, $i); my @indexs = @$possible_regions_Ref;

	for $i (0..@{$differences_Ref} - 1){
		push @positives_idxs, $indexs[$i] if((@$differences_Ref[$i] > $cutoffs[0]) and (@$e1s_Ref[$i] > $cutoffs[1]));
	}
	return (@positives_idxs);
}

# Boolean function of segments that overlap. 1 if segments shared nucleotides in range of $size_of_segment else 0.
# Parameters: 1) Array of positive indexs.
# Returns: 1) Array of booleans (1,0). 
# Notes: Uses global scalar $size_overlap.
sub check_for_segments_that_overlap
{
	my (@positives_idxs) = @_;
	return map {($positives_idxs[$_] + $size_overlap >= $positives_idxs[$_ + 1] - $size_overlap) ? 1 : 0} 0..(@positives_idxs - 2);
}

# Get starts and ends of segments.
# Parameters: 1) Reference array of overlaps; 2) Reference array of positive_idxs.
# Returns: 1) Reference array of starts. 2) Reference array of ends.
sub get_range_of_segments
{
	my ($overlaps, $positives_idxs) = @_;

	my (@starts, @ends, $start_idx, $i);
	my $start_segment = 0;

	for $i (0..@{$positives_idxs} - 1){

		if (@{$positives_idxs} == 1){
			push(@starts, ${$positives_idxs}[$i]);
			push(@ends, ${$positives_idxs}[$i]);
			next;
		}

		next if (${$overlaps}[$i] and $start_segment);
		if(${$overlaps}[$i]){
			$start_idx = ${$positives_idxs}[$i];
			$start_segment = 1;
		}
		else{
			if($start_segment > 0){
				push(@starts, $start_idx);
				push(@ends, ${$positives_idxs}[$i]);
				$start_segment = 0;
			}
			else{
				$start_segment = 0;
			}
		}
	}
	return(\@starts, \@ends);
}

# Check if the segments found are within a bounded region.
# Parameters: 1) Reference array of starts; 2) Reference array of ends;
# Returns: 1) Array with TRUE or FALSE string if the segment meets the condition.
sub check_true_false_segment
{
	my ($starts_Ref, $ends_Ref) = @_;
	my @starts = @{$starts_Ref};
	my @ends = @{$ends_Ref};
	my ($l, $u) = @delim_region_to_check;
	my @region = $region_of_seq[0]..$region_of_seq[1];

	return map {(($region[$starts[$_]] >= $l) and ($region[$ends[$_]] <= $u)) ? "TRUE" : "FALSE"} 0..@{$starts_Ref}-1;
}

# Join the elements of a vector with commas if it is not empty.
# Parameters: 1) Refenreces of Arrays to check.
# Returns: 1) NA if empty else string version of array.
sub na_if_empty_else_get_strings_versions
{
	my (@arrays_to_check) = @_;
	return map{@{$_} ? join(",", @{$_}) : "NA"} @arrays_to_check;
}