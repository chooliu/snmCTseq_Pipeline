
#!/usr/bin/perl -w
use strict;

# A04a_classify_mCT_reads_bismark.pl, v0.2 =====================================================
# original perl script written by Dr. Chongyuan Luo (@luogenomics)
# minor modifications by Choo Liu (@chooliu):
# - readability/documentation - minor bug fix (num total mCH over-counted) 
# inputs: - .sam file from Bismark (compatible with single-end or paired-end alignments)
# outputs: - "_annotations" .tsv recording each alignment's:
#             number of cytosines, mC/C fraction, and call (DNA, RNA, ambiguous)
#             this annotation file is subsequently appended to the .bam to keep DNA reads
# typical usage: perl A04a_classify_mCT_reads_bismark.pl alignments.sam
# =====================================================================================

my ($samfile)=($ARGV[0]);
my @name; my $name; my @sample; my $um; my $m;
my @row; my $xrtag;
my $fraction; my $totalch; my @call;

# default filtering criteria [*]
my $filter_num_CH=3;
my $filter_frac_mCHmin=0.5; # DNA (mCH/CH <0.5)
my $filter_frac_mCHmax=0.9; # RNA (mCH/CH >0.9)

# for each line in .sam file,
# count unmethylated and methylated cytosines in CHG and CHH context based on XR-tag 

# notes:
# - skip header header rows (starting with @)
# - assumes the XR-tag is in column 15 (bismark output), $row[14]
# - counting based on CHG and CHH, indicated with "h", "x", "H", "X"
#   (does not include "unknown context" CN or CHN, marked by "u" and "U")
# - default settings: keep reads with >=3 CHNs and mCHN/CHN fraction <=0.5 [*]
#   at present, we only retain call="DNA" and discard all others

open sam_in, "$samfile" or die $!;
open sam_annotations, ">$samfile\_annotations" or die $!;

while (<sam_in>)
{
	chop $_;
	if (substr($_,0,1) eq '@') { print sam_annotations "\n"; }
	else
	{
		@row=split(/\t/,$_);
		$xrtag=$row[14];
		$um=($xrtag=~tr/hx//); 
		$m=(($xrtag=~tr/HX//)-1); 
		$totalch=($um+$m); 
		if ($totalch==0) { 
			$fraction=999;
			@call="amb";
		} else {
			$fraction=(($m)/$totalch);
			if (($totalch>=$filter_num_CH) and ($fraction <= $filter_frac_mCHmin)) { @call="DNA"; }
				elsif (($totalch>=$filter_num_CH) and ($fraction>=$filter_frac_mCHmax)) { @call="RNA"; }
			else { @call="amb"; }
		}
		print sam_annotations "${totalch}\t${fraction}\t@call\n";
	}
}

close sam_annotations;
close sam_in;
