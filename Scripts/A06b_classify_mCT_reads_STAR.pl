
#!/usr/bin/perl -w
use strict;


# A06b_classify_mCT_reads_STAR.pl, v0.2 =====================================================
# based loosely on original perl script written by Dr. Chongyuan Luo (@luogenomics)
# modifications by Choo Liu (@chooliu):
# - readability/documentation
# - changes in logic for paired-end mapping (fwd/rev strand) & more efficient parsing of MD:Z flag
#   incl looking at only C/G positions in read vs. ref changes, calculating # cytosines, etc
# inputs: - .sam file from STAR (compatible with single-end or paired-end alignments)
# outputs: - "_annotations" .tsv recording each alignment's:
#             number of cytosines, mC/C fraction, and call (DNA, RNA, ambiguous)
#             this annotation file is subsequently appended to the. bam to keep RNA reads
# typical usage: perl A06b_classify_mCT_reads_STAR.pl alignments.sam
# =====================================================================================

# reading in .sam file
my ($samfile)=($ARGV[0]);
my @sample; my @samline; 
my $read; my @read; my @ref; my $ref; my @md; my $dir; my $unmch; 

# determining read classification
my $filter_num_CH=3;
my $filter_frac_mCHmin=0.5; # DNA (mCH/CH <0.5)
my $filter_frac_mCHmax=0.9; # RNA (mCH/CH >0.9)
my $mch_fraction; my $totalch; my $totcpg; my @call;




# load .sam file, loop through lines in .sam 
# exports info on each read to "(input .sam file name)_annotations"

# for each line in .sam file, --------------------------------------------------------
# count unmethylated and methylated cytosines in CHG and CHH context based on XR-tag 

# notes:
# - skip header header rows (starting with @)
# - assumes the XR-tag is in column 10 (STAR output), $samline[9]
# - currently ignores indels
# - default settings: keep reads with >=3 CHNs and mCHN/CHN fraction >=0.9 [*]
#   at present, we only retain call="RNA" and discard all others

open sam_in, "$samfile" or die $!;
open sam_annotations, ">$samfile\_annotations" or die $!;
while (<sam_in>)
{
  chop $_; 
  if (substr($_,0,1) eq '@') { print sam_annotations "\n"; }
  else
  {
  
    @samline = split(/\t/,$_);
    $read = $samline[9];
    @read = split(//,$read);
    @ref = @read;
    
    @md=split( /(\d+)/ , $samline[15]);
    
    # check mapping orientation relative to reference genome ----------------------------
    # dir=1 if fwd (99, 163, 0, 73), =0 if reverse (147, 83, 16, 89)
    # flaw in script is that these numbers manually specified--if see strange results,
    # check mapping output for other flags before filtering / add interpreter
    $dir = 0 + ( ($samline[1] eq 99) or ($samline[1] eq 163 ) or ($samline[1] eq 0) or ($samline[1] eq 73) ); 
    
    # compare read to reference genome --------------------------------------------------
    # split MD:Z: flag by numbers then examine resulting length
    # e.g., MD:Z:10A0B121 --> MD:Z: 10 A 0 B 121 has $num_md_features=5
    
    # assume fully methylated if perfect match to genomic ref
    # (if MD stores a single number indicating no changes, $#md = 1)
    my $unmch = 0;
    my $num_md_features=$#md;
    if ($num_md_features == 1) { }
    
    # <-- start ELSE for num_md_features != 1
    # otherwise, attempt to reconstruct reference sequence, 
    # looping through positions in read where the read != ref nucleotide
    else {
        my $pos=0;

        for (my $i=1; $i<=$num_md_features; $i += 2) {
        
            # loop through bases where there are changes
            $pos = $pos + @md[$i];
            
            # do nothing if starts with "^" (deletion)
            # modify @ref sequence if differences
             if ( (rindex("@md[$i+1]", "^", 0) == 0) or ($i+1 > $num_md_features)) { }  
             else {
                 @ref[$pos] = @md[$i + 1];
                 }
            
            # increment base by 1
            $pos += 1;
        }
  
    # tabulate # unmethylated cytosines --------------------------------------------------- 
    # again looping through each position where read != ref
    # (although same loop as above, re-loop because potential cases where adjacent bases differ btwn ref & read)
        my $pos=0;
        for (my $i=1; $i<=$num_md_features; $i += 2) {
        
           $pos = $pos + @md[$i];
           
           # if read maps to forward strand of genome
           # check if at a CH-site, and unmethylated cytosine converted to "T"
           if ( ($dir eq 1) and
                (@ref[$pos] eq "C") and (@ref[$pos + 1] ne "G") and (@read[$pos] eq "T") ) {
                    $unmch += 1;
                }
       
           # if read maps to reverse strand of genome
           # STAR rev compliments the $read sequence, so cytosines are represented by "G"
           # check if at a CH-site, and unmethylated cytosine converted to "A"
           # (note: mCH underestimated in niche case where dir=0 & 
           #        first base is cytosine, as @ref[$pos - 1] is undefined)
           if ( ($dir eq 0) and
                (@ref[$pos] eq "G") and (@ref[$pos - 1] ne "C") and (@read[$pos] eq "A") ) {
                $unmch += 1;
            }
            
            $pos += 1;
        }
    } # <-- end ELSE statement for num_md_features != 1

    $ref = join('', @ref);
    
    # count cytosines in CH-context --------------------------------------------------
    # note: since done via regex, faster to count [# C] and subtract [# CG] vs searching wild flag C[ACT]
    if ($dir eq 1) {
      my @totalch = $ref =~ /C/g;
      $totalch = scalar @totalch;
      my @totcpg = $ref =~ /CG/g;
      $totcpg = scalar @totcpg;
    }
    if ($dir eq 0) {
      my @totalch = $ref =~ /G/g;
      $totalch = scalar @totalch;
      my @totcpg = $ref =~ /CG/g;
      $totcpg = scalar @totcpg;
    }

    my $totalch = $totalch - $totcpg;
    
    # classify each read into modalities ----------------------------------------------
    if ($totalch==0) { # avoid div by zero error
          $mch_fraction=-999;
          @call="amb";
    } else {
          $mch_fraction = 1 - ($unmch/$totalch);
          if (($totalch>=$filter_num_CH) and ($mch_fraction < $filter_frac_mCHmin)) { @call="DNA"; }
          elsif (($totalch>=$filter_num_CH) and ($mch_fraction >= $filter_frac_mCHmax)) { @call="RNA"; }
          else { @call="amb"; } # exclude (low # CH or ambiguous mCH between 0.5-0.9)
    }

    print sam_annotations "${totalch}\t${mch_fraction}\t@call\n";

}
}

close bam_in;
