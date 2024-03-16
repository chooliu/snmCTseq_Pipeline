
#!/usr/bin/perl -w
use strict;

# A02a_demultiplex_fastq.pl, v0.2 ==============================================
# perl script originally written by Dr. Chongyuan Luo (@luogenomics)
# edits for readability, +1 Hamming CB mismatch by Choo Liu (@chooliu) 
# inputs: - R1 and R2 .fastq.gz file in 'fastq_raw' for whole plate;
#           containing all nuclei, first 8bp of Read 1 assumed to be 'cell barcode'
#         - index .fasta (> 384-well plate location \n cell barcode sequence)
# outputs: - R1 and R2 .fastq.gz for each well in 384-well plate
#          - summary .txt showing # of reads per cell barcode
# typical usage: perl A02a_demultiplex_fastq.pl r1.fq.gz r2.fq.gz index.fa
# caution: some hardcoded parsing of permitted cb index.fa list
# ==============================================================================

#!/usr/bin/perl -w
use strict;



# user parameters
my $enable_fuzzy=1; # <-- default=1 if enable 1bp hamming mismatch in CB
my $cblen=8; # <-- default=8, cell barcode length in bases


# load cell barcode --> well indices .fa
# & initialize %cb_counts (for summary .txts)

print "initializing perl demuliplexing script...\n";

my %hash_cb; my %index; my %cb_counts; 
my $hash_cb; my $index_str; my @index=();
open index_file, "$ARGV[2]" or die $!;
while (<index_file>)
{
  chop $_; $index_str=substr($_,1,length($_)-1);
  $_=<index_file>; 
  chop $_; $hash_cb{substr($_,1,length($_)-1)}=$index_str;
  $index{$index_str}=substr($_,1,length($_)-1);
  push (@index,$index_str);
}
close index_file;

my $cb; my $well;
while ( ($cb,$well) = each %hash_cb ) {
    $cb_counts{$well}=0;
}
$cb_counts{'undetermined'}=0;



# hamming dist subroutine: via perlmonks.org/?node_id=500235
# if enable_fuzzy=1, allow 1bp cell barcode mismatch
# note: our cell barcode designs are presently 2bp hamming dist

sub hd{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
my $numcorrectcb=0;



# prep output files

my @r1_list; my @r2_list; my @r1_name; my @r2_name;  my $r1_name; my $r2_name; 
my %r1_out; my %r2_out;
my $r1_in="$ARGV[0]";
my $r2_in="$ARGV[1]";

$r1_name=$r1_in;
$r1_name =~ s/fastq_raw/fastq_demultip/ig;
$r2_name=$r2_in;
$r2_name =~ s/fastq_raw/fastq_demultip/ig;

# 0: platename, 1: .fastq.gz suffix
 
@r1_name=split(/_R1/,$r1_name);
@r2_name=split(/_R2/,$r2_name);

my $cb; my $well;
while ( ($cb,$well) = each %hash_cb ) {
    local *FILE;
    open FILE, " | gzip -c > $r1_name[0]\_$well\_indexed_R1$r1_name[1]" or die $!;
    $r1_out{$well}=*FILE;
    
    local *FILE;
    open FILE, " | gzip -c > $r2_name[0]\_$well\_indexed_R2$r2_name[1]" or die $!;
    $r2_out{$well}=*FILE;
}



# loop through read-pairs

print "parsing cell barcodes...\n";

my $totreads; my $eight;
my @r1_tmp; my @r2_tmp; 
my $r1; my $r2;
my $foundfuzzy;

open r1_fastq_in, "gzip -cd $r1_in | " or die $!;
open r2_fastq_in, "gzip -cd $r2_in | " or die $!;
$totreads=0; @r1_tmp=(); @r2_tmp=();
while (<r1_fastq_in>)
{
    chop $_; $totreads++;
    $r1=$_; $r2=<r2_fastq_in>; chop $r2;
    push (@r1_tmp, $r1); push (@r2_tmp,$r2);
    if (!($totreads%4)) {
        $eight=substr($r1_tmp[1],0,$cblen);
        
        # if perfect cell barcode match found
        if (exists $hash_cb{$eight}) {
            $cb_counts{$hash_cb{$eight}}++;
                print {$r1_out{$hash_cb{$eight}}} "$r1_tmp[0]\n$r1_tmp[1]\n$r1_tmp[2]\n$r1_tmp[3]\n";
                print {$r2_out{$hash_cb{$eight}}} "$r2_tmp[0]\n$r2_tmp[1]\n$r2_tmp[2]\n$r2_tmp[3]\n";
        }
        
        # --> attempt to 1bp-mismatch correct ($myhd==1)
        # if observed eight-bp CB in read doesn't match a key in hash_cb
        else {
            if ($enable_fuzzy==1) {
                my $foundfuzzy=0;
                while ( (($cb,$well) = each %hash_cb)) {
                    my $myhd=hd($cb,$eight);
                    if ($myhd==1) {
                        $numcorrectcb++;
                        $foundfuzzy=1;
                        $cb_counts{$hash_cb{$cb}}++;
                        print {$r1_out{$hash_cb{$cb}}} "$r1_tmp[0]\n$r1_tmp[1]\n$r1_tmp[2]\n$r1_tmp[3]\n";
                        print {$r2_out{$hash_cb{$cb}}} "$r2_tmp[0]\n$r2_tmp[1]\n$r2_tmp[2]\n$r2_tmp[3]\n";
                        last;
                    }
                }
                # if no 1bp mismatch found after looping thru valid CBs
                if ($foundfuzzy == 0) {
                    $cb_counts{'undetermined'}++;
                }
            } # <-- attempt to correct
            
            # or, if enable_fuzzy != 1;
            # i.e., no attempt as mismatch correct
            else 
            {
              $cb_counts{'undetermined'}++;
            }
        }
        
    # go to next read in input fq
    @r1_tmp=(); @r2_tmp=();
}
}

close r1_fastq_in;
close r2_fastq_in;



# close open files &
# remove empty demultip.fq.gz (wells without any CB matches)
# results in these files being "missing" in logs but avoids some downstream issues

print "wrapping up & checking for empty .fqs...\n";

my $percreads; my $perccorr; my $command; my $return;
while ( ($cb,$well) = each %hash_cb ) {

    close $r1_out{$well};
    close $r2_out{$well};
    
    if ($cb_counts{$well} == 0) {
        $command="rm $r1_name[0]\_$well\_indexed_R1$r1_name[1]";
        print "$command\n"; system($command);
        $command="rm $r2_name[0]\_$well\_indexed_R2$r2_name[1]";
        print "$command\n"; system($command);
    }
    
}

# prep summary .txt files where $cbset = 1 or 2
# (careful with updated cb valid .fa sets, as below is hardcoded)

my $cbset = substr($ARGV[2], -4, 1);
open summary_out, ">$r1_name[0]\_summary\_$cbset.txt" or die $!;
$totreads=int($totreads/4);
$perccorr=sprintf('%.2f', 100*($numcorrectcb/$totreads)); $perccorr="$perccorr%";
print summary_out "# total read pairs\t$totreads\n";
print summary_out "# pairs 1bp mm corrected\t$numcorrectcb\t$perccorr\n";

for (sort {$a cmp $b} keys %cb_counts) 
{
    $percreads=sprintf('%.2f', 100*($cb_counts{$_}/$totreads)); $percreads="$percreads%";
    print summary_out "$_\t$index{$_}\t$cb_counts{$_}\t$percreads\n";
}
close summary_out;



print "completed perl script.\n";
