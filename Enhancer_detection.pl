#!/usr/bin/perl
use lib "/nfs/central/home/ma90/lib/perl5/site_perl/5.8.8/";
our $moduleid =  $ARGV[0];
our $windowsize = 500;
chomp($moduleid); chomp($windowsize);

use constant KVALUE => 6;  #k stands for kmer
use constant fraction => 0.8;
use constant W => 0.1;
use constant STEP => 50;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::BPbl2seq;


my @alphabet = ("A","T","G","C");
# generate random set of oligonucleotides
our @kmerarray = ("A","G","T","C");

for(my $times=0; $times<KVALUE-1;$times++) {
        my $size = scalar @kmerarray;
        for(my $i = 0; $i < 4 ; $i++) {
                for(my $j=0;$j<$size;$j++) {
                        push(@temparray, $kmerarray[$j].$alphabet[$i]);
                }
        }
        @kmerarray = @temparray;
        undef @temparray;
}

our @backarray1;
open(FREAD1,"Dmel_background");
my @lines1 = <FREAD1>;
foreach $ele1 (@lines1)
{
        chomp($ele1);
        my @array1 = split(/\t/,$ele1);
        $backarray1{"$array1[0]"} = $array1[1] * ($windowsize - KVALUE +1);
}      

our @backarray2;
open(FREAD2,"$ARGV[1]_background");
my @lines2 = <FREAD2>;
foreach $ele2 (@lines2)
{
        chomp($ele2);
        my @array2 = split(/\t/,$ele2);
        $backarray2{"$array2[0]"} = $array2[1] * ($windowsize - KVALUE +1);
}

sub fac
{
	return $_[0]>1?$_[0]*fac($_[0]-1):1;	
}

sub poi_dm
{
	local($x,$mer) = @_;
	my $exp_no1 = $backarray1{"$mer"};
	my $func_value1 = 0;
    for(my $ii=0;$ii<=$x;$ii++)
    {
        $func_value1 += ( (exp(-$exp_no1)) * ($exp_no1**$ii) / (&fac($ii)) );
        
	}
	return($func_value1);
}

sub poi_dp
{
	local($x,$mer) = @_;
	my $exp_no2 = $backarray2{"$mer"};
	my $func_value2 = 0;
	for(my $ii=0;$ii<=$x;$ii++)
	{
		$func_value2 += ( (exp(-$exp_no2)) * ($exp_no2**$ii) / (&fac($ii)) );
		
		}
		return($func_value2);
}

sub getmin
{
	local(%array) = @_;
	my $min = 100000;
	foreach $diff (values %array)
	{
		if($diff < $min)
		{$min = $diff;}
		else{next;}
	}
	return($min);

}
sub fillhash
{
    local($dm,$dp) = @_;
    my %mhash;my %dhash;
    foreach $mer(@kmerarray)
    {
        my $count_dm = () = $dm =~ m/(?=($mer))/g;
        my $count_dp = () = $dp =~ m/(?=($mer))/g;
        $mhash{$mer} = $count_dm; $dhash{$mer} = $count_dp;
    }
    return(\%mhash,\%dhash);
    
}

sub caldist
{
	local($ref_hashdm,$ref_dp,$ref_hashdp) = @_;
    my %hashdm = %{$ref_hashdm}; my %hashdp = %{$ref_hashdp}; my $dp = ${$ref_dp};
	my $sum_sim = 0,$sum_dissim = 0,$min,$prob_sim;
	foreach $mer(@kmerarray)
	{
		my $count_dm = $hashdm{$mer};
		my $count_dp = $hashdp{$mer};
        if($count_dm > $count_dp) {$min = $count_dp;}
		else{$min = $count_dm;}
		if($min == 0){$prob_sim = 1;}
		else
		{
			$prob_sim = ( (1 - &poi_dm($min - 1,$mer)) * (1- &poi_dp($min - 1,$mer)) );
		}
		
        my $p = (1 - $prob_sim);
		$sum_sim += $p;
        my $q = abs(&poi_dp($count_dp,$mer) - &poi_dm($count_dm,$mer));
        $sum_dissim += $q;
	} 
	my $sim =  $sum_sim/(4**KVALUE);
	my $dissim = $sum_dissim/(4**KVALUE);
    return ($sim - (W * $dissim));	
}

sub updatehash
{
    local($ref_hash,$subseq,$addseq) = @_;
    my %hash = %{$ref_hash};
    my $len = length($subseq);
    for(my $i=0;$i<$len - KVALUE+1;$i++)
    {
        my $word = substr($subseq,$i,KVALUE);
        $hash{$word}--;
    }
    my $len = length($addseq);
    for(my $i=0;$i<$len - KVALUE+1;$i++)
    {
        my $word = substr($addseq,$i,KVALUE);
        $hash{$word}++;
    }
    return(\%hash);
}

sub slidewindow
{

local($ref_dpseu_reg,$windowsize,$ref_dmel_reg) = @_;
my $dpseu_reg = ${$ref_dpseu_reg};  my $dmel_reg = ${$ref_dmel_reg};
my $lastindex_dpseu = length($dpseu_reg); my $lastindex_dmel = length($dmel_reg);
our %assocarray;
my $maxsim = -100;
my %dmelhash,%dpseuhash;
$index = 0; # index for checking if windowsize is greater than dmel_reg
if($lastindex_dmel <= $windowsize){$index = 1;}
my $dmel_enh;
for(my $k=0;($k< ($lastindex_dmel - $windowsize) || $index ==1) ;$k = $k+ STEP)
{
    if($index == 0){
    $dmel_enh = substr($dmel_reg,$k,$windowsize);
    }
    else{
	$dmel_enh = $dmel_reg;}
    my $flag = 0, $output = "", $percent = 10;
    my %assocarray,%dmelhash,%dpseuhash;
for(my $i = 0 ; $i < $lastindex_dpseu - $windowsize; $i = $i + STEP)
{
    my $winseq = substr($dpseu_reg,$i,$windowsize);
    if($i==0){($ref_dmelhash,$ref_dpseuhash) = fillhash($dmel_enh,$winseq);
                %dmelhash = %{$ref_dmelhash}; %dpseuhash = %{$ref_dpseuhash}; }
    my $sim = &caldist(\%dmelhash,\$winseq,\%dpseuhash);
    print OOO "$sim\t";
    my $f = $i,$l = $i+$windowsize;
    $assocarray{"$f-$l"} = $sim;
    if($sim > $maxsim) {$maxsim = $sim; $max_f = $f ; $max_l = $l; $flag = 1; }
    else {}
    
    $output .= "$i\t$sim\n";
    print OUTPUT $i."\t".$sim."\n";
    if($i > $percent * $lastindex_dpseu/ 100)
    {
        $percent += 10;
    }
    my $frontseq = substr($dpseu_reg,$i,STEP + KVALUE -1);
    my $backseq = substr($dpseu_reg,$i+$windowsize - KVALUE +1,STEP+ KVALUE-1);
    ($ref_dpseuhash) = updatehash(\%dpseuhash,$frontseq,$backseq);
    %dpseuhash = %{$ref_dpseuhash};
}

    my $frontseq_dmel = substr($dmel_reg,$k,STEP + KVALUE -1);
    my $backseq_dmel = substr($dmel_reg,$k+$windowsize - KVALUE +1,STEP+ KVALUE-1);
    ($ref_dmelhash) = updatehash(\%dmelhash,$frontseq_dmel,$backseq_dmel);
    %dmelhash = %{$ref_dmelhash};  
    if($flag == 1){ %maxassocarray = %assocarray; $maxoutput = $output; $flag =0;  }
    if($index == 1){last;}
}
    print MAXOUT $maxoutput;
    return(\%maxassocarray,$maxsim,$max_f,$max_l);
}

open(OOO,">simout");
# MAIN PROGRAM START
my $fname = $moduleid."_".$ARGV[1].".csv";
open(OUTPUT,">plots/$fname");
my $fmaxname = $moduleid."_".$ARGV[1]."_max.csv";
open(MAXOUT,">plots/$fmaxname");

my $file1 = $moduleid."_Dmel.fa";
my $file2 = $moduleid."_$ARGV[1].fa";

$seqio_obj_dm = Bio::SeqIO->new(-file => "$file1", -format => "fasta" );
 $seq_obj_dm = $seqio_obj_dm->next_seq;
 my $dmel_enh = $seq_obj_dm->seq;
 
$seqio_obj_dp = Bio::SeqIO->new(-file => "$file2", -format => "fasta" );
 $seq_obj_dp = $seqio_obj_dp->next_seq;
 my $dpseu_reg = $seq_obj_dp->seq;

$timestart = time();
my ($ref_assocarray1,$maxsim1,$max_f1,$max_l1) = &slidewindow(\$dpseu_reg,$windowsize,\$dmel_enh);
$timeend = time();
$comptime = $timeend - $timestart;

my $dpseu_reg_rev = $dpseu_reg;
$dpseu_reg_rev =~ tr/ATGC/TACG/;
$dpseu_reg_rev = reverse($dpseu_reg_rev);
$timestart = time();
my ($ref_assocarray2,$maxsim2,$max_f2,$max_l2) = &slidewindow(\$dpseu_reg_rev,$windowsize,\$dmel_enh);
$timeend = time();
$comptime = $timeend - $timestart;
my %assocarray,$maxsim,$max_f,$max_l;
if($maxsim1 > $maxsim2)
{
    %assocarray = %{$ref_assocarray1};
    $maxsim = $maxsim1;
    $max_f = $max_f1;
    $max_l = $max_l1;
    $dpseq_out = $dpseu_reg;
}
else
{
    %assocarray = %{$ref_assocarray2};
    $maxsim = $maxsim2;
    $max_f = $max_f2;
    $max_l = $max_l2;
    $dpseq_out = $dpseu_reg_rev;
}


my $cutoffsim = fraction * $maxsim;
$store_max_f = $max_f; $store_max_l = $max_l;
our $left = $max_f;
our $right = $max_l;


while()
{
	$test_f = $max_f - STEP; $test_l = $max_l - STEP;
	if($test_f < 0){last;}
    $key = "$test_f-$test_l";
	if(defined($assocarray{"$key"})){} else{last;}
	if( $assocarray{"$key"} > $cutoffsim)
	{
		$left = $test_f;
		$max_f = $test_f;
		$max_l = $test_l;
		next;
	}
	else {last;}
}
$max_f = $store_max_f; $max_l = $store_max_l; 
while()
{
        $test_f = $max_f + STEP; $test_l = $max_l + STEP;
        if($test_l > $lastindex) {last;}
        $key = "$test_f-$test_l";	
        if(defined($assocarray{"$key"})){} else{last;}
        if( $assocarray{"$key"} > $cutoffsim)
        {       
                $right = $test_l;
                $max_f = $test_f;
                $max_l = $test_l;
                next;
        }
        else {last;}
}
open(OT,">$moduleid-$ARGV[1].res");
print OT "score: $maxsim\n";
print OT "max: $store_max_f , $store_max_l\n";
print OT "ends: $left , $right\n";
$sq = substr($dpseq_out,$left,$right - $left)."\n";
`rm plots/$fname`;
