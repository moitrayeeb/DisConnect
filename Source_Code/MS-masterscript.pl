###############################################################################
 # Copyright (c) 2012, Moitrayee Bhattacharyya, Kallol Gupta, P. Balaram, Indian Institute of Science, Bangalore, India.
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification are permitted provided that the following conditions are met:
 # 
 # (1) Redistributions of source code must retain the above copyright notice,
 # this list of conditions and the following disclaimer.
 # (2) Redistributions in binary form must reproduce the above copyright notice,
 # this list of conditions and the following disclaimer in the documentation
 # and/or other materials provided with the distribution.
 # (3) Neither the name of the Indian Institute of Science, Bangalore, India 
 # nor the names of its contributors may be used to endorse or promote products 
 # derived from this software without specific prior written permission.
 #  
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 # BUSINESSINTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 # IN CONTRACT, STRICT LIABILITY,OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 # POSSIBILITY OF SUCH DAMAGE.
 ##############################################################################

### This script predicts proteolytic peptide fragments ###
#use strict;
#use warnings;

##### Initializations #####

my %hash=(); my %term=(); my @uniq; my $key; my $i1; my $j1; my $j; my $a_me; my $compare1; my $choice; my @frag1=(); my @seq=(); my $line; my $hit; my $comp; my @user1=(); my $ele1; my @mul_term=(); my $len_mul; my $t; my @var_mult_prot=(); my  $len_mult_prot; my @user_prot=(); my @user_term_stat=(); my $mb; my $len_term_stat; my $mult_prot;
my @sort_hold1 = ();

print "Enter name of protein: ";
$prot=<STDIN>;

open(Read,"prot_seq") or die $!;
my @hold_seq=<Read>;
foreach $line (@hold_seq)
{
        $line=~s/\n//g;
        push (@seq1,$line);
}
close(Read);
foreach $line (@seq1)
{
        $line=~s/\n//g;
        push (@seq_tmp,$line);
        $ele=join('',@seq_tmp);
        $seq[0]=$ele;
        $len=@seq;
}
if ($seq[0]=~m/X/g)
{
        @ele=split('X',$seq[0]);
        $seq[0]="";
        foreach $line1(@ele)
        {
                $line1=~s/\n//g;
                push(@seq,$line1);
        }
#push(@seq,"\*");
shift(@seq);
}
#print @seq,"\n";

if ((join('',@seq))!~m/\*/)
{
	print "Error in input: Protein sequence input format is Sequence followed by *\n";
}

##### Rescaling of user mass; taking care of charge states ######

readpipe("perl MS-rescale_peak.pl");

##### User input for usage of monoisotopic / average mass ######

print  "Enter 1/2 for monoisotopic/average mass\n";
$mass_choice=<STDIN> ;

if ($mass_choice == 1)
{
        readpipe("cp mass_list_monoisotopic-MS mass_list");
}
else
{
        readpipe("cp mass_list_average-MS mass_list");
}
print "Enter the C-term status of your protein (Amidated/Free) as 1/2\n";
my $cterm_s=<STDIN>;
&prot_digest();

######## Additional rounds of digestion ###############

print "Enter number of additional proteases to be added from the list: ";
$num_prot=<STDIN>;
for ($tmp=0;$tmp<$num_prot;$tmp++)
{
#print "Enter your choice for $tmp_prot th protease: ";	
#$add_prot=<STDIN>;
#if ($add_prot=~m/y/i)
#{
	&add_proteases();
#}
}

sub add_proteases
{
	readpipe("cp proteo_frag.out prot_seq");
	open(Read1,"prot_seq") or die $!;
	@seq=<Read1>;
#	print @seq;
	close(Read1);
        &prot_digest();
	return();
}
################################################################

sub prot_digest
{
open(W,">proteo_frag.out");
##### User Input for proteolytic fragment generation #####
print "Enter your choice of protease\n";
print "Enter 1 for Trypsin\n";
print "Enter 2 for Trypsin-P\n";
print "Enter 3 for V8(DE)\n";
print "Enter 4 for V8(D)\n";
print "Enter 5 for Asp-N\n";
print "Enter 6 for Lys-C\n";
print "Enter 7 for Arg-C\n";
print "Enter 8 for CNBr\n";
#print "Enter 9 for any combination of above proteases (multiple proteases)\n";
print "Enter 10 for user defined protease\n";
print "####################################################\n";
$choice=<STDIN>;
$choice=~s/\n//g;

##### Initializations of hash #####

push (@{$hash{1}},"K","R");
push (@{$hash{2}},"K","R");
push (@{$hash{3}},"D","E");
push (@{$hash{4}},"D");
push (@{$hash{5}},"D");
push (@{$hash{6}},"K");
push (@{$hash{7}},"R");
push (@{$hash{8}},"M");

%term=(1=>'C', 2=>'C', 3=>'C', 4=>'C', 5=>'N', 6=>'C', 7=>'C', 8=>'C');

############ Special conditions 9 and 10 and basic calculations ##############
if ($choice == 9)
{
	my $item;
	print "Enter the desired combination of proteases from the list above in the format 1+2+3+.... \n";
	$mult_prot=<STDIN>;
	$mult_prot=~s/\n//g;
	@var_mult_prot=split('\+',$mult_prot);
	$len_mult_prot=@var_mult_prot;
	for($i1=0;$i1<$len_mult_prot;$i1++)
	{
		my $a1=$var_mult_prot[$i1];		
        	push (@user1,@{$hash{$a1}});
	}
	my %seen = (); foreach $item (@user1) { push(@uniq, $item) unless $seen{$item}++; }
	$compare1=join('',@uniq);

		
}
if ($choice == 10)
{
	print "User entry of the specifics of the protease\n####################################################\n";
	print "Enter the cleavage site/s in the format X+Y+...., where X and Y are the amino acid residues: ";
	$compare1=<STDIN>;
	$compare1=~s/\+//g;
	$compare1=~s/\n//g;
	@user_prot=split('',$compare1);
	print "Enter C/N for cleavage for each site: ";
	my $termini_stat=<STDIN>;
	$termini_stat=~s/\n//g;
	@user_term_stat=split('',$termini_stat);
	$len_term_stat=@user_term_stat;
}
if ($choice<=8)
{
	$compare1=join('',@{$hash{$choice}});
}
############################################################################################################

my $count=0;
foreach $line_tmp (@seq)
{
$line_tmp=~s/\n//g;
if ($line_tmp=~m/\*/)
{
	$flag=1;
}
else
{
	$flag=0;
}
$line_tmp=~s/\*//g;
$line_tmp1=uc($line_tmp);
my @ele_seq=split('',$line_tmp1);
my $len_seq=@ele_seq;
#print "$line_tmp1	$len_seq\n";
####### Option 1 : Trypsin (no X-Pro bond cleavage #######
if ($choice == 1)
{
for($j=1;$j<$len_seq;$j++)
{

        if (($compare1=~m/$ele_seq[$j]/) || ($ele_seq[$j]=~m/[a-z]/) || ($ele_seq[$j]=~m/X/))
## && ($j ne $len_seq-2))
        {
                if ($ele_seq[$j+1] ne "P")
                {
                        if ($ele_seq[$j]=~m/X/i)
                        {
                        $ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        @frag1=@ele_seq[$count..$j-1];
                        }
                        else
                        {
                        $ele_seq[$j]=~tr/[A-Z]/[a-z]/;
                        @frag1=@ele_seq[$count..$j];
                        }
                        print W @frag1,"\n";
                        $count=$j+1;
                }
        }
}
$ele_seq[$len_seq-1]=~tr/[A-Z]/[a-z]/;
@frag1=@ele_seq[$count..$len_seq-1];
if ($flag == 1)
{
	print W @frag1,"*\n";
}
else
{
	print W @frag1,"\n";
}
}
####### Option 2-8 cleavage #######
if ($choice > 1 && $choice < 9)
{
for($j=1;$j<$len_seq;$j++)
{
        if (($compare1=~m/$ele_seq[$j]/) || ($ele_seq[$j]=~m/[a-z]/) || ($ele_seq[$j]=~m/X/))
        {
                if ($term{$choice} eq "C") 
                {
                        if ($ele_seq[$j]=~m/X/i)
                        {
                        $ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        @frag1=@ele_seq[$count..$j-1];
                        }
                        else
                        {
                        $ele_seq[$j]=~tr/[A-Z]/[a-z]/;
                        @frag1=@ele_seq[$count..$j];
                        }
                        print W @frag1,"\n";
                        $count=$j+1;
                }
                else
                {
#			if (($j ne 1) && ($j ne $len_seq-1)) Remark: Commented out to remove endoprotease selectivity
#			{
                        	if ($ele_seq[$j]=~m/X/i)
                        	{
                        		$ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        		@frag1=@ele_seq[$count..$j-1];
                        		$ele_seq[$j]=~s/[a-z]//i;
                        	}
                        	else
                        	{
                        		$ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        		@frag1=@ele_seq[$count..$j-1];
                        	}
                        	print W @frag1,"\n";
                        	$count=$j;
#			}
                }
        }

}
$ele_seq[$len_seq-1]=~tr/[A-Z]/[a-z]/;
@frag1=@ele_seq[$count..$len_seq-1];
if ($flag == 1)
{
	print W @frag1,"*\n";
}
else
{
	print W @frag1,"\n";	
}
}

####### Option 10 cleavage - when user inputs the specifics of desired protease/s not present in the given list #######
if ($choice == 10)
{
for($j=1;$j<$len_seq;$j++)
{
        if (($compare1=~m/$ele_seq[$j]/)|| ($ele_seq[$j]=~m/[a-z]/) || ($ele_seq[$j]=~m/X/))
        {
                $hit=$&;
                for($mb=0;$mb<$len_term_stat;$mb++)
                {
                        my $a2=$user_prot[$mb];
                #        print $a2,"\n";
                        if ($a2=~m/$hit/)
                        {
                                $term{$choice}=$user_term_stat[$mb];
                        }
                }
#                print "$term{$choice}";
                if ($term{$choice} eq "C") 
                {
                        if ($ele_seq[$j]=~m/X/i)
                        {
                        $ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        @frag1=@ele_seq[$count..$j-1];
                        }
                        else
                        {
                        $ele_seq[$j]=~tr/[A-Z]/[a-z]/;
                        @frag1=@ele_seq[$count..$j];
                        }
                        print W @frag1,"\n";
                        $count=$j+1;
                }
		else
                {
#                        if (($j ne 1) && ($j ne $len_seq-1))
#                        {
                        	if ($ele_seq[$j]=~m/X/i)
                        	{
                        		$ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        		@frag1=@ele_seq[$count..$j-1];
                        		$ele_seq[$j]=~s/[a-z]//i;
                        	}
                        	else
                        	{
                        		$ele_seq[$j-1]=~tr/[A-Z]/[a-z]/;
                        		@frag1=@ele_seq[$count..$j-1];
                        	}
                        	print W @frag1,"\n";
                        	$count=$j;
#			}
                }
        }

}
$ele_seq[$len_seq-1]=~tr/[A-Z]/[a-z]/;
@frag1=@ele_seq[$count..$len_seq-1];
print W @frag1,"*\n";
}
$count=0;
}
close(W);
$comp .= $compare1;
return();
}
###### Allowance for miscleavage ###########
open(W1,">mis_comb_proteo_frag.out");
open(Read,"proteo_frag.out");
@hold=<Read>;
close(Read);
$len=@hold;
print "Enter number of miscleavage allowed (no miscleavage enter 0): ";
$miscleav = <STDIN>;
if ($miscleav > 0)
{
for($i=0;$i<$len-1;$i++)
{
        $hold[$i]=~s/\n//g;
	print W1 $hold[$i],"\n";
        $hold1=uc($hold[$i]);
        $j=$i+1;
        if ($j < $len)
        {
                $hold[$j]=~s/\n//g;
                $new=join('',$hold1,$hold[$j]);
                print W1 $new,"\n";
                if ($miscleav > 1)
                {
                $k=$j+1;
                if ($k < $len)
                {
                        $hold[$k]=~s/\n//g;
                        $new_tmp=uc($new);
                        $new1=join('', $new_tmp, $hold[$k]);
                        print W1 "$new1    \n";
                        if ($miscleav > 2)
                        {
                        $l=$k+1;
                        if ($l <$len)
                        {
                                $hold[$l]=~s/\n//g;
                                $new_tmp1=uc($new1);
                                $new2=join('', $new_tmp1, $hold[$l]);
                                print W1 "$new2\n";
                        }
                        }
                }
                }
        }
}
close(W1);
}
else
{
	readpipe("cp proteo_frag.out mis_comb_proteo_frag.out");
}
#######################################################################################

##### Separating the proteolytic peptides with and without Cys #####
my $line2; my $nium_cys; my $r=0; my $r1; my %mass=(); my @ele_frag=(); my $r2; my $sum; my $cys_counter;

open(W,">pep_nocys.out");
open(W1,">pep_cys.out");
open(Read1,"mis_comb_proteo_frag.out") or die $!;
my @hold1=<Read1>;
close(Read1);
open(Read2,"mass_list") or die $!;
my @hold_mass=<Read2>;
close(Read2);

foreach $line2(@hold_mass)
{
        $line2=~s/\n//g;
        @mass_val=split('	',$line2);
        $mass{$mass_val[0]}=$mass_val[1];
}
foreach $line2 (@hold1)
{
	$sum=0;
	$line2=~s/\n//g;
	my $apo=$line2;
	$line2=~s/\*//g;
	$num_cys=$apo=~s/C/$r/gi;
	if ($num_cys eq "")
	{
		$num_cys=0;
		@ele_frag=split('',$line2);
		$r2=@ele_frag;
		for($r1=0;$r1<$r2;$r1++)
		{
			my $hi=$ele_frag[$r1];
			$sum=$sum+$mass{$hi};
		}		
		$sum=$sum+19;
		if (($cterm_s == 1) && ($apo=~m/\*/g))
		{
			$sum=$sum-1;
		}
		print W "$line2	$sum	$num_cys\n";
	}
	else
	{
#		$cys_counter++;                
		@ele_frag=split('',$line2);
                $r2=@ele_frag;
                for($r1=0;$r1<$r2;$r1++)
                {
			my $hi=$ele_frag[$r1];
                        $sum=$sum+$mass{$hi};
                }
		$sum=$sum+18;
		if (($cterm_s == 1) && ($apo=~m/\*/g))
                {
                        $sum=$sum-1;
                }
                print W1 "$line2	$sum	$num_cys\n";
	}
}
#print $cys_counter,"\n";
close(W);
close(W1);


########################################################################################

###### Estimating the comb length cut-off #######

open(Read,"pep_cys.out");
@hold = <Read>; # read data
close(Read);
foreach (@hold)
{
#       print $_,"\n";
        ($seq,$weight,$cys_no) = split; # get score
        $weight{$_} = $weight; # record it
}
@sort_hold = sort {
$weight{$a} <=> $weight{$b};
} @hold;

open(Read1,"peak_mass_ms");
@hold1=<Read1>;
close(Read1);
$len1=@hold1;
@sort_hold1=sort {$a <=> $b} @hold1;
#print $sort_hold1[$len1-1];

my $sum=0;
my $i=0;

foreach $line(@sort_hold)
{
        $i++;
        $line=~s/\n//g;
        @var=split('\t',$line);
        $sum=$sum+$var[1];
        if (($sum >= $sort_hold1[$len1-1]) || ($sum >= ($sort_hold1[$len1-1]-1)))
        {
		$cutoff_len=$i;
#                print $i,"\n";
                last;
        }
}
if (($miscleav > 0) && ($cutoff_len > 4))
{
	$cutoff_len=4;
}
#print $cutoff_len,"\n";

###### Filtering pep_cys.out based on maximum MS #########

open(W2,">fil_pep_cys.out");
foreach $line3 (@hold)
{
	$line3=~s/\n//g;
	@var_tmp=split('\t',$line3);
	if ($var_tmp[1] < $sort_hold1[$len1-1])
	{
		print W2 "$line3\n";
		$cys_counter++;
	}
}
close(W2);

###### Getting the combination of Cys containing proteolytic fragments to match with MS data ######

##### Generating the possible combinations of Cys containing proteolytic fragments ######
$tmp_cut = $sort_hold1[$len1-1];
$tmp_cut=~s/\n//g;
$tmp_cut=~s/\s+//g;
$tmp_cut=$tmp_cut+100;
readpipe ("perl resolv_conflict.pl -f fil_pep_cys.out -d $comp -n $cutoff_len -l $tmp_cut > pep_comb_MS");
readpipe ("cat pep_comb_MS pep_nocys.out > match_MS.inp");

##### Matching mass of fragments with MS profile #####

print "Enter the mass window range for match with MS data\n";
my $window_mass=<STDIN>;

readpipe("perl match_MS.pl $window_mass");
$prot=~s/\n//g;
$prot=~s/\s+//g;
readpipe("cp match_MS.out ../Results_MS/match_MS_$prot.out");
system("perl -p -i -e 's/ \\| /X/g' match_MS.out");
readpipe("mv match_MS.out ../Results_MS/inp_MSn_match_MS_$prot.out");
#readpipe("rm -f comb_MS match_MS.inp pep_*.out pep_comb_MS fil_pep_cys.out mis_comb_proteo_frag.out");
