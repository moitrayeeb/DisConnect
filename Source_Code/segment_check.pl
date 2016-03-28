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

use constant S_EMPTY    => '';

open(W,">remap1.out");
open(W2,">sing_seq");
open(Read,"match_mass");
@hold=<Read>;
close(Read);

my $seq1=@ARGV[0];
@ele1=split('',$seq1);
$len1=@ele1;

### Modification of parent sequence ###
for($i=0;$i<$len1;$i++)
{
   if( length($s) > 0 && $ele1[$i] !~ /^\d/ ){
        push @new_ele1, $s;
        $s = S_EMPTY;
   }
   $s .= $ele1[$i];
  }
  if( length($s) > 0 ){
   push @new_ele1, $s;
   $s = S_EMPTY;
  }
$len_p=@new_ele1;
print $len_p,"\n";
#######################################################
my @result ;
foreach $line (@hold)
{
#	print $line,"\n";
	$line =~s/\n//g;
	@var=split('	',$line);
	$var[4]=~s/\|/\-/g;
   
### Calculation begins on each fragment ###
	my $seq2=$var[0];
	my $mass_tot=$var[1];
	my $comb_cys=$var[4];
	$comb_cys=~s/\-/\|/g;
	@ele2=split('',$seq2);
	$len2=@ele2;
	$k=0;

### Modification of fragment sequence ###
	for($i=0;$i<$len2;$i++)
	{
		$ele2[$i]=~s/ //g;
   		if( length($s) > 0 && $ele2[$i] !~ /^\d/ ){
        	push @new_ele2, $s;
        	$s = S_EMPTY;
   	}
   	$s .= $ele2[$i];
  	}
  	if( length($s) > 0 ){
   	push @new_ele2, $s;
   	$s = S_EMPTY;
  	}
	$len_f=@new_ele2;
#	print "$len_f @new_ele2\n";

### Calculations of index of alignment ###
	my @position=();
	for ($i=0;$i<$len_f;$i++)
	{
		for ($j=0;$j<$len_p;$j++)
		{
			if ($new_ele2[$i] eq $new_ele1[$j])
			{
				$index_res=$i;
				$index_match=$j;
				$position[$i]->[$k]=$index_match;
				$k++;
			}
		}
	#	print "$position[$i]->[$k]\n";
		$k=0;
	}

### Getting all the comparisons ###
	@result = () ;
	foreach my $arr (@position)
	{
	#	print @{$arr},"\n";
		@result = map {@{conct($_)} }  @{$arr} ;
	}
	foreach (@result) 
	{
		@adj_match=split(' ',$_);
		$len_mb=@adj_match;
		for($chk=0;$chk<$len_mb;$chk++)
		{
			$tst_san=$adj_match[$chk+1]-$adj_match[$chk];
			if ($tst_san > 0)
			{	$count++;
			}
		}
		if ($count == $len_mb-1)
		{

### Checking for continuous and discontinuous segments ###
			for($t=0;$t<$len_mb-1;$t++)
			{
				$diff=$adj_match[$t+1]-$adj_match[$t];
				if ($diff == 1)
				{
					$y1=$adj_match[$t];
					push(@cat,$new_ele1[$y1]);
				}
				else
				{
					$y3=$adj_match[$t];
					push(@cat,$new_ele1[$y3]);
					$break=" | ";
					push(@cat,$break);
				}
			}
			push(@cat,$new_ele2[$len_f-1]);
			$out=join('',@cat);
#			print "$out	$mass_tot	$comb_cys\n";
			$line2=join('	',$out,$mass_tot,$comb_cys);
			$len_cat=@cat;
			splice(@cat,0,$len_cat);
#        splice(@count_cys,0,$len4);
#        }
#        $count=0;


############################################################

### Further calculations ###
	$line2=~s/\n//g;
	@kal1=split('\t',$line2);
	$gag[0]=$kal1[0];
	$kal1[0]=~s/ \| /\t/g;
	@frag=split('\t',$kal1[0]);
	$len1=@frag;
	$len=$gag[0]=~s/\|/$1/g;
	$len_c=$gag[0]=~s/C/$2/gi;
	for($r=0;$r<$len1;$r++)
	{
		if ($frag[$r] =~m/C/i)
		{
			$count1++;
		}
	}
	if ($len > 0)
	{
		if (($len_c ge (2*$len)) && ($count1 == $len1)) 
		{
			for($rt=0;$rt<$len1;$rt++)
			{
				$count_go[$rt]=$frag[$rt]=~s/C/$1/gi;
				push(@count_cys,$count_go[$rt]);
			}
			$new_cys_no=join('	',@count_cys);
			$len4=@count_cys;
			print  W "$line2	$new_cys_no\n";
			$k++;
		}
	}
	else
	{
		print W2 "$line2\n";
	}
	$count1=0;
	splice(@count_cys,0,$len4);
	}
	$count=0;
       }
        splice(@new_ele2,0,$len_f);

}
close(W);
close(W2);
system("cat remap1.out sing_seq > remap2.out");

### Subroutine 1 ###
sub conct 
{
#	print "$#result\n";
	if ($#result == -1 ) 
	{
		[$_[0]] ;
    }
    else 
	{
        [ map {$_." ".$_[0]} @result ] ;
    }
}
