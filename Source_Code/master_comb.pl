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
my @aa=();

$seq=@ARGV[0];
$seq=~s/\n//g;
@ele1=split('',$seq);
$len=@ele1;

### Modification of parent sequence ###
my @mapArray = ();
my $s = S_EMPTY;

for($i=0;$i<$len;$i++)
{
#        if (($ele1[$i] eq "C") || ($ele1[$i]=~m/[a-z]/))
#        {
#		$new=$ele1[$i].$ele1[$i+1];
#                push(@new_ele1,$new);
#                $i++;
#        }
#        else
#        {
#                push(@new_ele1,$ele1[$i]);
#        }


#  for( my $i=0; $i < scalar @flds; ++$i ){
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
#print @new_ele1,"\n";
#print $len_p,"\n";

####################################################################

for($j=0;$j<$len_p;$j++)
{
	if (($new_ele1[$j] =~m/C\d/i) || ($new_ele1[$j] =~m/X/))
	{
		$len_a=@aa;
		if ($len_a > 0)
		{
			$e++;
			$pep[$j]=join('',@aa);
			push(@a,$pep[$j]);
			push(@b,"P$e");
		}
		splice(@aa,0,$len);
        	push(@b,$new_ele1[$j]);
	}
	else
	{
		$new_ele1[$j]=~s/\n//g;
		push (@aa,$new_ele1[$j]);
	}
}
#print $new_ele1[$len_p-1],"\n";

if (($new_ele1[$len_p-1] !~ "C") && ($new_ele1[$len_p-1] !~ "c"))
{
	$pep[$j]=join('',@aa);
	push(@a,$pep[$j]);
	$e++;
	push(@b,"P$e");
}
$len_f=@a;

######## Generating combinations for peptides without cysteine ########

readpipe("rm -f out_comb_P");

for($t=0;$t<$len_f;$t++)
{
	$t1=$t+1;
	$len10=length($a[$t]);
	$tmp1=$a[$t];
	$len11=$tmp1=~s/\d/$1/g;
	$len_pep=$len10-$len11;
#	print  $len_pep,"\n";
	for (  $po = 1 ;  $po <= $len_pep;  $po++  )
	{
    		readpipe("./comb_MSn -s $a[$t] -n $po >> out_comb_P");
	}
}

######## Generating Cys containing sequences as blocks with peptides #########

open(W11,">ip");

readpipe("rm -f out_comb_cys");

$new=join('',@b);
print W11 $new;
#$new=~s/X//g;
$len_cys=@b;
print "$len_cys	$new\n";
for (  $ii = 1 ;  $ii <= $len_cys;  $ii++  )
{
#    print "$new $ii\n";
#    readpipe("./comb_cys $ii $len_cys $new >> out_comb_cys");
    readpipe("./comb_MSn -s $new -n $ii -c >> out_comb_cys");
}
readpipe("perl comb_with_cys_filter.pl $new");
readpipe("perl X_filter.pl > noX_comb_with_cys");

######## Peptide blocks to be combined with Cys-containing filtered combinations #######

open(W1,">pep_frag");

for($t=0;$t<$len_f;$t++)
{
	print W1 "$a[$t]\n";
}





