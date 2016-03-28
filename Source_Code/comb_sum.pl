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

# All possible permutations of four possible Cys residue masses are added for every Cys present in a given combination of sequence #

readpipe("rm -f cys_mass* tot_mass_final");
open(W,">>tot_mass_final");

open(INIS_mb,"<@ARGV");
@hold=split('',@ARGV[0]);
$meaow1=@ARGV[0]=~s/C/$5/gi;
$z=$meaow1/2;
my $x=0;
my $pl=0;
open(Read,"sum_out");
@var=<Read>;
close(Read);
foreach $line (@var)
{	$line=~s/\n//g;
	@ele=split('	',$line);
 	push (@a,$ele[2]);
  	push (@m,$ele[1]);
  	push(@seq,$ele[0]);
}
@sort=sort {$a <=> $b} @a;
$len=@sort;
$e=$len-1;
$max=$sort[$e];

for ($no1=0;$no1<=$max;$no1++)
{	$x=$no1;
	open(W2,">>cys_mass2");
	@mass_c=(69,135,103,101);
	my $wa=69.0; my $wb=135.0; my $wc=103.0; my $wd=101.0; my $a=0;
	for ($i=0;$i<=$x;$i++)
	{       $r=$x-$i;
        	for ($j=0;$j<=$r;$j++)
        	{       $y=$x-$i-$j;
                	for ($o=0;$o<=$y;$o++)
                	{       $e=$x-$i-$j-$o;
                        	for ($p=0;$p<=$e;$p++)
                        	{       $f=$i+$j+$o+$p;
                                	if (($f == $x) && ($i <= $z) && ($j <= $z) && ($o <= $z) && ($p <= $z))
                               		{       print W2 "$i   $j      $o      $p\n";
#                                        	$Wcys=($i*$wa)+($j*$wb)+($o*$wc)+($p*$wd);
						$Wcys[$no1]->[$k1]=($i*$wa)+($j*$wb)+($o*$wc)+($p*$wd);
                                      		$pos[$no1]->[$k1]=$i.$j.$o.$p;
#						print "$no1	$pos[$no1]\n" ;
						$k1++;
                                	}
                        	}	
                	}
        	}
	}
	$k1=0;
}

foreach $arr(@Wcys)
{
	$lo[$pl]=@{$arr};
	$pl++;
}

#print $lo[4];

for($rt=0;$rt<$len;$rt++)
{	$x=$a[$rt]; 
        $lo1=$lo[$x];
#	print "$x	$lo1\n";
	for($u=0;$u<$lo1;$u++)
	{#	print $Wcys[$x][$u],"\n";
  		$sum1=$m[$rt]+$Wcys[$x][$u];
   		@gh=split('',$pos[$x][$u]);
   		$pen1=$gh[0]*$mass_c[0];
   		$pen2=$gh[1]*$mass_c[1];
  		$pen3=$gh[2]*$mass_c[2];
   		$pen4=$gh[3]*$mass_c[3];
   		@kal=split('',$seq[$rt]);
   		$len_k=@kal;
		for ($hj=0;$hj<$len_k;$hj++)
		{
			if ($kal[$hj]=~m/[a-z]/g)
			{	$count_hj++;	 
			}
		}
		$sum2=$sum1+($count_hj*17);
		$count_hj=0;
		print W "$seq[$rt]\t$sum2\t$x\t$pos[$x][$u]\t$pen1|$pen2|$pen3|$pen4\n"; 
	}
   		$sum1=0;
   		$sum2=0;
}
