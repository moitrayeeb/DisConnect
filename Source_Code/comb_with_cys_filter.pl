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

open(W,">comb_with_cys");

my %hash=();

$seq=@ARGV[0];
chomp $seq;
#print $seq,"\n";
@ele=split('',$seq);
$len_inp=@ele;
for($j=0;$j<$len_inp;$j++)
{
	if ($ele[$j]=~m/P/)
	{	
		$var=join('',$ele[$j],$ele[$j+1]);
		if (($j-2 >= 0) && ($j-1 > 0) && ($ele[$j-1] ne "X"))
		{
			$var1=join('',$ele[$j-2],$ele[$j-1]);
		}
		if (($j+2 < $len_inp) && ($ele[$j+2] ne "X"))
		{
			$var2=join('',$ele[$j+2],$ele[$j+3]);
		}
		$neigh=join('',$var1,$var2);
		$hash{$var}=$neigh;
#		print $hash{$var},"\n";
	}
}

open(Read,"out_comb_cys");
@hold=<Read>;
close(Read);

foreach (@hold)
{
	chomp;
	$line=$_;
	$_=~s/\d//g;
	$line=~s/^\d.*\. //g;
	if ($line=~m/P/)
	{
		@ele1=split('',$line);
		$len1=@ele1;
		for($i=0;$i<$len1;$i++)
		{
			if ($ele1[$i]=~m/P/)	
			{
				$count_p++;
				$var_p=join('',$ele1[$i],$ele1[$i+1]);
	                	$var1_p=join('',$ele1[$i-2],$ele1[$i-1]);
	                	$var2_p=join('',$ele1[$i+2],$ele1[$i+3]);
				$neigh_p=join('',$var1_p,$var2_p);
				if (($hash{$var_p}=~m/$var1_p/) || ($hash{$var_p}=~m/$var2_p/) || ($neigh_p=~m/$hash{$var_p}/))
				{
					$counter++;
				}
			}
		}
		if ($counter == $count_p)
		{	
			print W $line,"\n";
		}
		$counter=0;
		$count_p=0;
	}
	else
		{
			print W $line,"\n";
		}
}
#readpipe("rm -f out_comb_cys");	
