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

# For every possible combination, this code calculates the sum of residue masses #

open(W,">sum_out");
open(W1,">cys_num");
$mb=0;

my %mass;

open(Read1,"peak_mass");
@var_match=<Read1>;
close(Read1);
$len_match=@var_match;
@sort_var_match=sort{$a<=>$b}@var_match;
$mb_match=$sort_var_match[$len_match-1];
#print $mb_match;

open(Read,"<mass_list");
@var=<Read>;
close(Read);

open(R,"out_comb");
@hold=<R>;
close(R);

foreach $line1(@var)
{
	$line1=~s/\n//g;
	@mass_val=split('	',$line1);
	$mass{$mass_val[0]}=$mass_val[1];
#	print "$mass{$mass_val[0]} $mass_val[1]\n";
}


foreach $line (@hold)
{
	$line=~s/\n//g;
    	@ele1=split('\d',$line);
    	$len=@ele1;
#	print "@ele1 $len\n";
	for ($i=0;$i<$len;$i++)
	{
		$hi=$ele1[$i];
#		print $hi,"\n";
		#print "$hi $mass{$hi}","\n";
		$sum=$sum+$mass{$hi};		
		if ($ele1[$i]=~m/C/i)
		{
   			$mb++;
		}	
	}
	if ($sum <= $mb_match+5)
	{
		print W "$line	$sum	$mb\n";
	}
#print W1 "$mb\n";
$sum=0;
$mb=0;
}
