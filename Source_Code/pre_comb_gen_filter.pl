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

open(Read1,"peak_mass");
@var_match=<Read1>;
close(Read1);
$len_match=@var_match;
@sort_var_match=sort{$a<=>$b}@var_match;
$match_max=$sort_var_match[$len_match-1];
$match_min=$sort_var_match[0];

my %mass=();
my $sum=0;

$seq=@ARGV[0];
$seq=~s/\d//g;
$seq=~s/X//g;
my @ele=split('',$seq);
my $len=@ele;

open(Read,"<mass_list_1st");
@var=<Read>;
close(Read);

foreach $line1(@var)
{
        $line1=~s/\n//g;
        @mass_val=split('\t',$line1);
        $mass{$mass_val[0]}=$mass_val[1];
}

#print $seq,"\n";

for($i=0;$i<$len;$i++)
{
	$r=$ele[$i];
	push(@a,$mass{$r});
}
@sort_a=sort{$a<=>$b}@a;

for($i=0;$i<$len;$i++)
{
	$sum=$sum+$sort_a[$i];
	if ($sum >= $match_max)
	{
		print $i+1,"\n";
		last;
	}
}
if ($sum <= $match_max)
{
	print $i+1,"\n";
}
$sum=0;

