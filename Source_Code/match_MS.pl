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

#/use/bin/perl

#print "Enter the mass window range for match with MS data (default 0.2)\n";
#my $win=<STDIN>;

my $win=$ARGV[0];

my $j=0;
my $i;
open(W,">match_MS.out");
print W "MS_(M+H)+\tPeptide_frag\tNum_of_Cys\n";
open(Read,"match_MS.inp") or die $!;
my @mass=<Read>;
close(Read);
open(Read1,"peak_mass_ms") or die $!;
my @hold=<Read1>;
close(Read1);
my $len2=@hold;

foreach $line (@mass)
{
	$line=~s/\n//g;
	my @ele=split('\t',$line);
#	print $ele[1],"\n";
	my $match1=$ele[1]+$win;
	my $match2=$ele[1]-$win;
#	print "$match1	$match2\n";
	for($i=0;$i<$len2;$i++)
	{
		$hold[$i]=~s/\n//g;
#		print "$match1 $hold[$i] $match2\n";
		if (($hold[$i] <= $match1) && ($hold[$i] >= $match2))
		{
			$j++;
		}
	}
	if ($j > 0)
	{
		$ele1_tmp = uc ($ele[0]);
		print W "$ele[1]\t$ele1_tmp\t$ele[2]\n";
	}
	$j=0;
}
readpipe("sort -n match_MS.out > a");
readpipe("mv a match_MS.out");

