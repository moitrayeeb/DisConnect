###############################################################################
 # Copyright (c) 2012, Moitrayee Bhattacharyya, Kallol Gupta, P. Balaram, Indian Institute of Science, Bangalore, India
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

$str=$ARGV[0];
@tmp=split('',$str);

for($j=0;$j<scalar(@tmp);$j++)
{
	if ($tmp[$j] =~m/C/)
	{
		$k++;
		push(@tmp_arr,$tmp[$j],$k);
	}
	else
	{
		push(@tmp_arr,$tmp[$j]);
	}
}
$str1=join('',@tmp_arr);
#print $str1;

open(Read,"rev_out");
@hold=<Read>;
close(Read);

foreach $line1 (@hold)
{
$line1=~s/\n//g;
@str2=split('\t',$line1);
#$str2[0]=~s/\d//g;
#print $str2[0],"\n";
@ele=split(' \| ',$str2[0]);
foreach $line (@ele)
{
	$line=~s/\n//g;
	if ($str1=~m/$line/) 
	{
#		print $line,"\n";
		$a=$`;
		$b=$';
		push (@arr,$a,$b);
	}
}
#print @arr,"\n";
=commentfor($i=0;$i<scalar(@arr);$i++)
{
	print $arr[$i],"\n";
}
=cut

if ($arr[0] ne "")
{
	push(@newarr,$arr[0]);
	if ($arr[0]=~m/C/gi)
	{
		$sum++;
	}
}

### Longest substring problem ###

use Data::Dumper;
for($i=1;$i<scalar(@arr)-1;$i++)
{
$tmp1=substr($arr[$i+1],length($arr[$i+1])-1,1);

my $needle = $arr[$i];
my $haystack = $arr[$i+1];
my @matches;

for my $start (0..length $needle) {
    for my $len ($start+1 .. length $needle) {
        my $substr = substr($needle, $start, $len);
	my $tmp=substr($substr,length($substr)-1,1);
        push @matches, $haystack =~ m[($substr)]g if ($tmp =~m/$tmp1/g);
    }
}
my ($len, $longest) = 0;
length > $len and ($longest, $len) = ($_, length) for @matches;
#print $longest,"\t";
if ($longest=~m/^X$/)
{
}
if ($longest=~m/^.*X.*/)
{
	@ele_long=split('X',$longest);
	if ($ele_long[0] ne "")
	{
		push(@newarr,$ele_long[0]);
	}
	if ($ele_long[1] ne "")
        {
                push(@newarr,$ele_long[1]);
        }
	if ($ele_long[0]=~m/C/gi)
	{
		$sum++;
	}
	if ($ele_long[1]=~m/C/gi)
	{
		$sum++;
	}
}
else
{
	push(@newarr,$longest);
	if ($longest=~m/C/gi)
	{
		$sum++;
	}	
}
$i++;
}
if ($arr[scalar(@arr)-1] ne "")
{
	push(@newarr,$arr[scalar(@arr)-1]);
	if ($arr[scalar(@arr)-1]=~m/C/gi)
	{
		$sum++;
	}
}
#$sum=$cys_count1+$cys_count2+$cys_count3+$cys_count4+$cys_count5;
#print @newarr,"\t",scalar(@newarr),"\t",$sum,"\n";
if (scalar(@newarr) le $sum+1)
{
	print $line1,"\n";
}
@newarr=();
splice(@arr,0,scalar(@arr));
$sum=0;
}
