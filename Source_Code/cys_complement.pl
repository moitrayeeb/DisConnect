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

# This code retains the complementarity in the fragment ions considering both the present and the lost Cys residues #

open(W1,">final_list1.out");
open(Read,"final_list.out");
open(INIS,"<@ARGV");
@hold=split('',@ARGV[0]);
$len=@hold;
for ($j=0;$j<$len;$j++)
{
  if ($hold[$j]=~m/C/)
{ $too++;}}
$len1=$len-$too; $seq=@ARGV[0];
if (@ARGV[0]=~m/c/)
{ $too1=$too+1;}
else
{$too1=$too;}
while (chomp (my $line = <Read>))
{ @var=split('\t',$line);
  if ($var[2] ne "")
{  $gag1=$var[0];
   $var[0]=~s/ \| //g;
   $gag=$var[0]; $count_cys=$gag=~s/C/$1/gi; $x=$too1-$count_cys;
   @ele=split('\|',$var[2]);
   $sum_cys=$ele[0]+$ele[1]+$ele[2]+$ele[3];
   $count_pipe=$gag1=~s/ \| /$2/g;
   $sum_disul=$count_pipe*204;
   $sum_tot= $sum_cys+$sum_disul;
#   print $var[2],"\n";
@mass_c=(69,135,103,101);
my $wa=69; my $wb=135; my $wc=103; my $wd=101; my $a=0;
for ($i=0;$i<=$x;$i++)
{  $r=$x-$i;
    for ($j=0;$j<=$r;$j++)
{   $y=$x-$i-$j;
for ($o=0;$o<=$y;$o++)
{   $e=$x-$i-$j-$o;
for ($p=0;$p<=$e;$p++)
{   $f=$i+$j+$o+$p;
if ($f >= $x)
{   $Wcys[$a]=($i*$wa)+($j*$wb)+($o*$wc)+($p*$wd);
    $a++;
}}}}}
for($u=0;$u<$a;$u++)
{  print $Wcys[$u],"\n";
   $sum=$sum_tot+$Wcys[$u];
   $n_sum=$too1*102;
if ($sum eq $n_sum)
{ 
print W1 $line,"\n"; 
}}}
else
{ 
print W1 $line,"\n";
}
}


