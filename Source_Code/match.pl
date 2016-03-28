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

# User defined m/z values for the MS/MS fragment ions are incorporated (peak_mass) and these values are then compared with that of every possible fragment #

open(W,">match_mass");

open(Read,"peak_mass");
@var=<Read>;
close(Read);

open(R1,"tot_mass_final");
@hold=<R1>;
close(R1);

foreach $line (@var)
{ $line=~s/\n//g;
  $line=~m/\./g;
  $m_line=$`;
#  print $m_line,"\n";
  foreach $line1(@hold)
{  $line1=~s/\n//g;
   @map=split('\t',$line1);
if ($map[1]=~m/\./g)
{
   $m_map[1]=$`;
}
else
{
   $m_map[1]=$map[1];
}
#   print "$m_map[1]\n";
   $map1[1]=$m_map[1]+2;
   $map2[1]=$m_map[1]+3;
   $map4[1]=$m_map[1]+4;
   $map3[1]=$m_map[1]+1;
#   print "$map1[1] $map2[1] $map4[1] $map3[1] $m_line\n";
   if (($m_map[1]=~m/^$m_line$/) || ($map1[1]=~m/^$m_line$/) || ($map2[1]=~m/^$m_line$/) || ($map3[1]=~m/^$m_line$/) || ($map4[1]=~m/^$m_line$/))
{ $line1=~s/^ //g;
  print W $line1,"\n"; }
}
}

