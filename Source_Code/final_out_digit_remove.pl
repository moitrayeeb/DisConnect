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

$stage=$ARGV[0];

open (W1,">strict_out_ms$stage_tmp");
open (W2,">out_ms$stage_tmp");

open(Read1,"strict_out_ms$stage");
@hold1=<Read1>;
close(Read1);

open(Read2,"out_ms$stage");
@hold2=<Read2>;
close(Read2);

foreach $line(@hold1)
{
$line=~s/\n//g;
@var=split('\t',$line);
@ele=split('',$var[0]);

  my @flds = @ele;
#  print scalar (@flds);
  my @mapArray = ();
  my $s = S_EMPTY;
  for( my $i=0; $i < scalar @flds; ++$i ){
  if( length($s) > 0 && $flds[$i] !~ /^\d/ ){
	push @mapArray, $s;
        $s = S_EMPTY;
   }
       $s .= $flds[$i];
  }
  if( length($s) > 0 ){
   push @mapArray, $s;
   $s = S_EMPTY;
  }
#print @mapArray,"\n";

for($j=0;$j<scalar @mapArray;$j++)
{
	$mapArray[$j]=~s/\n//g;
	if ($mapArray[$j]!~m/c/i)
	{
		$mapArray[$j]=~s/\d+//g;
		push(@newArray,$mapArray[$j]);
	}
	else
	{
		push(@newArray,$mapArray[$j]);
        }
}
$new_var=uc(join('',@newArray));
print W1 "$new_var	$var[1]	$var[2]\n";
$len_tmp=@newArray;
splice(@newArray,0,$len_tmp);
}
readpipe("mv strict_out_ms$stage_tmp strict_out_ms$stage");
close(W1);
########################################################

foreach $line1(@hold2)
{
$line1=~s/\n//g;
@var1=split('\t',$line1);
@ele1=split('',$var1[0]);

  my @flds1 = @ele1;
#  print scalar (@flds);
  my @mapArray1 = ();
  my $s = S_EMPTY;
  for( my $i=0; $i < scalar @flds1; ++$i ){
  if( length($s) > 0 && $flds1[$i] !~ /^\d/ ){
        push @mapArray1, $s;
        $s = S_EMPTY;
   }
       $s .= $flds1[$i];
  }
  if( length($s) > 0 ){
   push @mapArray1, $s;
   $s = S_EMPTY;
  }
#print @mapArray,"\n";

for($j=0;$j<scalar @mapArray1;$j++)
{
        $mapArray1[$j]=~s/\n//g;
        if ($mapArray1[$j]!~m/c/i)
        {
                $mapArray1[$j]=~s/\d+//g;
                push(@newArray1,$mapArray1[$j]);
        }
        else
        {
                push(@newArray1,$mapArray1[$j]);
        }
}
$new_var1=uc(join('',@newArray1));
print W2 "$new_var1      $var1[1] $var1[2]\n";
$len_tmp1=@newArray1;
splice(@newArray1,0,$len_tmp1);
}
readpipe("mv out_ms$stage_tmp out_ms$stage");
