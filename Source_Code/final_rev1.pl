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

$inp_seq=@ARGV[0];
$cap1=$inp_seq=~s/C/$1/gi;
#$stage=@ARGV[1];
$stage=n;

open(Read1,"out_ms$stage");
@hold1=<Read1>;
close(Read1);
open(W1,">rev_out_ms$stage");

open(Read2,"strict_out_ms$stage");
@hold2=<Read2>;
close(Read2);
open(W2,">rev_strict_out_ms$stage");

foreach $line1 (@hold1)
{
	$line1=~s/\n//g;
	@ele=split('\t',$line1);
	$cap=$ele[0]=~s/C/$1/gi;
	if ($cap == $cap1)
	{
		print W1 "$line1\tAll Cys connected\n";
	}
	else
	{
		print W1 "$line1\n";
	}
}

foreach $line2 (@hold2)
{
        $line2=~s/\n//g;
        @ele1=split('\t',$line2);
        $cap2=$ele1[0]=~s/C/$1/gi;
        if ($cap2 == $cap1)
        {
                print W2 "$line2\tAll Cys connected\n";
        }
        else
        {
                print W2 "$line2\n";
        }
}

