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

#! /usr/bin/perl
#use strict;
#use warnings;
use Getopt::Std;

open(W1,">cys_comb.out");

use constant B_TRUE  => (1 == 1);
use constant B_FALSE => (1 == 0);
use constant MAP_PREFIX => 'P';
use constant S_EMPTY    => '';

#my %globalMap_termP = ();
my %globalMap_midP = ();
my $o;

open(Read11,"ip");
my @u=<Read11>;
close(Read11);

my $seq=$u[0];

my %hash=();
my $me=0;

my $seq1=$ARGV[4];
readpipe("perl pre_comb_gen_filter.pl $seq1 > range.out");

open(Read1,"range.out");
my @var=<Read1>;
close(Read1);
$var[0]=~s/\n//g;

sub trim
{
   my @localcpy = @_;
   my @trimcpy = ();
   foreach (@localcpy){
	next if(! defined );
	s/^\s+//g;  s/\s+$//g;
	push @trimcpy, $_ if( defined $_ and $_ ne '');
   }
   return (wantarray)?@trimcpy:$trimcpy[0];
}

sub isArray($){
 my ($ref) = @_;
 eval{ my @arr = @$ref; };  
 return (! $@ );
}

sub isInteger($){
  my ( $ref ) = @_;
  return ( defined $ref and $ref =~ /^\d+$/ );
}

sub combGen($$$$);
my $combi;
sub combGen($$$$)
{
  my ( $pstrArr, $posArr,$dir,$pos) = @_;
  die "Error: invalid array pointer !\n" if( 
					     ! isArray($pstrArr) or 
					     ! isArray($posArr) or 
                                             ! isArray($dir) or
					     ! isInteger($pos) or 
					     (
                                              ( $pos < scalar @$posArr ) and 
			#		      ! defined $globalMap_termP{uc $pstrArr->[$posArr->[$pos]]} and
                                              ! defined $globalMap_midP{uc $pstrArr->[$posArr->[$pos]]}
                                             )
					   );
                                         
if( $pos == scalar @$posArr ){
        my $new=join('',@$pstrArr);
        my $len_new=length($new);
        my $dig=$new=~s/\d/$me/g;
        my $len_new1=$len_new-$dig;
        if ($len_new1 <= $var[0]) {
        print W1 join('',@$pstrArr)," \n";
   }
  }else{  
 if (@$dir[$pos] ne 0)
  {
   $combi=join('',uc $pstrArr->[$posArr->[$pos]],@$dir[$pos]);
  }
  else
  {
   $combi=uc $pstrArr->[$posArr->[$pos]];
  }
#   print $combi,"\n";
   foreach my $mapVal (@{$globalMap_midP{$combi}})
   {
  	my @strArr = @$pstrArr;
  	$strArr[$posArr->[$pos]] = $mapVal;
	combGen(\@strArr,$posArr,$dir,$pos+1);
   }
  }
}

sub usage
{
   print STDERR <<_EOM_
   ---------------------------------------------------
   Usage: $^X $0 -p <peptide file> -m <map files> [-h]
   ---------------------------------------------------
_EOM_
;
   exit(1);
}

sub extractArray
{
  my @flds = @_;
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
  return @mapArray;
}

##############################
##	Main Programme
##############################
my %opt = ();
getopts('m:p:h', \%opt);

die usage if(defined $opt{h});

############### P neighbours #########################


my $var1; my $var2; my $p_handle;
$seq=~s/\n//g;
my @ele=split('',$seq);
my $len_inp=@ele;
for(my $j=0;$j<$len_inp;$j++)
{
        if ($ele[$j]=~m/P/)
        {
		my $var=join('',$ele[$j],$ele[$j+1]);
                if (($j-2 >= 0) && ($j-1 > 0) && ($ele[$j-1] ne "X"))
                {
                       $var1=join('',$ele[$j-2],$ele[$j-1]);
                }
                else
                {
                        $var1="";
                }
                if (($j+2 < $len_inp) && ($ele[$j+2] ne "X"))
		{
                       $var2=join('',$ele[$j+2],$ele[$j+3]);
                }
                else
                {
                        $var2="";
                }
#		print $j,$j-2 , $j-1 , $j+2, $j+3, "$var	$var1	$var2\n";
                my $neigh=join('.',$var1,$var2);
		my $len_hash=length($neigh);
#		print "$var	$neigh $len_hash\n";
		if ($len_hash == 5)
		{
			$p_handle=$p_handle.$var;
#			print "$p_handle\n";
		}
#		if (($len_hash == 2) && 
#		{
#			$p_handle_r=$p_handle_r.$var;
#		}
		else
		{
			$p_handle=$p_handle."P";
		}
                $hash{$var}=$neigh;
#                print "$hash{$var} $var\n";
	}
}

########################################################
my $num;
my ($mapfile, $peptidefile) = undef;

$mapfile = $opt{m} if( defined $opt{m} && -f $opt{m});
$peptidefile = $opt{p} if(defined $opt{p} && -f $opt{p});

die usage if( ! defined $mapfile || ! defined $peptidefile);

my $fh;
my $counter = 1;

open $fh, '<' , "$mapfile" or die "Error: opening file :$mapfile : $! \n";
while( my $line = <$fh> ){
 chomp $line;
 my @flds = trim split(//, $line );
 next if( scalar @flds == 0 );
 my $ptag = MAP_PREFIX . "$counter";
 my $ptag1=MAP_PREFIX . "$counter"."_l";
 my $ptag2=MAP_PREFIX . "$counter"."_r";
# print "$ptag	$ptag1	$ptag2\n";
 $counter++;
 my @mapArray = extractArray(@flds);
 my @reverseMapArray = reverse @mapArray;
 my @intArray;
# print "0 $mapArray[0] $reverseMapArray[0]\n";
 for( my $i=1; $i < scalar @mapArray; ++$i ){
   $mapArray[$i] = $mapArray[$i - 1].$mapArray[$i];
   $reverseMapArray[$i] .= $reverseMapArray[$i - 1];
  }

if ($p_handle =~m/$ptag/)
{
# print $ptag,"\n";
 for( my $i=0; $i < scalar @mapArray-1; ++$i ){
   for( my $j=0; $j < scalar @mapArray-1; ++$j ){
   my $newline=join('',$mapArray[$i],$reverseMapArray[$j]);
    if ($newline eq $line){
     last;
    }
    else {
  #   print $i,"\n";
      push (@intArray,$newline) ;
    }
   }   
 }
   $globalMap_midP{$ptag} = [ @mapArray, @reverseMapArray, @intArray];
   $globalMap_midP{$ptag1} = [ @mapArray];
   $globalMap_midP{$ptag2} = [ @reverseMapArray];
   my $len_int=@intArray;
   splice(@intArray,0,$len_int);
}
else
{
   $globalMap_midP{$ptag} = [ @mapArray, @reverseMapArray];
   $globalMap_midP{$ptag1} = [ @mapArray];
   $globalMap_midP{$ptag2} = [ @reverseMapArray];
}
}
close $fh;


my $count_p=0; my $var_p=""; my $var1_p=""; my $var2_p=""; my $neigh_p="";  my @hash_var=();

open $fh, '<', "$peptidefile" or die "Error: opening file :$peptidefile : $! \n"; 

while( my $line = <$fh> ){
 $num++;
 chomp $line;
# print "$line\n";
 my @flds = trim split(//, $line );
 next if( scalar @flds == 0 );
 print "Processing : $line $num\n";
 my @mapArray = extractArray(@flds); 
 my @positions = ();
 my @dir_counter=();
 my $pattern = MAP_PREFIX . '\d+'; 
 for( my $i = 0; $i < scalar @mapArray; ++$i ){
#	print $mapArray[$i],"	",$i,"\n";
   if( $mapArray[$i] =~ /$pattern/i )
    {
	push @positions, $i;
        $var_p=join('',$mapArray[$i]);
#	print "$var_p	$i\n";
        if (($i-1 >= 0) && ($i+1 < scalar @mapArray))
        {
         $neigh_p=join('',$mapArray[$i-1],$mapArray[$i+1]);
#	 print "$neigh_p $hash{$var_p}\n";
	 my $tmp1=$hash{$var_p};
	 $tmp1=~s/\.//g;
	 if ($neigh_p eq $tmp1)
	{
	 $count_p++;
	 push(@dir_counter,0); 
#	 print "hit1\n";
	}
	}
        @hash_var=split('\.',$hash{$var_p});
#	print $hash{$var_p},"\n";
#	print $i-1," ","$mapArray[$i-1]"," ",$hash_var[0],"\n";
	if ((($i-1 < 0) || ($mapArray[$i-1]=~ /$pattern/i)) || ($mapArray[$i-1] ne $hash_var[0]))
	{
	  push(@dir_counter,"_r");
#	  print "hit2\n";
	}
#	print $i+1," ","$mapArray[$i+1]"," ",$hash_var[1],"\n";
        if ((($i+1 >= scalar @mapArray) || ($mapArray[$i+1]=~ /$pattern/i)) || ($mapArray[$i+1] ne $hash_var[1]))
	{
	  push (@dir_counter,"_l");
#	  print "hit3\n";
	}
#	else
#	{
#	  push(@dir_counter,0);
#	  print "hit4\n";
#	}
   #  print MAP_PREFIX.$i,"\n";
    }     
#    else
#        {
#          push(@dir_counter,0);
#          print "hit4\n";
#        }
 $count_p=0;
 }
#    print scalar @dir_counter, " ","@dir_counter"," ",@positions;
#    print "\n";
 combGen(\@mapArray, \@positions,\@dir_counter,0); 
} 
close $fh;
close(W1);

readpipe("cat out_comb_P cys_comb.out > out_comb");
