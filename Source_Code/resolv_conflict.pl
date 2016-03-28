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
use strict;
use warnings;
use Getopt::Std;

use constant TAG => 'A';
use constant B_TRUE  => ( 1 == 1 );
use constant B_FALSE => ( 1 == 0 );
use constant COMBINATION_EXE => './comb';

sub trim
{
  my @localcpy = @_;
	my @trimcpy = ();
	foreach (@localcpy){
		next if( ! defined );
		s/^\s+//;  s/\s+$//;
    push @trimcpy , $_ if( defined and $_ ne '');
	}
	return (wantarray)?@trimcpy:$trimcpy[0];
}

sub usage
{
	print STDERR <<_EOM_
	------------------------------------------------
  $^X $0 -f <input file>  -d <cut delimitor> -n <max tuple size> -l <mass limit> [-h]
	------------------------------------------------
_EOM_
;
 exit 1;
}
my ($pep_cnt_copy , $cys_cnt_copy) = undef;

sub recursion($$$$$$$$$$)
{
	my ($map_hash,$peptide_hash,$comb_ele,$array_count,$str,$pep_mass,$pep_cnt,$limit,$cys_num,$cys_cnt)=@_;	
	if ($pep_cnt < $limit)
	{
		if (scalar @$comb_ele == $array_count)
		{
			if (($cys_cnt_copy % 2) == 0)	
			{
				print $str,"\t",$pep_cnt_copy+1,"\t","$cys_cnt_copy\n";
			}
			return;	
		}
		my @local_array = (${$comb_ele}[$array_count]);
		push @local_array , @{$map_hash->{${$comb_ele}[$array_count]}};
		foreach my $line (@local_array)
		{
			my $local_copy = $str;
			$local_copy = $local_copy." | ".$peptide_hash->{$line};
			$local_copy =~s/^ \| //g;
			$pep_cnt_copy = $pep_cnt + $pep_mass->{$line};
			$cys_cnt_copy = $cys_cnt + $cys_num->{$line};
			&recursion ($map_hash,$peptide_hash,$comb_ele,$array_count + 1,$local_copy,$pep_mass,$pep_cnt_copy,$limit,$cys_num,$cys_cnt_copy);
		}
	}
}

my %opts = ();
getopts('f:d:n:l:h', \%opts);

die usage if(defined $opts{h});

my ($filename , $delim , $new_delim, $tuple , $limit) = undef;
my @ele_delim = ();
my @peptides = ();
my @pep_mass = ();
my @cys_num = ();
my %peptide_m = ();
my %cys_count = ();
my %peptide_map = ();
my %conflict_list = ();
my %conflict_map = ();

$filename = $opts{f} if(defined $opts{f} && -f $opts{f});
$delim = $opts{d} if(defined $opts{d});
@ele_delim=split('',$delim);
my %hashTemp = map { $_ => 1 } @ele_delim;
my @array_out = sort keys %hashTemp;
$new_delim=join('\|',@array_out);
#print $new_delim,"\n";
$tuple = $opts{n} if (defined $opts{n} && $opts{n} =~ /^\d+$/ && $opts{n} > 1);
$limit = $opts{l};

die usage if( ! defined $filename || ! defined $delim); 


my $fh;
open $fh, '<' , "$filename" || die "Error: opening file $filename : $! \n";

while( my $line = <$fh>)
{
	chomp $line;
	my @flds = trim split(/\s+/, $line);
	push @peptides , uc $flds[0] if( @flds && scalar @flds > 1 && $flds[0] =~ /^[A-Z]+$/i);
	push @pep_mass, $flds[1] if ( @flds && scalar @flds > 1 );
	push @cys_num, $flds[2] if ( @flds && scalar @flds > 1 );
}
close $fh;
foreach my $i ( 0 .. $#peptides )
{
	my $tag  = TAG . "$i";
	$peptide_map{$tag} = $peptides[$i]; 
	$peptide_m{$tag} = $pep_mass[$i];
	$cys_count{$tag} = $cys_num[$i]; 
#	print "$peptide_m{$tag}\n";
	$conflict_list{$tag} = {};
	my @localsplits = split( /$new_delim/ , $peptides[$i]);
	foreach (@localsplits)
	{
#		print $_,"\n";
		$conflict_list{$tag}{$_} = 1;
	}
}

foreach my $tag_i ( keys %peptide_map )
{
	$conflict_map{$tag_i} = {};
	foreach my $tag_j ( keys %peptide_map )
	{
		next if( "$tag_i" eq "$tag_j");
		foreach my $str ( keys %{ $conflict_list{$tag_i} } )
		{
			if( $peptide_map{$tag_j} =~ /^$str/ )
			{
				$conflict_map{$tag_i}{$tag_j} = 1;
			}
		}
	}
}

my $sstr = '';
my @bag = ();
my %map = ();
foreach my $tag ( sort { substr($a,1) <=> substr($b,1)} keys %peptide_map )
{
#	print $tag,"\n";
	my $conflict = B_FALSE;
	foreach my $bin (@bag)
	{
#		print $bin,"\n";
		if( defined $conflict_map{$bin}{$tag} or defined $conflict_map{$tag}{$bin})
		{
			$conflict = B_TRUE; 
			$map{$bin} = [] if(! defined $map{$bin});
			push @{$map{$bin}}, $tag;
		}
	}
	if( ! $conflict )
	{
		$map{$tag} = [];
		push @bag, $tag;
		$sstr .= $tag;
	}
}
my $ii;my @result_str1 = (); 
#print "$sstr\n";

for ($ii=1;$ii<=$tuple;$ii++)
{
	my $cmd = "" . COMBINATION_EXE . " -s $sstr -n $ii -c";
	my @result_str = qx/$cmd/;
	chomp @result_str;
	@result_str = trim @result_str;
	push (@result_str1,@result_str);
}

foreach my $str(@result_str1)
	{
#		print $str,"\n"; 
		my $t = TAG;
   		my @elems = map { TAG . $_ } trim split( /$t/ ,$str);
		recursion(\%map,\%peptide_map,\@elems,0,'',\%peptide_m,0,$limit,\%cys_count,0);
	}

=foreach my $tag (@bag){
	print "$tag : [",join(' ',@{$map{$tag}} ),"] \n";
}
=cut



