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

# This is the masterscript to run all codes #


readpipe("rm -f tot_mass_final peak_mass* out_comb* sum_out cys_num cys_mass* match_mass cys_count_finale remap*.out sing_seq post_fil_sum_out final_list*.out out* map_mass tmp1 comb_with_cys pep_frag cys_comb.out mass_range range.out ip noX_comb_with_cys");

print "Enter 1 for monoisotopic mass\n";
print "Enter 2 for average mass\n";
$mass_choice=<STDIN>;

if ($mass_choice == 1)
{
	readpipe("cp mass_list_monoisotopic mass_list");
}
else
{
	readpipe("cp mass_list_average mass_list");
}

print "Enter 1 if you are calculating on whole sequence\n";
print "Enter 2 if you are calculating on daughter fragments\n";
$num1=<STDIN>;

print "Enter 1 for smart calculation\n";
print "Enter 2 for rigorous calculation\n";
$num_rigor=<STDIN>;

open(W,">mass_range");	
print "Enter the error tolerance range (+/- Da)\n";
$window=<STDIN>;
print W $window;
close(W);

#############CALCULATION ON MS2 FRAGMENT IONS############

if ($num1 == 1)
{	print "Enter the sequence: ";
  	$name_tmp=<STDIN>;
	my $counter=1;
	my $count_c=1;
	$seq=$name_tmp;
	$seq=~s/\n//g;
	@ele=split('',$seq);

	for($i=0;$i<scalar @ele;$i++)
	{
        	if ($ele[$i]=~m/C/i)
        	{
                	$new=$ele[$i].$count_c;
                	push @array,$new;
                	$count_c++;
        	}
        	else
        	{
                	$new=$ele[$i].$counter;
                	push (@array,$new);
                	$counter++;
        	}
	}
        $len_tmp_array = @array;
        for($lta = 0; $lta <$len_tmp_array;$lta++)
        {
                $array[$lta]=~s/\n//g;
                if ($array[$lta] =~m/X/)
                {
                        $array[$lta-1] = lc($array[$lta-1]);
                }
        }
        $tmp = lc(pop(@array));
        push (@array,$tmp);
	$name=join('',@array);
	#print $name,"\n";
	$namet=$name;
	$namet=~s/\d//g;
  	$name_p=$name;
  	$num_cterm=$name_p=~s/[a-z]/$1/g;
  	readpipe("cp ../peak_mass_ms2 peak_mass");
	readpipe("perl rescale_peak.pl");
  	readpipe("perl master_comb.pl $name");
	readpipe("perl replace_list.pl -p noX_comb_with_cys -m pep_frag $name");
	readpipe("perl duplicate_rem1.pl");	
  	readpipe("perl sum_out.pl");
	if ($num_rigor == 1)
        {
                readpipe("perl comb_sum_fast.pl $name");
        }
        else
        {
                readpipe("perl comb_sum.pl $name");
        }
	print "\nEnter 1 for normal peptide fragments\n";
  	print "Enter 2 for peptide fragments with CO loss\n";
  	print "Enter 3 for normal peptide with H2O loss\n";
 	print "Enter 4 for normal peptide with NH3 loss\n";
  	print "Enter 5 for normal peptide with NH3+CO loss\n";
  	$user_in=<STDIN>;
  	open(Read,"peak_mass");
  	@hold=<Read>;
  	$len=@hold;
  	open(W,">peak_mass1");
		if ($user_in == 1)
		{
   			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+0;
   				print W $line1,"\n";
			}
		}
		if ($user_in == 2)
		{
   			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+28;
   				print W $line1,"\n";
			}   
		}
		if ($user_in == 3)
		{
     			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+18;
   				print W $line1,"\n"; 
			}
		}
		if ($user_in == 4)
		{
			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+17;
   				print W $line1,"\n";
			}
		}		
		if ($user_in == 5)
		{
       			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+45;
   				print W $line1,"\n";
			}
		}


	readpipe("cp peak_mass1 peak_mass");
	readpipe("perl match.pl");
	readpipe("perl segment_check.pl $name");
	readpipe("perl duplicate_rem2.pl");
	for($kgm=1;$kgm<=$num_cterm;$kgm++)
	{
  		print "Enter mass to be added for C-term $kgm status\n";
	}
	
	readpipe("perl post_filter.pl $name");
	readpipe("perl duplicate_rem3.pl");
	readpipe("perl cys_complement.pl $name");
	readpipe("perl duplicate_rem4.pl");
	readpipe("cp final_list2.out out_ms2");
	readpipe("perl post_fil_disul_frag_filter.pl out_ms2 > strict_out_ms2");
	readpipe("cp strict_out_ms2 inp_MSN_strict_out_ms2");
	readpipe("cp out_ms2 inp_MSN_out_ms2");
        system("perl -p -i -e 's/ \\| //g' inp_MSN_strict_out_ms2");
        system("perl -p -i -e 's/ \\| //g' inp_MSN_out_ms2");
	readpipe("perl final_out_digit_remove.pl 2");
#	readpipe("final_rev1.pl $name 2");

########################Final Filter###################################################

$inp_seq=$name;
$cap1=$inp_seq=~s/C/$1/gi;
$stage_t=2;

open(Read1,"out_ms$stage_t");
@hold1=<Read1>;
close(Read1);
open(W1,">rev_out");

open(Read2,"strict_out_ms$stage_t");
@hold2=<Read2>;
close(Read2);
#open(W2,">rev_strict_out_ms$stage_t");

if ($num_rigor == 2)
{
foreach $line1 (@hold1)
{
        $line1=~s/\n//g;
        @ele_t=split('\t',$line1);
        $cap=$ele_t[0]=~s/C/$1/gi;
        if ($cap == $cap1)
        {
                print W1 "$line1\t\/\tAll Cys connected\n";
        }
        else
        {
                print W1 "$line1\n";
        }
}
}
if ($num_rigor == 1)
{
foreach $line2 (@hold2)
{
        $line2=~s/\n//g;
        @ele1_t=split('\t',$line2);
        $cap2=$ele1_t[0]=~s/C/$1/gi;
        if ($cap2 == $cap1)
        {
                print W1 "$line2\t\/\tAll Cys connected\n";
        }
        else
        {
                print W1 "$line2\n";
        }
}
}
###########################################################################################################################
close(W1);
#close(W2);	

print "Enter 1 for Ion trap MSMS filtering\n";
print "Enter 2 for no Ion trap MSMS filtering\n";
$ion_trap=<STDIN>;
$name_tmp=~s/\n//g;

	if ($ion_trap == 1)
	{
		readpipe("perl ion_trap_filter.pl $name_tmp > results.out");
	}
	else
	{
		readpipe("mv rev_out results.out");
	}


	if ($num_rigor==2)
	{
		readpipe("mkdir -p ../Results_MSn/Result_MS2_rigorous_$namet");
		readpipe("mv results.out ../Results_MSn/Result_MS2_rigorous_$namet");
		readpipe("mv inp_MSN_out_ms2 ../Results_MSn/Result_MS2_rigorous_$namet/inp_Result_MS2_rigorous_$namet");
		readpipe("rm -f tot_mass_final peak_mass* out_comb sum_out cys_num cys_mass* match_mass cys_count_finale remap*.out sing_seq post_fil_sum_out final_list*.out out* map_mass tmp1 mass_range range.out ip");
	}
        if ($num_rigor==1)
        {
                readpipe("mkdir -p ../Results_MSn/Result_MS2_smart_$namet");
                readpipe("mv results.out ../Results_MSn/Result_MS2_smart_$namet");
                readpipe("mv inp_MSN_strict_out_ms2 ../Results_MSn/Result_MS2_smart_$namet/inp_Result_MS2_smart_$namet");
                readpipe("rm -f tot_mass_final peak_mass* out_comb sum_out cys_num cys_mass* match_mass cys_count_finale remap*.out sing_seq post_fil_sum_out final_list*.out out* map_mass tmp1 mass_range range.out  ip");
        }
}
close(W);

#############CALCULATION ON MSN FRAGMENT IONS############

if ($num1 == 2)
{	print "Enter the complete sequence: ";
  	$name1_tmp=<STDIN>;
        my $counter=1;
        my $count_c=1;
        $seq=$name1_tmp;
        $seq=~s/\n//g;
        @ele=split('',$seq);
        for($i=0;$i<scalar @ele;$i++)
        {
                if ($ele[$i]=~m/C/i)
                {
                        $new=$ele[$i].$count_c;
                        push @array,$new;
                        $count_c++;
                }
                else
                {
                        $new=$ele[$i].$counter;
                        push (@array,$new);
                        $counter++;
                }
        }
	$len_tmp_array = @array;
        for($lta = 0; $lta <$len_tmp_array;$lta++)
        {
                $array[$lta]=~s/\n//g;
                if ($array[$lta] =~m/X/)
                {
                        $array[$lta-1] = lc($array[$lta-1]);
                }
        }
        $tmp = lc(pop(@array));
        push (@array,$tmp);
        $name1=join('',@array);
  	$name_p=$name1;
  	$num_cterm=$name_p=~s/[a-z]/$1/g;
	#print $num_cterm,"\n";
 	print "Enter the fragment sequence: ";
  	$name2=<STDIN>;
	$namet=$name2;
        $namet=~s/\d//g;
	$namet=~s/\n//g;
	
  	print "Enter the stage of MS/MS: ";
  	$stage=<STDIN>;
  	readpipe("cp ../peak_mass_msn peak_mass");
	readpipe("perl rescale_peak.pl");
        readpipe("perl master_comb.pl $name2 ");
        readpipe("perl replace_list.pl -p noX_comb_with_cys -m pep_frag $name2 ");
        readpipe("perl duplicate_rem1.pl");
  	readpipe("perl sum_out.pl");
	if ($num_rigor == 1)
        {
                readpipe("perl comb_sum_fast.pl $name1");
        }
        else
        {
                readpipe("perl comb_sum.pl $name1");
        }
	print "\nEnter 1 for normal peptide fragments\n";
  	print "Enter 2 for peptide fragments with CO loss\n";
  	print "Enter 3 for normal peptide with H2O loss\n";
  	print "Enter 4 for normal peptide with NH3 loss\n";
  	print "Enter 5 for normal peptide with NH3+CO loss\n";
  	$user_in=<STDIN>;
  	open(Read,"peak_mass");
  	@hold=<Read>;
  	$len=@hold;
  	open(W,">peak_mass1");
		if ($user_in == 1)
		{
   			foreach $line (@hold)
			{
  				$line=~s/\n//g;
   				$line1=$line+0;
   				print W $line1,"\n";
			}
		}
		if ($user_in == 2)
		{
   			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+28;
   				print W $line1,"\n";
			}
		}
		if ($user_in == 3)
		{
     			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+18;
   				print W $line1,"\n";
			}
		}	
		if ($user_in == 4)
		{
       			foreach $line (@hold)
			{
   				$line=~s/\n//g;
   				$line1=$line+17;
   				print W $line1,"\n";
			}
		}
		if ($user_in == 5)
		{
       			foreach $line (@hold)
			{
  	 			$line=~s/\n//g;
   				$line1=$line+45;
   				print W $line1,"\n";
			}
		}
	readpipe("cp peak_mass1 peak_mass");
	readpipe("perl match.pl");
	readpipe("perl segment_check.pl $name1");
	for($kgm=1;$kgm<=$num_cterm;$kgm++)
	{
  		print "Enter mass to be added for C-term $kgm status\n";
	}
	readpipe("perl post_filter.pl $name1");
	readpipe("perl duplicate_rem3.pl");
	readpipe("perl cys_complement.pl $name1");
	readpipe("perl duplicate_rem4.pl");
	readpipe("cp final_list2.out out_msn");
	readpipe("perl post_fil_disul_frag_filter.pl out_msn > strict_out_msn");
        readpipe("cp strict_out_msn inp_MSN_strict_out_msn");
        readpipe("cp out_msn inp_MSN_out_msn");
	system("perl -p -i -e 's/ \\| //g' inp_MSN_strict_out_msn");
	system("perl -p -i -e 's/ \\| //g' inp_MSN_out_msn");
	readpipe("perl final_out_digit_remove.pl n");
#	readpipe("perl final_rev1.pl $name1 n");

##########################Final Filter##############################################################
$inp_seq=$name1_tmp;
#print $name1_tmp,"\n";
$cap1=$inp_seq=~s/C/$1/gi;
$stage_t="n";

open(W1,">rev_out");
open(Read1,"out_ms$stage_t");
@hold1=<Read1>;
close(Read1);
open(Read2,"strict_out_ms$stage_t");
@hold2=<Read2>;
close(Read2);

if ($num_rigor == 2)
{
foreach $line1 (@hold1)
{
        $line1=~s/\n//g;
        @ele_t=split('\t',$line1);
        $cap=$ele_t[0]=~s/C/$1/gi;
        if ($cap == $cap1)
        {
                print W1 "$line1\t\/\tAll Cys connected\n";
        }
        else
        {
                print W1 "$line1\n";
        }
}
}
if ($num_rigor == 1)
{
foreach $line2 (@hold2)
{
	#print $line2,"\n";
        $line2=~s/\n//g;
        @ele1_t=split('\t',$line2);
        $cap2=$ele1_t[0]=~s/C/$1/gi;
        if ($cap2 == $cap1)
        {
                print W1 "$line2\t\/\tAll Cys connected\n";
        }
        else
        {
                print W1 "$line2\n";
        }
}
}
#####################################################################################################
	close(W1);
#	close(W2);

print "Enter 1 for Ion trap MSMS filtering\n";
print "Enter 2 for no Ion trap MSMS filtering\n";
$ion_trap=<STDIN>;
$name_tmp=~s/\n//g;

print "Enter the MS3 parent ion sequence: ";
$name_ion=<STDIN>;
$name_ion=~s/ \| /X/g;
$name_ion=~s/\n//g;

        if ($ion_trap == 1)
        {
                readpipe("perl ion_trap_filter_num.pl $name_ion > results.out");
        }
        else
        {
                readpipe("cp rev_out results.out");
        }


        if ($num_rigor==2)
        {
                readpipe("mkdir -p ../Results_MSn/Result_MSN_rigorous_$namet");
                readpipe("mv results.out ../Results_MSn/Result_MSN_rigorous_$namet");
		readpipe("mv inp_MSN_out_msn ../Results_MSn/Result_MSN_rigorous_$namet/inp_Result_MSN_rigorous_$namet");
                readpipe("rm -f tot_mass_final peak_mass* out_comb sum_out cys_num cys_mass* match_mass cys_count_finale remap*.out sing_seq post_fil_sum_out final_list*.out out* map_mass tmp1 mass_range range.out  ip");
        }
        if ($num_rigor==1)
        {
                readpipe("mkdir -p ../Results_MSn/Result_MSN_smart_$namet");
                readpipe("mv results.out ../Results_MSn/Result_MSN_smart_$namet");
		readpipe("mv inp_MSN_strict_out_msn ../Results_MSn/Result_MSN_smart_$namet/inp_Result_MSN_smart_$namet");
#                readpipe("rm -f tot_mass_final peak_mass* out_comb sum_out cys_num cys_mass* match_mass cys_count_finale remap*.out sing_seq post_fil_sum_out final_list*.out out* map_mass tmp1 mass_range range.out  ip");
        }
}

close(W);
