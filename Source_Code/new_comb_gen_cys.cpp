/*******************************************************************************
 * Copyright (c) 2012, Moitrayee Bhattacharyya, Kallol Gupta, P. Balaram, Indian Institute of Science, Bangalore, India.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * (2) Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * (3) Neither the name of the Indian Institute of Science, Bangalore, India 
 * nor the names of its contributors may be used to endorse or promote products 
 * derived from this software without specific prior written permission.
 *  
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESSINTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY,OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

# include<iostream>
# include<cstdlib>
# include<cstdio>
# include<cstring>
# include<string>
# include<cassert>
# include<vector>
# include<unistd.h>

typedef std::vector<std::string>  ExtendedSequenceRep;

std::vector<int> global_marker;

ExtendedSequenceRep convertToExtended(const char* str, const int len){
	int ilen = strlen(str);
	ExtendedSequenceRep sequence;
	char const* ptr = str;
	assert(len > 0 && len <= ilen);
	for( int i=0; i<len; ){
		std::string tmpStore;
		do{
			tmpStore += *ptr;
			ptr++;
			i++;
		}while( i<len && isdigit(*ptr) );
		sequence.push_back(tmpStore);
	}
	return sequence;
}

std::ostream& operator << (std::ostream& os, ExtendedSequenceRep const& sequence){
	for(int i=0; i<sequence.size(); ++i) os<<sequence[i];
	return os;
}

void enumerate_combination(
							ExtendedSequenceRep const& sequence, 
							const int substr_sz, 
							const int comb_size, 
							const int n_i, 
							const int j_i )
{
	assert(substr_sz <= sequence.size() && comb_size <= substr_sz);
	if(comb_size == 0){
		for(int i=0; i<global_marker.size(); ++i)
			std::cout<<sequence[global_marker[i]];
		std::cout<<std::endl;
	}else
		for(int i=n_i, j = j_i ; i < substr_sz && comb_size <= (substr_sz - n_i) ; ++i,++j){
			global_marker.push_back(j);
			enumerate_combination(sequence,substr_sz, comb_size - 1, i+1,j+1);
			global_marker.erase( global_marker.begin() + global_marker.size() - 1);
		}
}

void enumerate_contiguous( 
						  ExtendedSequenceRep const& sequence, 
						  const int substr_sz, 
						  const int comb_size, 
						  const int strt )
{
	assert( substr_sz <= sequence.size() && comb_size <= substr_sz);
	for(int i=strt; i <= substr_sz - comb_size ; ++i ){
		for(int j=0; j<comb_size; ++j )
			std::cout<<sequence[i+j];
		std::cout<<std::endl;
	}
}

void usage(char const* prog)
{
	std::cerr<<"Usage: "<<prog<<" -s <input string> -n <substring size> [-c ]"<<std::endl;
	std::cerr<<"\tc : for switching to all combination output."<<std::endl;
	exit(1);
}

int main(int argc, char** argv){
	int c, prev_ind;
	char *str = NULL;
	char *prog = argv[0];
	int substr_size = 0;
	bool comb = false;

	while(prev_ind = optind, ( c= getopt(argc,argv,":s:n:c")) != EOF ){
		if ( optind == prev_ind + 2 && *optarg == '-' ) { c = ':'; --optind; }

		switch(c){
			case 's': 
					str = optarg; 
					break;
			case 'n': 
					substr_size = atoi(optarg); 
					break;
			case 'c':
					comb = true;
					break;
			default: 
					std::cerr<<"Error: Invalid Argument"<<optarg<<std::endl;
					usage(prog);
		}
	}
	if( str == NULL || substr_size == 0 ) usage(prog);

	ExtendedSequenceRep sequence = convertToExtended(str,strlen(str));
	if( comb )
		enumerate_combination(sequence,sequence.size(),substr_size,0,0);
	else
		enumerate_contiguous(sequence,sequence.size(),substr_size,0);
}


