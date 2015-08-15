#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "bed_read.h"
#include "gzstream.h"

using namespace std;

int bin = 100;

void usage()
{
   cout << "./bed_read -b [bins] <fa.fai> <out.bed>" << endl;
	cout << "  -b <int>  bin-size for a density point, default = " << bin << endl;
	cout << "  -h        get help information"   << endl;
   exit (0);
}

int main(int argc, char *argv[])
{
	int c;
	while ( (c=getopt(argc,argv,"b:h")) != -1 )
	{
        switch(c)
        {
            case 'b' : bin = atoi(optarg);break;
            case 'h' : usage();break;
            default : usage();
        }
	}
	if (argc < 5) usage();

	string fa_fai  = argv[optind++];
	string out_bed = argv[optind++];
	
	map< string,uint64_t >      ChrLen;
	
	map< string,vector<int>   > ChrBin_site ;
   
	//	Question1: How many chromosomes to consider and how long is it for each chromosome?
	load_chrlen(fa_fai,ChrLen);
	/*
		data structure for ChrLen:
		ChrLen[lineVec[0]] = length;
		Chromosomes:  
			keys   for map< string,uint64_t >ChrLen
		Length for each Chr:
			values for map< string,uint64_t >ChrLen
	*/	
	
	
	// Question2: Given a chr with length ChrLen[chr], we scan the chr with bin-size of bin, how we initiate two linear-table structure(Vector) for this chromosome for CpG sites covered.

	//Initate a hash table point to a vector with length of (length/bin) . Key for hash table is chrom.
	load_bins(ChrLen,ChrBin_site,bin);
	/*
		chr->bins, CpG sites covered:   ChrBin_site  e.g. chr11->[  0,   0,   0 ,  0...] , length = int(ChrLen["chr11"]/bin )
	*/	
	// Read SingleC file to write values for ChrBin_site
	
	cout << "Step2" << endl;
	report_bin(ChrBin_site,ChrLen,out_bed,bin);
}