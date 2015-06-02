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
#include "get_bin.h"

using namespace std;

int len_upstream = 5000;
int len_downstream = 5000;
int len_TSS_downstream = 5000;
int bin_size = 50;
int body_bin = 100;

void usage()
{
   cout << "./get_bin <sample_peak_hg19_ref.bdg> [outfile_TSS_up_down_Bp-region] [outfile_geneBody_Bin_up_down_Bp-region]" << endl;
   cout << "  -U <int>  length extends from upstream,       default = " << len_upstream << endl;
   cout << "  -T <int>  length extends from TSS downstream, default = " << len_upstream << endl;
   cout << "  -D <int>  length extends from downstream,     default = " << len_downstream << endl;
   cout << "  -b <int>  bin-size for a density point,       default = " << bin_size << endl;
	cout << "  -B <int>  divided bins number in gene-body,   default = " << body_bin << endl;
   cout << "  -h        get help information"   << endl;
   cout << " U \% b or D \% b should better to be zero to avoid mistake edge-caculation in TSS or TES" << endl;
   exit (0);
}

int main(int argc, char *argv[])
{
   int c;
   while ( (c=getopt(argc,argv,"U:T:D:b:B:h")) != -1 )
   {
      switch(c)
      {
         case 'h' : usage();break;
         case 'U' : len_upstream       = atoi(optarg);break;
         case 'T' : len_TSS_downstream = atoi(optarg);break;
         case 'D' : len_downstream     = atoi(optarg);break;
         case 'b' : bin_size           = atoi(optarg);break;
			case 'B' : body_bin           = atoi(optarg);break;
         default : usage();
      }
   }
   if (argc < 2) usage();
   string f_bdg_file     = argv[optind++];
   string f_out_bp_file  = argv[optind++];
   string f_out_bin_file = argv[optind++];

   int bin_Total = (len_upstream+len_downstream ) / bin_size + body_bin;
   vector<double> V_tot_bin_depth(bin_Total,0);
   boost::unordered_map< string,vector<double> > M_Gen_bin_depth;
   
   int bp_Total  = (len_upstream+len_TSS_downstream) / bin_size ;
   vector<double> V_tot_bp_depth(bp_Total,0);
   boost::unordered_map< string,vector<double> > M_Gen_bp_depth;

   int total = load_posdep(M_Gen_bin_depth,M_Gen_bp_depth,f_bdg_file,len_upstream,len_TSS_downstream,len_downstream,V_tot_bin_depth,V_tot_bp_depth,body_bin,bin_size);

   ofstream outFile1;
   ofstream outFile2;

   outFile1.open(f_out_bp_file.c_str());
   outFile2.open(f_out_bin_file.c_str());

   outFile1 << total << endl;
   outFile1 << "gene\ttranscritps\ttotal";
   for (int i=0;i<bin_Total;i++){
      outFile1 << "\t" << V_tot_bin_depth[i];
   }
   outFile1 << endl;
   boost::unordered_map< string,vector<double> >::iterator I_it;
   for (I_it=M_Gen_bin_depth.begin(); I_it!=M_Gen_bin_depth.end();I_it++ ){
      outFile1 << I_it->first ;
      for (int i=0;i<bin_Total;i++){
         outFile1 << "\t" << I_it->second[i];
      }
      outFile1 << endl;
   }
   
   outFile2 << total << endl;
   outFile2 << "gene\ttranscritps\ttotal";
   for (int i=0;i<bp_Total;i++){
      outFile2 << "\t" << V_tot_bp_depth[i];
   }
   outFile2 << endl;
   for (I_it=M_Gen_bp_depth.begin();  I_it!=M_Gen_bp_depth.end(); I_it++ ){
      outFile2 << I_it->first ;
      for (int i=0;i<bp_Total;i++){
         outFile2 << "\t" << I_it->second[i];
      }
      outFile2 << endl;
   }

   outFile1.close();
   outFile2.close();
}
