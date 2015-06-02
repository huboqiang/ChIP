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

int len_TSS_upstream   = 1000;
int len_TSS_downstream = 500;
int bin_size = 1;

void usage()
{
   cout << "./get_bin <sample_peak_hg19_ref.bdg> [outfile_TSS_up_down_Bp-region] " << endl;
   cout << "  -U <int>  length extends from upstream,       default = " << len_TSS_upstream   << endl;
   cout << "  -D <int>  length extends from downstream,     default = " << len_TSS_downstream << endl;
   cout << "  -b <int>  bin-size for a density point,       default = " << bin_size << endl;
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
         case 'U' : len_TSS_upstream  	= atoi(optarg);break;
         case 'D' : len_TSS_downstream	= atoi(optarg);break;
         case 'b' : bin_size           = atoi(optarg);break;
         default : usage();
      }
   }
   if (argc < 2) usage();
   string f_bdg_file     = argv[optind++];
   string f_out_bp_file  = argv[optind++];
   
   int bp_Total  = (len_TSS_upstream+len_TSS_downstream) / bin_size ;
   vector<double> V_tot_bp_depth_s(bp_Total,0);
   vector<double> V_tot_bp_depth_a(bp_Total,0);
   vector<double> V_tot_bp_methy_s(bp_Total,0);
   vector<double> V_tot_bp_methy_a(bp_Total,0);
   boost::unordered_map< string,vector<double> > M_Gen_bp_depth;
   boost::unordered_map< string,vector<double> > M_Gen_bp_methy;

   int total = load_posdep(M_Gen_bp_depth,M_Gen_bp_methy,f_bdg_file,len_TSS_upstream,len_TSS_downstream,V_tot_bp_depth_s,V_tot_bp_methy_s,V_tot_bp_depth_a,V_tot_bp_methy_a,bin_size);

   ofstream outFile1;

   outFile1.open(f_out_bp_file.c_str());

   outFile1 << total << endl;
   outFile1 << "gene\ttranscritps\tstrand\tpromoter_type\ttotal(sense)";
   for (int i=0;i<bp_Total;i++){
		double meth_level = V_tot_bp_methy_s[i] / boost::lexical_cast<double>(V_tot_bp_depth_s[i]);
      outFile1 << "\t" << meth_level;
   }
   outFile1 << endl;
   outFile1 << "gene\ttranscritps\tstrand\tpromoter_type\ttotal(anti-sense)";
   for (int i=0;i<bp_Total;i++){
		double meth_level = V_tot_bp_methy_a[i] / boost::lexical_cast<double>(V_tot_bp_depth_a[i]);
      outFile1 << "\t" << meth_level;
   }
   outFile1 << endl;
	
	
	boost::unordered_map< string,vector<double> >::iterator I_it;
   for (I_it=M_Gen_bp_depth.begin();  I_it!=M_Gen_bp_depth.end(); I_it++ ){
      outFile1 << I_it->first ;
      for (int i=0;i<bp_Total;i++){
			outFile1 << "\t" <<  M_Gen_bp_methy[I_it->first][i] << "," <<  M_Gen_bp_depth[I_it->first][i];
      }
      outFile1 << endl;
   }
   outFile1.close();
}
