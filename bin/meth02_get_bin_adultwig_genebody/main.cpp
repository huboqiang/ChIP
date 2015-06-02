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
int bin_size = 50;
int body_bin = 100;

void usage()
{
   cout << "./get_bin <sample_peak_hg19_ref.bdg> [outfile_TSS_up_down_Bp-region] [outfile_geneBody_Bin_up_down_Bp-region]" << endl;
   cout << "  -U <int>  length extends from upstream,       default = " << len_upstream << endl;
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
         case 'D' : len_downstream     = atoi(optarg);break;
         case 'b' : bin_size           = atoi(optarg);break;
			case 'B' : body_bin           = atoi(optarg);break;
         default : usage();
      }
   }
   if (argc < 2) usage();
   string f_bdg_file     = argv[optind++];
   string f_out_bin_file = argv[optind++];

   int bin_Total = (len_upstream+len_downstream ) / bin_size + body_bin;
   vector<double> V_tot_bin_depth_s(bin_Total,0);
   vector<double> V_tot_bin_depth_a(bin_Total,0);
   vector<double> V_tot_bin_methy_s(bin_Total,0);
   vector<double> V_tot_bin_methy_a(bin_Total,0);
   boost::unordered_map< string,vector<double> > M_Gen_bin_depth;
   boost::unordered_map< string,vector<double> > M_Gen_bin_methy;

   int total = load_posdep(M_Gen_bin_depth,M_Gen_bin_methy,f_bdg_file,len_upstream,len_downstream,V_tot_bin_depth_s,V_tot_bin_methy_s,V_tot_bin_depth_a,V_tot_bin_methy_a,body_bin,bin_size);

   ofstream outFile1;

   outFile1.open(f_out_bin_file.c_str());
   outFile1 << total << endl;
   outFile1 << "gene\ttranscritps\tpromoter_type\tstrand\ttotal(sense)";
   for (int i=0;i<bin_Total;i++){
		double meth_level = V_tot_bin_methy_s[i] / boost::lexical_cast<double>(V_tot_bin_depth_s[i]);
		outFile1 << "\t" << meth_level;
   }
   outFile1 << endl;
   outFile1 << "gene\ttranscritps\tpromoter_type\tstrand\ttotal(anti-sense)";
   for (int i=0;i<bin_Total;i++){
		double meth_level = V_tot_bin_methy_a[i] / boost::lexical_cast<double>(V_tot_bin_depth_a[i]);
		outFile1 << "\t" << meth_level;
   }
   outFile1 << endl;
	
   boost::unordered_map< string,vector<double> >::iterator I_it;
	
   for (I_it=M_Gen_bin_depth.begin(); I_it!=M_Gen_bin_depth.end();I_it++ ){
		outFile1 << I_it->first ;
      for (int i=0;i<bin_Total;i++){
			outFile1 << "\t" <<  M_Gen_bin_methy[I_it->first][i] << "," <<  M_Gen_bin_depth[I_it->first][i];
      }
      outFile1 << endl;
   }
   outFile1.close();
}
