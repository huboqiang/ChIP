#include "bed_read.h"
#include "gzstream.h"

using namespace std;

void load_chrlen(string fa_fai,map< string,uint64_t > &ChrLen){
   ifstream infile;
   string file_name = fa_fai;
   infile.open(file_name.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << file_name << endl;
      exit(0);
   }
   string lineStr;
   while (getline(infile,lineStr,'\n')){
      if (lineStr[0] == ' ' || lineStr[0] == '\n'){
         continue;
      }
      vector<string> lineVec;
      boost::split(lineVec,lineStr,boost::is_any_of(":, \t\n"), boost::token_compress_on);
      uint64_t length =  boost::lexical_cast<uint64_t>(lineVec[1]);
      ChrLen[lineVec[0]] = length;
   }
   infile.close();
}

void load_bins(map< string,uint64_t > &ChrLen, map< string,vector<int> > &ChrBin_site,int bin){
   map< string,uint64_t >::iterator I_it;
   for ( I_it=ChrLen.begin(); I_it!=ChrLen.end();I_it++ ){
      string chr = I_it->first;
      uint64_t    len = I_it->second;
      uint64_t    bin_len = len/bin + 1;
      vector<int>   len_zeros_site( bin_len,0);
		vector<float> len_zeros_ratio(bin_len,0.0);
      ChrBin_site.insert(  map< string,vector<int>   >::value_type(chr,len_zeros_site ) );
   }
}

void report_bin( map< string,vector<int> > &ChrBin_site, map< string,uint64_t > &ChrLen, string out_bed,int bin){
	//ogzstream outFile;
	ofstream outFile;
   outFile.open(out_bed.c_str());
   if ( ! outFile )
   {
      cerr << "fail to open input file" << outFile << endl;
      exit(0);
   }
   map< string,vector<int>  >::iterator I_it;
   for ( I_it=ChrBin_site.begin(); I_it!=ChrBin_site.end();I_it++ ){
      string chr = I_it->first;
      for (uint64_t i=0;i<ChrBin_site[chr].size();i++){
         uint64_t begin = i*bin + 1;
         uint64_t end   = (i+1)*bin; 
         if (end > ChrLen[chr]){
             end = ChrLen[chr];
         }
         outFile << chr << "\t" << begin << "\t" << end << endl;
      }
   }
   outFile.close();
}