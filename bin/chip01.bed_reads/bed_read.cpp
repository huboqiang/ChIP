#include "bed_read.h"
#include "gzstream.h"

using namespace std;

void load_chrlen(string fa_fai,map< string,unsigned > &ChrLen){
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
      unsigned length =  boost::lexical_cast<unsigned>(lineVec[1]);
      ChrLen[lineVec[0]] = length;
   }
   infile.close();
}

void load_bins(map< string,unsigned > &ChrLen, map< string,vector<int> > &ChrBin_site,int bin){
   map< string,unsigned >::iterator I_it;
   for ( I_it=ChrLen.begin(); I_it!=ChrLen.end();I_it++ ){
      string chr = I_it->first;
      unsigned    len = I_it->second;
      unsigned    bin_len = len/bin + 1;
      vector<int>   len_zeros_site( bin_len,0);
		vector<float> len_zeros_ratio(bin_len,0.0);
      ChrBin_site.insert(  map< string,vector<int>   >::value_type(chr,len_zeros_site ) );
   }
}

void bed_bin(string in_bed,map< string,vector<int> > &ChrBin_site,int bin){
   ifstream infile;
   infile.open(in_bed.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << in_bed << endl;
      exit(0);
   }
   string lineStr;
   while (getline(infile,lineStr,'\n')){
      if (lineStr[0] == ' ' || lineStr[0] == '\n'){
         continue;
      }
      vector<string> lineVec;
      boost::split(lineVec,lineStr,boost::is_any_of(":, \t\n"), boost::token_compress_on);
		// my @f = split(/\t/,$lineStr);
		
      string chr = lineVec[0];
      unsigned begin = boost::lexical_cast<unsigned>(lineVec[1]);
		unsigned end   = boost::lexical_cast<unsigned>(lineVec[2]);
		unsigned bin_n = (begin + end) / (2*bin);		
		if ( bin_n >= ChrBin_site[chr].size() or bin_n < 0 ){
			continue;
		}

      ChrBin_site[chr][bin_n]  += 1;
   }
   infile.close();
}

void report_bin( map< string,vector<int> > &ChrBin_site, string out_bed,int bin){
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
      for (unsigned i=0;i<ChrBin_site[chr].size();i++){
         unsigned begin = i*bin + 1;
         unsigned end   = (i+1)*bin; 
			int  sum_site  = ChrBin_site[chr][i];
//         outFile << chr << "\t" << begin << "\t" << end << "\t" << sum_site << endl;
         outFile << sum_site << endl;
      }
   }
   outFile.close();
}