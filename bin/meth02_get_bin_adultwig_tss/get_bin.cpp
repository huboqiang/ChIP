#include "get_bin.h"

using namespace std;

int  load_posdep( 
	boost::unordered_map< std::string,std::vector<double> > &M_Gen_bp_depth, 	
	boost::unordered_map< std::string,std::vector<double> > &M_Gen_bp_methy,
	std::string &bdg_file, 
	int len_TSS_upstream, 
	int len_TSS_downstream,
	std::vector<double> &V_tot_bp_depth_s,  
	std::vector<double> &V_tot_bp_methy_s, 
	std::vector<double> &V_tot_bp_depth_a, 
	std::vector<double> &V_tot_bp_methy_a,
	int bin_size){
   
   boost::unordered_map< std::string,std::vector<double> >::iterator I_it;

   ifstream infile;
   infile.open(bdg_file.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << bdg_file << endl;
      exit(0);
   }
   int total = 0;
   
	int bp_Total  = (len_TSS_upstream+len_TSS_downstream)/bin_size;
   
	vector<double> len_zeros(bp_Total,0.0);
	
   string lineStr;
   while (getline(infile,lineStr,'\n')){
      if (lineStr[0] == ' ' || lineStr[0] == '\n'){
         continue;
      }
      vector<string> LineVec;
      boost::split(LineVec,lineStr,boost::is_any_of("\t\n"),boost::token_compress_on);
      string tag = LineVec[11] + "\t" + LineVec[9] + "\t" + LineVec[12];
		if (LineVec[8][0] == LineVec[3][0]){
			tag += "\tsense\t"     + LineVec[10];
		}
		else{
			tag += "\tantisense\t" + LineVec[10];
		}
      int    t_beg = boost::lexical_cast<int>(LineVec[6]);  // t for gene tss - tes
      int    t_end = boost::lexical_cast<int>(LineVec[7]);
      int    p_beg = boost::lexical_cast<int>(LineVec[1]);   // p for peak
      int    p_end = boost::lexical_cast<int>(LineVec[2]);
		double methy = boost::lexical_cast<double>(LineVec[4]); 
      double depth = 1.0;
		
		I_it = M_Gen_bp_depth.find(tag);
      if ( I_it==M_Gen_bp_depth.end()){
			M_Gen_bp_methy.insert( boost::unordered_map< std::string,std::vector<double> >::value_type(tag,len_zeros));
         M_Gen_bp_depth.insert( boost::unordered_map< std::string,std::vector<double> >::value_type(tag,len_zeros));
      }
		
      for (int i=p_beg;i<p_end+1;++i){

			/* bp level */

			if (LineVec[8][0] == '+'){
				int delta_tss = i - t_beg ;
				
           	int pos_num = ( i - t_beg )/bin_size;
				if (pos_num >= bp_Total){
					pos_num = bp_Total - 1;
				}
				if (LineVec[8][0] == LineVec[3][0]){
					V_tot_bp_methy_s[pos_num]      += methy;
					V_tot_bp_depth_s[pos_num]      += depth;
				}
				else{
					V_tot_bp_methy_a[pos_num]      += methy;
					V_tot_bp_depth_a[pos_num]      += depth;
				}
				
				M_Gen_bp_methy[tag][pos_num] += methy;
				M_Gen_bp_depth[tag][pos_num] += depth;
			}
			if (LineVec[8][0] == '-'){

				int delta_tss = t_end - i;

            int pos_num = ( t_end - i )/bin_size;
				if (pos_num >= bp_Total){
					pos_num = bp_Total - 1;
				}
				if (LineVec[8][0] == LineVec[3][0]){
					V_tot_bp_methy_s[pos_num]      += methy;
					V_tot_bp_depth_s[pos_num]      += depth;
				}
				else{
					V_tot_bp_methy_a[pos_num]      += methy;
					V_tot_bp_depth_a[pos_num]      += depth;
				}
				M_Gen_bp_methy[tag][pos_num] += methy;
				M_Gen_bp_depth[tag][pos_num] += depth;
			}
      }
   }
	infile.close();
   return total;
}
