#include "get_bin.h"

using namespace std;

int  load_posdep( 
	boost::unordered_map< std::string,std::vector<double> > &M_Gen_bin_depth, 	
	boost::unordered_map< std::string,std::vector<double> > &M_Gen_bin_methy, 
	std::string &bdg_file, 
	int len_upstream, 
	int len_downstream,
	std::vector<double> &V_tot_bin_depth_s, 
	std::vector<double> &V_tot_bin_methy_s, 
	std::vector<double> &V_tot_bin_depth_a, 
	std::vector<double> &V_tot_bin_methy_a, 
	int body_bin,
	int bin_size){

   boost::unordered_map< std::string,std::vector<double> >::iterator I_it;

   ifstream infile;
   infile.open(bdg_file.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << bdg_file << endl;
      exit(0);
   }
   int total = 0;
   
   int bin_Total = (len_upstream+len_downstream)/bin_size + body_bin;
   
	vector<double> len_zero(bin_Total,0.0);
	
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
      
		I_it = M_Gen_bin_depth.find(tag);
      if ( I_it==M_Gen_bin_depth.end()){
			M_Gen_bin_methy.insert( boost::unordered_map< std::string,std::vector<double> >::value_type(tag,len_zero));
         M_Gen_bin_depth.insert( boost::unordered_map< std::string,std::vector<double> >::value_type(tag,len_zero));
      }

      for (int i=p_beg;i<p_end+1;++i){


			/* gene-bin level */			
			           
         if (LineVec[8][0] == '+'){
			
				int delta_tss = i - (t_beg + len_upstream);
				int delta_tes = i - (t_end - len_downstream);
				
            if (delta_tss < -1*len_upstream){
               continue;
            }
            if (delta_tes > len_downstream){
               break;
            }
            if (delta_tss < 0){
               int pos_num  = (len_upstream + delta_tss)/bin_size;
					if (pos_num >= bin_Total){
						pos_num = bin_Total - 1;
					}
               total += depth;
					if (LineVec[8][0] == LineVec[3][0]){
						V_tot_bin_methy_s[pos_num]      += methy;
						V_tot_bin_depth_s[pos_num]      += depth;
					}
					else{
						V_tot_bin_methy_a[pos_num]      += methy;
						V_tot_bin_depth_a[pos_num]      += depth;
					}
               M_Gen_bin_methy[tag][pos_num] += methy;
               M_Gen_bin_depth[tag][pos_num] += depth;
            }
            else{
               if (delta_tes > 0){
                  int pos_num  = len_upstream/bin_size + body_bin + delta_tes/bin_size;
						if (pos_num >= bin_Total){
							pos_num = bin_Total - 1;
						}
						total += depth;
						if (LineVec[8][0] == LineVec[3][0]){
							V_tot_bin_methy_s[pos_num]      += methy;
							V_tot_bin_depth_s[pos_num]      += depth;
						}
						else{
							V_tot_bin_methy_a[pos_num]      += methy;
							V_tot_bin_depth_a[pos_num]      += depth;
						}
                  M_Gen_bin_methy[tag][pos_num] += methy;
                  M_Gen_bin_depth[tag][pos_num] += depth;
               }
               else{
                  int tran_len = (t_end - t_beg) - (len_upstream + len_downstream);
                  int tran_bin_size = tran_len / body_bin;
                  
                  if (tran_bin_size == 0 ){
                     continue;
                  }

                  int pos_num  = len_upstream/bin_size + delta_tss/tran_bin_size;
						if (pos_num >= bin_Total){
							pos_num = bin_Total - 1;
						}
                  total += depth;
						if (LineVec[8][0] == LineVec[3][0]){
							V_tot_bin_methy_s[pos_num]      += methy;
							V_tot_bin_depth_s[pos_num]      += depth;
						}
						else{
							V_tot_bin_methy_a[pos_num]      += methy;
							V_tot_bin_depth_a[pos_num]      += depth;
						}
                  M_Gen_bin_methy[tag][pos_num] += methy;
                  M_Gen_bin_depth[tag][pos_num] += depth;
					}
            }
         }
         if (LineVec[8][0] == '-'){

				int delta_tss = i - (t_beg + len_upstream);
				int delta_tes = i - (t_end - len_downstream);
            
				if (delta_tss < -1*len_upstream){
               continue;
            }
            if (delta_tes > len_downstream){
               break;
            }
            if (delta_tss < 0){
               int pos_num  =  -1*delta_tss/bin_size + body_bin + len_downstream/bin_size;
					if (pos_num >= bin_Total){
						pos_num = bin_Total - 1;
					}
               total += depth;
					if (LineVec[8][0] == LineVec[3][0]){
						V_tot_bin_methy_s[pos_num]      += methy;
						V_tot_bin_depth_s[pos_num]      += depth;
					}
					else{
						V_tot_bin_methy_a[pos_num]      += methy;
						V_tot_bin_depth_a[pos_num]      += depth;
					}
               M_Gen_bin_methy[tag][pos_num] += methy;
               M_Gen_bin_depth[tag][pos_num] += depth;
            }
            else{
               if (delta_tes > 0){
                  int pos_num  = (len_downstream-delta_tes)/bin_size;
						if (pos_num >= bin_Total){
							pos_num = bin_Total - 1;	
						}
	               total += depth;
						if (LineVec[8][0] == LineVec[3][0]){
							V_tot_bin_methy_s[pos_num]      += methy;
							V_tot_bin_depth_s[pos_num]      += depth;
						}
						else{
							V_tot_bin_methy_a[pos_num]      += methy;
							V_tot_bin_depth_a[pos_num]      += depth;
						}
	               M_Gen_bin_methy[tag][pos_num] += methy;
	               M_Gen_bin_depth[tag][pos_num] += depth;
					}
               else{
                  int tran_len = (t_end - t_beg) - (len_upstream + len_downstream);
                  int tran_bin_size = tran_len / body_bin;
                  
                  if (tran_bin_size == 0 ){
                     continue;
                  }

                  int pos_num  = -1*delta_tes/tran_bin_size + len_downstream/bin_size;
						if (pos_num >= bin_Total){
							pos_num = bin_Total - 1;
						}
	               total += depth;
						if (LineVec[8][0] == LineVec[3][0]){
							V_tot_bin_methy_s[pos_num]      += methy;
							V_tot_bin_depth_s[pos_num]      += depth;
						}
						else{
							V_tot_bin_methy_a[pos_num]      += methy;
							V_tot_bin_depth_a[pos_num]      += depth;
						}
	               M_Gen_bin_methy[tag][pos_num] += methy;
	               M_Gen_bin_depth[tag][pos_num] += depth;
					}
            }
         }
      }
   }
	infile.close();
   return total;
}
