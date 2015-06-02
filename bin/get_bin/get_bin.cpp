#include "get_bin.h"

using namespace std;

int  load_posdep( boost::unordered_map< std::string,std::vector<double> > &M_Gen_bin_depth, boost::unordered_map< std::string,std::vector<double> > &M_Gen_bp_depth, std::string &bdg_file, int len_upstream, int len_TSS_downstream,int len_downstream, std::vector<double> &V_tot_bin_depth, std::vector<double> &V_tot_bp_depth, int body_bin,int bin_size){
   //( M_Gen_bin_depth,M_Gen_bp_depth,f_bdg_file,len_upstream,len_TSS_downstream,len_downstream,V_tot_bin_depth,V_tot_bp_depth,body_bin,bin_size )
   boost::unordered_map< std::string,std::vector<double> >::iterator I_it;

   ifstream infile;
   infile.open(bdg_file.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << bdg_file << endl;
      exit(0);
   }
   int total = 0;
   
	int bp_Total  = (len_upstream+len_TSS_downstream)/bin_size;
   int bin_Total = (len_upstream+len_downstream)/bin_size + body_bin;
   
	vector<double> len_zero(bin_Total,0.0);
	vector<double> len_zeros(bp_Total,0.0);
	
   string lineStr;
   while (getline(infile,lineStr,'\n')){
      if (lineStr[0] == ' ' || lineStr[0] == '\n'){
         continue;
      }
      vector<string> LineVec;
      boost::split(LineVec,lineStr,boost::is_any_of("\t\n"),boost::token_compress_on);
      string tag = LineVec[16] + "\t" + LineVec[14] + "\t" + LineVec[15];
      int t_beg = boost::lexical_cast<int>(LineVec[11]);  // t for gene tss - tes
      int t_end = boost::lexical_cast<int>(LineVec[12]);
      int p_beg = boost::lexical_cast<int>(LineVec[6]);   // p for peak
      int p_end = boost::lexical_cast<int>(LineVec[7]) - 1;
      int depth = boost::lexical_cast<double>(LineVec[8]);
      
		I_it = M_Gen_bin_depth.find(tag);
      if ( I_it==M_Gen_bin_depth.end()){
         M_Gen_bin_depth.insert( boost::unordered_map< std::string,std::vector<double> >::value_type(tag,len_zero));
      }
		
		I_it = M_Gen_bp_depth.find(tag);
      if ( I_it==M_Gen_bp_depth.end()){
         M_Gen_bp_depth.insert( boost::unordered_map< std::string,std::vector<double> >::value_type(tag,len_zeros));
      }
		
      for (int i=p_beg;i<p_end+1;++i){
			int delta_tss = i - (t_beg + len_upstream);
         int delta_tes = i - (t_end - len_downstream);
			if (LineVec[13][0] == '+'){
				if (delta_tss >= -1*len_upstream and delta_tss <= len_TSS_downstream){
            	int pos_num = ( i - t_beg )/bin_size;
//					cerr << lineStr << '\t' << i << '\t' << pos_num << '\t'<< delta_tss << '\t'<< delta_tes << '\t' << bp_Total << endl;
					if (pos_num >= bp_Total){
						pos_num = bp_Total - 1;
					}
					V_tot_bp_depth[pos_num] += depth/boost::lexical_cast<double>(bin_size);
					M_Gen_bp_depth[tag][pos_num] += depth/boost::lexical_cast<double>(bin_size);
            }
			}
			if (LineVec[13][0] == '-'){
            if (delta_tes >= -1*len_TSS_downstream and delta_tes <= len_upstream){
            	int pos_num = ( t_end - i )/bin_size;
					if (pos_num >= bp_Total){
						pos_num = bp_Total - 1;
					}
//					cerr << lineStr << '\t' << i << '\t' << pos_num << '\t' << delta_tss << '\t' << delta_tes << '\t' << bp_Total << endl;
					V_tot_bp_depth[pos_num] += depth/boost::lexical_cast<double>(bin_size);
					M_Gen_bp_depth[tag][pos_num] += depth/boost::lexical_cast<double>(bin_size);
            }
			}
               
         if (LineVec[13][0] == '+'){
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
               V_tot_bin_depth[pos_num] += depth/boost::lexical_cast<double>(bin_size);
               total += depth;
               M_Gen_bin_depth[tag][pos_num] += depth/boost::lexical_cast<double>(bin_size);
//					cout << i << "\t" << pos_num << "\t" << lineStr << "\t" << endl; 
            }
            else{
               if (delta_tes > 0){
                  int pos_num  = len_upstream/bin_size + body_bin + delta_tes/bin_size;
						if (pos_num >= bin_Total){
							pos_num = bin_Total - 1;
						}
                  V_tot_bin_depth[pos_num] += depth/boost::lexical_cast<double>(bin_size);
                  total += depth;
                  M_Gen_bin_depth[tag][pos_num] += depth/boost::lexical_cast<double>(bin_size);
//						cout << i << "\t" << pos_num << "\t" << lineStr << "\t" << endl; 
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
                  V_tot_bin_depth[pos_num] += depth/boost::lexical_cast<double>(tran_bin_size);
                  total += depth;
                  M_Gen_bin_depth[tag][pos_num] += depth/boost::lexical_cast<double>(bin_size);
//						cout << i << "\t" << pos_num << "\t" << lineStr << "\t" << endl; 
					}
            }
         }
         if (LineVec[13][0] == '-'){
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
               V_tot_bin_depth[pos_num] += depth / boost::lexical_cast<double>(bin_size);
               total += depth;
               M_Gen_bin_depth[tag][pos_num] += depth / boost::lexical_cast<double>(bin_size);
//					cout << i << "\t" << pos_num << "\t" << lineStr << "\t" << endl; 
            }
            else{
               if (delta_tes > 0){
                  int pos_num  = (len_downstream-delta_tes)/bin_size;
						if (pos_num >= bin_Total){
							pos_num = bin_Total - 1;	
						}
                  V_tot_bin_depth[pos_num] += depth / boost::lexical_cast<double>(bin_size);
                  total += depth;
                  M_Gen_bin_depth[tag][pos_num] += depth / boost::lexical_cast<double>(bin_size);
//						cout << i << "\t" << pos_num << "\t" << lineStr << "\t" << endl; 
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
                  V_tot_bin_depth[pos_num] += depth / boost::lexical_cast<double>(tran_bin_size);
                  total += depth;
                  M_Gen_bin_depth[tag][pos_num] += depth / boost::lexical_cast<double>(bin_size);
//     				cout << i << "\t" << pos_num << "\t" << lineStr << "\t" << endl; 
					}
            }
         }

      }
   }
	infile.close();
   return total;
}
