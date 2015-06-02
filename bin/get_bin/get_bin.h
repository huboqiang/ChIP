#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

int  load_posdep( boost::unordered_map< std::string,std::vector<double> > &M_Gen_bin_depth, boost::unordered_map< std::string,std::vector<double> > &M_Gen_bp_depth, std::string &bdg_file, int len_upstream, int len_TSS_downstream,int len_downstream, std::vector<double> &V_tot_bin_depth, std::vector<double> &V_tot_bp_depth, int body_bin,int bin_size);
