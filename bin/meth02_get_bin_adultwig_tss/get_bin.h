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
	int bin_size);