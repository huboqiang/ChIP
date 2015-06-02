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
	int bin_size);