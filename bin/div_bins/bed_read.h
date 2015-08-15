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
#include "gzstream.h"

void load_chrlen(std::string fa_fai,std::map< std::string,uint64_t > &ChrLen);
void load_bins(std::map< std::string,uint64_t > &ChrLen,std::map< std::string,std::vector<int> > &ChrBin_site,int bin);
void report_bin(std::map< std::string,std::vector<int> > &ChrBin_site, std::map< std::string,uint64_t > &ChrLen, std::string out_bed,int bin);
