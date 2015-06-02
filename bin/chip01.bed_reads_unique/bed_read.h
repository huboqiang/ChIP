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

void load_chrlen(std::string fa_fai,std::map< std::string,unsigned > &ChrLen);
void load_bins(std::map< std::string,unsigned > &ChrLen,std::map< std::string,std::vector<int> > &ChrBin_site,int bin);
void bed_bin(std::string in_bed,std::map< std::string,std::vector<int> > &ChrBin_site,int bin);
void report_bin(std::map< std::string,std::vector<int> > &ChrBin_site,std::string out_bed,int bin);
