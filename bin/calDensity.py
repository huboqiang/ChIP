from __future__ import division
import re,sys,os,gzip
import subprocess
import time
import numpy   as np
import scipy.stats
from optparse   import OptionParser

def put_enhancer(enhancer,line,samp_cnt):
   f = line.split()
   
   np_peak_cnt = np.array( f[          3: 3+samp_cnt  ], dtype=int   )
   np_reg_dens = np.array( f[ 6+samp_cnt: 6+samp_cnt*2], dtype=float )
            

   enhancer['dens'] += np_reg_dens
   enhancer['cnt']    += 1
   if enhancer['chr'] == "":
      enhancer['chr']     =      f[ 0 ]
      enhancer['beg']     = int( f[ 1 ] )
      enhancer['end']     = int( f[ 2 ] )
      
         
      
def report_enhancer(enhancer,ext_len,samp_cnt):
   chrom  = enhancer['chr']
   center = int( (enhancer['beg']+enhancer['end'])/2 )
   
   beg = center - ext_len
   end = center + ext_len
   
   
   np_dens = enhancer['dens'] / enhancer['cnt']
   l_dens    = list( np.array( np_dens,dtype='string' ) )
   dens_info = "\t".join(l_dens)
   
   out = "%s\t%d\t%d\t%s" % (chrom,  beg,  end,  dens_info )
   
   return out


def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   bedtools intersect -wo -sorted -a /datd/huboqiang/ChIP_human/Week12/analysis/Merge_2kb_Peaks.bed -b Merge_RPKM_bin100.mrg.rev.bed | python %s /dev/stdin

   """ % (sys.argv[0],sys.argv[0])

   description = " Select core-region in a given enhancer. "

   optparser = OptionParser(version="%s v0.2 20150124" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-e", "--ext_len" , default=1000  , help="\nPair end or Single end data [default: %default]")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   


def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      in_file =     args[0]
      ext_len = int( options.ext_len )
   
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
   
   f_infile    = open(in_file    ,"r")
   
   enhancer = {}
   enhancer['chr'] = ""
   enhancer['beg'] = -1
   enhancer['end'] = -1
   enhancer['dens'] = np.zeros(1)
   enhancer['cnt']  = 0
   samp_cnt = 0
   pre_inf = ""

   for line in f_infile:
      line = line.strip('\n')
      f    = line.split()
      
      if samp_cnt == 0:
         samp_cnt = int( (len(f)-7)/2 )
         enhancer['dens'] = np.zeros(samp_cnt)
      
      pos_inf = '_'.join(f[0:3])
      if pos_inf == pre_inf:
         put_enhancer(enhancer,line,samp_cnt)
      else:
         if pre_inf != "":
            out = report_enhancer(enhancer,ext_len,samp_cnt)
            print  out
            enhancer['chr']   = ""
            enhancer['beg']   = -1
            enhancer['end']   = -1
            enhancer['dens']  = np.zeros( samp_cnt )
            enhancer['cnt']   = 0
            
         put_enhancer(enhancer,line,samp_cnt)
      pre_inf = pos_inf
   
   out = report_enhancer(enhancer,ext_len,samp_cnt)
   print  out
   
   f_infile.close()   

if __name__ == '__main__':
   main()
