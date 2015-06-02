from __future__ import division
import re,sys,os,gzip
import subprocess
import time
import numpy   as np
import scipy.stats
from optparse   import OptionParser

def put_enhancer(enhancer,line,samp_cnt):
   f = line.split()
   
   np_peak_cnt = np.array( f[           8: 8+samp_cnt  ], dtype=int   )
   np_reg_dens = np.array( f[ 12+samp_cnt:12+samp_cnt*2], dtype=float )
      
   np_sum_dens = 0
   
   for i,val in enumerate( np_peak_cnt ):
      if val > 0:
         np_sum_dens += np_reg_dens[ i ]
      
   enhancer['peak'] += np_peak_cnt

   if enhancer['chr'] == "":
      enhancer['chr']     =      f[ 9+samp_cnt]
      enhancer['beg']     = int( f[10+samp_cnt] )
      enhancer['end']     = int( f[11+samp_cnt] )
      enhancer['sum_den'] = float(np_sum_dens)

   elif np_sum_dens >  enhancer['sum_den']:
      enhancer['chr']     =      f[ 9+samp_cnt]
      enhancer['beg']     = int( f[10+samp_cnt] )
      enhancer['end']     = int( f[11+samp_cnt] )
      enhancer['sum_den'] = float(np_sum_dens)
   
   elif np_sum_dens == enhancer['sum_den']:
      enhancer['end']     = int( f[11+samp_cnt] )
         
      
def report_enhancer(enhancer,ext_len,samp_cnt):
   chrom  = enhancer['chr']
   center = int( (enhancer['beg']+enhancer['end'])/2 )
   
   beg = center - ext_len
   end = center + ext_len
   
   np_peak = enhancer['peak']
   np_peak[ np_peak > 0 ] = 1
   l_peak    = list( np.array( np_peak,dtype='string' ) )
   peak_info = "\t".join(l_peak)
   
   out = "%s\t%d\t%d\t%s" % (chrom,  beg,  end,  peak_info )
   
   return out


def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   python %s -e 1000 /datd/huboqiang/ChIP_human/Week12/analysis/Merge_RPKM_overlap.500bp_merge.multiinter.100bp_density.bed

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
   enhancer['sum_den'] = 0.0
   
   samp_cnt = 0
   enhancer['peak']    = np.zeros( 7,dtype=int )
   
   
   pre_inf = ""

   for line in f_infile:
      line = line.strip('\n')
      f    = line.split()
      
      if samp_cnt == 0:
         samp_cnt = int( (len(f)-13)/2 )
         enhancer['peak']    = np.zeros( samp_cnt,dtype=int )
         
      pos_inf = '_'.join(f[0:3])
      if pos_inf == pre_inf:
         put_enhancer(enhancer,line,samp_cnt)
      else:
         if pre_inf != "":
            out = report_enhancer(enhancer,ext_len,samp_cnt)
            print  out
            enhancer['chr'] = ""
            enhancer['beg'] = -1
            enhancer['end'] = -1
            enhancer['sum_den'] = 0.0
            enhancer['peak']    = np.zeros( samp_cnt,dtype=int )
            
         put_enhancer(enhancer,line,samp_cnt)
      pre_inf = pos_inf
   
   out = report_enhancer(enhancer,ext_len,samp_cnt)
   print  out
   
   f_infile.close()   

if __name__ == '__main__':
   main()
