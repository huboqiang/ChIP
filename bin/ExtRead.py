from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
from optparse   import OptionParser
import pysam
def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s -l 300 -p 2 test.dup.bam test
   
   """ % (sys.argv[0],sys.argv[0])

   description = " Select reads while converting bam to bed "

   optparser = OptionParser(version="%s v0.2 20150124" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-l", "--ext_len" , default=300, help="\nLength for extending in 3 if using Single End data [default: %default]")
   optparser.add_option("-p", "--pair_cnt", default=2  , help="\nPair end or Single end data [default: %default]")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser

def cigar_len(cigar):
   l0 = 0
   for pairs in cigar:
      ltype, leng = pairs
      if ltype == 0 or ltype == 3:
         l0 += leng
   return l0

def get_str(flag):
   strand = "+"
   bin_code = bin(flag)[2:].zfill(11)
   if int( bin_code[-5] ) != 1:
      strand = "-"
   return strand

def main():
   prepare_optparser()
   ( options,args ) = prepare_optparser().parse_args()
   try:
      bam_file = args[0]
      out_prefix= args[1]
      ext_len = int( options.ext_len  )
      pair_cnt= int( options.pair_cnt )
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
   
   try:
      f_bam = pysam.Samfile( bam_file,'rb' )
   except:
      sys.exit('Please make sure %s is a BAM file!' % bam_file )
   
   out_unique = "%s.unique.bed" % (out_prefix)
   out_multi  = "%s.multi.bed"  % (out_prefix)
   
   f_out_unique = open( out_unique,"w" )
   f_out_multi  = open( out_multi ,"w" )
   for read in f_bam:
      chrom = f_bam.header['SQ'][read.rname]['SN']
      if chrom == "*":
         continue
      begin = read.pos + 1
      cigar = read.cigar
      end   = read.pos + 1 + cigar_len(cigar)
      qname = read.qname
      mapq  = read.mapq
      strand= get_str( read.flag )
      
      if pair_cnt == 1:
         if strand == "+":
            end = begin + ext_len
         else:
            beg = end - ext_len
      
      out = "%s\t%d\t%d\t%s\t%d\t%s" % (chrom,begin,end, qname,mapq,strand)
      if mapq >=5:
         print >>f_out_unique, out
      else:
         print >>f_out_multi , out
      
   f_out_unique.close()
   f_out_multi.close()      
      
if __name__ == '__main__':
   main()