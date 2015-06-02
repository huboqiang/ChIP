from __future__ import division
import re,sys,os,gzip
import cPickle as pickle
import numpy   as np

import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib as mpl

from optparse   import OptionParser

class LoadSamp(object):
   def __init__(self,sam_file_chip):
      self.infile = sam_file_chip
   
   def load_samp(self):
      self.mrg_ctrl   = {} 
      f_sam_file_chip = open(self.infile,"r")
      in_h = f_sam_file_chip.readline()
      
      for line in f_sam_file_chip:
         line = line.strip('\n')
         f    = line.split()
         sam    = f[0]
         stage  = f[1]
         ltype  = f[2]
         tissue = f[3]
         brief  = f[4]
         merge_name = f[5]
         end_type   = f[6]
         control    = f[7]
         
         if ltype == "input":
            continue
            
         if merge_name not in self.mrg_ctrl:
            self.mrg_ctrl[ merge_name ] = control
      
      f_sam_file_chip.close()

class SampCol(LoadSamp):
   def __init__(self,chip_sam_file ):
      LoadSamp.__init__(self,chip_sam_file)
      
      self.load_samp()
   
   def load_Matrix(self,mat_file):
      self.f_mat_file = gzip.open( mat_file,"rb" )
      self.__load_header()
      self.__get_ChIP_idx()
      self.__parse_line()
   
   def __load_header(self):
      h_line   = self.f_mat_file.readline()
      self.l_mrgsam = h_line.split()[3:]
      
   def __get_ChIP_idx(self):
      """ 1.For each sample, is input? """
      self.l_mrgsam_out = []
      for mrgsam in self.l_mrgsam:
         if mrgsam in self.mrg_ctrl:
            self.l_mrgsam_out.append( mrgsam )

      """ 2.Get the sub-index for sample. """
      self.l_mrgsam_out_idx = []
      for mrgsam in self.l_mrgsam_out:
         idx = self.l_mrgsam.index( mrgsam )
         self.l_mrgsam_out_idx.append( idx+3 )

      """ 3.Get the sub-index for sample's control. """
      self.l_mrgsam_ctrl_out_idx = []
      for mrgsam in self.l_mrgsam_out:
         mrgsam_ctrl = self.mrg_ctrl[ mrgsam ]
         idx = self.l_mrgsam.index( mrgsam_ctrl )
         self.l_mrgsam_ctrl_out_idx.append( idx+3 )
         
   def __parse_line(self):
      out_samp = "#chr\tbeg\tend\t%s" % ( "\t".join(self.l_mrgsam_out) )
      print out_samp
      
      for line in self.f_mat_file:
         line = line.strip('\n')
         f    = line.split()
         pos  = "\t".join( f[0:3] )
         l_val = [  str( max(float(f[self.l_mrgsam_out_idx[ i ]]) - float(f[self.l_mrgsam_ctrl_out_idx[ i ]]),0) ) for i in xrange(0,len(self.l_mrgsam_out))  ]
         print "%s\t%s" % ( pos,"\t".join(l_val) )
      self.f_mat_file.close()

def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   python %s sample.xls Merge_RPKM_bin100.mrg.bed.gz

   """ % (sys.argv[0],sys.argv[0])

   description = " Select reads while converting bam to bed "

   optparser = OptionParser(version="%s v0.2 20150124" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)

   optparser.add_option("-r", "--region",   help="\nInput intervals.")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser

def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      chip_samp = args[0]
      RPM_matrix= args[1]

   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)

   m_revRPM = SampCol(chip_samp)
   m_revRPM.load_Matrix( RPM_matrix )
   
if __name__ == '__main__':
   main()