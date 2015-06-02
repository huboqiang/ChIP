from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
from optparse   import OptionParser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import scipy.stats
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import seaborn as sns
import pandas  as pd
sys.path.append('/datc/huboqiang/test_chipPipe')
import module_Matrix as m_mat

class MatrixCorr(object):
   def __init__(self,infile,method="pearson"):
      self.infile = infile
      self.method = method
   
   def load_matrix(self):
      m_matrix = m_mat.Matrix_info( self.infile, 3,"float",1 )
      m_matrix.load_mat()
      self.mat = m_matrix
      self.mat.matrix = np.log10( self.mat.matrix+1 )
   
   def get_cor_matrix( self,method="pearson" ):
      self.method = method
      
      out_cor_file   = "%s.corMat.%s.pdf" % ( ".".join( self.infile.split(".")[:-2] ), self.method )
      
      pd_mat = pd.DataFrame( self.mat.matrix )
      pd_mat.columns = self.mat.colname
      pd_mat.index   = self.mat.rowname
      self.cor_mat   = pd_mat.corr( self.method ).values
      
      sns.set(style="darkgrid")
      f, ax = plt.subplots(figsize=(9, 9))
      cmap = sns.diverging_palette(220, 10, as_cmap=True)
      sns.corrplot(pd_mat, annot=False, sig_stars=False, diag_names=False, cmap=cmap, ax=ax, cmap_range=(0.0, 1.0),method=self.method  )
      f.savefig( out_cor_file,format="pdf" )
      
      
      
      
   def output_cor_matrix(self):
      out_cor_file   = "%s.corMat.%s.xls" % ( ".".join( self.infile.split(".")[:-2] ), self.method )
      f_out_cor_file = open(out_cor_file,"w")
      out = "Sample\t%s" % ( "\t".join( self.mat.colname ) )
      print    >>f_out_cor_file, out
      for i,samp in enumerate( self.mat.colname ):
         out = "%s\t%s" % ( samp, "\t".join( np.array(self.cor_mat[i,:] ,dtype="string") ) )
         print >>f_out_cor_file, out
      f_out_cor_file.close()

def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   python %s /datc/huboqiang/test_chipPipe/02.RPM_Merged/Merge_RPKM_bin1kb.mrg.bed.gz

   """ % (sys.argv[0],sys.argv[0])

   description = "Printing ChIP seq correlation using read density (RPKM)."

   optparser = OptionParser(version="%s v0.1 20150131" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      infile = args[0]
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
   
   m_cor = MatrixCorr( infile )
   m_cor.load_matrix()
   m_cor.get_cor_matrix( "pearson" )
   m_cor.output_cor_matrix()
   m_cor.get_cor_matrix( "spearman" )
   m_cor.output_cor_matrix()

if __name__ == '__main__':
   main()