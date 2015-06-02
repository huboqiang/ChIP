from __future__ import division
import re,sys,os
import scipy
import numpy as np
import gzip

class Matrix_info(object):
   def __init__(self,infile,inf_column=1,in_dtype="int",header=1):
      self.infile     = infile
      self.inf_column = inf_column
      self.in_dtype   = in_dtype
      self.header     = header
      self.gzipped    = True
      self.__check_gzipped()
      
   def load_mat(self):
      self.colname= []
      l_matrix    = []
      self.rowname= []
      
      f_infile = ""
      if not self.gzipped:
         f_infile = open( self.infile,"r" )
      else:
         f_infile = gzip.open( self.infile,"rb" )
      if self.header == 1:
         head = f_infile.readline()
         head = head.strip('\n')
         f_h  = head.split()
         self.colname = f_h[ self.inf_column: ]
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         gene = "\t".join( f[0:self.inf_column] )
         self.rowname.append( gene )
         l_matrix.append( np.array( f[self.inf_column:],dtype=self.in_dtype ) )
      f_infile.close()
      self.matrix = np.array( l_matrix,dtype=self.in_dtype )

   def __check_gzipped(self):
      f_infile = gzip.open(self.infile)
      try:
         f_infile.read( 10 )
      except IOError:
         self.gzipped = False
      f_infile.close()