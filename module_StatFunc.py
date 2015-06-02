from __future__ import division
import re,sys,os
import subprocess,time
import numpy as np

class QcStat(object):
   def __init__(self,infile):
      self.infile = infile
      self.__init_statInfo()
   
   def __init_statInfo(self):
      self.raw_reads = 0
      self.raw_bases = 0
      self.cln_reads = 0
      self.cln_bases = 0
      self.ER_l_rate = 0.0
      self.ER_r_rate = 0.0
      self.Q20_l_rate= 0.0
      self.Q20_r_rate= 0.0
      self.Q30_l_rate= 0.0
      self.Q30_r_rate= 0.0
      self.GC_l_rate = 0.0
      self.GC_r_rate = 0.0
      self.N_remove  = 0
      self.Q_remove  = 0
      self.A_remove  = 0
      self.A_trimed  = 0
      
   def read_infile(self):
      if os.path.isfile( self.infile ):
         f_infile = open( self.infile,"r" )
         l_line = f_infile.readlines()
         f_1    = l_line[1].split('\t')
         self.raw_reads = int(f_1[0])
         self.raw_bases = int(f_1[1])
         self.cln_reads = int(f_1[2])
         self.cln_bases = int(f_1[3])
                  
         if len( f_1[5].split(';') ) > 1:
            self.ER_l_rate, self.ER_r_rate  = [ float(i) for i in f_1[5].split(';') ]
            self.Q20_l_rate,self.Q20_r_rate = [ float(i) for i in f_1[6].split(';') ]
            self.Q30_l_rate,self.Q30_r_rate = [ float(i) for i in f_1[7].split(';') ]
            self.GC_l_rate, self.GC_r_rate  = [ float(i) for i in f_1[8].split(';') ]
         else:
            self.ER_l_rate  = [ float(i) for i in f_1[5].split(';') ]
            self.Q20_l_rate = [ float(i) for i in f_1[6].split(';') ]
            self.Q30_l_rate = [ float(i) for i in f_1[7].split(';') ]
            self.GC_l_rate  = [ float(i) for i in f_1[8].split(';') ]
            
         self.N_remove  = l_line[2].split()[-1]
         self.Q_remove  = l_line[3].split()[-1]
         self.A_remove  = l_line[4].split()[2]
         self.A_trimed  = l_line[4].split()[-1]
         f_infile.close()

class BwaPicard(object):
   def __init__(self):
      self.dup    = 0
      self.unique = 0
      self.multi  = 0
      
   def read_picard_file( self,picard_info ):
      if os.path.isfile( picard_info ):
         f_infile = open( picard_info,"r" )
         l_line = f_infile.readlines()
         f      = l_line[7].split("\t")
         f_num  = []
         for num in f[1:]:
            if len(num.split())>0:
               f_num.append( float(num) )
         
         np_line= np.array( f_num,dtype='float' )
         self.dup += np_line[4] + np_line[5]*2
         f_infile.close()
   
   def read_unique_bed( self,unique_bed ):
      self.unique = self.__get_line_cnt( unique_bed )
      
   def read_multi_bed(  self,multi_bed  ):
      self.multi  = self.__get_line_cnt( multi_bed  )
   
   def __get_line_cnt( self,infile_bed ):
      shell_work = 'wc -l %s' % ( infile_bed )
      p       = subprocess.Popen(shell_work,shell='True',stdout=subprocess.PIPE)
      l_lines = p.stdout.readlines()
      lincnt  = int( l_lines[0].split()[0] )
      p.stdout.close()
      return lincnt
      
#class PeakInfo(object):
#   def __init__(self):
   
   