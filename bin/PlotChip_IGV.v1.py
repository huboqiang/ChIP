from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy   as np
import pandas  as pd

import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib as mpl

from optparse   import OptionParser
import pysam
import tabix




class FetchRegion(object):
   def __init__(self,l_bdg_file,region):
      self.l_bdg_file = l_bdg_file
      self.region     = region
   
      self.__parse_region()
      self.__search_pos()
   
   def plot(self):
      out_pdf = "%s.pdf" % ( self.region  )
      
      sns.set_style("white")
      fig = plt.figure(figsize=(9,5))
      
      for i,bdg_file in enumerate( self.l_bdg_file ):
         out_xls = "%s.%s.xls" % ( ".".join( bdg_file.split(".")[:-1] ), self.region  )
         self.pd_frame[bdg_file].to_csv( out_xls,sep="\t" )
         ax = plt.subplot( len(self.l_bdg_file), 1, i+1 )
      
         bin_size = int( (self.end-self.beg)/10 )
      
         ax.stackplot( self.pd_frame[bdg_file]['pos'],self.pd_frame[bdg_file]['con'],color=None )
      
         ax.spines['right'].set_visible(False)
         ax.spines['top'].set_visible(False)
         if i != len( self.l_bdg_file )-1:
            ax.spines['bottom'].set_visible(False)
            ax.get_xaxis().set_ticks( [] )
            ax.get_xaxis().set_ticklabels( [] )
         else:
            ax.get_xaxis().set_ticks( self.l_xpos )
            ax.get_xaxis().set_ticklabels( self.l_xticks,rotation=70 )

      fig.savefig( out_pdf,format="pdf")
      
      
   
   def __parse_region(self):
      self.chrom =     re.split(":|-",self.region)[0]
      self.beg   = int(re.split(":|-",self.region)[1])-1
      self.end   = int(re.split(":|-",self.region)[2])+1
   
   def __search_pos(self):
      
      self.pd_frame = {}
      for i,bdg_file in enumerate(self.l_bdg_file):
         tb = tabix.open( bdg_file )
         record = tb.query( self.chrom, self.beg, self.end )
         l_pos  = []
         l_cons = []
         
         l_xticks = [ self.beg+1 ]
         bin_size = int( (self.end-self.beg)/10 )
         bin_size = 10**int(np.log10(bin_size))
         
         
         pre_pos  = 0
               
         for rec in record:
            for pos in xrange( int(rec[1]),int(rec[2]) ):
               cons = float(rec[3])
               
               if pre_pos == 0:
                  pre_pos = int(rec[1])
                  
               """ Only consider the given region. """
               if self.__is_intersect(pos):
                  
                  """ If bedGraph has gaps, """
                  if pos > pre_pos+1:
                     """ Using zero to fill bedGraph gaps """
                     for p in xrange( pre_pos+1,pos ):
                        l_pos.append( p )
                        l_cons.append( 0.0 )
                  
                  l_pos.append( pos  )
                  l_cons.append(cons )
                  
                  if pos % bin_size == 0:
                     l_xticks.append( pos )
                  
               pre_pos  = pos
         
         l_xticks.append( self.end )
         data = { 'pos':l_pos, 'con':l_cons }
         self.pd_frame[ bdg_file ] = pd.DataFrame( data )
         
         if i == 0:
            self.l_xpos   = l_xticks
            self.l_xticks = [ str(tick) for tick in l_xticks ]
            
            
   
   def __is_intersect(self,pos):
      idx = 0
      if pos >= self.beg+1  and pos < self.end-1:
         idx = 1
      return idx
   


def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   python %s -r chr10:121410319-121411974 /datd/huboqiang/ChIP_human/Week12/03.2.Peak_mrg/W12_brain_H3K27ac/W12_brain_H3K27ac_VS_Input_treat_pileup.sort.bdg.gz

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
      l_in_bdg_file = args
      region        = options.region

   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)

   m_reg = FetchRegion( l_in_bdg_file,region )
   m_reg.plot()
   
if __name__ == '__main__':
   main()
