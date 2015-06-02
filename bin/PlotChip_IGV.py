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

import svgwrite
      

class FetchRegion(object):
   def __init__(self,l_bdg_file,region):
      self.l_bdg_file = l_bdg_file
      self.region     = region
   
      self.parse_region()
      self.__search_pos()

   def parse_region(self):
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
   
   




class SvgIGV(FetchRegion):
   def __init__(self, out_file, l_bdg_file, region):
      FetchRegion.__init__(self,l_bdg_file, region)
      
      self.left_edge_pxi   = 40
      self.right_edge_pxi  = 40
      self.top_edge_pxi    = 30
      self.bottom_edge_pxi = 60 # 30 for gene track + 30 for blank      
      
      self.frame_width     = 600
      self.frame_height    = 600
      self.figure_width    = self.left_edge_pxi + self.right_edge_pxi  + self.frame_width
      self.figure_height   = self.top_edge_pxi  + self.bottom_edge_pxi + self.frame_height
      
      self.out_file = out_file
      self.region   = region
      
      self.parse_region()
      
   def draw_box(self):
      self.dwg =svgwrite.Drawing( self.out_file,width=self.figure_width, height=self.figure_height,profile="full",debug=True )
      self.dwg.viewbox( width=self.figure_width+100,height=self.frame_height+100 )
      self.dwg.add( self.dwg.rect( (self.left_edge_pxi,self.top_edge_pxi),(self.frame_width,self.frame_height),fill="white",stroke="black" ) )
      
   def make_refGene_track(self,bed_tabix="/datd/huboqiang/ChIP_human/Week12/Database/refGene.sort.bed.gz"):
      '''
cut -f 2- refGene.txt | awk '{OFS="\t";print $1"__"$12,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13,$14,$15}' | /data/Analysis/huboqiang/software/UCSC/genePredToBed /dev/stdin /dev/stdout | bedtools sort -i /dev/stdin >refGene.sort.bed
bgzip -fc  refGene.sort.bed >refGene.sort.bed.gz
tabix -p bed -s 1 -b 2 -e 3  refGene.sort.bed.gz
      '''
      tb = tabix.open( bed_tabix )
      record = tb.query( self.chrom, self.beg, self.end )
      
      self.geneGraph = {}
      for rec in record:
         tran,gene = rec[3].split("__")
         if gene not in self.geneGraph:
            self.geneGraph[ gene ] = {}
         if tran not in self.geneGraph[gene]:
            self.geneGraph[ gene ][ tran ] = { 'beg':int(rec[1]),'end':int(rec[2]), 'exon_beg':[],'exon_ext':[],'cds_beg':[],'cds_ext':[], 'strand':rec[5]  }
         
         l_beg = [ int(beg)+int(rec[1]) for beg in rec[11].split(",")[:-1] ]
         l_ext = [ int(ext)             for ext in rec[10].split(",")[:-1] ]
         
         self.geneGraph[ gene ][ tran ][ 'exon_beg' ] = l_beg
         self.geneGraph[ gene ][ tran ][ 'exon_ext' ] = l_ext
         
         cds_beg = int(rec[6])
         cds_end = int(rec[7])
         exon_cnt= int(rec[9])
         
         if cds_beg == cds_end:
            continue
         for i in xrange( 0,exon_cnt ):
            beg = l_beg[i]
            end = l_beg[i] + l_ext[i]
            
            if cds_end < beg or cds_beg > end:
               continue
            
            self.geneGraph[ gene ][ tran ][ 'cds_beg' ].append( max(cds_beg,beg) )
            self.geneGraph[ gene ][ tran ][ 'cds_ext' ].append( min(cds_end,end)-max(cds_beg,beg) )
      
      self.__Only_Region()
      self.__Only_Longest_tran()
      
      
   def plot_Peak(self,y_lab_max=0):
      minus_idx = self.beg+1
      total_len = self.end-self.beg
      total_sam = len(self.l_bdg_file)

      y_lab_max_used = y_lab_max

      for i,bdg_file in enumerate(self.l_bdg_file):
         track_bottom_height = self.top_edge_pxi + (i+1)/total_sam*self.frame_height
         
         pd_track = self.pd_frame[ bdg_file ]
         
         l_str_path = [ "M %f,%f" % (  self.left_edge_pxi,track_bottom_height ) ]
         
         for pos,con in enumerate(self.pd_frame[ bdg_file ]['con']):
            
            if y_lab_max == 0:
               y_lab_max_used = max( self.pd_frame[ bdg_file ]['con'] )
               y_lab_max_used = int( (y_lab_max_used)*1.25 )
               
            pos_x = pos/total_len*self.frame_width + self.left_edge_pxi
            pos_y = (i+1-con/y_lab_max_used)*(self.frame_height/total_sam)+self.top_edge_pxi
            if pos % 4 == 0:
               l_str_path.append( "L %f,%f" % ( pos_x,pos_y ) )
         l_str_path.append( " L %f,%f Z" % (self.left_edge_pxi+self.frame_width,track_bottom_height) )
#         print " ".join(l_str_path)
         self.dwg.add( self.dwg.path( d=" ".join(l_str_path), fill="blue",stroke="none" ) )
      
   
   def plot_gene(self):
      """
         First, a axis and a hash table for beg-end to 0-len pix.
         This pix were then normalized into a 0-frame_width range.
         
         Line plot for gene's begin-end.
         Rect plot for exons' begin-end.
         Rect plot for cds'   begin-end.
      """
      minus_idx = self.beg+1
      total_len = self.end-self.beg
      
      gene_track_height = self.figure_height-30   # 30 for bottom blank
      
      for gene in self.geneGraphLongest:
         for tran in self.geneGraphLongest[ gene ]:
            begin = (self.geneGraphLongest[ gene ][ 'beg_new' ] - minus_idx)/total_len*self.frame_width + self.left_edge_pxi
            endin = (self.geneGraphLongest[ gene ][ 'end_new' ] - minus_idx)/total_len*self.frame_width + self.left_edge_pxi
            
            self.dwg.add( self.dwg.line( start=(begin,gene_track_height),end=(endin,gene_track_height),stroke="blue" ) )
            self.dwg.add( self.dwg.text( "%s,%s" % (gene,self.geneGraphLongest[ gene ][ 'strand' ]), ( (begin+endin)/2,gene_track_height-10 ),fill="black" ) )
            
            for i in xrange( 0,len( self.geneGraphLongest[ gene ]['exon_beg_new'] ) ):
               exon_beg = (self.geneGraphLongest[ gene ][ 'exon_beg_new' ][ i ] - minus_idx)/total_len*self.frame_width + self.left_edge_pxi
               exon_ext = (self.geneGraphLongest[ gene ][ 'exon_ext_new' ][ i ]            )/total_len*self.frame_width + self.left_edge_pxi
               self.dwg.add( self.dwg.rect( (exon_beg,gene_track_height-3),(exon_ext,6),stroke="blue",fill="blue" ) )
               
            for i in xrange( 0,len( self.geneGraphLongest[ gene ]['cds_beg_new'] ) ):
               exon_beg = (self.geneGraphLongest[ gene ][ 'cds_beg_new' ][ i ] - minus_idx)/total_len*self.frame_width + self.left_edge_pxi
               exon_ext = (self.geneGraphLongest[ gene ][ 'cds_ext_new' ][ i ]            )/total_len*self.frame_width + self.left_edge_pxi
               self.dwg.add( self.dwg.rect( (exon_beg,gene_track_height-6),(exon_ext,12),stroke="blue",fill="blue" ) )
         
      self.dwg.save()
      
   
      
   def __Only_Region(self):
      for gene in self.geneGraph:
         for tran in self.geneGraph[ gene ]:
            self.geneGraph[ gene ][ tran ][ 'beg_new' ]    = self.geneGraph[ gene ][ tran ][ 'beg' ]
            self.geneGraph[ gene ][ tran ][ 'end_new' ]    = self.geneGraph[ gene ][ tran ][ 'end' ]
            if self.geneGraph[ gene ][ tran ][ 'beg' ]     < self.beg:
               self.geneGraph[ gene ][ tran ][ 'beg_new' ] = self.beg
            if self.geneGraph[ gene ][ tran ][ 'end' ]     > self.end:
               self.geneGraph[ gene ][ tran ][ 'end_new' ] = self.end
            
            self.geneGraph[ gene ][ tran ][ 'exon_beg_new'] = []
            self.geneGraph[ gene ][ tran ][ 'exon_ext_new'] = []
            self.geneGraph[ gene ][ tran ][ 'cds_beg_new' ] = []
            self.geneGraph[ gene ][ tran ][ 'cds_ext_new' ] = []
            for exon_i in xrange( len(self.geneGraph[ gene ][ tran ][ 'exon_beg' ]) ):
               exon_beg = self.geneGraph[ gene ][ tran ][ 'exon_beg' ][ exon_i ]
               exon_end = self.geneGraph[ gene ][ tran ][ 'exon_beg' ][ exon_i ] + self.geneGraph[ gene ][ tran ][ 'exon_ext' ][ exon_i ]
               if exon_beg < self.beg:
                  exon_beg = self.beg
               if exon_end > self.end:
                  exon_end = self.end
               
               if exon_beg < exon_end:               
                  self.geneGraph[ gene ][ tran ][ 'exon_beg_new' ].append( exon_beg )
                  self.geneGraph[ gene ][ tran ][ 'exon_ext_new' ].append( exon_end-exon_beg )
               
            for cds_i in xrange( len(self.geneGraph[ gene ][ tran ][ 'cds_beg' ]) ):
               cds_beg = self.geneGraph[ gene ][ tran ][ 'cds_beg' ][ cds_i ]
               cds_end = self.geneGraph[ gene ][ tran ][ 'cds_beg' ][ cds_i ] + self.geneGraph[ gene ][ tran ][ 'cds_ext' ][ cds_i ]
               if cds_beg  < self.beg:
                  cds_beg  = self.beg
               if cds_end  > self.end:
                  cds_end  = self.end
               
               if cds_beg < cds_end:               
                  self.geneGraph[ gene ][ tran ][ 'cds_beg_new' ].append( cds_beg )
                  self.geneGraph[ gene ][ tran ][ 'cds_ext_new' ].append( cds_end-cds_beg )
      
   def __Only_Longest_tran(self):
      self.geneGraphLongest = {}
      for gene in self.geneGraph:
         self.geneGraphLongest[ gene ] = {}
         longest_tran = ""
         longest_leng = 0
         for tran in self.geneGraph[ gene ]:
            leng = self.geneGraph[gene][tran]['end'] - self.geneGraph[gene][tran]['beg']
            if leng > longest_leng:
               longest_leng = leng
               longest_tran = tran
               self.geneGraphLongest[ gene ] = self.geneGraph[ gene ][ tran ]      
               
      
      
      


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

   m_refGene = SvgIGV("test.out.svg",l_in_bdg_file,region)
   m_refGene.draw_box()
   m_refGene.make_refGene_track()
   m_refGene.plot_Peak(  )
   m_refGene.plot_gene()
   
if __name__ == '__main__':
   main()
