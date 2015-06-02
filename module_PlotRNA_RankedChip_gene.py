from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import scipy.stats
import time
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from scipy import stats

import module_Matrix as m_mat

def color_map( ltype ):
   cdictTest = {'red':( (0.00, 0.07, 0.07),
                     (0.25, 0.09, 0.09),
                     (0.50, 1.00, 1.00),
                     (0.75, 0.98, 0.98),
                     (1.00, 0.91, 0.91)),
            'green':((0.00, 0.05, 0.05),
                     (0.25, 0.36, 0.36),
                     (0.50, 1.00, 1.00),
                     (0.75, 0.71, 0.71),
                     (1.00, 0.37, 0.37)),
            'blue':( (0.00, 0.36, 0.36),
                     (0.25, 0.76, 0.76),
                     (0.50, 1.00, 1.00),
                     (0.75, 0.12, 0.12),
                     (1.00, 0.02, 0.02))}
   cdictRKG= {'red':((0.0, 1.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 1.0, 0.0)),
            'green':((0.0, 0.0, 1.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 1.0)),  
            'blue': ((0.0, 0.0, 0.5),
                     (0.5, 0.0, 0.0),
                     (1.0, 0.0, 0.0))}
   
   cdictRKB= {'red':((0.0, 1.00, 0.00),
                     (0.5, 0.69, 0.69),
                     (1.0, 1.00, 0.00)),
            'green':((0.0, 0.00, 0.50),
                     (0.5, 0.69, 0.69),
                     (1.0, 0.00, 0.00)),
            'blue': ((0.0, 0.00, 1.00),
                     (0.5, 0.84, 0.84),
                     (1.0, 0.00, 1.00))}
   
                     
   cdictRG = {'red':((0.0, 1.0, 1.0),
                     (0.2, 0.8, 0.8),
                     (0.4, 0.5, 0.5),
                     (1.0, 0.0, 0.0)),
            'green':((0.0, 1.0, 1.0),
                     (0.2, 0.9, 0.9),
                     (0.4, 0.9, 0.9),
                     (1.0, 0.8, 0.8)),
            'blue':( (0.0, 1.0, 1.0),
                     (0.2, 0.8, 0.8),
                     (0.4, 0.5, 0.5),
                     (1.0, 0.0, 0.0))} 
   
   cdictRR = {'red':((0.0, 1.0, 1.0),
                     (0.2, 0.9, 0.9),
                     (0.4, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
            'green':((0.0, 1.0, 1.0),
                     (0.2, 0.8, 0.8),
                     (0.4, 0.5, 0.5),
                     (1.0, 0.0, 0.0)),
            'blue':( (0.0, 1.0, 1.0),
                     (0.2, 0.8, 0.8),
                     (0.4, 0.5, 0.5),
                     (1.0, 0.0, 0.0))}
   
   cdictRP = {'red':((0.0, 1.0, 1.0),
                     (0.2, 0.8, 0.8),
                     (0.4, 0.5, 0.5),
                     (1.0, 0.52, 0.52)),
            'green':((0.0, 1.0, 1.0),
                     (0.2, 0.8, 0.8),
                     (0.4, 0.5, 0.5),
                     (1.0, 0.11, 0.11)),
            'blue':( (0.0, 1.0, 1.0),
                     (0.2, 0.9, 0.9),
                     (0.4, 0.9, 0.9),
                     (1.0, 0.79, 0.79))}
   
   cdictRTu={'red':( (0.0, 1.0, 1.0),
                     (0.2, 0.92, 0.92),
                     (0.4, 0.84, 0.84),
                     (1.0, 0.75, 0.75)),
            'green':((0.0, 1.0, 1.0),
                     (0.2, 0.84, 0.84),
                     (0.4, 0.68, 0.68),
                     (1.0, 0.51, 0.51)),
            'blue':( (0.0, 1.0, 1.0),
                     (0.2, 0.68, 0.68),
                     (0.4, 0.36, 0.36),
                     (1.0, 0.05, 0.05))}

   my_cmap1 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictTest,256)
   my_cmap2 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRKG ,256)
   my_cmap22= matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRKB ,256)
   my_cmap3 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRG  ,256)       
   my_cmap4 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRR  ,256)
   my_cmap5 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRP  ,256)
   my_cmap6 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRTu ,256)
   out = my_cmap1
   if   ltype == "RNA":
      out = my_cmap1
   elif ltype == "k27ac":
      out = my_cmap3
   elif ltype == "k9me3":
      out = my_cmap4
   elif ltype == "k4me3":
      out = my_cmap5
   elif ltype == "k27me3":
      out = my_cmap6
   
   return out
   
   

class PlotRNA_Ranked_ChIP_Heatmap(object):
   def __init__( self, samp_chip, dir_name, stat_Info, l_merge_name, merge_name, TSS_promoter_up=5000,TSS_promoter_down=5000,width=500,bodybin=100 ):
      self.samp_chip    = samp_chip
      self.dir_name     = dir_name
      self.l_merge_name = l_merge_name
      self.stat_Info    = stat_Info
      
      self.TSS_promoter_up   = TSS_promoter_up
      self.TSS_promoter_down = TSS_promoter_down
      self.width        = width
      self.bodybin      = bodybin

      tissue= self.samp_chip['merge_tissue'][ merge_name ]
      stage = self.samp_chip['merge_stage'][  merge_name ]
      self.name = "%s_%s" % (stage,tissue)

      self.MrgSam_mat   = {}
      self.GeneList     = []
      self.GeneRPKM     = []

   def load_samp_gene(self, read_len=101):
      l_geneList = []
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      for i,merge_name in enumerate(self.l_merge_name):
         
         infile   = "%s/%s.tss.up%d_down%d.RPKM_ranked_ChIP.allGene.xls" % ( Peak_mrg_TSS_Cor , merge_name, self.TSS_promoter_up, self.TSS_promoter_down )
         
         mat_ChIP = m_mat.Matrix_info( infile,2,'float',0 )              #   infile,  inf_column=1,  in_dtype="int",  header=1
         mat_ChIP.load_mat()
         mapped_reads = np.sum( self.stat_Info.StatInfo[ merge_name ]['unique'] )
         
         self.MrgSam_mat[ merge_name ] = mat_ChIP.matrix * self.width * 1000000 /  (read_len * mapped_reads)
         if i == 0:
            l_geneList = mat_ChIP.rowname

      for geneList in l_geneList:
         self.GeneRPKM.append( float(geneList.split()[0]) )
         self.GeneList.append(       geneList.split()[1]  )

      self.GeneList = np.array( self.GeneList,dtype="string" )
      self.GeneRPKM = np.array( self.GeneRPKM,dtype="float"  )
   
   def generate_matrix(self):
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      def plot_matrix( fig, M_matrix_info ):
         
         main_mat    = M_matrix_info['main']['matrix']
         main_mat_pos= M_matrix_info['main']['pos']
         main_cmap   = M_matrix_info['main']['cmap']
   
         color_pos      = M_matrix_info['cmap']['pos']
         color_tick     = M_matrix_info['cmap']['tick']
         color_label    = M_matrix_info['cmap']['label']
         
         print main_mat_pos
         if color_pos is None:
            upper_space      = 1 - main_mat_pos[3] - main_mat_pos[1]
            color_tag_pos_x  = main_mat_pos[0] + main_mat_pos[2]/2 - 0.05
            color_pos        = [ color_tag_pos_x, 0.05, 0.07,      0.02  ]
                  
         main_matrix = fig.add_axes( main_mat_pos )
         main_matrix.grid(True,color='black',linestyle='-',linewidth=3)
         cax01 = main_matrix.imshow(main_mat, aspect='auto', cmap=main_cmap,interpolation='nearest')
         main_matrix.get_xaxis().set_ticks([])
         main_matrix.get_yaxis().set_ticks([])
         main_matrix.get_yaxis().set_ticklabels([])
         main_matrix.get_xaxis().set_ticklabels([])
         main_matrix.tick_params(colors="white")
#         for sp in main_matrix.spines.values():
#            sp.set_visible(False)
         
         cbaxes = fig.add_axes(color_pos)
         
         cbar = plt.colorbar( cax01,cax=cbaxes, orientation='horizontal',ticks=color_tick )
         cbar.ax.yaxis.set_ticks_position('left')
         cbar.set_label( color_label, size=12)
         cbar.ax.tick_params(labelsize=10)

      fig = plt.figure(figsize=(56,20))
      l_width_beg,l_width_end = self.__get_plot_edge_x()
      
      dataMatrixOrdered = np.array( [ self.GeneRPKM ] ).T
      dataMatrixOrdered = np.log10( dataMatrixOrdered +1 )
      M_matrix_info_rpkm = { 'main': {                                                                            \
                           'matrix' :dataMatrixOrdered,                                                           \
                           'pos'    :[0.07, 0.10, 0.03, 0.80 ],                                                   \
                           'cmap'   :color_map("RNA"),                                                            \
                        },                                                                                        \
                        'cmap':{                                                                                  \
                           'pos'    :None,                                                                        \
                           'tick'   :[ 0,np.max(dataMatrixOrdered)                                  ],            \
                           'label'  :"log10(RKPM) for %s \n ChIP Signal" % self.name ,                            \
                        },                                                                                        \
      }
      plot_matrix( fig,M_matrix_info_rpkm )
      
      for i,merge_name in enumerate(self.l_merge_name):
         ltype = self.samp_chip['merge_type'][ merge_name ]
         dataMatrixOrdered1 = self.MrgSam_mat[ merge_name ]
         dataMatrixOrdered1 = np.log10( dataMatrixOrdered1 +1 )
         dataMatrixOrdered1[ dataMatrixOrdered1>0.1 ] = 0.1
         M_matrix_info_rpkm = { 'main': {                                                                            \
                              'matrix' :dataMatrixOrdered1,                                                          \
                              'pos'    :[l_width_beg[i], 0.10, l_width_end[i]-l_width_beg[i], 0.80 ],                \
                              'cmap'   :color_map(ltype),                                                            \
                           },                                                                                        \
                           'cmap':{                                                                                  \
                              'pos'    :[l_width_beg[i], 0.05, (l_width_end[i]-l_width_beg[i])*0.8, 0.02 ],          \
                              'tick'   :[ 0,np.max(dataMatrixOrdered1)                             ],                \
                              'label'  :"log10(RPM) for \n%s \n ChIP Signal" % merge_name,                           \
                           },                                                                                        \
         }
         plot_matrix( fig,M_matrix_info_rpkm )
      
      out_prefix = "%s/%s" % ( Peak_mrg_TSS_Cor,self.name )
      plt.savefig('%s.pdf' % (out_prefix),format='pdf')


   def __get_plot_edge_x(self):
      idx_beg = 0
      idx_end = 0
      pre_len = 0
      
      l_beg = []
      l_end = []
      for i,merge_name in enumerate(self.l_merge_name):
         idx_beg += pre_len
         idx_end = idx_beg+ 1
         l_beg.append( idx_beg )
         l_end.append( idx_end )
         pre_len = 1+0.1
      
      l_beg_pos = []
      l_end_pos = []
      y_end = max(l_end)*1.2
      for i in range( 0,len(l_beg) ):
         l_beg_pos.append( 0.85 - l_beg[i]*(0.85)/float(y_end) )
         l_end_pos.append( 0.85 - l_end[i]*(0.85)/float(y_end) )
   
      return l_end_pos[::-1],l_beg_pos[::-1]