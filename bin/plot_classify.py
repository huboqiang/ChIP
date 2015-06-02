from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import scipy.stats
from optparse   import OptionParser

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class heatmap(dict):
   def __init__(self,in_rpkm_file,left_infor_col):
      self['infile'] = {'rpkm':in_rpkm_file}
      self['pos_info']   = {}
      self['left_infor_col']= left_infor_col
      self['pos_gene']   = {'nearest':{},'revised':{} }
   
   def load_matrix(self):
      self['pos_info']['l_gene'] = []
      self['pos_info']['tag_tissue'] = []
      self['pos_info']['tag_class']  = []
      f_infile = open( self['infile']['rpkm'],"r" )
      in_h = f_infile.readline()
      in_h = in_h.strip('\n')
      f_h  = in_h.split()
      self['sample'] = f_h[ self['left_infor_col']+1: ]
      l_mat = []
      l_tag = []
      
      index = -1
      M_tag_index = {}
      
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         gene       = f[0]
         tag_tissue = f[1]
         tag_class  = f[1]
         self['pos_info']['l_gene'].append(gene)
         self['pos_info']['tag_tissue'].append( tag_tissue )
         self['pos_info']['tag_class'].append(  tag_class  )
         
         tag = tag_tissue
         if tag not in l_tag :
            index += 1
         M_tag_index[tag] = index
         
         l_tag.append( tag )   
         
         line_vec = np.array(f[1+self['left_infor_col']: ],dtype='string')
         if np.sum(line_vec=="NA") > 0:
            continue

         l_mat.append(  np.array( line_vec,dtype=float)  )
#         l_mat.append(  np.array(f[2:-1],dtype=float)  )
         
      f_infile.close()         
      
#      self['raw_matrix'] = np.log10(np.array( l_mat,dtype=float )+1)      
      self['raw_matrix'] = np.array( l_mat,dtype=float )
      self['pos_info']['tag_tissue_less']  = []
      for tag_tissue in self['pos_info']['tag_tissue']:
         if tag_tissue not in self['pos_info']['tag_tissue_less']:
            self['pos_info']['tag_tissue_less'].append( tag_tissue )

      self['pos_info']['tag_class_less']  = []
      for tag_class in self['pos_info']['tag_class']:
         if tag_class not in self['pos_info']['tag_class_less']:
            self['pos_info']['tag_class_less'].append( tag_class )
            
      self['pos_info']['tag_index']       = []
      self['pos_info']['tag_index_less']  = []
      for tag in l_tag:
         self['pos_info']['tag_index'].append( M_tag_index[tag] )

      self['pos_info']['tag_index'] = np.array( self['pos_info']['tag_index'],dtype=int )

      #### First, we should get a "less" list for all tags
      self['pos_info']['tag_index_less'] = []
      for tag in self['pos_info']['tag_index']:
         if tag not in self['pos_info']['tag_index_less']:
            self['pos_info']['tag_index_less'].append( tag )

   def load_samp_tag(self):
#      self['tag'] = { 'l_stage':{},'l_tissue':{},'l_ltype':{},'l_ltype_less':{}, 'stage_idx':{},'tissue_idx':{},'ltype_idx':{}, 'lab_p':{},'y_pos':{},'y_lab':{} }

      self['tag'] = {}
      
      self['tag']['l_stage'] = []
      self['tag']['l_tissue']= []
      self['tag']['l_ltype'] = []
      
      dict_stage = {}
      dict_tissue = {}
      dict_ltype = {}
      
      id_stage = 0
      id_tissue= 0
      id_ltype = 0
      
      for samp in self['sample']:
         fs = samp.split('_')
         tissue = fs[0]
         stage  = fs[0]
         datatype  = "RNA"
         
         if stage  not in self['tag']['l_stage']:
            dict_stage[  stage ]  = id_stage
            id_stage             += 1
         if tissue not in self['tag']['l_tissue']:
            dict_tissue[ tissue]  = id_tissue
            id_tissue            += 1
         if datatype  not in self['tag']['l_ltype']:
            dict_ltype[ datatype  ]  = id_ltype
            id_ltype             += 1
      
         self['tag']['l_stage'].append(  stage  )
         self['tag']['l_tissue'].append( tissue )
         self['tag']['l_ltype'].append(  datatype  )
         
      ## Initiate index for each tissue 
      
      l_stage_idx  = []
      l_tissue_idx = []
      l_ltype_idx  = []
      
      y_pos = []
      lab_p = []
      y_lab = []
      
      for stage  in self['tag']['l_stage']:
         l_stage_idx.append(  dict_stage[  stage  ] )
      for tissue in self['tag']['l_tissue']:
         l_tissue_idx.append( dict_tissue[ tissue ] )   
      for i,datatype  in enumerate( self['tag']['l_ltype'] ):
         if dict_ltype[  datatype  ] not in l_ltype_idx:
            y_pos.append( i-0.5 )
            y_lab.append( datatype )
         l_ltype_idx.append(  dict_ltype[  datatype  ] )
      
      np_ltype = np.array( self['tag']['l_ltype'],dtype='string' )
      
      self['tag']['l_ltype_less'] = []
      for datatype in self['tag']['l_ltype']:
         if datatype not in self['tag']['l_ltype_less']:
            self['tag']['l_ltype_less'].append( datatype )
      
      for i,datatype  in enumerate(self['tag']['l_ltype_less']):
         lab_p.append( y_pos[i]-0.5+len( np_ltype[ np_ltype == datatype ] )/2 )
      
      self['tag']['stage_idx']  = np.array(l_stage_idx ,dtype=float)
      self['tag']['tissue_idx'] = np.array(l_tissue_idx,dtype=float)
      self['tag']['ltype_idx']  = np.array(l_ltype_idx ,dtype=float)
      
      self['tag']['lab_p'] = lab_p
      self['tag']['y_pos'] = y_pos
      self['tag']['y_lab'] = y_lab
   
   def plot_matrix(self):

      y_pos_tag,y_pos_tag_cnt,y_pos_tag_cntpos = self.__load_tag_idx()
      
      y_pos_tag_cnt     = np.array( y_pos_tag_cnt,dtype=int )
      y_pos_tag_cntpos  = np.array( y_pos_tag_cntpos,dtype=int )
      
      fig = plt.figure(figsize=(14,20))
   
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

      my_cmap3 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictTest,256)
   

      cxmatrix = fig.add_axes([ 0.15, 0.3, 0.80, 0.6])
      cxmatrix.grid(True,color='black',linestyle='-',linewidth=3)
      dataMatrixOrdered = np.log10(self['raw_matrix']+1)
      dataMatrixOrdered[ dataMatrixOrdered > 10 ] = 10
      
#      data_new = dataMatrixOrdered[:, [0,1, 2,3,4,4, 5,6,7,7, 8,9,10,10] ]
      cax02 = cxmatrix.imshow(dataMatrixOrdered, aspect='auto', cmap=my_cmap3,interpolation='nearest')
      cxmatrix.get_xaxis().set_ticks([])
      cxmatrix.get_yaxis().set_ticks(y_pos_tag)
      cxmatrix.get_yaxis().set_ticklabels([])
      cxmatrix.get_xaxis().set_ticklabels([])
      for sp in cxmatrix.spines.values():
         sp.set_visible(False)

      ax1 = fig.add_axes([ 0.13, 0.3  ,0.01, 0.6 ]  )
      ax1.get_xaxis().set_ticks( [] )
      ax1.get_yaxis().set_ticks( y_pos_tag )
      ax1.get_yaxis().set_ticklabels( [] )
      ax1.get_xaxis().set_ticklabels( [] )
      ax1.tick_params(colors="white")
      ax1.grid(True,color='white',linestyle='-',linewidth=3)
      np_show = np.array( [self['pos_info']['tag_index']],dtype=float  )
      cax1 = ax1.imshow( np_show.T, aspect='auto', cmap='jet',interpolation='nearest')
      for i,cnt in enumerate( y_pos_tag_cntpos[ y_pos_tag_cnt > 100 ] ):
         ax1.text( -1,cnt+1,'%d' % ( y_pos_tag_cnt[ y_pos_tag_cnt > 100 ][i]),ha = "right",size=15)


      ax2 = fig.add_axes([ 0.15, 0.275, 0.80, 0.005])
      
      tick_idx = [0] + [ len(self['tag']['tissue_idx'][ self['tag']['tissue_idx']<=c ]) for c in sorted(set(self['tag']['tissue_idx'])) ]
      tick_idx = np.array( tick_idx ) - 0.5
      print tick_idx
      ax2.get_xaxis().set_ticks( tick_idx )
      ax2.get_yaxis().set_ticks( [] )
      ax2.get_yaxis().set_ticklabels( [] )
      ax2.get_xaxis().set_ticklabels( [] )
      ax2.tick_params(colors="white")
      ax2.grid(True,color='white',linestyle='-',linewidth=3)   
      
      #[0, 0, 0,  0, 0,  0, 0,   0, 0,  0, 0,   0, 0, 0,   1, 1])
      
      np_stage_idx = np.array( [ self['tag']['tissue_idx'] ] )
      cax2 = ax2.imshow( np_stage_idx, aspect='auto', cmap="jet",interpolation='nearest')
      l_fpkm_sample = np.array(self['sample'],dtype='string')
      for i in  range( 0,len(l_fpkm_sample) ) :
         ax2.text( i+0.5,3,'%s' % ( l_fpkm_sample[i]),rotation=270,ha = "right",size=15)
      for sp in ax2.spines.values():
         sp.set_visible(False)
      
      cbaxes = fig.add_axes([0.25,0.95,0.10,0.02])
      cbar = plt.colorbar(cax02,cax=cbaxes, orientation='horizontal',ticks=[-2,0,2,4])
      cbar.ax.yaxis.set_ticks_position('left')

      cbar.set_label("Z-score of different \nexpressed gene-FPKM", size=12)
      cbar.ax.tick_params(labelsize=10)

      file_prefix = ".".join( self['infile']['rpkm'].split('.')[:-1] )

      plt.savefig('%s.rev.pdf' % (file_prefix),format='pdf')
      
         
   def __load_tag_idx(self):
      l_pos_ind = [-0.5]
      l_pos_cnt = []
      l_lab_p   = []     
      pos_ind   = 0
      self['pos_info']['tag_index']       = np.array( self['pos_info']['tag_index']     ,dtype='int' )
      self['pos_info']['tag_index_less']  = np.array( self['pos_info']['tag_index_less'],dtype='int' )
      
      for i,cnt in enumerate(self['pos_info']['tag_index_less']) :
         pos_ind += len(self['pos_info']['tag_index'][  self['pos_info']['tag_index']==cnt ])
         pos_ind_1 = pos_ind - 0.5
         l_pos_ind.append(pos_ind_1)
         l_pos_cnt.append(           len( self['pos_info']['tag_index'][  self['pos_info']['tag_index']==cnt ] )   )
         l_lab_p.append( pos_ind_1 - len( self['pos_info']['tag_index'][  self['pos_info']['tag_index']==cnt ] )/2 )
         
      return l_pos_ind,l_pos_cnt,l_lab_p


def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s -l 1 merge.dexseq_clean.gene.no_Mol.result.gene.embryo_zscore.xls
   
   """ % (sys.argv[0],sys.argv[0])

   description = " Plot heatmap using 0/1 correlation method. Large data-site could be hard to calculate scipy.stats.pearsonR, so try to speed up."

   optparser = OptionParser(version="%s v0.1 20141119" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-l", "--left_infor_col",default=1,help="\nInformation in the left l column. The rest column in this line is the data for plotting [default: %default]")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      fpkm_file      = args[0]
      left_infor_col = int(options.left_infor_col)
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
      
   heatmap_info = heatmap( fpkm_file,left_infor_col )
   heatmap_info.load_matrix()
   heatmap_info.load_samp_tag()
   heatmap_info.plot_matrix()
   
if __name__ == '__main__':
   main()
      
