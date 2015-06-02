from __future__ import division
import re,sys,os,gzip,string
import subprocess
import time
import cPickle as pickle
import numpy   as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import module_refGene      as m_refG
import module_running_jobs as m_jobs
import module_PlotChIP_RNA as m_plot
import module_PlotRNA_RankedChip_gene as m_heatmap

class ChIP_RNAcor(object):
   def __init__( self,  samp_info, dir_name,stat_Info ):
      self.samp_info = samp_info
      self.samp_chip = samp_info['samp']['chip']
      self.samp_RNA  = samp_info['samp']['RNA']
      self.dir_name  = dir_name
      self.stat_Info = stat_Info
      
      self.samp_gene_ChIP = {}
      self.samp_gene_RPKM = {}
      self.samp_gene_Class= {}
      self.samp_gene_RNA_Rank  = {}
      self.samp_gene_ChIP_Rank = {}
      
      self.samp_RPKM_rank = {}
      self.samp_gene_rank = {}
      self.samp_ChIP_rank = {}
      self.samp_Classify  = {}
      
      self.all_gene = []
      self.all_RPKM = []
   
   def load_data(self, TSS_promoter_up=5000,TSS_promoter_down=5000,width=500,bodybin=100):
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      if not os.path.isdir( Peak_mrg_TSS_Cor ):
         os.mkdir( Peak_mrg_TSS_Cor )
      for tissue in self.samp_chip['tissue_stage_merge']:
         for stage in self.samp_chip['tissue_stage_merge'][tissue]:
            l_merge_name = self.samp_chip['tissue_stage_merge'][tissue][stage]
            for merge_name in l_merge_name:
               self.__load_ChIP_TSS( merge_name,TSS_promoter_up,TSS_promoter_down )
               stage_tis = "%s_%s" % ( stage,tissue )
               l_RNA_samples = self.samp_RNA['stage_tis_sam'][stage_tis]
               self.__load_RNA_RPKM( merge_name, l_RNA_samples )
               self.__commonGenes( merge_name )
            self.__allGenes( l_merge_name, TSS_promoter_up,TSS_promoter_down,width )
   
   def output_RPKM_rank_ChIP(self,TSS_promoter_up=5000,TSS_promoter_down=5000):
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      
      if not os.path.isdir( Peak_mrg_TSS_Cor ):
         os.mkdir( Peak_mrg_TSS_Cor )
         
      for merge_name in self.samp_Classify:
         
         ChIP_mean_TSS   = "%s/%s.tss.up%d_down%d.Ranked_ChIP_Mean.xls" % ( Peak_mrg_TSS_Cor , merge_name, TSS_promoter_up, TSS_promoter_down )
         f_ChIP_mean_TSS = open( ChIP_mean_TSS,"w" )
         for i in [ 2,1,0 ]:
            idx         = ( self.samp_Classify[ merge_name ]==i )
            np_ChIP     = self.samp_ChIP_rank[ merge_name ][ idx, ]
            np_ChIP_mean= np.mean( np_ChIP,axis=0 )
            out_mean    = "\t".join( np.array(np_ChIP_mean,dtype="string") )
            out = "%d\t%s" % (i,out_mean)
            print >>f_ChIP_mean_TSS, out
         f_ChIP_mean_TSS.close()
   
   def plot_RPKM_rank_ChIP(self,TSS_promoter_up=5000,TSS_promoter_down=5000,width=500,bodybin=100 ):
      PlotChIP = m_plot.PlotChIPRNA( self.dir_name, TSS_promoter_up, TSS_promoter_down, width, bodybin )
      for merge_name in self.samp_chip['merge']:
         mapped_reads = np.sum( self.stat_Info.StatInfo[ merge_name ]['unique'] )
         PlotChIP.load_samp_class( merge_name )
         PlotChIP.plot_class( merge_name,mapped_reads )
         PlotChIP.clean_class(merge_name )
      
   def plot_RPKM_ChIP_heatmap(self,TSS_promoter_up=5000,TSS_promoter_down=5000,width=500,bodybin=100):
#      for merge_name in ['W12_brain_H3K4me1','W12_heart_H3K4me1','W12_liver_H3K4me1']:
      for merge_name in self.samp_chip['merge']:
         PlotHeatmap = m_heatmap.PlotRNA_Ranked_ChIP_Heatmap( self.samp_chip, self.dir_name, self.stat_Info, self.samp_chip['merge'],merge_name, TSS_promoter_up, TSS_promoter_down, width, bodybin )
         PlotHeatmap.load_samp_gene()
         PlotHeatmap.generate_matrix()
      
   
   def __load_ChIP_TSS(self,merge_name, TSS_promoter_up,TSS_promoter_down):
      self.samp_gene_ChIP[ merge_name ] = {}
      Peak_mrg_TSS = self.dir_name.Peak_mrg_TSS
      ChIP_TSS     = "%s/%s/%s.tss.up%d_down%d.xls" % ( Peak_mrg_TSS , merge_name,merge_name, TSS_promoter_up, TSS_promoter_down )
      f_ChIP_TSS   = open( ChIP_TSS,"r" )
      f_ChIP_TSS.readline()
      f_ChIP_TSS.readline()
      for line in f_ChIP_TSS:
         line   = line.strip('\n')
         f      = line.split()
         gene   = f[0]
         np_dens= np.array( f[3:],dtype='float' )
         self.samp_gene_ChIP[ merge_name ][ gene ] = np_dens
      f_ChIP_TSS.close()
      
   def __load_RNA_RPKM(self,merge_name, l_RNA_samples ):
      self.samp_gene_RPKM[merge_name] = {}
      for sam in l_RNA_samples:
         samp_RPKM    = "/datc/huboqiang/human_embryo_sequencing/RNA_new/07.Div_RPKM_Group/%s.fpkm_group.xls.sort" % ( sam )
         f_samp_RPKM  = open( samp_RPKM,"r" )
         for line in f_samp_RPKM:
            line   = line.strip('\n')
            f      = line.split()
            gene   = f[1]
            np_RPKM= np.array( f[2],dtype='float' )
            if gene not in self.samp_gene_RPKM[ merge_name ]:
               self.samp_gene_RPKM[ merge_name ][ gene ] = []
            self.samp_gene_RPKM[ merge_name ][ gene ].append( np_RPKM )
         f_samp_RPKM.close()
      
   def __commonGenes(self,merge_name):
      
      l_gene = []
      l_RPKM = []
      l_ChIP = []
      
      l_all_gene = []
      l_all_RPKM = []
      
      for gene in self.samp_gene_RPKM[ merge_name ]:
         RPKM_val = sum(self.samp_gene_RPKM[ merge_name ][ gene ])/len(self.samp_gene_RPKM[ merge_name ][ gene ])
         l_all_gene.append( gene )
         l_all_RPKM.append( RPKM_val )
         if gene not in self.samp_gene_ChIP[ merge_name ]:
            continue
         ChIP_val = self.samp_gene_ChIP[ merge_name ][ gene ]
         l_gene.append( gene     )
         l_RPKM.append( RPKM_val )
         l_ChIP.append( ChIP_val )
         
      np_gene       = np.array( l_gene,dtype="string" )
      np_RPKM       = np.array( l_RPKM,dtype="float"  )
      np_ChIP       = np.array( l_ChIP,dtype="float"  )
      np_all_gene   = np.array( l_all_gene,dtype="string" )
      np_all_RPKM   = np.array( l_all_RPKM,dtype="float"  )

      self.all_gene = np_all_gene[ np_all_RPKM.argsort()[::-1] ]
      self.all_RPKM = np_all_RPKM[ np_all_RPKM.argsort()[::-1] ]

      self.samp_RPKM_rank[ merge_name ] = np_RPKM[ np_RPKM.argsort()[::-1] ]
      self.samp_gene_rank[ merge_name ] = np_gene[ np_RPKM.argsort()[::-1] ]
      self.samp_ChIP_rank[ merge_name ] = np_ChIP[ np_RPKM.argsort()[::-1] ]
      
      self.samp_Classify[  merge_name ]  = np.zeros( len(np_RPKM),dtype="int" )
      self.samp_Classify[  merge_name ][ self.samp_RPKM_rank[ merge_name ] > 1 ] = 1
      self.samp_Classify[  merge_name ][ self.samp_RPKM_rank[ merge_name ] > 10] = 2
      

   def __allGenes( self,l_merge_name,   TSS_promoter_up=5000,TSS_promoter_down=5000,width=500 ):
      bin_cnt= (TSS_promoter_up+TSS_promoter_down)/width
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      
      for i,merge_name in enumerate(l_merge_name):
         outFile = "%s/%s.tss.up%d_down%d.RPKM_ranked_ChIP.allGene.xls" % ( Peak_mrg_TSS_Cor , merge_name, TSS_promoter_up, TSS_promoter_down )
         f_file  = open( outFile,"w" )
      
         for j,gene in enumerate(self.all_gene):
            RPKM_val = self.all_RPKM[j]
            
            if gene not in self.samp_gene_ChIP[ merge_name ]:
               self.samp_gene_ChIP[ merge_name ][ gene ] = np.zeros( int(bin_cnt),dtype="int" )
            ChIP_val = "\t".join( np.array(self.samp_gene_ChIP[ merge_name ][ gene ],dtype="string") )
            print >>f_file, "%f\t%s\t%s" % ( RPKM_val,gene,ChIP_val )
      
         f_file.close()