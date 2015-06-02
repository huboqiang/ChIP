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

import module_StatFunc as m_Stat

class StatInfo(object):
   def __init__( self, samp_info,dir_name,TSS_promoter_up=5000,TSS_promoter_down=5000,width=500,bodybin=100 ):
      self.samp_chip  = samp_info['samp']['chip']
      self.dir_name   = dir_name
      self.group_dens = {}
      self.StatInfo = {}
      
   def Basic_Stat(self):
      '''
      Stat for QC, bwa mapping and peak calling.
      '''
      cln_dir     = self.dir_name.clean_data
      bam_dir     = self.dir_name.bam
      bed_rep_dir = self.dir_name.bed_rep
      
      stat_dir    = self.dir_name.stat
      f_out = open( "%s/01.Basic_info.xls" % (stat_dir),"w" )
      
      for merge_name in self.samp_chip['merge']:
         if merge_name not in self.StatInfo:
            self.StatInfo[ merge_name ] = {}
            
         l_totRaw = []
         l_totCln = []
         l_totDup = []
         l_totUnique = []
         l_totMulti  = []
         
         for sam in self.samp_chip['merge_sam'][merge_name]:
            QC_file = "%s/%s/log" % ( cln_dir,sam )
            QC_info = m_Stat.QcStat( QC_file )
            QC_info.read_infile()
            l_totRaw.append( QC_info.raw_reads )
            l_totCln.append( QC_info.cln_reads )

            brief_name = self.samp_chip['samp_brief'][sam]
            map_dup_info = m_Stat.BwaPicard()
            rmdup_file = "%s/%s/%s.picard_info.txt" % ( bam_dir,brief_name,brief_name )
            bed_unique = "%s/%s.unique.bed" % ( bed_rep_dir,brief_name )
            bed_multi  = "%s/%s.multi.bed" % ( bed_rep_dir,brief_name )
            
            map_dup_info.read_picard_file( rmdup_file )
            map_dup_info.read_unique_bed(  bed_unique )
            map_dup_info.read_multi_bed(   bed_multi  )
            
            l_totDup.append(    map_dup_info.dup )
            l_totUnique.append( map_dup_info.unique )
            l_totMulti.append(  map_dup_info.multi  )
            
         self.StatInfo[ merge_name ]['raw']    = np.array( l_totRaw    )
         self.StatInfo[ merge_name ]['cln']    = np.array( l_totCln    )
         self.StatInfo[ merge_name ]['dup']    = np.array( l_totDup    )
         self.StatInfo[ merge_name ]['unique'] = np.array( l_totUnique )
         self.StatInfo[ merge_name ]['multi']  = np.array( l_totMulti  )
         
         
         out = "%s\t%s\t%d\t%d\t%d\t%d\t%d" % ( merge_name, ",".join(self.samp_chip['merge_sam'][merge_name]),  np.mean(self.StatInfo[ merge_name ]['raw']),np.mean(self.StatInfo[ merge_name ]['cln']),np.mean(self.StatInfo[ merge_name ]['dup']),np.mean(self.StatInfo[ merge_name ]['unique']),np.mean(self.StatInfo[ merge_name ]['multi']) )
         print >>f_out, out

      
      f_out.close()