#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import numpy   as np
import cPickle as pickle

class LoadSamp(dict):
   def __init__(self,sam_file_chip):
      self['infile'] = sam_file_chip
      self['samp'] = { 'chip':{},'RNA':{} }
   
   def load_samp(self):
      f_sam_file_chip = open(self['infile'],"r")
      in_h = f_sam_file_chip.readline()
      self['samp']['chip']['list'] = []
      self['samp']['chip']['tag'] = []
      self['samp']['chip']['merge'] = []
      
      self['samp']['chip']['sam_type'] = {}
      self['samp']['chip']['sam_stage'] = {}
      self['samp']['chip']['sam_tissue'] = {}
      self['samp']['chip']['samp_brief'] = {}
      self['samp']['chip']['sam_end']    = {}
      self['samp']['chip']['sam_ctrl']   = {}
      self['samp']['chip']['sam_merge']  = {}
      
      self['samp']['chip']['merge_type']   = {}
      self['samp']['chip']['merge_stage']  = {}
      self['samp']['chip']['merge_tissue'] = {}
      self['samp']['chip']['merge_ctrl']   = {}
      self['samp']['chip']['merge_sam']    = {}      
      
      self['samp']['chip']['type_sam'] = {}
      self['samp']['chip']['type_stage_sam'] = {}
      self['samp']['chip']['stage_sam'] = {}
      self['samp']['chip']['tissue_sam'] = {}
      self['samp']['chip']['tissue_type_sam'] = {}
      self['samp']['chip']['tissue_type_stage_sam'] = {}
      self['samp']['chip']['tissue_stage_type_sam'] = {}

      self['samp']['chip']['tissue_stage_merge']  =  {}
      self['samp']['chip']['type_tissue_stage_merge'] = {}
      self['samp']['chip']['tissue_stage_type_merge'] = {}
      
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
         
         tag = "%s_%s" % (tissue,stage)
         
         if tag        not in self['samp']['chip']['tag']:
            self['samp']['chip']['tag'].append(tag)

         if merge_name not in self['samp']['chip']['merge']:
            self['samp']['chip']['merge'].append(merge_name)
            
         self['samp']['chip']['list'].append( sam )
         self['samp']['chip']['sam_type'][sam]   = ltype
         self['samp']['chip']['sam_stage'][sam]  = stage
         self['samp']['chip']['sam_tissue'][sam] = tissue
         self['samp']['chip']['samp_brief'][sam] = brief
         self['samp']['chip']['sam_end'][sam]    = end_type
         self['samp']['chip']['sam_ctrl'][sam]   = control
         self['samp']['chip']['sam_merge'][sam]  = merge_name
         
         self['samp']['chip']['merge_type'][merge_name]   = ltype
         self['samp']['chip']['merge_stage'][merge_name]  = stage
         self['samp']['chip']['merge_tissue'][merge_name] = tissue
         self['samp']['chip']['merge_ctrl'][merge_name]   = control

         if merge_name not in self['samp']['chip']['merge_sam']:
            self['samp']['chip']['merge_sam'][merge_name] = []
         self['samp']['chip']['merge_sam'][merge_name].append(sam)
         
         if ltype  not in self['samp']['chip']['type_sam']:
            self['samp']['chip']['type_sam'][ltype] = []

         if ltype  not in self['samp']['chip']['type_stage_sam']:
            self['samp']['chip']['type_stage_sam'][ltype] = {}
         if stage  not in self['samp']['chip']['type_stage_sam'][ltype]:
            self['samp']['chip']['type_stage_sam'][ltype][stage] = []


         if stage  not in self['samp']['chip']['stage_sam']:
            self['samp']['chip']['stage_sam'][stage] = []
         if tissue not in self['samp']['chip']['tissue_sam']:
            self['samp']['chip']['tissue_sam'][tissue] = []

         if tissue not in self['samp']['chip']['tissue_type_sam']:
            self['samp']['chip']['tissue_type_sam'][tissue] = {}
         if ltype  not in self['samp']['chip']['tissue_type_sam'][tissue]:
            self['samp']['chip']['tissue_type_sam'][tissue][ltype] = []
            
         if tissue not in self['samp']['chip']['tissue_type_stage_sam']:
            self['samp']['chip']['tissue_type_stage_sam'][tissue] = {}
         if ltype  not in self['samp']['chip']['tissue_type_stage_sam'][tissue]:
            self['samp']['chip']['tissue_type_stage_sam'][tissue][ltype] = {}
         if stage  not in self['samp']['chip']['tissue_type_stage_sam'][tissue][ltype]:
            self['samp']['chip']['tissue_type_stage_sam'][tissue][ltype][stage] = []

         if tissue not in self['samp']['chip']['tissue_stage_type_sam']:
            self['samp']['chip']['tissue_stage_type_sam'][tissue] = {}
         if stage  not in self['samp']['chip']['tissue_stage_type_sam'][tissue]:
            self['samp']['chip']['tissue_stage_type_sam'][tissue][stage] = {}
         if ltype  not in self['samp']['chip']['tissue_stage_type_sam'][tissue][stage]:
            self['samp']['chip']['tissue_stage_type_sam'][tissue][stage][ltype] = []
            
         if tissue not in self['samp']['chip']['tissue_stage_merge']:
            self['samp']['chip']['tissue_stage_merge'][tissue] = {}
         if stage  not in self['samp']['chip']['tissue_stage_merge'][tissue]:
            self['samp']['chip']['tissue_stage_merge'][tissue][stage] = []
            
         if ltype  not in self['samp']['chip']['type_tissue_stage_merge']:
            self['samp']['chip']['type_tissue_stage_merge'][ltype] = {}
         if tissue not in self['samp']['chip']['type_tissue_stage_merge'][ltype]:
            self['samp']['chip']['type_tissue_stage_merge'][ltype][tissue] = {}

         if tissue not in self['samp']['chip']['tissue_stage_type_merge']:
            self['samp']['chip']['tissue_stage_type_merge'][tissue] = {}
         if stage  not in self['samp']['chip']['tissue_stage_type_merge'][tissue]:
            self['samp']['chip']['tissue_stage_type_merge'][tissue][stage] = {}

         self['samp']['chip']['type_sam'][ltype].append(sam)
         self['samp']['chip']['type_stage_sam'][ltype][stage].append(sam)
         self['samp']['chip']['stage_sam'][stage].append(sam)
         self['samp']['chip']['tissue_sam'][tissue].append(sam)
         self['samp']['chip']['tissue_type_sam'][tissue][ltype].append(sam)
         self['samp']['chip']['tissue_type_stage_sam'][tissue][ltype][stage].append(sam)
         self['samp']['chip']['tissue_stage_type_sam'][tissue][stage][ltype].append(sam)
         self['samp']['chip']['type_tissue_stage_merge'][ltype][tissue][stage] = merge_name         
         self['samp']['chip']['tissue_stage_merge'][tissue][stage].append( merge_name )
         self['samp']['chip']['tissue_stage_type_merge'][tissue][stage][ltype] = merge_name
      
      f_sam_file_chip.close()
      
   def load_RNA(self,infile):
      f_sam_file_RPKM = open( infile,"r")
      in_h = f_sam_file_RPKM.readline()
      self['samp']['RNA']['samp_brief'] = {}
      self['samp']['RNA']['type']     = {}
      self['samp']['RNA']['stage']    = {}
      self['samp']['RNA']['dilute']   = {}
      self['samp']['RNA']['stage_tis_sam']    = {}
      self['samp']['RNA']['RFP_mols'] = {}
      self['samp']['RNA']['GFP_mols'] = {}
      self['samp']['RNA']['CRE_mols'] = {}
      
      self['samp']['RNA']['stage_name']= []
      self['samp']['RNA']['sample']    = []
      for line in f_sam_file_RPKM:
         line = line.strip('\n')
         f = line.split()
         samp       = f[0]
         brief_name = f[1]
         stage_tis  = f[2]
         ltype      = f[3]
         ERCC_dilute= float( f[4] )
         RFP_mols   = float( f[5] )
         GFP_mols   = float( f[6] )
         CRE_mols   = float( f[7] )

         self['samp']['RNA']['sample'].append( samp )
         self['samp']['RNA']['samp_brief'][ samp ] = brief_name
         self['samp']['RNA']['type'][  samp ]      = ltype
         self['samp']['RNA']['stage'][ samp ]      = stage_tis
         self['samp']['RNA']['dilute'][samp ]      = ERCC_dilute
         self['samp']['RNA']['RFP_mols'][ samp ]   = RFP_mols
         self['samp']['RNA']['GFP_mols'][ samp ]   = GFP_mols
         self['samp']['RNA']['CRE_mols'][ samp ]   = CRE_mols
         
         if stage_tis not in self['samp']['RNA']['stage_name']:
            self['samp']['RNA']['stage_name'].append( stage_tis )
            
         if stage_tis not in self['samp']['RNA']['stage_tis_sam']:
            self['samp']['RNA']['stage_tis_sam'][stage_tis] = []
         self['samp']['RNA']['stage_tis_sam'][stage_tis].append( samp )
         
      f_sam_file_RPKM.close()
