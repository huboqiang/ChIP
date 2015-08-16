#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import gzip
import string
import subprocess
import time
import cPickle as pickle

import numpy   as np

import ChIP.utils.module_refGene            as m_refG
import ChIP.utils.module_running_jobs       as m_jobs
import ChIP.settings.scripts                as m_scpt
import ChIP.frame.module01_mapping_from_raw as m01
import ChIP.frame.module00_StatInfo         as m00


class GeneBasedInfo(m_scpt.Scripts):

    def __init__(self, sam_ChIPinfo, ref, is_debug = 1):
        super(GeneBasedInfo,self).__init__()
        self.is_debug = is_debug
        self.sam_ChIPinfo= sam_ChIPinfo
        s_idx = ".".join(sam_ChIPinfo.split("/")[-1].split(".")[:-1])
        self.load_ChIP_samInfo( sam_ChIPinfo )
        self.define_scripts(s_idx)
        self.define_files(ref)
        
        self.__load_StatInfo()
                

    def __load_StatInfo(self):
        prefix  = ".".join(self.sam_ChIPinfo.split(".")[:-1])
        data_db = "%s.dat" % (prefix)
        try:
            self.stat_Info = pickle.load(open(data_db))
        except:
            self.stat_Info = m00.StatInfo(self.sam_ChIPinfo)
            self.stat_Info.Basic_Stat()
            pickle.dump(self.stat_Info,open(data_db,"wb"),True)
      
    def extend_gene_region(self,
            TSS_genebody_up,TSS_genebody_down,
            TSS_promoter_up,TSS_promoter_down):
      
        data_db = "%s.dat" % (self.refGeneTxt)
        try:
            refG_info = pickle.load(open(data_db))
        except:
            refG_info = m_refG.RefGeneTxt(self.refGeneTxt,self.genome_ref)
            refG_info.refGene2bed(0,0)
            refG_info.refGene2bed(
                TSS_genebody_up, TSS_genebody_down, ext_type="genebody.")

            refG_info.refGene2bed(
                TSS_promoter_up, TSS_promoter_down, ext_type="promoter.")

            refG_info.refGeneInfo(TSS_promoter_up, TSS_promoter_down )
            refG_info.Only_LongestTid_Bed(
                TSS_genebody_up, TSS_genebody_down, ext_type="genebody.")

            refG_info.Only_LongestTid_Bed(
                TSS_promoter_up, TSS_promoter_down, ext_type="promoter.")

            pickle.dump( refG_info,open(data_db,"wb"),True )


    def run_anno_peak(self,  
                        TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,
                        TSS_promoter_down,ext_binlen=50,body_bincnt=100,
                        tss_binlen=1):
      
        sh_file       = "%s/s12.PeakGeneRegion.sh"      % (self.scripts)
        sh_work_file  = "%s/s12.PeakGeneRegion_work.sh" % (self.scripts)

        l_sh_info     = self.s12_PeakGeneRegion(
                        TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,
                        TSS_promoter_down,ext_binlen=50,body_bincnt=100,
                        tss_binlen=1)
                        
        l_sh_work     = []

        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            m01.make_dir([ self.dir_Peak_mrg_TSS,  merge_name ])
            m01.make_dir([ self.dir_Peak_mrg_Gene, merge_name ])
            l_sh_work.append("sh %s %s" % ( sh_file, merge_name ))
            
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#        my_job.running_SGE(vf="400m",maxjob=100,is_debug = self.is_debug)

    def div_bed_to_bins_unique_rep(self):
        sh_file      = "%s/s13.RPM_density_rep.sh"      % (self.scripts)
        sh_work_file = "%s/s13.RPM_density_rep_work.sh" % (self.scripts)

        l_sh_info = self.s13_RPM_density_rep()
        l_sh_work = []
        
        for brief_name in self.samInfo_pd_ChIP['brief_name']:
            m01.make_dir([ self.dir_RPM_bins_rep, brief_name])
            
            idx = (self.samInfo_pd_ChIP['brief_name'] == brief_name)
            merge_name = self.samInfo_pd_ChIP[ idx ]['merge_name'].values[0]
            
            l_brief = self.stat_Info.StatInfo[merge_name]['l_brief']
            idx2 = l_brief.index(brief_name)
            mapped_reads = self.stat_Info.StatInfo[merge_name]['q30'][idx2]
            l_sh_work.append("sh %s %s %d" % (sh_file, brief_name, mapped_reads))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#        my_job.running_SGE(vf="400m",maxjob=100,is_debug = self.is_debug)

    def div_bed_to_bins_unique_mrg(self):
        sh_file      = "%s/s14.RPM_density_rep.sh"      % (self.scripts)
        sh_work_file = "%s/s14.RPM_density_rep_work.sh" % (self.scripts)

        l_sh_info = self.s14_RPM_density_mrg()
        l_sh_work = []
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            m01.make_dir([ self.dir_RPM_bins_mrg, merge_name])
            
            mapped_reads = np.sum(self.stat_Info.StatInfo[merge_name]['q30'])
            l_sh_work.append("sh %s %s %d" % (sh_file,merge_name, mapped_reads))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#        my_job.running_SGE(vf="400m",maxjob=100,is_debug = self.is_debug)


    def merge_RPKM_uniq( self ):

        m01.make_dir([ self.dir_RPM_mrg ])
        sh_file      = "%s/s15.MergeRPKM.sh"        % (self.scripts)
        sh_work_file = "%s/s15.1.MergeRPKM_work.sh" % (self.scripts)
        
        l_brief = self.samInfo_pd_ChIP['brief_name']
        l_merge = set(self.samInfo_pd_ChIP['merge_name'])
        
        l_sh_info = self.s15_merge_RPKM()
        l_sh_work = []
        
        for window in [ "100","1kb" ]:
            for ltype in ["rep","mrg"]:
                l_sam = l_brief
                RPKM_dir = self.dir_RPM_bins_rep
                if ltype == "mrg":
                    l_sam = l_merge
                    RPKM_dir = self.dir_RPM_bins_mrg
                
                header = "\"#chr\\tbeg\\tend\\t%s\"" % ("\\t".join(l_sam))
                l_RPKM_file = [ 
                    "%s/%s/%s.RPKM.uniq.%s"                                 %\
                    (RPKM_dir,sam,sam,window) for sam in l_sam
                ]
                l_sh_work.append(
                    "sh %s  %s %s %s %s"                                    %\
                    (sh_file, header, window, ltype, " ".join(l_RPKM_file))
                )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#        my_job.running_SGE(vf="400m",maxjob=100,is_debug = self.is_debug)
