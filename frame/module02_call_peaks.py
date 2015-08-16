#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import time

import numpy   as np

import ChIP.settings.scripts                as m_scpt
import ChIP.frame.module00_StatInfo         as m00
import ChIP.frame.module01_mapping_from_raw as m01
import ChIP.utils.module_running_jobs       as m_jobs

def get_shiftSize(infile):
    f_infile = open(infile,"r")
    l_line   = f_infile.readlines()
    line     = l_line[0]
    line     = line.strip()
    f        = line.split()
    f_infile.close()
    l_val = np.array(f[2].split(','),dtype=int)
    return int(np.mean(l_val)/2)
    


class Macs2Peaks(m_scpt.Scripts):
   
    def __init__(self, sam_ChIPinfo, ref, is_debug = 1):
        super(Macs2Peaks,self).__init__()
        
        self.sam_ChIPinfo= sam_ChIPinfo
        
        s_idx = ".".join(sam_ChIPinfo.split("/")[-1].split(".")[:-1])
        self.load_ChIP_samInfo( sam_ChIPinfo )
        self.define_scripts(s_idx)
        self.define_files(ref)
        
        self.is_debug    = is_debug
       
    def get_shift_size_rep(self):
        sh_file       = "%s/s05.1.spp_rep_shiftSize.sh"      % (self.scripts)
        sh_work_file  = "%s/s05.1.spp_rep_shiftSize_work.sh" % (self.scripts)

        l_sh_info     = self.s05_1_spp_rep_shiftSize()
        l_sh_work     = []
       
        for brief_name in self.samInfo_pd_ChIP['brief_name']:
            idx   =(self.samInfo_pd_ChIP['brief_name'] == brief_name)
            if self.__is_input(idx):
                continue
                
            m01.make_dir([self.dir_spp_rep_shiftSize, brief_name, "test"])
            m01.make_dir([self.dir_spp_rep_shiftSize, brief_name, "out" ])
            l_sh_work.append("sh %s %s" % ( sh_file, brief_name ))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
#        my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)

    
    def get_shift_size_mrg(self):
        sh_file      = "%s/s05.2.spp_mrg_shiftSize.sh"      % (self.scripts)
        sh_work_file = "%s/s05.2.spp_mrg_shiftSize_work.sh" % (self.scripts)
        
        l_sh_info    = self.s05_2_spp_mrg_shiftSize()
        l_sh_work    = []
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx   =(self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue

            m01.make_dir([self.dir_spp_mrg_shiftSize, merge_name, "test"])
            m01.make_dir([self.dir_spp_mrg_shiftSize, merge_name, "out" ])       
            l_sh_work.append("sh %s  %s" % ( sh_file, merge_name ))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug=self.is_debug)
#        my_job.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)
    
    
    def run_macs_rep( self, pvalue=0.001, ref="hs" ):
        sh_file      = "%s/s06.1.macs2PeakRep.sh"      % (self.scripts)
        sh_work_file = "%s/s06.1.macs2PeakRep_work.sh" % (self.scripts)
        
        l_sh_info    = self.s06_1_macs2PeakRep( ref )
        l_sh_work    = []
    
        for brief_name in self.samInfo_pd_ChIP['brief_name']:
            idx   =(self.samInfo_pd_ChIP['brief_name'] == brief_name)
            if self.__is_input(idx):
                continue
            
            m01.make_dir([self.dir_Peak_rep, brief_name])
            ctrl_name  = self.samInfo_pd_ChIP[ idx ]['control'].values[0]
            
            shift_size = 300
            f_shiftSize = "%s/%s/out/out.tab"                               %\
                (self.dir_spp_rep_shiftSize, brief_name )
                
            if os.path.isfile( f_shiftSize ):
                val = get_shiftSize( f_shiftSize )
                if val > 0:
                    shift_size  = val
            
            l_sh_work.append("sh %s  %s %s %f %d"                           %\
                (sh_file, brief_name, ctrl_name, pvalue, shift_size))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug=self.is_debug)
#        my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)


    def run_macs_mrg( self, pvalue=0.001, ref="hs" ):
        sh_file      = "%s/s06.2.macs2PeakMrg.sh"      % (self.scripts)
        sh_work_file = "%s/s06.2.macs2PeakMrg_work.sh" % (self.scripts)
        
        l_sh_info    = self.s06_2_macs2PeakMrg( ref )
        l_sh_work    = []
    
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx = (self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue
            
            m01.make_dir([self.dir_Peak_mrg, merge_name])
            ctrl_name  =  self.samInfo_pd_ChIP[idx]['control'].values[0]
            
            shift_size = 300
            f_shiftSize = "%s/%s/out/out.tab"                               %\
                (self.dir_spp_mrg_shiftSize, merge_name)
                
            if os.path.isfile( f_shiftSize ):
                val = get_shiftSize( f_shiftSize )
                if val > 0:
                    shift_size  = val
            
            l_sh_work.append("sh %s  %s %s %f %d"                           %\
                (sh_file, merge_name, ctrl_name, pvalue, shift_size))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug=self.is_debug)
#        my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)

    def prepare_idr_input(self, top_peak = 100000):
        sh_file       = "%s/s07.1.IDR_prepare.sh"      % (self.scripts)
        sh_work_file  = "%s/s07.1.IDR_prepare_work.sh" % (self.scripts)

        l_sh_info     = self.s07_1_IDR_prepare()
        l_sh_work     = []
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx   =(self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue
            
            m01.make_dir( [ self.dir_Peak_idr, merge_name ] )
            l_brief   = list(self.samInfo_pd_ChIP[ idx ]['brief_name'])
            list_brief = " ".join(l_brief)
            l_sh_work.append(
                "sh %s  %s %d %s" % (sh_file, merge_name, top_peak,list_brief)
            )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug=self.is_debug)
#        my_job.running_SGE(vf="200m", maxjob=100, is_debug=self.is_debug)


    def run_idr(self):
        sh_file            =                                                 \
            "%s/s07.2.IDR_usingMacs2Peak.sh"               % (self.scripts)
            
        sh_work_file_rep   =                                                 \
            "%s/s07.2.IDR_usingMacs2Peak_rep_work.sh"      % (self.scripts)
            
        sh_work_file_self  =                                                 \
            "%s/s07.2.IDR_usingMacs2Peak_selfReps_work.sh" % (self.scripts)
            
        sh_work_file_pool  =                                                 \
            "%s/s07.2.IDR_usingMacs2Peak_poolReps_work.sh" % (self.scripts)
            
        
        l_sh_info      = self.s07_2_usingMacs2Peak()
        l_sh_work_rep  = []
        l_sh_work_self = []
        l_sh_work_pool = []
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx   =(self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue
            
            m01.make_dir([self.dir_Peak_idr, merge_name])
            l_brief = list(self.samInfo_pd_ChIP[idx]['brief_name'])
#            print merge_name,l_brief,list(l_brief)
            ### IDR ANALYSIS ON ORIGINAL REPLICATES
            subtype = "rep"
            m01.make_dir([self.dir_Peak_idr, merge_name, subtype])
            if len(l_brief) > 1:
                for i    in range( 0,len(l_brief)-1 ):
                    for j in range( i+1,len(l_brief) ):
                        in_sam1 = "%s" % (l_brief[i])
                        in_sam2 = "%s" % (l_brief[j])
                        l_sh_work_rep.append("sh %s  %s %s %s %s"           %\
                            (sh_file,in_sam1,in_sam2,merge_name,subtype))

            ### IDR ANALYSIS ON SELF-PSEUDOREPLICATES
            subtype = "selfPseudoReps"
            m01.make_dir([self.dir_Peak_idr, merge_name, subtype])
            for sam in l_brief:
                in_sam1 = "%s.pr1" % (sam)
                in_sam2 = "%s.pr2" % (sam)
                l_sh_work_self.append("sh %s   %s %s %s %s"                 %\
                     (sh_file, in_sam1, in_sam2, merge_name, subtype))
                
            
            ### IDR ANALYSIS ON POOLED-PSEUDOREPLICATES
            subtype = "pooledPseudoReps"
            m01.make_dir([ self.dir_Peak_idr, merge_name, subtype ])
            
            in_sam1 = "%s.pr1" % (merge_name)
            in_sam2 = "%s.pr2" % (merge_name)
            l_sh_work_pool.append("sh %s   %s %s %s %s"                     %\
                 (sh_file, in_sam1, in_sam2, merge_name, subtype))
        
        my_job_rep        = m_jobs.run_jobs(                                 \
            sh_file, sh_work_file_rep,  l_sh_info, l_sh_work_rep)

        my_job_self = m_jobs.run_jobs(                                       \
            sh_file, sh_work_file_self, l_sh_info, l_sh_work_self)

        my_job_pool = m_jobs.run_jobs(                                       \
            sh_file, sh_work_file_pool, l_sh_info, l_sh_work_pool)
            
        my_job_rep.running_multi(cpu=8, is_debug=self.is_debug)
        my_job_self.running_multi(cpu=8,is_debug=self.is_debug)
        my_job_pool.running_multi(cpu=8,is_debug=self.is_debug)
#       my_job_self.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)
#       my_job_rep.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)
#       my_job_pool.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)


    def get_idr_stat(self):
        mod_Stat = m00.StatInfo(self.sam_ChIPinfo)
        mod_Stat.IDR_Stat()

    def get_idr_Peak(self):
        prefix      = ".".join( self.sam_ChIPinfo.split(".")[:-1] )
        file_idr_out= "%s/IDR_result.%s.xls" % (self.dir_StatInfo, prefix)
        f_idr_out   = open(file_idr_out,"r")
        
        sh_file       = "%s/s07.3.IDR_pass_Peaks.sh"      % (self.scripts)
        sh_work_file  = "%s/s07.3.IDR_pass_Peaks_work.sh" % (self.scripts)
        
        l_sh_info     = self.s07_3_IDR_passPeaks()
        l_sh_work     = []
        
        h = f_idr_out.readline()
        for line in f_idr_out:
            line      = line.strip()
            f         = line.split()
            merge_name= f[0]
            peak_cnt  = int(f[1])
            
            l_sh_work.append("sh %s  %s %d" %  (sh_file, merge_name,peak_cnt))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi( cpu=8, is_debug = self.is_debug )
#        my_job.running_SGE( vf="200m",maxjob=100,is_debug = self.is_debug )


    def run_macs_rep_broad(self, pvalue=0.05, ref="hs"):
        sh_file       = "%s/s08.macs2BroadPeakRep.sh"      % (self.scripts)
        sh_work_file  = "%s/s08.macs2BroadPeakRep_work.sh" % (self.scripts)
        
        l_sh_info    = self.s08_macs2BroadPeakRep( ref )
        l_sh_work    = []
    
        for brief_name in self.samInfo_pd_ChIP['brief_name']:
            idx   =(self.samInfo_pd_ChIP['brief_name'] == brief_name)
            if self.__is_input(idx):
                continue

            m01.make_dir([self.dir_BroadPeak_rep, brief_name])
            ctrl_name  = self.samInfo_pd_ChIP[idx]['control'].values[0]

            l_sh_work.append("sh %s  %s %s %f"                              %\
                ( sh_file, brief_name, ctrl_name, pvalue ))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8,is_debug=self.is_debug)
#        my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)



    def run_macs_mrg_broad(self, pvalue=0.05, ref="hs"):
        sh_file       = "%s/s09.macs2BroadPeakMrg.sh"      % (self.scripts)
        sh_work_file  = "%s/s09.macs2BroadPeakMrg_work.sh" % (self.scripts)
        
        l_sh_info = self.s09_macs2BroadPeakMrg(ref)
        l_sh_work = []
    
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx = (self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue

            m01.make_dir([self.dir_Peak_mrg, merge_name])
            ctrl_name = self.samInfo_pd_ChIP[idx]['control'].values[0]

            l_sh_work.append("sh %s  %s %s %f" % \
                (sh_file, merge_name, ctrl_name, pvalue))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
            
        my_job.running_multi(cpu=8, is_debug=self.is_debug)
#        my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)
        
    def sort_bdg(self):
        sh_file                                                             =\
            "%s/s10.sortbdg.sh"                     % (self.scripts)
            
        sh_work_repPeak_file                                                =\
            "%s/s10.1.sortbdg.repPeak_work.sh"      % (self.scripts)
            
        sh_work_mrgPeak_file                                                =\
            "%s/s10.2.sortbdg.mrgPeak_work.sh"      % (self.scripts)
            
        sh_work_repBroadPeak_file                                           =\
            "%s/s10.3.sortbdg.repBroadPeak_work.sh" % (self.scripts)
            
        sh_work_mrgBroadPeak_file                                           =\
            "%s/s10.4.sortbdg.mrgBroadPeak_work.sh" % (self.scripts)
        
        l_sh_info = self.s10_sortbdg()
        l_sh_work_repPeak      = []
        l_sh_work_mrgPeak      = []
        l_sh_work_repBroadPeak = []
        l_sh_work_mrgBroadPeak = []
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx = (self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue

            l_brief = list(self.samInfo_pd_ChIP[idx]['brief_name'])
            
            for brief in l_brief:
                l_sh_work_repBroadPeak.append("sh %s  %s %s"                %\
                    (sh_file, brief, self.dir_BroadPeak_rep))
         
            l_sh_work_mrgBroadPeak.append("sh %s  %s %s"                    %\
                 (sh_file, merge_name, self.dir_BroadPeak_mrg))

            for brief in l_brief:
                l_sh_work_repPeak.append("sh %s  %s %s"                     %\
                    (sh_file, brief, self.dir_Peak_rep))
         
            l_sh_work_mrgPeak.append("sh %s  %s %s"                         %\
                (sh_file, merge_name, self.dir_Peak_mrg))


        my_job_repPeak = m_jobs.run_jobs(sh_file,                            \
            sh_work_repPeak_file, l_sh_info, l_sh_work_repPeak)

        my_job_mrgPeak = m_jobs.run_jobs(sh_file,                            \
            sh_work_mrgPeak_file, l_sh_info, l_sh_work_mrgPeak)

        my_job_repBroadPeak = m_jobs.run_jobs(sh_file,                       \
            sh_work_repBroadPeak_file, l_sh_info, l_sh_work_repBroadPeak)

        my_job_mrgBroadPeak = m_jobs.run_jobs(sh_file,                       \
            sh_work_mrgBroadPeak_file, l_sh_info, l_sh_work_mrgBroadPeak)
        
        
        my_job_repPeak.running_multi(cpu=8, is_debug=self.is_debug)
        my_job_mrgPeak.running_multi(cpu=8, is_debug=self.is_debug)
        my_job_repBroadPeak.running_multi(cpu=8, is_debug=self.is_debug)
        my_job_mrgBroadPeak.running_multi(cpu=8, is_debug=self.is_debug)
                
#        my_job_repBroadPeak.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)
#        my_job_mrgPeak.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)
#        my_job_repPeak.running_SGE(vf="200m",maxjob=100,is_debug=self.is_debug)
#        my_job_mrgBroadPeak.running_SGE( vf="200m",maxjob=100,is_debug = self.is_debug )


    def make_igv_IDR(self):
        sh_file = "%s/s11.makeIGV.sh" % (self.scripts)
        sh_work_repPeak_file = "%s/s11.1.makeIGV.repPeak_work.sh"           %\
            (self.scripts)
        sh_work_mrgPeak_file = "%s/s11.2.makeIGV.mrgPeak_work.sh"           %\
            (self.scripts)
        
        l_sh_info = self.s11_makeIGV()
        l_sh_work_rep = []
        l_sh_work_mrg = []
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx = (self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue
                
            m01.make_dir([self.dir_Peak_TDF, merge_name])
            l_brief = list(self.samInfo_pd_ChIP[ idx ]['brief_name'])
            for brief in l_brief:
                l_sh_work_rep.append("sh %s  %s %s %s"                      %\
                    (sh_file, brief, merge_name, self.dir_Peak_rep))
            
            l_sh_work_mrg.append("sh %s  %s %s %s"                          %\
                (sh_file, merge_name, merge_name, self.dir_Peak_mrg))

        my_job_rep = m_jobs.run_jobs(sh_file, sh_work_repPeak_file,          \
            l_sh_info, l_sh_work_rep)
        my_job_mrg = m_jobs.run_jobs(sh_file, sh_work_mrgPeak_file,          \
            l_sh_info, l_sh_work_mrg)

        my_job_rep.running_multi(cpu=8, is_debug=self.is_debug)
        my_job_mrg.running_multi(cpu=8, is_debug=self.is_debug)
        
#       my_job_rep.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)
#       my_job_mrg.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)


    def make_igv_broad(self):
        sh_file              = "%s/s11.makeIGV_broad.sh"                    %\
            (self.scripts)
        sh_work_repPeak_file = "%s/s11.3.makeIGV.repBroadPeak_work.sh"      %\
            (self.scripts)
        sh_work_mrgPeak_file = "%s/s11.4.makeIGV.mrgBroadPeak_work.sh"      %\
            (self.scripts)
        
        l_sh_info = self.s11_makeIGV_broad()
        l_sh_work_rep = []
        l_sh_work_mrg = []
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx     = (self.samInfo_pd_ChIP['merge_name'] == merge_name)
            if self.__is_input(idx):
                continue
            
            m01.make_dir([self.dir_BroadPeak_TDF, merge_name])            
            l_brief = list(self.samInfo_pd_ChIP[idx]['brief_name'])
            for brief in l_brief:
                l_sh_work_rep.append("sh %s   %s %s %s"                     %\
                     (sh_file, brief, merge_name, self.dir_BroadPeak_rep))

            l_sh_work_mrg.append("sh %s  %s %s %s"                          %\
                 (sh_file, merge_name, merge_name, self.dir_BroadPeak_mrg))

        my_job_rep = m_jobs.run_jobs(sh_file, sh_work_repPeak_file,          \
            l_sh_info, l_sh_work_rep)

        my_job_mrg = m_jobs.run_jobs(sh_file, sh_work_mrgPeak_file,          \
            l_sh_info, l_sh_work_mrg)
        
        my_job_rep.running_multi(cpu=8, is_debug=self.is_debug)
        my_job_mrg.running_multi(cpu=8, is_debug=self.is_debug)
        
#       my_job_rep.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)
#       my_job_mrg.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)

    def __is_input(self,idx):
        ltype = self.samInfo_pd_ChIP[idx]['type'].values[0]
        is_input = 0
        if ltype.upper() == "input".upper():
            is_input = 1

        return is_input

