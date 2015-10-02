from __future__ import division
import re,sys,os,gzip,string
import subprocess
import time
import cPickle as pickle
import numpy   as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import ChIP.settings.projpath         as m_proj
import ChIP.utils.module_StatFunc          as m_Stat

def get_PassPeak(infile, FDR):
    cnt = 0
    if os.path.isfile(infile):
        f_infile = open(infile,"r")
        f_infile.readline()
        for line in f_infile:
            line = line.strip('\n')
            f    = line.split()
            if float(f[10]) <= FDR:
                cnt += 1
        f_infile.close()

    return cnt

class StatInfo(m_proj.ProjInfo):
    def __init__(self, sam_ChIPinfo):
        self.sam_ChIPinfo = sam_ChIPinfo
        super(StatInfo,self).__init__()
        self.load_ChIP_samInfo(self.sam_ChIPinfo)
        self.StatInfo = {}
        
    def Basic_Stat(self):
        '''
        Stat for QC, bwa mapping and peak calling.
        '''
        prefix = ".".join( self.sam_ChIPinfo.split(".")[:-1] )
        f_out  = open("%s/01.Basic_info.%s.xls" % (self.dir_StatInfo, prefix),"w")
        
        out  = "Merge_name\tSamples_inclueds\tRaw_reads\tClean_reads\t"
        out += "Duplicate_reads\tUniqueMapped_reads\tMultipleMapped_Reads\t"
        out += "Q30Mapped_Reads"
        print >>f_out, out
        
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            if merge_name not in self.StatInfo:
                self.StatInfo[ merge_name ] = {}

            l_totRaw = []
            l_totCln = []
            l_totDup = []
            l_totUnique = []
            l_totMulti  = []
            l_totQ30    = []

            idx       =(self.samInfo_pd_ChIP['merge_name'] == merge_name)
            l_sample  = self.samInfo_pd_ChIP[ idx ]['sample'] 
            l_brief   = list(self.samInfo_pd_ChIP[ idx ]['brief_name'] )
            for i,sam in enumerate(l_sample):
                QC_file = "%s/%s/log" % (self.dir_clean_data,sam)
                QC_info = m_Stat.QcStat(QC_file)
                QC_info.read_infile()
                l_totRaw.append(QC_info.raw_reads)
                l_totCln.append(QC_info.cln_reads)

                brief_name = l_brief[i]
                map_dup_info = m_Stat.BwaPicard()
                rmdup_file = "%s/%s/%s.picard_info.txt"                     %\
                    (self.dir_bam,     brief_name, brief_name)
                bed_unique = "%s/%s.unique.bed"                             %\
                    (self.dir_bed_rep, brief_name)
                bed_multi  = "%s/%s.multi.bed"                              %\
                    (self.dir_bed_rep, brief_name)
                bed_q30    = "%s/%s.sort.tagAlign.gz"                       %\
                    (self.dir_bed_rep, brief_name)

                map_dup_info.read_picard_file(rmdup_file)
                map_dup_info.read_unique_bed(bed_unique)
                map_dup_info.read_multi_bed(bed_multi)
                map_dup_info.read_q30_bed(bed_q30)

                l_totDup.append(map_dup_info.dup)
                l_totUnique.append(map_dup_info.unique)
                l_totMulti.append(map_dup_info.multi)
                l_totQ30.append(map_dup_info.q30)
            
            self.StatInfo[ merge_name ]['l_brief']= l_brief
            self.StatInfo[ merge_name ]['raw']    = np.array(l_totRaw)
            self.StatInfo[ merge_name ]['cln']    = np.array(l_totCln)
            self.StatInfo[ merge_name ]['dup']    = np.array(l_totDup)
            self.StatInfo[ merge_name ]['unique'] = np.array(l_totUnique)
            self.StatInfo[ merge_name ]['multi']  = np.array(l_totMulti)
            self.StatInfo[ merge_name ]['q30']    = np.array(l_totQ30)

            out = "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d"                          %\
                ( merge_name, 
                  ",".join(l_brief),  
                  np.sum(self.StatInfo[ merge_name ]['raw']),
                  np.sum(self.StatInfo[ merge_name ]['cln']),
                  np.sum(self.StatInfo[ merge_name ]['dup']),
                  np.sum(self.StatInfo[ merge_name ]['unique']),
                  np.sum(self.StatInfo[ merge_name ]['multi']),
                  np.sum(self.StatInfo[ merge_name ]['q30']) )
            print >>f_out, out

        f_out.close()

    def IDR_Stat(self):
        prefix      = ".".join(self.sam_ChIPinfo.split(".")[:-1])
        file_idr_out= "%s/IDR_result.%s.xls" % (self.dir_StatInfo, prefix)
        f_idr_out   = open(file_idr_out,"w")

        out  = "Merge_sam\tPeak\trepCmp_peak\tpseudoRep_peak\t"
        out += "Range_rep\tRange_cmp\trepCmp_sam\tpseudoRep_sam"
        print >>f_idr_out, out
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            print merge_name
            idx   =(self.samInfo_pd_ChIP['merge_name'] == merge_name)
            ltype = self.samInfo_pd_ChIP[ idx ]['type'].values[0]
            if ltype.upper() == "input".upper():
                continue
            
            idx       =(self.samInfo_pd_ChIP['merge_name'] == merge_name)
            l_sam     = self.samInfo_pd_ChIP[ idx ]['sample'] 
            l_brief   = list(self.samInfo_pd_ChIP[ idx ]['brief_name'])

            l_repCmp         = []
            l_repCmpPeak     = []
            l_selfPseCmp     = []
            l_selfPseCmpPeak = []

            if len(l_brief) > 1:
                ### IDR RESULTS ON ORIGINAL REPLICATES         
                for i    in range( 0,len(l_brief)-1 ):
                    for j in range( i+1,len(l_brief) ):
                        sam1 = l_brief[i]
                        sam2 = l_brief[j]
                        in_sam1 = "%s/%s/%s_VS_Input_peaks.regionPeak.gz"   %\
                            (self.dir_Peak_idr, merge_name, sam1)
                            
                        in_sam2 = "%s/%s/%s_VS_Input_peaks.regionPeak.gz"   %\
                            (self.dir_Peak_idr, merge_name, sam2)
                            
                        rep_cmp = "%s/%s/rep/%s_VS_%s-overlapped-peaks.txt" %\
                            (self.dir_Peak_idr, merge_name, sam1, sam2)
                            
                        peak_cnt= get_PassPeak(rep_cmp, 0.03)
                        l_repCmp.append("%s_%s" % (sam1,sam2))
                        l_repCmpPeak.append(peak_cnt)

            ### IDR RESULTS ON SELF-PSEUDOREPLICATES
            for sam in l_brief:
                in_sam_pr1 = "%s/%s/%s.pr1_VS_Input_peaks.regionPeak.gz"%\
                    (self.dir_Peak_idr, merge_name, sam)
                    
                in_sam_pr2 = "%s/%s/%s.pr2_VS_Input_peaks.regionPeak.gz"%\
                    (self.dir_Peak_idr, merge_name, sam)
                    
                selfPse_cmp = "%s/%s/selfPseudoReps/"                   %\
                    (self.dir_Peak_idr, merge_name)
                selfPse_cmp+= "%s.pr1_VS_%s.pr2-overlapped-peaks.txt"   %\
                    (sam, sam)
                    
                peak_cnt= get_PassPeak(selfPse_cmp, 0.02 )
                l_selfPseCmp.append(sam)
                l_selfPseCmpPeak.append(peak_cnt)


            ### IDR RESULTS ON POOLED-PSEUDOREPLICATES
            in_mrg_pr1 = "%s/%s/%s.pr1_VS_Input_peaks.regionPeak.gz"        %\
                (self.dir_Peak_idr, merge_name, merge_name)
                
            in_mrg_pr2 = "%s/%s/%s.pr2_VS_Input_peaks.regionPeak.gz"        %\
                (self.dir_Peak_idr, merge_name, merge_name)
                
            poolPse_cmp = "%s/%s/pooledPseudoReps/%s.pr1_VS_%s.pr2"         %\
                (self.dir_Peak_idr, merge_name, merge_name, merge_name)
            poolPse_cmp+= "-overlapped-peaks.txt"

            PeakCnt    = get_PassPeak(poolPse_cmp, 0.01)
            repCmp         = "-"
            repCmpPeak     = "0"
            selfPseCmp     = ",".join(l_selfPseCmp)
            selfPseCmpPeak = ",".join(["%s" % i for i in l_selfPseCmpPeak])
            range_rep      = np.nan
            range_selfPse  = np.nan

            if len( l_brief ) > 1:         
                PeakCnt       = max(l_repCmpPeak)
                repCmp        = ",".join(l_repCmp)
                repCmpPeak    = ",".join(["%s" % i for i in l_repCmpPeak])
                selfPseCmp    = ",".join(l_selfPseCmp)
                selfPseCmpPeak= ",".join(["%s" % i for i in l_selfPseCmpPeak])
                div_rep       = (min(l_repCmpPeak    )+0.001)
                div_selfPse   = (min(l_selfPseCmpPeak)+0.001)
                range_rep     = max(l_repCmpPeak    )/div_rep
                range_selfPse = max(l_selfPseCmpPeak)/div_selfPse

            print >>f_idr_out,  "%s\t%d\t%s\t%s\t%f\t%f\t%s\t%s"            %\
                (merge_name, PeakCnt, repCmpPeak, selfPseCmpPeak, range_rep, \
                 range_selfPse,  repCmp,selfPseCmp )

        f_idr_out.close()
