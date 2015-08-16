#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import numpy   as np
import cPickle as pickle
import pandas  as pd

class LoadSamp(object):
    def __init__(self):
        super( LoadSamp,self ).__init__()
        
    def load_ChIP_samInfo(self,ChIP_infile):
        self.samInfo_ChIP_infile = ChIP_infile
        self.samInfo_pd_ChIP = pd.read_csv(self.samInfo_ChIP_infile, sep="\t")
        
    def load_RNA_samInfo(self,RNA_infile):
        self.samInfo_RNA_infile  = sam_file_RNA
        self.samInfo_pd_RNA  = pd.read_csv(self.samInfo_RNA_infile , sep="\t")
        
      
class DirSystem(object):
    def __init__(self):
        super( DirSystem,self ).__init__()
        self.home_dir               = os.path.abspath('./')
        path_file                   = os.path.realpath(__file__)
        self.dir_raw_data           = "%s/00.0.raw_data"                    % (self.home_dir)
        self.dir_clean_data         = "%s/00.1.clean_data"                  % (self.home_dir)
        self.dir_bam                = "%s/01.bam"                           % (self.home_dir)
        self.dir_bed_rep            = "%s/01.1.ext300_bed_rep"              % (self.home_dir)
        self.dir_bed_mrg            = "%s/01.2.ext300_bed_mrg"              % (self.home_dir)
        self.dir_spp_rep_shiftSize  = "%s/01.3.ext300_bed_rep_shiftSize"    % (self.home_dir)
        self.dir_spp_mrg_shiftSize  = "%s/01.4.ext300_bed_mrg_shiftSize"    % (self.home_dir)
        self.dir_RPM_mrg            = "%s/02.RPM_Merged"                    % (self.home_dir)
        self.dir_RPM_bins_rep       = "%s/02.1.RPM_bins_rep"                % (self.home_dir)
        self.dir_RPM_bins_mrg       = "%s/02.2.RPM_bins_mrg"                % (self.home_dir)
        self.dir_Peak_rep           = "%s/03.1.Peak_rep"                    % (self.home_dir)
        self.dir_Peak_mrg           = "%s/03.2.Peak_mrg"                    % (self.home_dir)
        self.dir_Peak_idr           = "%s/03.3.Peak_idr"                    % (self.home_dir)
        self.dir_Peak_TDF           = "%s/03.4.Peak_tdf"                    % (self.home_dir)
        self.dir_BroadPeak_TDF      = "%s/03.5.BroadPeak_tdf"               % (self.home_dir)
        self.dir_BroadPeak_rep      = "%s/04.1.BroadPeak_rep"               % (self.home_dir)
        self.dir_BroadPeak_mrg      = "%s/04.2.BroadPeak_mrg"               % (self.home_dir)
        self.dir_Peak_mrg_TSS       = "%s/05.1.Peak_mrg.density.TSS"        % (self.home_dir)
        self.dir_Peak_mrg_Gene      = "%s/05.2.Peak_mrg.density.Genebody"   % (self.home_dir)
        self.dir_Peak_mrg_TSS_Cor   = "%s/06.Peak_mrg.TSS.RPKM_Cor"         % (self.home_dir)
        self.dir_StatInfo           = "%s/StatInfo"                         % (self.home_dir)
        self.Database               = "%s/Database"                         % (self.home_dir)
        self.path                   = "%s"                                  % ("/".join(path_file.split('/')[:-2]))
        self.bin                    = "%s/bin"                              % ("/".join(path_file.split('/')[:-2]))

        
        if not os.path.exists( self.dir_StatInfo ):
            os.mkdir( self.dir_StatInfo )
            
        #####################################
        ######## Revise this path!!! ########
        #####################################
        self.Database               = "/data/Analysis/huboqiang/Database_ChIP_v2"
        
    def define_scripts(self, s_idx):
        dir_script = "%s/scripts"      % (self.home_dir)
        self.scripts = "%s/scripts/%s" % (self.home_dir, s_idx)

        if not os.path.exists(dir_script):
            os.mkdir(dir_script)
        
        if not os.path.exists(self.scripts):
            os.mkdir(self.scripts)


class UsedSoftware(object):
    def __init__(self):
        super( UsedSoftware,self ).__init__()
        self.sftw_py             = "/data/Analysis/huboqiang/software/anaconda/bin/python"
        self.sftw_pl             = "/data/Analysis/huboqiang/lib/local_perl/bin/perl"
        self.sftw_java           = "/usr/bin/java"
        self.sftw_R              = "/data/Analysis/huboqiang/bin/Rscript"
        self.sftw_MarkDup        = "/data/Analysis/huboqiang/software/picard-tools-1.119/MarkDuplicates.jar"
        self.sftw_bwa            = "/data/Analysis/huboqiang/bin/bwa"
        self.sftw_samtools       = "/data/Analysis/huboqiang/software/samtools-0.1.18/samtools"    
        self.sftw_macs2          = "/data/Analysis/huboqiang/software/MACS/bin/macs2"
        self.sftw_bedtools       = "/data/Analysis/huboqiang/software/bedtools-2.17.0/bin/bedtools"
        self.sftw_bgzip          = "/data/Analysis/huboqiang/software/tabix/bgzip"
        self.sftw_tabix          = "/data/Analysis/huboqiang/software/tabix/tabix"
        self.sftw_igvtools       = "/data/Analysis/huboqiang/bin/igvtools"
        self.sftw_batchIDR       = "/data/Analysis/huboqiang/software/study_broadMacsScripts/idrCode/batch-consistency-analysis.r"
        self.sftw_spp            = "/data/Analysis/huboqiang/software/spp_rev/run_spp.R"
        self.sftw_get_psudoCount = "/data/Analysis/huboqiang/project/human_sc/chIp_analysis/bin/run_psudoCount.sh"
        self.sftw_sort_bdg       = "/data/Analysis/huboqiang/project/human_sc/chIp_analysis/bin/sort_bdg.sh"
        self.sftw_ucsc_dir       = "/data/Analysis/huboqiang/software/UCSC"


class ProjInfo(LoadSamp,DirSystem,UsedSoftware):
    def __init__(self):
        super( ProjInfo,self ).__init__()

