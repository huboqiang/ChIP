#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import numpy   as np
from scipy import stats
import scipy   as sp
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd

import ChIP.utils.module_running_jobs as m_jobs
import ChIP.settings.projpath         as m_proj
import ChIP.settings.scripts          as m_scpt


def make_dir( l_args ):
    """docstring for make_dir"""
    if len(l_args) == 1:
        """ only dictionary """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )

    elif len(l_args) == 2:
        """ dictionary/sample """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
        if not os.path.isdir( "%s/%s" % (l_args[0],l_args[1]) ):
            os.mkdir(         "%s/%s" % (l_args[0],l_args[1]) )

    elif len(l_args) == 3:
        """ dictionary/sample/subname """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
        if not os.path.isdir( "%s/%s"    % (l_args[0],l_args[1]) ):
            os.mkdir(         "%s/%s"    % (l_args[0],l_args[1]) )
        if not os.path.isdir( "%s/%s/%s" % (l_args[0],l_args[1],l_args[2]) ):
            os.mkdir(         "%s/%s/%s" % (l_args[0],l_args[1],l_args[2]) )


class Map_From_raw( m_scpt.Scripts):

    def __init__(self, sam_ChIPinfo, ref, is_debug=1 ):
        super( Map_From_raw,self ).__init__()

        s_idx = ".".join(sam_ChIPinfo.split("/")[-1].split(".")[:-1])

        self.load_ChIP_samInfo( sam_ChIPinfo )

        self.define_scripts(s_idx)
        self.define_files(ref)

        self.M_CvtEnd = {"PE":2,"SE":1}
        self.is_debug = is_debug

    def run_QC(self, core_num=4):
        sh_file       = "%s/s01.QC.sh"              % (self.scripts)
        sh_work_file  = "%s/s01.QC_work.sh"         % (self.scripts)

        l_sh_info = self.s01_QC()
        l_sh_work = []
        for samp in self.samInfo_pd_ChIP['sample']:
            make_dir(  [ self.dir_clean_data, samp ] )
            idx       = (self.samInfo_pd_ChIP['sample'] == samp)
            end       = self.samInfo_pd_ChIP[ idx ]['end_type'].values[0]
            data_dype = self.M_CvtEnd[ end ]
            l_sh_work.append("sh %s %s %d" % (sh_file, samp, data_dype) )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)


    def run_bwa(self, core_num=2):
        sh_file      = "%s/s02.bwa.sh"      % (self.scripts)
        sh_work_file = "%s/s02.bwa_work.sh" % (self.scripts)

        l_sh_info = self.s02_bwa()
        l_sh_work = []
        for samp in self.samInfo_pd_ChIP['sample']:
            idx        =(self.samInfo_pd_ChIP['sample'] == samp)
            brief_name = self.samInfo_pd_ChIP[ idx ]['brief_name'].values[0]
            end       = self.samInfo_pd_ChIP[ idx ]['end_type'].values[0]
            data_dype = self.M_CvtEnd[ end ]
            make_dir( [ self.dir_bam, brief_name ] )
            l_sh_work.append(
                "sh %s %s %s %d" % (sh_file, samp, brief_name, data_dype)
            )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=core_num, is_debug = self.is_debug)
#       my_job.running_SGE( vf="7g",maxjob=100,is_debug = self.is_debug )


    def bam2repbed(self, ext_len=300, core_num=4):
        sh_file      = "%s/s03.bam2bedrep.sh"      % (self.scripts)
        sh_work_file = "%s/s03.bam2bedrep_work.sh" % (self.scripts)

        l_sh_info = self.s03_bam2bedrep()
        l_sh_work = []

        make_dir( [ self.dir_bed_rep ] )
        for samp in self.samInfo_pd_ChIP['sample']:
            idx        =(self.samInfo_pd_ChIP['sample'] == samp)
            brief_name = self.samInfo_pd_ChIP[ idx ]['brief_name'].values[0]
            end       = self.samInfo_pd_ChIP[ idx ]['end_type'].values[0]
            data_dype = self.M_CvtEnd[ end ]

            l_sh_work.append("sh %s  %s %s  %d" % \
                (sh_file, brief_name, data_dype, ext_len))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)

    def mrgbed(self, core_num=4):
        sh_file      = "%s/s04.1.bedmrg.sh"      % (self.scripts)
        sh_work_file = "%s/s04.1.bedmrg_work.sh" % (self.scripts)

        l_sh_info = self.s04_bedmrg()
        l_sh_work = []

        make_dir( [ self.dir_bed_mrg ] )
        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx       =(self.samInfo_pd_ChIP['merge_name'] == merge_name)

            merged    = "%s/%s" % (self.dir_bed_mrg, merge_name)
            l_brief   = self.samInfo_pd_ChIP[ idx ]['brief_name']
            l_bed_rep = [ "%s/%s" % (self.dir_bed_rep,brief_name )           \
                 for brief_name in l_brief ]

            l_sh_work.append(
                "sh %s %s %s" % (sh_file, merged, " ".join(l_bed_rep))
            )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)

    def mrgbed_multi(self, core_num=4):
        sh_file      = "%s/s04.2.bedmrg.multi.sh"      % (self.scripts)
        sh_work_file = "%s/s04.2.bedmrg.multi_work.sh" % (self.scripts)

        l_sh_info = self.s04_2_bedmrg_multi()
        l_sh_work = []

        for merge_name in set(self.samInfo_pd_ChIP['merge_name']):
            idx       =(self.samInfo_pd_ChIP['merge_name'] == merge_name)

            merged    = "%s/%s" % (self.dir_bed_mrg, merge_name)
            l_brief   = self.samInfo_pd_ChIP[ idx ]['brief_name']
            l_bed_rep = [ "%s/%s"        % (self.dir_bed_rep,brief_name)     \
                for brief_name in l_brief ]

            l_sh_work.append(
                "sh %s  %s %s" % ( sh_file, merged, " ".join(l_bed_rep) )
            )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=core_num, is_debug = self.is_debug)
#       my_job.running_SGE(vf="200m",maxjob=100,is_debug = self.is_debug)
