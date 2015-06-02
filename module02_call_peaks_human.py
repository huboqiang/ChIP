#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import numpy   as np
import cPickle as pickle
import scipy   as sp
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats 
import module_running_jobs as m_jobs

def get_shiftSize(infile):
   f_infile = open( infile,"r" )
   l_line   = f_infile.readlines()
   line     = l_line[0]
   line     = line.strip()
   f        = line.split()
   f_infile.close()
   l_val = np.array( f[2].split(','),dtype=int )
   return int( np.mean(l_val)/2 )

def get_PassPeak(infile,FDR):
   print infile
   f_infile = open( infile,"r" )
   cnt = 0
   f_infile.readline()
   for line in f_infile:
      line = line.strip('\n')
      f    = line.split()
      if float(f[10]) <= FDR:
         cnt += 1
   f_infile.close()
   return cnt

class Macs2Peaks(dict):
   
   def __init__(self,samp_info, genome_file,anno_file, dir_name, sftw_name ):
      self['sam_info'] = samp_info['samp']['chip']
      self['infile']   = { 'genome_file':genome_file,'anno_file':anno_file }
      self['stage']    = { 'name':[] }
      self['dir_name']  = dir_name
      self['sftw_name'] = sftw_name
      
   def get_shift_size_rep(self):
      home_dir    = os.path.abspath('./')
      bed_rep_dir = self['dir_name'].bed_rep
      spp_rep_dir = self['dir_name'].spp_rep_shiftSize
      
      if not os.path.isdir( spp_rep_dir ):
         os.mkdir( spp_rep_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      R_exe        = self['sftw_name'].R
      spp_exe      = self['sftw_name'].spp
      sh_file       = "%s/scripts/s05.1.spp_rep_shiftSize.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s05.1.spp_rep_shiftSize_work.sh" % (home_dir)
      
      sh_info = """
R_exe=$1
spp_exe=$2
bed_rep_dir=$3
spp_rep_dir=$4
samp=$5

rm $spp_rep_dir/$samp/out/out.tab
$R_exe $spp_exe -c=$bed_rep_dir/$samp.sort.tagAlign.gz -odir=$spp_rep_dir/$samp/test  -savp -rf -out=$spp_rep_dir/$samp/out/out.tab
      """
      
      sh_work = ""
      for samp in self['sam_info']['list']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         if not os.path.isdir( "%s/%s"      % ( spp_rep_dir,brief_name ) ):
            os.mkdir(          "%s/%s"      % ( spp_rep_dir,brief_name ) )
            os.mkdir(          "%s/%s/test" % ( spp_rep_dir,brief_name ) )
            os.mkdir(          "%s/%s/out"  % ( spp_rep_dir,brief_name ) )
                     
         sh_work  += " sh %s   %s %s %s %s %s\n" % ( sh_file, R_exe,  spp_exe, bed_rep_dir, spp_rep_dir, brief_name )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )

   def get_shift_size_mrg(self):
      home_dir    = os.path.abspath('./')
      bed_mrg_dir = self['dir_name'].bed_mrg
      spp_mrg_dir = self['dir_name'].spp_mrg_shiftSize
      
      if not os.path.isdir( spp_mrg_dir ):
         os.mkdir( spp_mrg_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      R_exe        = self['sftw_name'].R
      spp_exe      = self['sftw_name'].spp
      sh_file       = "%s/scripts/s05.2.spp_mrg_shiftSize.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s05.2.spp_mrg_shiftSize_work.sh" % (home_dir)
      
      sh_info = """
R_exe=$1
spp_exe=$2
bed_mrg_dir=$3
spp_mrg_dir=$4
samp=$5

rm $spp_mrg_dir/$samp/out/out.tab
$R_exe $spp_exe -c=$bed_mrg_dir/$samp.sort.tagAlign.gz -odir=$spp_mrg_dir/$samp/test  -savp -rf -out=$spp_mrg_dir/$samp/out/out.tab
      """
      
      sh_work = ""
      for merge_name in self['sam_info']['merge']:
         
         if not os.path.isdir( "%s/%s"      % ( spp_mrg_dir,merge_name ) ):
            os.mkdir(          "%s/%s"      % ( spp_mrg_dir,merge_name ) )
            os.mkdir(          "%s/%s/test" % ( spp_mrg_dir,merge_name ) )
            os.mkdir(          "%s/%s/out"  % ( spp_mrg_dir,merge_name ) )
                     
         sh_work  += " sh %s  %s %s %s %s %s\n" % ( sh_file, R_exe,  spp_exe, bed_mrg_dir, spp_mrg_dir, merge_name )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
   
   

   def run_macs_rep(self, pvalue=0.001):
      home_dir    = os.path.abspath('./')
      bed_rep_dir = self['dir_name'].bed_rep
      bed_mrg_dir = self['dir_name'].bed_mrg
      macs_rep_dir= self['dir_name'].Peak_rep
      spp_rep_dir = self['dir_name'].spp_rep_shiftSize
      
      if not os.path.isdir( macs_rep_dir ):
         os.mkdir( macs_rep_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      py_exe       = self['sftw_name'].py
      macs2_exe    = self['sftw_name'].macs2
      get_psudoCount=self['sftw_name'].get_psudoCount
      sh_file       = "%s/scripts/s06.1.macs2PeakRep.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s06.1.macs2PeakRep_work.sh" % (home_dir)
      
      sh_info = """
py_exe=$1
macs2_exe=$2
treat_bed=$3
ctrl_bed=$4
name_prefix=$5
pvalue=$6
shift_size=$7
get_psudoCount=$8

$py_exe $macs2_exe callpeak -t $treat_bed                 -c $ctrl_bed -f BED -n ${name_prefix}_VS_Input     -g hs -p $pvalue --to-large --nomodel --extsize $shift_size -B
sh $get_psudoCount $treat_bed
$py_exe $macs2_exe callpeak -t $treat_bed.pr1.tagAlign.gz -c $ctrl_bed -f BED -n ${name_prefix}.pr1_VS_Input -g hs -p $pvalue --to-large --nomodel --extsize $shift_size -B
$py_exe $macs2_exe callpeak -t $treat_bed.pr2.tagAlign.gz -c $ctrl_bed -f BED -n ${name_prefix}.pr2_VS_Input -g hs -p $pvalue --to-large --nomodel --extsize $shift_size -B
      """
      
      sh_work = ""
      for samp in self['sam_info']['list']:
         brief_name = self['sam_info']['samp_brief'][samp]
         ctrl_name  = self['sam_info']['sam_ctrl'][samp]
         if not os.path.isdir( "%s/%s" % (macs_rep_dir,brief_name) ):
            os.mkdir(          "%s/%s" % (macs_rep_dir,brief_name) )

         name_prefix = "%s/%s/%s"                 % (macs_rep_dir,brief_name,brief_name)
         treat_bed   = "%s/%s.sort.tagAlign.gz"   % (bed_rep_dir ,brief_name)
         ctrl_bed    = "%s/%s.sort.tagAlign.gz"   % (bed_mrg_dir ,ctrl_name)
         
         shift_size = 300
         f_shiftSize = "%s/%s/out/out.tab"        % (spp_rep_dir, brief_name )
         if os.path.isfile( f_shiftSize ):
            shift_size  = get_shiftSize( f_shiftSize )
         
         sh_work  += " sh %s  %s %s %s %s  %s %f %d %s \n" % ( sh_file, py_exe,  macs2_exe, treat_bed, ctrl_bed, name_prefix, pvalue,shift_size, get_psudoCount )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )
      

   def run_macs_mrg(self, pvalue=0.001):
      home_dir    = os.path.abspath('./')
      bed_mrg_dir = self['dir_name'].bed_mrg
      macs_mrg_dir= self['dir_name'].Peak_mrg
      spp_mrg_dir = self['dir_name'].spp_mrg_shiftSize
      
      if not os.path.isdir( macs_mrg_dir ):
         os.mkdir( macs_mrg_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      py_exe       = self['sftw_name'].py
      macs2_exe    = self['sftw_name'].macs2
      get_psudoCount=self['sftw_name'].get_psudoCount

      sh_file       = "%s/scripts/s06.2.macs2PeakMrg.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s06.2.macs2PeakMrg_work.sh" % (home_dir)
      
      sh_info = """
py_exe=$1
macs2_exe=$2
treat_bed=$3
ctrl_bed=$4
name_prefix=$5
pvalue=$6
shift_size=$7
get_psudoCount=$8

$py_exe $macs2_exe callpeak -t $treat_bed                 -c $ctrl_bed -f BED -n ${name_prefix}_VS_Input     -g hs -p $pvalue --to-large --nomodel --extsize $shift_size -B
sh $get_psudoCount $treat_bed
$py_exe $macs2_exe callpeak -t $treat_bed.pr1.tagAlign.gz -c $ctrl_bed -f BED -n ${name_prefix}.pr1_VS_Input -g hs -p $pvalue --to-large --nomodel --extsize $shift_size -B
$py_exe $macs2_exe callpeak -t $treat_bed.pr2.tagAlign.gz -c $ctrl_bed -f BED -n ${name_prefix}.pr2_VS_Input -g hs -p $pvalue --to-large --nomodel --extsize $shift_size -B
      """
      
      sh_work = ""
      for merge_name in self['sam_info']['merge']:
         ctrl_name  = self['sam_info']['merge_ctrl'][merge_name]
         if not os.path.isdir( "%s/%s" % (macs_mrg_dir,merge_name) ):
            os.mkdir(          "%s/%s" % (macs_mrg_dir,merge_name) )
            
         name_prefix = "%s/%s/%s"                 % (macs_mrg_dir,merge_name,merge_name)
         treat_bed   = "%s/%s.sort.tagAlign.gz"   % (bed_mrg_dir ,merge_name)
         ctrl_bed    = "%s/%s.sort.tagAlign.gz"   % (bed_mrg_dir ,ctrl_name)
         
         shift_size = 300
         f_shiftSize    = "%s/%s/out/out.tab"        % (spp_mrg_dir, merge_name )
         if os.path.isfile( f_shiftSize ):
            shift_size  = get_shiftSize( f_shiftSize )
         
         sh_work  += " sh %s  %s %s %s %s  %s %f %d %s \n" % ( sh_file, py_exe,  macs2_exe, treat_bed, ctrl_bed, name_prefix, pvalue,shift_size,  get_psudoCount )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )
   

   def prepare_idr_input(self,top_peak=100000):
      home_dir    = os.path.abspath('./')
      Peak_rep_dir= self['dir_name'].Peak_rep
      Peak_mrg_dir= self['dir_name'].Peak_mrg
      Peak_idr_dir= self['dir_name'].Peak_idr
      
      if not os.path.isdir( Peak_idr_dir ):
         os.mkdir( Peak_idr_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      R_exe        = self['sftw_name'].R
      batchIDR_exe = self['sftw_name'].batchIDR
      
      bgzip_exe = self['sftw_name'].bgzip
      tabix_exe = self['sftw_name'].tabix

      sh_file       = "%s/scripts/s07.1.IDR_prepare.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s07.1.IDR_prepare_work.sh" % (home_dir)
      
      sh_info = """      
mrg_peak=$1
top_peak=$2

bgzip_exe=$3
tabix_exe=$4
Peak_rep_dir=$5
Peak_mrg_dir=$6
Peak_idr_dir=$7

shift
shift
shift
shift
shift
shift
shift

for i in $@
do 
   sort -k 8nr,8nr  $Peak_rep_dir/${i}/${i}_VS_Input_peaks.narrowPeak | head -n $top_peak              | \
   $bgzip_exe -cf > $Peak_rep_dir/${i}/${i}_VS_Input_peaks.regionPeak.gz                              && \
   sort -k 8nr,8nr  $Peak_rep_dir/${i}/${i}.pr1_VS_Input_peaks.narrowPeak | head -n $top_peak          | \
   $bgzip_exe -cf > $Peak_rep_dir/${i}/${i}.pr1_VS_Input_peaks.regionPeak.gz                          && \
   sort -k 8nr,8nr  $Peak_rep_dir/${i}/${i}.pr2_VS_Input_peaks.narrowPeak | head -n $top_peak          | \
   $bgzip_exe -cf > $Peak_rep_dir/${i}/${i}.pr2_VS_Input_peaks.regionPeak.gz                          && \
   mv    $Peak_rep_dir/${i}/${i}_VS_Input_peaks.regionPeak.gz*        $Peak_idr_dir/${mrg_peak}       && \
   mv    $Peak_rep_dir/${i}/${i}.pr1_VS_Input_peaks.regionPeak.gz*    $Peak_idr_dir/${mrg_peak}       && \
   mv    $Peak_rep_dir/${i}/${i}.pr2_VS_Input_peaks.regionPeak.gz*    $Peak_idr_dir/${mrg_peak}
done

sort -k 8nr,8nr  $Peak_mrg_dir/${mrg_peak}/${mrg_peak}_VS_Input_peaks.narrowPeak     | head -n $top_peak          | \
$bgzip_exe -cf > $Peak_mrg_dir/${mrg_peak}/${mrg_peak}_VS_Input_peaks.regionPeak.gz                              && \
sort -k 8nr,8nr  $Peak_mrg_dir/${mrg_peak}/${mrg_peak}.pr1_VS_Input_peaks.narrowPeak | head -n $top_peak          | \
$bgzip_exe -cf > $Peak_mrg_dir/${mrg_peak}/${mrg_peak}.pr1_VS_Input_peaks.regionPeak.gz                          && \
sort -k 8nr,8nr  $Peak_mrg_dir/${mrg_peak}/${mrg_peak}.pr2_VS_Input_peaks.narrowPeak | head -n $top_peak          | \
$bgzip_exe -cf > $Peak_mrg_dir/${mrg_peak}/${mrg_peak}.pr2_VS_Input_peaks.regionPeak.gz                          && \
mv $Peak_mrg_dir/${mrg_peak}/${mrg_peak}_VS_Input_peaks.regionPeak.gz      $Peak_idr_dir/${mrg_peak}             && \
mv $Peak_mrg_dir/${mrg_peak}/${mrg_peak}.pr1_VS_Input_peaks.regionPeak.gz  $Peak_idr_dir/${mrg_peak}             && \
mv $Peak_mrg_dir/${mrg_peak}/${mrg_peak}.pr2_VS_Input_peaks.regionPeak.gz  $Peak_idr_dir/${mrg_peak}
      """
      
      sh_work = ""
      for merge_name in self['sam_info']['merge']:
         ctrl_name  = self['sam_info']['merge_ctrl'][merge_name]
         if not os.path.isdir( "%s/%s" % (Peak_idr_dir,merge_name) ):
            os.mkdir(          "%s/%s" % (Peak_idr_dir,merge_name) )
            
         l_sam     = self['sam_info']['merge_sam'][merge_name]
         l_brief   = [ "%s" % (self['sam_info']['samp_brief'][sam]) for sam in l_sam ]
                  
         sh_work  += " sh %s %s %d   %s %s   %s %s %s   %s\n" % (  sh_file,merge_name,top_peak,   bgzip_exe,tabix_exe,   Peak_rep_dir,Peak_mrg_dir,Peak_idr_dir ," ".join(l_brief) )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )



   def run_idr(self,top_peak=100000):
      home_dir    = os.path.abspath('./')
      Peak_idr_dir= self['dir_name'].Peak_idr
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      R_exe        = self['sftw_name'].R
      batchIDR_exe = self['sftw_name'].batchIDR
      
      bgzip_exe = self['sftw_name'].bgzip
      tabix_exe = self['sftw_name'].tabix

      sh_file                        = "%s/scripts/s07.2.IDR_usingMacs2Peak.sh"                       % (home_dir)
      sh_work_file_rep               = "%s/scripts/s07.2.IDR_usingMacs2Peak_rep_work.sh"              % (home_dir)
      sh_work_file_selfPseudoReps    = "%s/scripts/s07.2.IDR_usingMacs2Peak_selfPseudoReps_work.sh"   % (home_dir)
      sh_work_file_pooledPseudoReps  = "%s/scripts/s07.2.IDR_usingMacs2Peak_pooledPseudoReps_work.sh" % (home_dir)
      
      sh_info = """
R_exe=$1
batchIDR_exe=$2
in_sam1=$3
in_sam2=$4
out=$5
genome_table=$6

$R_exe $batchIDR_exe $in_sam1 $in_sam2  -1  $out 0 F p.value $genome_table
      """

      sh_work_rep              = ""
      sh_work_selfPseudoReps   = ""
      sh_work_pooledPseudoReps = ""
      for merge_name in self['sam_info']['merge']:   
         l_sam     = self['sam_info']['merge_sam'][merge_name]
         l_brief   = [ "%s" % (self['sam_info']['samp_brief'][sam]) for sam in l_sam ]
         
         ### IDR ANALYSIS ON ORIGINAL REPLICATES
         if not os.path.isdir( "%s/%s/rep"              % (Peak_idr_dir,merge_name) ):
            os.mkdir(          "%s/%s/rep"              % (Peak_idr_dir,merge_name) )
         
         if len( l_brief ) > 1:
            for i    in range( 0,len(l_brief)-1 ):
               for j in range( i+1,len(l_brief) ):
                  sam1 = l_brief[i]
                  sam2 = l_brief[j]
                  in_sam1 = "%s/%s/%s_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, sam1       )
                  in_sam2 = "%s/%s/%s_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, sam2       )
                  out     = "%s/%s/rep/%s_VS_%s"                    % ( Peak_idr_dir, merge_name, sam1, sam2 )
                  sh_work_rep  += " sh %s   %s %s   %s %s   %s %s\n" % (  sh_file,  R_exe,batchIDR_exe,   in_sam1,in_sam2, out,"%s.len"%(self['infile']['genome_file'])  )
         
         ### IDR ANALYSIS ON SELF-PSEUDOREPLICATES
         if not os.path.isdir( "%s/%s/selfPseudoReps"   % (Peak_idr_dir,merge_name) ):
            os.mkdir(          "%s/%s/selfPseudoReps"   % (Peak_idr_dir,merge_name) )
         for sam in l_brief:
            in_sam_pr1 = "%s/%s/%s.pr1_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, sam )
            in_sam_pr2 = "%s/%s/%s.pr2_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, sam )
            out        = "%s/%s/selfPseudoReps/%s.pr1_VS_%s.pr2"     % ( Peak_idr_dir, merge_name, sam, sam )
            sh_work_selfPseudoReps  += " sh %s   %s %s   %s %s   %s %s\n" % (  sh_file,  R_exe,batchIDR_exe,   in_sam_pr1,in_sam_pr2, out,"%s.len"%(self['infile']['genome_file'])  )
            
         
         ### IDR ANALYSIS ON POOLED-PSEUDOREPLICATES
         if not os.path.isdir( "%s/%s/pooledPseudoReps"   % (Peak_idr_dir,merge_name) ):
            os.mkdir(          "%s/%s/pooledPseudoReps"   % (Peak_idr_dir,merge_name) )
         in_mrg_pr1 = "%s/%s/%s.pr1_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, merge_name )
         in_mrg_pr2 = "%s/%s/%s.pr2_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, merge_name )
         out        = "%s/%s/pooledPseudoReps/%s.pr1_VS_%s.pr2"     % ( Peak_idr_dir, merge_name, merge_name, merge_name )
         sh_work_pooledPseudoReps  += " sh %s   %s %s   %s %s   %s %s\n" % (  sh_file,  R_exe,batchIDR_exe,   in_mrg_pr1,in_mrg_pr2, out,"%s.len"%(self['infile']['genome_file'])  )
      
      my_job_rep = m_jobs.running_jobs(sh_file,sh_work_file_rep)
      my_job_rep.load_sh_file(      sh_info )
      my_job_rep.load_sh_work_file( sh_work_rep )
      my_job_rep.running_multi( cpu=8 )

      my_job_selfPseudoReps = m_jobs.running_jobs(sh_file,sh_work_file_selfPseudoReps)
      my_job_selfPseudoReps.load_sh_file(      sh_info )
      my_job_selfPseudoReps.load_sh_work_file( sh_work_selfPseudoReps )
      my_job_selfPseudoReps.running_multi( cpu=8 )

      my_job_pooledPseudoReps = m_jobs.running_jobs(sh_file,sh_work_file_pooledPseudoReps)
      my_job_pooledPseudoReps.load_sh_file(      sh_info )
      my_job_pooledPseudoReps.load_sh_work_file( sh_work_pooledPseudoReps )
      my_job_pooledPseudoReps.running_multi( cpu=8 )


   def get_idr_stat( self ):
      home_dir    = os.path.abspath('./')
      Peak_idr_dir= self['dir_name'].Peak_idr
      bin_dir     = self['dir_name'].bin
      script_dir  = self['dir_name'].script
      stat_dir    = self['dir_name'].stat
      
      if not os.path.isdir( stat_dir ):
         os.mkdir( stat_dir )
      
      file_idr_out= "%s/IDR_result.xls" % ( stat_dir )
      f_idr_out   = open( file_idr_out,"w" ) 
      
      
      print >>f_idr_out, "Merge_sam\tPeak\trepCmp_peak\tpseudoRep_peak\tRange_rep\tRange_cmp\trepCmp_sam\tpseudoRep_sam"
      for merge_name in self['sam_info']['merge']:   
         l_sam     = self['sam_info']['merge_sam'][merge_name]
         l_brief   = [ "%s" % (self['sam_info']['samp_brief'][sam]) for sam in l_sam ]
         
         l_repCmp         = []
         l_repCmpPeak     = []
         l_selfPseCmp     = []
         l_selfPseCmpPeak = []
         
         if len( l_brief ) > 1:
            ### IDR RESULTS ON ORIGINAL REPLICATES         
            for i    in range( 0,len(l_brief)-1 ):
               for j in range( i+1,len(l_brief) ):
                  sam1 = l_brief[i]
                  sam2 = l_brief[j]
                  in_sam1 = "%s/%s/%s_VS_Input_peaks.regionPeak.gz"   % ( Peak_idr_dir, merge_name, sam1       )
                  in_sam2 = "%s/%s/%s_VS_Input_peaks.regionPeak.gz"   % ( Peak_idr_dir, merge_name, sam2       )
                  rep_cmp = "%s/%s/rep/%s_VS_%s-overlapped-peaks.txt" % ( Peak_idr_dir, merge_name, sam1, sam2 )
                  peak_cnt= get_PassPeak( rep_cmp, 0.03 )
                  l_repCmp.append(     "%s_%s" % (sam1,sam2) )
                  l_repCmpPeak.append( peak_cnt )
         
            ### IDR RESULTS ON SELF-PSEUDOREPLICATES
            for sam in l_brief:
               in_sam_pr1 = "%s/%s/%s.pr1_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, sam )
               in_sam_pr2 = "%s/%s/%s.pr2_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, sam )
               selfPse_cmp= "%s/%s/selfPseudoReps/%s.pr1_VS_%s.pr2-overlapped-peaks.txt"     % ( Peak_idr_dir, merge_name, sam, sam )
               peak_cnt= get_PassPeak( selfPse_cmp, 0.02 )
               l_selfPseCmp.append(     sam )
               l_selfPseCmpPeak.append( peak_cnt )

         
         ### IDR RESULTS ON POOLED-PSEUDOREPLICATES
         in_mrg_pr1 = "%s/%s/%s.pr1_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, merge_name )
         in_mrg_pr2 = "%s/%s/%s.pr2_VS_Input_peaks.regionPeak.gz" % ( Peak_idr_dir, merge_name, merge_name )
         poolPse_cmp= "%s/%s/pooledPseudoReps/%s.pr1_VS_%s.pr2-overlapped-peaks.txt"     % ( Peak_idr_dir, merge_name, merge_name, merge_name )
         
         PeakCnt    = get_PassPeak( poolPse_cmp, 0.006 )
         repCmp         = "-"
         repCmpPeak     = "0"
         selfPseCmp     = "-"
         selfPseCmpPeak = "0"
         range_rep      = np.nan
         range_selfPse  = np.nan
         
         if len( l_brief ) > 1:
            PeakCnt        = max( l_repCmpPeak )
            repCmp         = ",".join( l_repCmp )
            repCmpPeak     = ",".join( [ "%s" % i for i in l_repCmpPeak     ] )
            selfPseCmp     = ",".join( l_selfPseCmp )
            selfPseCmpPeak = ",".join( [ "%s" % i for i in l_selfPseCmpPeak ] )
            range_rep      = max( l_repCmpPeak     )/( min( l_repCmpPeak     )+0.001)
            range_selfPse  = max( l_selfPseCmpPeak )/( min( l_selfPseCmpPeak )+0.001)

         print >>f_idr_out,  "%s\t%d\t%s\t%s\t%f\t%f\t%s\t%s" % ( merge_name,  PeakCnt,  repCmpPeak,selfPseCmpPeak,  range_rep,range_selfPse,  repCmp,selfPseCmp )

      f_idr_out.close()


   def get_idr_Peak( self ):
      home_dir    = os.path.abspath('./')
      Peak_idr_dir= self['dir_name'].Peak_idr
      bin_dir     = self['dir_name'].bin
      script_dir  = self['dir_name'].script
      stat_dir    = self['dir_name'].stat
      
      file_idr_out= "%s/IDR_result.xls" % ( stat_dir )
      f_idr_out   = open( file_idr_out,"r" ) 
            
      sh_file       = "%s/s07.3.IDR_pass_Peaks.sh"      % ( script_dir )
      sh_work_file  = "%s/s07.3.IDR_pass_Peaks_work.sh" % ( script_dir )
      
      bgzip_exe     = self['sftw_name'].bgzip
      tabix_exe     = self['sftw_name'].tabix
      bedtools_exe  = self['sftw_name'].bedtools
      
      sh_info = """
infile=$1
outfile=$2
peak=$3
bedtools_exe=$4
bgzip_exe=$5
tabix_exe=$6

zcat $infile | sort -k7nr,7nr | head -n $peak | $bedtools_exe sort -i /dev/stdin | $bgzip_exe -cf >$outfile && $tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 0 $outfile
      """
      
      sh_work = ""
      
      h = f_idr_out.readline()
      for line in f_idr_out:
         line     = line.strip()
         f        = line.split()
         merge_name=     f[0]

         infile    = "%s/%s/%s_VS_Input_peaks.regionPeak.gz"              % ( Peak_idr_dir, merge_name, merge_name )
         outfile   = "%s/%s/%s_VS_Input_peaks.conservative.regionPeak.gz" % ( Peak_idr_dir, merge_name, merge_name )
         peak_cnt  = int(f[1])
         
         sh_work  += "sh %s  %s %s %d   %s %s %s\n" %  (  sh_file,  infile,outfile,peak_cnt,   bedtools_exe,bgzip_exe,tabix_exe  )
         
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )



   def run_macs_rep_broad(self, pvalue=0.05):
      home_dir    = os.path.abspath('./')
      bed_rep_dir = self['dir_name'].bed_rep
      bed_mrg_dir = self['dir_name'].bed_mrg
      macs_rep_dir= self['dir_name'].BroadPeak_rep
      
      if not os.path.isdir( macs_rep_dir ):
         os.mkdir( macs_rep_dir )
      
      script_dir   = self['dir_name'].script
      bin_dir      = self['dir_name'].bin
      stat_dir     = self['dir_name'].stat

      py_exe       = self['sftw_name'].py
      macs2_exe    = self['sftw_name'].macs2
      
      sh_file       = "%s/scripts/s08.macs2BroadPeakRep.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s08.macs2BroadPeakRep_work.sh" % (home_dir)
      
      sh_info = """
py_exe=$1
macs2_exe=$2
treat_bed=$3
ctrl_bed=$4
out_dir=$5
name_prefix=$6
pvalue=$7

$py_exe $macs2_exe callpeak --broad -t $treat_bed -c $ctrl_bed -f BED -g hs --outdir $out_dir -n $name_prefix -p $pvalue -B
      """
      
      sh_work = ""
      for samp in self['sam_info']['list']:
         brief_name = self['sam_info']['samp_brief'][samp]
         ctrl_name  = self['sam_info']['sam_ctrl'][samp]
         if not os.path.isdir( "%s/%s" % (macs_rep_dir,brief_name) ):
            os.mkdir(          "%s/%s" % (macs_rep_dir,brief_name) )

         name_prefix = brief_name
         out_dir     = "%s/%s"                    % (macs_rep_dir,brief_name)
         treat_bed   = "%s/%s.sort.tagAlign.gz"   % (bed_rep_dir ,brief_name)
         ctrl_bed    = "%s/%s.sort.tagAlign.gz"   % (bed_mrg_dir ,ctrl_name)
         
         sh_work  += " sh %s  %s %s %s %s %s %s %f\n" % ( sh_file, py_exe,  macs2_exe, treat_bed, ctrl_bed, out_dir, name_prefix, pvalue )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )
      

   def run_macs_mrg_broad(self, pvalue=0.05):
      home_dir    = os.path.abspath('./')
      bed_mrg_dir = self['dir_name'].bed_mrg
      macs_mrg_dir= self['dir_name'].BroadPeak_mrg
      
      if not os.path.isdir( macs_mrg_dir ):
         os.mkdir( macs_mrg_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      py_exe       = self['sftw_name'].py
      macs2_exe    = self['sftw_name'].macs2
      
      sh_file       = "%s/scripts/s09.macs2BroadPeakMrg.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s09.macs2BroadPeakMrg_work.sh" % (home_dir)
      
      sh_info = """
py_exe=$1
macs2_exe=$2
treat_bed=$3
ctrl_bed=$4
out_dir=$5
name_prefix=$6
pvalue=$7

$py_exe $macs2_exe callpeak --broad -t $treat_bed -c $ctrl_bed -f BED -g hs --outdir $out_dir -n $name_prefix -p $pvalue -B
      """
      
      sh_work = ""
      for merge_name in self['sam_info']['merge']:
         ctrl_name  = self['sam_info']['merge_ctrl'][merge_name]
         if not os.path.isdir( "%s/%s" % (macs_mrg_dir,merge_name) ):
            os.mkdir(          "%s/%s" % (macs_mrg_dir,merge_name) )

         name_prefix = merge_name
         out_dir     = "%s/%s"                    % (macs_mrg_dir,merge_name)
         treat_bed   = "%s/%s.sort.tagAlign.gz"        % (bed_mrg_dir ,merge_name)
         ctrl_bed    = "%s/%s.sort.tagAlign.gz"        % (bed_mrg_dir ,ctrl_name)
         
         sh_work  += " sh %s  %s %s %s %s %s %s %f\n" % ( sh_file, py_exe,  macs2_exe, treat_bed, ctrl_bed, out_dir, name_prefix, pvalue )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )


   def sort_bdg(self):
      bed_mrg_dir = self['dir_name'].bed_mrg
      
      macs_Rep_dir     = self['dir_name'].Peak_rep
      macs_Mrg_dir     = self['dir_name'].Peak_mrg
      macs_broadRep_dir= self['dir_name'].BroadPeak_rep
      macs_broadMrg_dir= self['dir_name'].BroadPeak_mrg
            
      sort_bdg_exe               = self['sftw_name'].sort_bdg
      script_dir                 = self['dir_name'].script
      bin_dir                    = self['dir_name'].bin
      sh_file                    = "%s/s10.sortbdg.sh"                     % ( script_dir )
      sh_work_repPeak_file       = "%s/s10.1.sortbdg.repPeak_work.sh"      % ( script_dir )
      sh_work_mrgPeak_file       = "%s/s10.2.sortbdg.mrgPeak_work.sh"      % ( script_dir )
      sh_work_repBroadPeak_file  = "%s/s10.3.sortbdg.repBroadPeak_work.sh" % ( script_dir )
      sh_work_mrgBroadPeak_file  = "%s/s10.4.sortbdg.mrgBroadPeak_work.sh" % ( script_dir )
      
      sh_info = """
sam=$1
dir=$2
sort_bdg_exe=$3

sh $sort_bdg_exe $dir/$sam/${sam}_VS_Input_treat_pileup.bdg   $dir/$sam/${sam}_VS_Input_treat_pileup.sort.bdg && \
sh $sort_bdg_exe $dir/$sam/${sam}_VS_Input_control_lambda.bdg $dir/$sam/${sam}_VS_Input_control_lambda.sort.bdg
      """
      
      sh_work_repPeak       = ""
      sh_work_mrgPeak       = ""
      sh_work_repBroadPeak  = ""
      sh_work_mrgBroadPeak  = ""
      
      for merge_name in self['sam_info']['merge']:         
         l_sam     = self['sam_info']['merge_sam'][merge_name]
         l_brief   = [ "%s" % (self['sam_info']['samp_brief'][sam]) for sam in l_sam ]
         for sam  in l_brief:
            insam_rb  = sam
            indir_rb  = macs_broadRep_dir
            sh_work_repBroadPeak += "sh %s  %s %s %s\n" % ( sh_file, insam_rb, indir_rb, sort_bdg_exe )
         
         insam_mb     = merge_name
         indir_mb     = macs_broadMrg_dir
         sh_work_mrgBroadPeak    += "sh %s  %s %s %s\n" % ( sh_file, insam_mb,  indir_mb, sort_bdg_exe )


         for sam  in l_brief:
            insam_r  = sam
            indir_r  = macs_Rep_dir
            sh_work_repPeak += "sh %s  %s %s %s\n" % ( sh_file, insam_r,   indir_r,  sort_bdg_exe )
         
         insam_m     = merge_name
         indir_m     = macs_Mrg_dir
         sh_work_mrgPeak    += " sh %s  %s %s %s\n"  % ( sh_file, insam_m,   indir_m,  sort_bdg_exe )

         
      
      my_job_r = m_jobs.running_jobs(sh_file,sh_work_repPeak_file)
      my_job_r.load_sh_file(      sh_info )
      my_job_r.load_sh_work_file( sh_work_repPeak )
      my_job_r.running_multi( cpu=8 )

      my_job_m = m_jobs.running_jobs(sh_file,sh_work_mrgPeak_file)
      my_job_m.load_sh_file(      sh_info )
      my_job_m.load_sh_work_file( sh_work_mrgPeak )
      my_job_m.running_multi( cpu=8 )

      my_job_rb = m_jobs.running_jobs(sh_file,sh_work_repBroadPeak_file)
      my_job_rb.load_sh_file(      sh_info )
      my_job_rb.load_sh_work_file( sh_work_repBroadPeak )
      my_job_rb.running_multi( cpu=8 )

      my_job_mb = m_jobs.running_jobs(sh_file,sh_work_mrgBroadPeak_file)
      my_job_mb.load_sh_file(      sh_info )
      my_job_mb.load_sh_work_file( sh_work_mrgBroadPeak )
      my_job_mb.running_multi( cpu=8 )
      


   def make_igv_IDR(self):
      bed_mrg_dir = self['dir_name'].bed_mrg
      TDF_mrg_dir = self['dir_name'].Peak_TDF
      macs_Rep_dir     = self['dir_name'].Peak_rep
      macs_Mrg_dir     = self['dir_name'].Peak_mrg
      macs_IDR_dir     = self['dir_name'].Peak_idr
      
      igvtools_exe               = self['sftw_name'].igvtools
      script_dir                 = self['dir_name'].script
      bin_dir                    = self['dir_name'].bin
      sh_file_rep                = "%s/s11.1.makeIGV.idrRep.sh"            % ( script_dir )
      sh_file_mrg                = "%s/s11.2.makeIGV.idrMrg.sh"            % ( script_dir )
      sh_work_repPeak_file       = "%s/s11.1.makeIGV.repPeak_work.sh"      % ( script_dir )
      sh_work_mrgPeak_file       = "%s/s11.2.makeIGV.mrgPeak_work.sh"      % ( script_dir )
      
      sh_info_rep = """
sam_in=$1
sam_out=$2
dir_in=$3
dir_out=$4
igvtools_exe=$5
genome_fa=$6

ln -s $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bdg   $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph
ln -s $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bdg $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph
$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph    $dir_out/${sam_out}/${sam_in}_treat_pileup.tdf   $genome_fa && \
$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph  $dir_out/${sam_out}/${sam_in}_control_lambda.tdf $genome_fa
rm  $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph
      """

      sh_info_mrg = """
sam_in=$1
sam_out=$2
dir_in=$3
dir_out=$4
igvtools_exe=$5
genome_fa=$6
dir_idrPeak=$7


ln -s $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bdg   $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph
ln -s $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bdg $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph
$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph    $dir_out/${sam_out}/${sam_in}_VS_Input_treat_pileup.tdf   $genome_fa && \
$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph  $dir_out/${sam_out}/${sam_in}_VS_Input_control_lambda.tdf $genome_fa
zcat $dir_idrPeak/$sam_in/${sam_in}_VS_Input_peaks.conservative.regionPeak.gz >$dir_idrPeak/${sam_in}/${sam_in}_VS_Input_peaks.conservative.regionPeak.bed
$igvtools_exe count $dir_idrPeak/${sam_in}/${sam_in}_VS_Input_peaks.conservative.regionPeak.bed $dir_out/${sam_out}/${sam_in}_peaks.IDR.tdf $genome_fa && \
rm $dir_idrPeak/${sam_in}/${sam_in}_VS_Input_peaks.conservative.regionPeak.bed $dir_in/$sam_in/${sam_in}_VS_Input_treat_pileup.bedGraph $dir_in/$sam_in/${sam_in}_VS_Input_control_lambda.bedGraph
      """
      
      sh_work_repPeak  = ""
      sh_work_mrgPeak  = ""
      
      if not os.path.isdir( TDF_mrg_dir ):
         os.mkdir(          TDF_mrg_dir )
      for merge_name in self['sam_info']['merge']:  
         if not os.path.isdir( "%s/%s" % ( TDF_mrg_dir,merge_name ) ):
            os.mkdir(          "%s/%s" % ( TDF_mrg_dir,merge_name ) )
                
         l_sam     = self['sam_info']['merge_sam'][merge_name]
         l_brief   = [ "%s" % (self['sam_info']['samp_brief'][sam]) for sam in l_sam ]
         for sam  in l_brief:
            insam_r  = sam
            outsam_r = merge_name
            indir_r  = macs_Rep_dir
            outdir_r = TDF_mrg_dir
            
            sh_work_repPeak += "sh %s  %s %s %s  %s %s %s\n" % ( sh_file_rep, insam_r, outsam_r,  indir_r, outdir_r, igvtools_exe, self['infile']['genome_file'] )
         
         insam_m  = merge_name
         outsam_m = merge_name
         indir_m  = macs_Mrg_dir
         outdir_m = TDF_mrg_dir
         
         sh_work_mrgPeak += "sh %s  %s %s %s  %s %s %s  %s\n" % ( sh_file_mrg, insam_m, outsam_m, indir_m, outdir_m, igvtools_exe,self['infile']['genome_file'], macs_IDR_dir )

      my_job_r = m_jobs.running_jobs(sh_file_rep,sh_work_repPeak_file)
      my_job_r.load_sh_file(      sh_info_rep )
      my_job_r.load_sh_work_file( sh_work_repPeak )
      my_job_r.running_multi( cpu=8 )

      my_job_m = m_jobs.running_jobs(sh_file_mrg,sh_work_mrgPeak_file)
      my_job_m.load_sh_file(      sh_info_mrg )
      my_job_m.load_sh_work_file( sh_work_mrgPeak )
      my_job_m.running_multi( cpu=8 )


      
   def make_igv_broad(self):
      bed_mrg_dir = self['dir_name'].bed_mrg
      TDF_mrg_dir = self['dir_name'].BroadPeak_TDF
      macs_Rep_dir     = self['dir_name'].Peak_rep
      macs_Mrg_dir     = self['dir_name'].Peak_mrg
      macs_broadRep_dir= self['dir_name'].BroadPeak_rep
      macs_broadMrg_dir= self['dir_name'].BroadPeak_mrg
            
      igvtools_exe               = self['sftw_name'].igvtools
      script_dir                 = self['dir_name'].script
      bin_dir                    = self['dir_name'].bin
      sh_file                    = "%s/s11.3.makeIGV_broad.sh"             % ( script_dir )

      sh_work_repBroadPeak_file  = "%s/s11.3.makeIGV.repBroadPeak_work.sh" % ( script_dir )
      sh_work_mrgBroadPeak_file  = "%s/s11.4.makeIGV.mrgBroadPeak_work.sh" % ( script_dir )
      
      sh_info = """
sam_in=$1
sam_out=$2
dir_in=$3
dir_out=$4
igvtools_exe=$5
genome_fa=$6

ln -s $dir_in/$sam_in/${sam_in}_treat_pileup.bdg   $dir_in/$sam_in/${sam_in}_treat_pileup.bedGraph
ln -s $dir_in/$sam_in/${sam_in}_control_lambda.bdg $dir_in/$sam_in/${sam_in}_control_lambda.bedGraph
$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_treat_pileup.bedGraph    $dir_out/${sam_out}/${sam_in}_treat_pileup.broadPeaks.tdf   $genome_fa && \
$igvtools_exe toTDF $dir_in/$sam_in/${sam_in}_control_lambda.bedGraph  $dir_out/${sam_out}/${sam_in}_control_lambda.broadPeaks.tdf $genome_fa
ln -s $dir_in/$sam_in/${sam_in}_peaks.broadPeak $dir_in/$sam_in/${sam_in}_peaks.broadPeak.bed 
$igvtools_exe count $dir_in/$sam_in/${sam_in}_peaks.broadPeak.bed $dir_out/${sam_out}/${sam_in}_peaks.broadPeak.tdf $genome_fa && \
rm $dir_in/$sam_in/${sam_in}_treat_pileup.bedGraph $dir_in/$sam_in/${sam_in}_control_lambda.bedGraph $dir_in/$sam_in/${sam_in}_peaks.broadPeak.bed
      """

      sh_work_repBroadPeak  = ""
      sh_work_mrgBroadPeak  = ""
      
      if not os.path.isdir( TDF_mrg_dir ):
         os.mkdir(          TDF_mrg_dir )
      for merge_name in self['sam_info']['merge']:  
         if not os.path.isdir( "%s/%s" % ( TDF_mrg_dir,merge_name ) ):
            os.mkdir(          "%s/%s" % ( TDF_mrg_dir,merge_name ) )
                
         l_sam     = self['sam_info']['merge_sam'][merge_name]
         l_brief   = [ "%s" % (self['sam_info']['samp_brief'][sam]) for sam in l_sam ]
         for sam  in l_brief:
            insam_rb  = sam
            outsam_rb = merge_name
            indir_rb  = macs_broadRep_dir
            outdir_rb = TDF_mrg_dir
            
            sh_work_repBroadPeak += "sh %s  %s %s %s  %s %s %s\n" % ( sh_file, insam_rb, outsam_rb,  indir_rb,  outdir_rb, igvtools_exe, self['infile']['genome_file'] )
         
         insam_mb  = merge_name
         outsam_mb = merge_name
         indir_mb  = macs_broadMrg_dir
         outdir_mb = TDF_mrg_dir
         
         sh_work_mrgBroadPeak += "sh %s  %s %s %s  %s %s %s\n" % ( sh_file, insam_mb, outsam_mb, indir_mb, outdir_mb, igvtools_exe,self['infile']['genome_file'] )

      my_job_rb = m_jobs.running_jobs(sh_file,sh_work_repBroadPeak_file)
      my_job_rb.load_sh_file(      sh_info )
      my_job_rb.load_sh_work_file( sh_work_repBroadPeak )
      my_job_rb.running_multi( cpu=8 )

      my_job_mb = m_jobs.running_jobs(sh_file,sh_work_mrgBroadPeak_file)
      my_job_mb.load_sh_file(      sh_info )
      my_job_mb.load_sh_work_file( sh_work_mrgBroadPeak )
      my_job_mb.running_multi( cpu=8 )