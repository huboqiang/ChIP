#-*- coding:utf-8 -*-
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

import module_refGene as m_refG
import module_running_jobs as m_jobs

class GeneBasedInfo(object):

   def __init__(self, refGeneTxt, genome_ref, samp_info, stat_Info, dir_name, sftw_name ):
      self.refGeneTxt  = refGeneTxt
      self.genomeFa    = genome_ref
      self.stat_Info   = stat_Info
      self.samp_info   = samp_info['samp']['chip']
      self.dir_name    = dir_name
      self.sftw_name   = sftw_name

   def extend_gene_region(self,TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,TSS_promoter_down):
      
      data_db = "%s.dat" % (self.refGeneTxt)
      try:
         refG_info = pickle.load( open(data_db) )
      except:
         refG_info = m_refG.RefGeneTxt( self.refGeneTxt,self.genomeFa  )
         refG_info.refGene2bed(  0,0 )
         refG_info.refGene2bed(  TSS_genebody_up,  TSS_genebody_down,  ext_type="genebody.")
         refG_info.refGene2bed(  TSS_promoter_up,  TSS_promoter_down,  ext_type="promoter.")
         refG_info.refGeneInfo(  TSS_promoter_up,  TSS_promoter_down )
         refG_info.Only_LongestTid_Bed( TSS_genebody_up,  TSS_genebody_down,  ext_type="genebody." )
         refG_info.Only_LongestTid_Bed( TSS_promoter_up,  TSS_promoter_down,  ext_type="promoter." )
         pickle.dump( refG_info,open(data_db,"wb"),True )
         

   def run_anno_peak(self,  TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,TSS_promoter_down, ext_binlen=50,body_bincnt=100,tss_binlen=1 ):
      home_dir = os.path.abspath('./')
      
      Peak_mrg_dir = self.dir_name.Peak_mrg
      Peak_idr_dir = self.dir_name.Peak_idr
      Peak_mrg_TSS = self.dir_name.Peak_mrg_TSS 
      Peak_mrg_Gene= self.dir_name.Peak_mrg_Gene
      
      script_dir  = self.dir_name.script
      bin_dir     = self.dir_name.bin
      
      if not os.path.exists( Peak_mrg_TSS ):
         os.mkdir( Peak_mrg_TSS )
      if not os.path.exists( Peak_mrg_Gene ):
         os.mkdir( Peak_mrg_Gene )
         
      sh_file        = "%s/s12.PeakGeneRegion.sh"       % ( script_dir )
      sh_work_file   = "%s/s12.PeakGeneRegion_work.sh"  % ( script_dir )
      
      bedtools_exe       = self.sftw_name.bedtools
      getbin_genebody_cpp="%s/get_bin_ignore_peak/get_bin"          % ( bin_dir )
      
      prefix = ".".join( self.refGeneTxt.split(".")[:-1] )
      genebody_region   = "%s.up%d_down%d.%sBsorted.longestTid.bed" % ( prefix, TSS_genebody_up, TSS_genebody_down, "genebody." )
      
      sh_info = """
bedtools_exe=$1
count_bdg=$2
peak_bed=$3
genebody_region=$4
getbin_genebody_cpp=$5
TSS_genebody_up=${6}
TSS_genebody_down=${7}
TSS_promoter_up=${8}
TSS_promoter_down=${9}
ext_binlen=${10}
body_bincnt=${11}
tss_binlen=${12}
stat_file1=${13}
stat_file2=${14}

zcat $count_bdg | tail -n +2 | $bedtools_exe intersect -sorted -wo -a /dev/stdin -b $genebody_region | $getbin_genebody_cpp -U $TSS_genebody_up -D $TSS_genebody_down -T $TSS_promoter_up -b $ext_binlen -B $body_bincnt /dev/stdin $stat_file1 $stat_file2
      """
      sh_work = ""
      
      for merge_name in self.samp_info['merge']:
         if not os.path.exists(  "%s/%s" % (Peak_mrg_TSS ,merge_name) ):
            os.mkdir(            "%s/%s" % (Peak_mrg_TSS ,merge_name) )            
         if not os.path.exists(  "%s/%s" % (Peak_mrg_Gene,merge_name) ):
            os.mkdir(            "%s/%s" % (Peak_mrg_Gene,merge_name) )            

         count_bdg = "%s/%s/%s_VS_Input_treat_pileup.sort.bdg.gz"         % ( Peak_mrg_dir,merge_name,merge_name )
         peak_bed  = "%s/%s/%s_VS_Input_peaks.conservative.regionPeak.gz" % ( Peak_idr_dir,merge_name,merge_name )
         
         stat_file1         = "%s/%s/%s.genebody.up%d_down%d.xls"      % ( Peak_mrg_Gene, merge_name,merge_name, TSS_genebody_up, TSS_genebody_down )
         stat_file2         = "%s/%s/%s.tss.up%d_down%d.xls"           % ( Peak_mrg_TSS , merge_name,merge_name, TSS_promoter_up, TSS_promoter_down )
      
         sh_work += "sh %s   %s %s %s %s %s   %d %d %d %d   %d %d %d   %s %s \n" % ( sh_file,  bedtools_exe,count_bdg,peak_bed,genebody_region,getbin_genebody_cpp,\
              TSS_genebody_up,TSS_genebody_down,TSS_promoter_up,TSS_promoter_down,   ext_binlen,body_bincnt,tss_binlen,   stat_file1,stat_file2 )
         
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )

   def div_bed_to_bins_rep(self):
      home_dir = os.path.abspath('./')
      
      Bed_rep     = self.dir_name.bed_rep
      RPM_bins    = self.dir_name.RPM_bins_rep 
      
      script_dir  = self.dir_name.script
      bin_dir     = self.dir_name.bin
      
      if not os.path.exists( RPM_bins ):
         os.mkdir( RPM_bins )
         
      sh_file        = "%s/s13.RPM_density_rep.sh"       % ( script_dir )
      sh_work_file   = "%s/s13.RPM_density_rep_work.sh"  % ( script_dir )
      
      bedtools_exe       = self.sftw_name.bedtools
      getbin_readcnt_cpp = "%s/chip01.bed_reads/bed_read"% ( bin_dir )
      
      genome_fai     = "%s.fai" % (self.genomeFa)
      
      sh_info = """
getbin_readcnt_cpp=$1
genome_fai=$2
unique_bed=$3
multi_bed=$4
out_cnt=$5
mapped_reads=$6
out_RPKM=$7

$getbin_readcnt_cpp -b 100  $genome_fai $unique_bed $multi_bed $out_cnt.100 
awk '{ print $1*1000000*1000/(100 *'$mapped_reads') }' $out_cnt.100 >$out_RPKM.100
$getbin_readcnt_cpp -b 1000 $genome_fai $unique_bed $multi_bed $out_cnt.1kb 
awk '{ print $1*1000000*1000/(1000*'$mapped_reads') }' $out_cnt.1kb >$out_RPKM.1kb
      """
      sh_work = ""
      
      for sam in self.samp_info['list']:
         brief_name   = self.samp_info['samp_brief'][sam]
         
         if not os.path.exists(  "%s/%s" % (RPM_bins ,brief_name) ):
            os.mkdir(            "%s/%s" % (RPM_bins ,brief_name) )            

         unique_bed = "%s/%s.unique.bed" % ( Bed_rep ,brief_name )
         multi_bed  = "%s/%s.multi.bed"  % ( Bed_rep ,brief_name )
         out_cnt    = "%s/%s/%s.cnt"     % ( RPM_bins,brief_name,brief_name )
         out_RPKM   = "%s/%s/%s.RPKM"    % ( RPM_bins,brief_name,brief_name )
         
         merge_sam  = self.samp_info['sam_merge'][ sam ]
         idx        = self.samp_info['merge_sam'][merge_sam].index( sam )
         mapped_reads = self.stat_Info.StatInfo[ merge_sam ]['unique'][idx] + self.stat_Info.StatInfo[ merge_sam ]['multi'][idx]
#         print  merge_sam, self.samp_info['merge_sam'][merge_sam], sam, mapped_reads
         sh_work += "sh %s   %s %s %s %s %s   %d  %s \n" % ( sh_file,  getbin_readcnt_cpp,genome_fai, unique_bed,multi_bed,out_cnt,mapped_reads,out_RPKM )
         
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )

   def div_bed_to_bins_mrg(self):
      home_dir = os.path.abspath('./')
      
      Bed_mrg     = self.dir_name.bed_mrg
      RPM_bins    = self.dir_name.RPM_bins_mrg
      
      script_dir  = self.dir_name.script
      bin_dir     = self.dir_name.bin
      
      if not os.path.exists( RPM_bins ):
         os.mkdir( RPM_bins )
         
      sh_file        = "%s/s14.RPM_density_mrg.sh"       % ( script_dir )
      sh_work_file   = "%s/s14.RPM_density_mrg_work.sh"  % ( script_dir )
      
      bedtools_exe       = self.sftw_name.bedtools
      getbin_readcnt_cpp = "%s/chip01.bed_reads/bed_read"% ( bin_dir )
      
      genome_fai     = "%s.fai" % (self.genomeFa)
      
      sh_info = """
getbin_readcnt_cpp=$1
genome_fai=$2
unique_bed=$3
multi_bed=$4
out_cnt=$5
mapped_reads=$6
out_RPKM=$7

$getbin_readcnt_cpp -b 100  $genome_fai $unique_bed $multi_bed $out_cnt.100 
awk '{ print $1*1000000*1000/(100 *'$mapped_reads') }' $out_cnt.100 >$out_RPKM.100
$getbin_readcnt_cpp -b 1000 $genome_fai $unique_bed $multi_bed $out_cnt.1kb 
awk '{ print $1*1000000*1000/(1000*'$mapped_reads') }' $out_cnt.1kb >$out_RPKM.1kb
      """
      sh_work = ""
      
      for merge_name in self.samp_info['merge']:
         if not os.path.exists(  "%s/%s" % (RPM_bins ,merge_name) ):
            os.mkdir(            "%s/%s" % (RPM_bins ,merge_name) )            

         unique_bed = "%s/%s.unique.bed" % ( Bed_mrg ,merge_name )
         multi_bed  = "%s/%s.multi.bed"  % ( Bed_mrg ,merge_name )
         out_cnt    = "%s/%s/%s.cnt"     % ( RPM_bins,merge_name,merge_name )
         out_RPKM   = "%s/%s/%s.RPKM"    % ( RPM_bins,merge_name,merge_name )
         mapped_reads = np.sum( self.stat_Info.StatInfo[ merge_name ]['unique'] ) + np.sum( self.stat_Info.StatInfo[ merge_name ]['multi'] )
         sh_work += "sh %s   %s %s %s %s %s   %d  %s \n" % ( sh_file,  getbin_readcnt_cpp,genome_fai, unique_bed,multi_bed,out_cnt,mapped_reads,out_RPKM )
         
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )   


   def merge_RPKM( self ):
      RPM_bins_rep = self.dir_name.RPM_bins_rep
      RPM_bins_mrg = self.dir_name.RPM_bins_mrg
      RPM_mrg      = self.dir_name.RPM_mrg
      
      script_dir   = self.dir_name.script
      bin_dir      = self.dir_name.bin
      data_dir     = self.dir_name.database
      
      if not os.path.exists( RPM_mrg ):
         os.mkdir( RPM_mrg )
      
      col_100        = "%s/column_100.bin.bed"               % ( data_dir )
      col_1kb        = "%s/column_1kb.bin.bed"               % ( data_dir )
      bgzip_exe      = self.sftw_name.bgzip
      tabix_exe      = self.sftw_name.tabix
      
      mrg_rep_100    = "%s/Merge_RPKM_bin100.rep"            % ( RPM_mrg )
      mrg_rep_1kb    = "%s/Merge_RPKM_bin1kb.rep"            % ( RPM_mrg )
      mrg_mrg_100    = "%s/Merge_RPKM_bin100.mrg"            % ( RPM_mrg )
      mrg_mrg_1kb    = "%s/Merge_RPKM_bin1kb.mrg"            % ( RPM_mrg )
      
      header         = "\"#chr\\tbeg\\tend\t%s\""            % ( "\\t".join( self.samp_info['list'] ) )
      sh_file        = "%s/s15.1.MergeRPKM_100_rep.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.1.MergeRPKM_100_rep_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_100=$2
mrg_rep_100=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e "$header"   >$mrg_rep_100.tmp                  && \\
paste $col_100 $@  >>$mrg_rep_100.tmp                  && \\
$bgzip_exe -c -f $mrg_rep_100.tmp >$mrg_rep_100.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_rep_100.bed.gz
rm $mrg_rep_100.tmp
      """
      l_brief     = [ "%s"                 % ( self.samp_info['samp_brief'][sam]) for sam   in self.samp_info['list'] ]
      l_RPKM_file = [ "%s/%s/%s.RPKM.100"  % ( RPM_bins_rep,brief,brief )  for brief in l_brief ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_100,mrg_rep_100,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file) )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
      sh_file        = "%s/s15.2.MergeRPKM_1kb_rep.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.2.MergeRPKM_1kb_rep_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_1kb=$2
mrg_rep_1kb=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e $header     >$mrg_rep_1kb.tmp                  && \\
paste $col_1kb $@  >>$mrg_rep_1kb.tmp                  && \\
$bgzip_exe -c -f $mrg_rep_1kb.tmp >$mrg_rep_1kb.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_rep_1kb.bed.gz
rm $mrg_rep_1kb.tmp
      """
      l_RPKM_file = [ "%s/%s/%s.RPKM.1kb"  % ( RPM_bins_rep,brief,brief )  for brief in l_brief ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_1kb,mrg_rep_1kb,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file)  )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
      header         = "\"#chr\\tbeg\\tend\t%s\"" % ( "\\t".join( self.samp_info['merge'] ) )
      sh_file        = "%s/s15.3.MergeRPKM_100_mrg.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.3.MergeRPKM_100_mrg_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_100=$2
mrg_mrg_100=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e "$header"   >$mrg_mrg_100.tmp                  && \\
paste $col_100 $@  >>$mrg_mrg_100.tmp                  && \\
$bgzip_exe -c -f $mrg_mrg_100.tmp >$mrg_mrg_100.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_mrg_100.bed.gz
rm $mrg_mrg_100.tmp
      """
      l_RPKM_file = [ "%s/%s/%s.RPKM.100"  % ( RPM_bins_mrg,merge_samp,merge_samp ) for merge_samp in self.samp_info['merge'] ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_100,mrg_mrg_100,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file) )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
      sh_file        = "%s/s15.4.MergeRPKM_1kb_mrg.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.4.MergeRPKM_1kb_mrg_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_1kb=$2
mrg_mrg_1kb=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e $header     >$mrg_mrg_1kb.tmp                  && \\
paste $col_1kb $@  >>$mrg_mrg_1kb.tmp                  && \\
$bgzip_exe -c -f $mrg_mrg_1kb.tmp >$mrg_mrg_1kb.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_mrg_1kb.bed.gz
rm $mrg_mrg_1kb.tmp
      """
      l_RPKM_file = [ "%s/%s/%s.RPKM.1kb"  % ( RPM_bins_mrg,merge_samp,merge_samp ) for merge_samp in self.samp_info['merge'] ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_1kb,mrg_mrg_1kb,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file)  )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )


   def div_bed_to_bins_unique_rep(self):
      home_dir = os.path.abspath('./')
      
      Bed_rep     = self.dir_name.bed_rep
      RPM_bins    = self.dir_name.RPM_bins_rep 
      
      script_dir  = self.dir_name.script
      bin_dir     = self.dir_name.bin
      
      if not os.path.exists( RPM_bins ):
         os.mkdir( RPM_bins )
         
      sh_file        = "%s/s13.1.RPM_density_rep_unique.sh"       % ( script_dir )
      sh_work_file   = "%s/s13.1.RPM_density_rep_unique_work.sh"  % ( script_dir )
      
      bedtools_exe       = self.sftw_name.bedtools
      getbin_readcnt_cpp = "%s/chip01.bed_reads_unique/bed_read"% ( bin_dir )
      
      genome_fai     = "%s.fai" % (self.genomeFa)
      
      sh_info = """
getbin_readcnt_cpp=$1
genome_fai=$2
unique_bed=$3
out_cnt=$4
mapped_reads=$5
out_RPKM=$6

zcat $unique_bed | $getbin_readcnt_cpp -b 100  $genome_fai /dev/stdin $out_cnt.100 
awk '{ print $1*1000000*1000/(100 *'$mapped_reads') }' $out_cnt.100  >$out_RPKM.100
zcat $unique_bed | $getbin_readcnt_cpp -b 1000 $genome_fai /dev/stdin $out_cnt.1kb 
awk '{ print $1*1000000*1000/(1000*'$mapped_reads') }' $out_cnt.1kb  >$out_RPKM.1kb
      """
      sh_work = ""
      
      for sam in self.samp_info['list']:
         brief_name   = self.samp_info['samp_brief'][sam]
         
         if not os.path.exists(  "%s/%s" % (RPM_bins ,brief_name) ):
            os.mkdir(            "%s/%s" % (RPM_bins ,brief_name) )            

         unique_bed = "%s/%s.sort.tagAlign.gz"   % ( Bed_rep ,brief_name )
         out_cnt    = "%s/%s/%s.cnt.uniq"  % ( RPM_bins,brief_name,brief_name )
         out_RPKM   = "%s/%s/%s.RPKM.uniq" % ( RPM_bins,brief_name,brief_name )
         
         merge_sam  = self.samp_info['sam_merge'][ sam ]
         idx        = self.samp_info['merge_sam'][merge_sam].index( sam )
         mapped_reads = self.stat_Info.StatInfo[ merge_sam ]['unique'][idx]
#         print  merge_sam, self.samp_info['merge_sam'][merge_sam], sam, mapped_reads
         sh_work += "sh %s   %s %s %s %s   %d  %s \n" % ( sh_file,  getbin_readcnt_cpp,genome_fai, unique_bed,out_cnt,mapped_reads,out_RPKM )
         
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )

   def div_bed_to_bins_unique_mrg(self):
      home_dir = os.path.abspath('./')
      
      Bed_mrg     = self.dir_name.bed_mrg
      RPM_bins    = self.dir_name.RPM_bins_mrg
      
      script_dir  = self.dir_name.script
      bin_dir     = self.dir_name.bin
      
      if not os.path.exists( RPM_bins ):
         os.mkdir( RPM_bins )
         
      sh_file        = "%s/s14.1.RPM_density_mrg_unique.sh"       % ( script_dir )
      sh_work_file   = "%s/s14.1.RPM_density_mrg_unique_work.sh"  % ( script_dir )
      
      bedtools_exe       = self.sftw_name.bedtools
      getbin_readcnt_cpp = "%s/chip01.bed_reads_unique/bed_read"% ( bin_dir )
      
      genome_fai     = "%s.fai" % (self.genomeFa)
      
      sh_info = """
getbin_readcnt_cpp=$1
genome_fai=$2
unique_bed=$3
out_cnt=$4
mapped_reads=$5
out_RPKM=$6

zcat $unique_bed | $getbin_readcnt_cpp -b 100  $genome_fai /dev/stdin $out_cnt.100 
awk '{ print $1*1000000*1000/(100 *'$mapped_reads') }' $out_cnt.100  >$out_RPKM.100
zcat $unique_bed | $getbin_readcnt_cpp -b 1000 $genome_fai /dev/stdin $out_cnt.1kb 
awk '{ print $1*1000000*1000/(1000*'$mapped_reads') }' $out_cnt.1kb  >$out_RPKM.1kb
      """
      sh_work = ""
      
      for merge_name in self.samp_info['merge']:
         if not os.path.exists(  "%s/%s" % (RPM_bins ,merge_name) ):
            os.mkdir(            "%s/%s" % (RPM_bins ,merge_name) )            

         unique_bed = "%s/%s.sort.tagAlign.gz"   % ( Bed_mrg ,merge_name )
         out_cnt    = "%s/%s/%s.cnt.uniq"  % ( RPM_bins,merge_name,merge_name )
         out_RPKM   = "%s/%s/%s.RPKM.uniq" % ( RPM_bins,merge_name,merge_name )
         mapped_reads = np.sum( self.stat_Info.StatInfo[ merge_name ]['unique'] )
         sh_work += "sh %s   %s %s %s %s   %d  %s \n" % ( sh_file,  getbin_readcnt_cpp,genome_fai, unique_bed,out_cnt,mapped_reads,out_RPKM )
         
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )   



   def merge_RPKM_uniq( self ):
      RPM_bins_rep = self.dir_name.RPM_bins_rep
      RPM_bins_mrg = self.dir_name.RPM_bins_mrg
      RPM_mrg      = self.dir_name.RPM_mrg
      
      script_dir   = self.dir_name.script
      bin_dir      = self.dir_name.bin
      data_dir     = self.dir_name.database
      
      if not os.path.exists( RPM_mrg ):
         os.mkdir( RPM_mrg )
      
      col_100        = "%s/column_100.bin.bed"               % ( data_dir )
      col_1kb        = "%s/column_1kb.bin.bed"               % ( data_dir )
      bgzip_exe      = self.sftw_name.bgzip
      tabix_exe      = self.sftw_name.tabix
      
      mrg_rep_100    = "%s/Merge_RPKM_bin100.uniq.rep"       % ( RPM_mrg )
      mrg_rep_1kb    = "%s/Merge_RPKM_bin1kb.uniq.rep"       % ( RPM_mrg )
      mrg_mrg_100    = "%s/Merge_RPKM_bin100.uniq.mrg"       % ( RPM_mrg )
      mrg_mrg_1kb    = "%s/Merge_RPKM_bin1kb.uniq.mrg"       % ( RPM_mrg )
      
      header         = "\"#chr\\tbeg\\tend\t%s\""            % ( "\\t".join( self.samp_info['list'] ) )
      sh_file        = "%s/s15.5.MergeRPKM_100_rep_uniq.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.5.MergeRPKM_100_rep_uniq_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_100=$2
mrg_rep_100=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e "$header"   >$mrg_rep_100.tmp                  && \\
paste $col_100 $@  >>$mrg_rep_100.tmp                  && \\
$bgzip_exe -c -f $mrg_rep_100.tmp >$mrg_rep_100.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_rep_100.bed.gz
rm $mrg_rep_100.tmp
      """
      l_brief     = [ "%s"                     % ( self.samp_info['samp_brief'][sam]) for sam   in self.samp_info['list'] ]
      l_RPKM_file = [ "%s/%s/%s.RPKM.uniq.100" % ( RPM_bins_rep,brief,brief )  for brief in l_brief ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_100,mrg_rep_100,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file) )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
      sh_file        = "%s/s15.6.MergeRPKM_1kb_rep_uniq.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.6.MergeRPKM_1kb_rep_uniq_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_1kb=$2
mrg_rep_1kb=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e $header     >$mrg_rep_1kb.tmp                  && \\
paste $col_1kb $@  >>$mrg_rep_1kb.tmp                  && \\
$bgzip_exe -c -f $mrg_rep_1kb.tmp >$mrg_rep_1kb.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_rep_1kb.bed.gz
rm $mrg_rep_1kb.tmp
      """
      l_RPKM_file = [ "%s/%s/%s.RPKM.uniq.1kb"  % ( RPM_bins_rep,brief,brief )  for brief in l_brief ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_1kb,mrg_rep_1kb,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file)  )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
      header         = "\"#chr\\tbeg\\tend\t%s\"" % ( "\\t".join( self.samp_info['merge'] ) )
      sh_file        = "%s/s15.7.MergeRPKM_100_mrg_uniq.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.7.MergeRPKM_100_mrg_uniq_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_100=$2
mrg_mrg_100=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e "$header"   >$mrg_mrg_100.tmp                  && \\
paste $col_100 $@  >>$mrg_mrg_100.tmp                  && \\
$bgzip_exe -c -f $mrg_mrg_100.tmp >$mrg_mrg_100.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_mrg_100.bed.gz
rm $mrg_mrg_100.tmp
      """
      l_RPKM_file = [ "%s/%s/%s.RPKM.uniq.100"  % ( RPM_bins_mrg,merge_samp,merge_samp ) for merge_samp in self.samp_info['merge'] ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_100,mrg_mrg_100,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file) )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
      
      
      sh_file        = "%s/s15.8.MergeRPKM_1kb_mrg_uniq.sh"       % ( script_dir )
      sh_work_file   = "%s/s15.8.MergeRPKM_1kb_mrg_uniq_work.sh"  % ( script_dir )
      sh_info = """
header=$1
col_1kb=$2
mrg_mrg_1kb=$3
bgzip_exe=$4
tabix_exe=$5
shift
shift
shift
shift
shift

echo -e $header     >$mrg_mrg_1kb.tmp                  && \\
paste $col_1kb $@  >>$mrg_mrg_1kb.tmp                  && \\
$bgzip_exe -c -f $mrg_mrg_1kb.tmp >$mrg_mrg_1kb.bed.gz && \\
$tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 1 $mrg_mrg_1kb.bed.gz
rm $mrg_mrg_1kb.tmp
      """
      l_RPKM_file = [ "%s/%s/%s.RPKM.uniq.1kb"  % ( RPM_bins_mrg,merge_samp,merge_samp ) for merge_samp in self.samp_info['merge'] ]
      sh_work     = " sh %s   %s  %s %s  %s %s   %s" % ( sh_file,  header,  col_1kb,mrg_mrg_1kb,  bgzip_exe,tabix_exe, " ".join(l_RPKM_file)  )
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
