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

class Map_From_raw(dict):
   
   def __init__(self,samp_info, genome_file,anno_file, dir_name, sftw_name ):
      self['sam_info'] = samp_info['samp']['chip']
      self['infile']   = { 'genome_file':genome_file,'anno_file':anno_file }
      self['stage']    = { 'name':[] }
      self['dir_name']  = dir_name
      self['sftw_name'] = sftw_name

   def run_QC(self):
      home_dir    = os.path.abspath('./')
      raw_dir     = self['dir_name'].raw_data
      cln_dir     = self['dir_name'].clean_data
      
      if not os.path.isdir( cln_dir ):
         os.mkdir( cln_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      pl_QC       = "%s/bin/QC.pl" % ( home_dir )
      
      sh_file       = "%s/scripts/s01.QC.sh"      % (home_dir)
      sh_work_file  = "%s/scripts/s01.QC_work.sh" % (home_dir)
      
      sh_info = """
pl_exe=$1
pl_QC=$2
in_dir=$3
out_dir=$4
samp=$5
data_type=$6

$pl_exe $pl_QC --indir $in_dir --outdir $out_dir --sample $samp --end $data_type
      """
      
      sh_work = ""
      for samp in self['sam_info']['list']:
         if not os.path.isdir( "%s/%s" % (cln_dir,samp) ):
            os.mkdir(          "%s/%s" % (cln_dir,samp) )
         pl_exe    = self['sftw_name'].pl
         in_dir    = raw_dir
         out_dir   = cln_dir
         data_type = 2
         if self['sam_info']['sam_end'][samp ] == "SE":
            data_type = 1
         sh_work  += " sh %s  %s %s %s %s %s %d\n" % ( sh_file, pl_exe,  pl_QC, in_dir,out_dir,samp,data_type )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="400m",maxjob=100 )
      
      
   
   def run_bwa(self):
      home_dir     = os.path.abspath('./')
      
      cln_dir   = self['dir_name'].clean_data
      bam_dir   = self['dir_name'].bam
      
      if not os.path.isdir(  bam_dir):
         os.mkdir( bam_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      sh_file      = "%s/s02.bwa.sh"      % (script_dir)
      sh_work_file = "%s/s02.bwa_work.sh" % (script_dir)
      
      sh_info = """
bwa_exe=$1
end=$2
genome=$3
cln_dir=$4
bam_dir=$5
sam_name=$6
samtools_exe=$7

if [ $end == "1" ]
   then $bwa_exe aln -i 15 -q 10 -t 4 $genome $cln_dir/1.cln.fq.gz >$bam_dir/1.sai         && \\
        $bwa_exe samse $genome $bam_dir/1.sai $cln_dir/1.cln.fq.gz >$bam_dir/$sam_name.sam && \\
        $samtools_exe view -u -b -S -t $genome.fai $bam_dir/$sam_name.sam                   | \\
        $samtools_exe sort -m 200000000 - $bam_dir/$sam_name                               && \\
        rm $bam_dir/1.sai && rm $bam_dir/$sam_name.sam
fi
if [ $end == "2" ]
   then $bwa_exe aln -i 15 -q 10 -t 4 $genome $cln_dir/1.cln.fq.gz >$bam_dir/1.sai                                             && \\
        $bwa_exe aln -i 15 -q 10 -t 4 $genome $cln_dir/2.cln.fq.gz >$bam_dir/2.sai                                             && \\
        $bwa_exe sampe $genome $bam_dir/1.sai $bam_dir/2.sai $cln_dir/1.cln.fq.gz $cln_dir/2.cln.fq.gz >$bam_dir/$sam_name.sam && \\
        $samtools_exe view -u -b -S -t $genome.fai $bam_dir/$sam_name.sam                                                   | \\
        $samtools_exe sort -m 200000000 - $bam_dir/$sam_name                                               
        rm $bam_dir/1.sai $bam_dir/2.sai $bam_dir/$sam_name.sam
fi
      """ 
      sh_work = ""
      for samp in self['sam_info']['list']:
         bwa_exe      = self['sftw_name'].bwa
         samtools_exe = self['sftw_name'].samtools
         brief_name = self['sam_info']['samp_brief'][samp]
         
         sam_cln_dir = "%s/%s" % (cln_dir,samp)
         sam_bam_dir = "%s/%s" % (bam_dir,brief_name)
         if not os.path.isdir( sam_bam_dir ):
            os.mkdir(          sam_bam_dir )

         end = 2
         if self['sam_info']['sam_end'][samp ] == "SE":
            end = 1
         
         sh_work += "sh %s  %s %d  %s %s %s  %s %s \n" % ( sh_file, bwa_exe,end, self['infile']['genome_file'],sam_cln_dir,sam_bam_dir,  brief_name,  samtools_exe  )
               
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=3 )
#      my_job.running_SGE( vf="7g",maxjob=100 )

   def bam2repbed(self,ext_len=300):
      home_dir     = os.path.abspath('./')
      
      bam_dir    = self['dir_name'].bam
      bed_rep_dir= self['dir_name'].bed_rep
      bed_mrg_dir= self['dir_name'].bed_mrg
      
      java_exe   = self['sftw_name'].java
      python_exe = self['sftw_name'].py
      markDup_jar= self['sftw_name'].MarkDup
      
      samtools_exe= self['sftw_name'].samtools
      bedtools_exe= self['sftw_name'].bedtools
      
      if not os.path.isdir(  bed_rep_dir ):
         os.mkdir( bed_rep_dir )

      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      py_ExtRead   = "%s/ExtRead.py"      % (bin_dir)
      
      sh_file      = "%s/s03.bam2bedrep.sh"      % (script_dir)
      sh_work_file = "%s/s03.bam2bedrep_work.sh" % (script_dir)
      
      sh_info = """
java_exe=$1
markDup_jar=$2
python_exe=$3
py_ExtRead=$4
bam=$5
bed_rep=$6
ext_len=$7
end=$8
samtools_exe=${9}
bedtools_exe=${10}

$java_exe -Xmx1g -jar $markDup_jar INPUT=$bam.bam OUTPUT=$bam.dedup.bam METRICS_FILE=$bam.picard_info.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

$python_exe $py_ExtRead -l $ext_len -p $end $bam.dedup.bam $bed_rep

$samtools_exe view -b -F 1548 -q 30 $bam.dedup.bam | $bedtools_exe bamtobed -i - | awk 'BEGIN{FS="\\t";OFS="\\t"}{$4="N"; print $0}'  >$bed_rep.tagAlign
#rm $bam.bam
      """ 
      sh_work = ""
      for samp in self['sam_info']['list']:
         end = 2
         if self['sam_info']['sam_end'][samp ] == "SE":
            end = 1
         
         
         brief_name   = self['sam_info']['samp_brief'][samp]
         bam     = "%s/%s/%s" % (bam_dir    ,brief_name,brief_name)
         bed_rep = "%s/%s"    % (bed_rep_dir,brief_name )
         sh_work += "sh %s  %s %s  %s %s  %s %s  %d %d  %s %s \n" % ( sh_file, java_exe,markDup_jar, python_exe,py_ExtRead,  bam,bed_rep,  ext_len,end,  samtools_exe,bedtools_exe )
                                                                                                                                            
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="200m",maxjob=100 )


   def mrgbed(self):
      home_dir     = os.path.abspath('./')
      
      bed_rep_dir= self['dir_name'].bed_rep
      bed_mrg_dir= self['dir_name'].bed_mrg

      samtools_exe= self['sftw_name'].samtools
      bedtools_exe= self['sftw_name'].bedtools

      if not os.path.isdir(  bed_mrg_dir ):
         os.mkdir( bed_mrg_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      sh_file      = "%s/s04.1.bedmrg.sh"      % (script_dir)
      sh_work_file = "%s/s04.1.bedmrg_work.sh" % (script_dir)
      
      sh_sortbed   = "%s/sort_bed.sh"        % (bin_dir)
      bgzip_exe = self['sftw_name'].bgzip
      tabix_exe = self['sftw_name'].tabix
      
      sh_info = """
sh_sortbed=$1
bed_mrg=$2
bgzip_exe=$3
tabix_exe=$4
samtools_exe=$5
bedtools_exe=$6

shift
shift
shift

shift
shift
shift

for i in $@
do 
   cat $i.unique.bed
done >$bed_mrg.unique.bed

for i in $@
do
   cat $i.tagAlign
done >$bed_mrg.tagAlign

for i in $bed_mrg $@
do 
   sh $sh_sortbed $i.unique.bed $i.unique.sort.bed $bgzip_exe $tabix_exe
done

for i in $bed_mrg $@
do 
   sh $sh_sortbed $i.tagAlign $i.sort.tagAlign $bgzip_exe $tabix_exe
done

      """ 
      sh_work = ""
      for merge_sam in self['sam_info']['merge']:
         bed_mrg   = "%s/%s"    % (bed_mrg_dir,merge_sam )
         
         l_sam     = self['sam_info']['merge_sam'][merge_sam]
         l_brief   = [ "%s"           % (self['sam_info']['samp_brief'][sam]) for sam        in l_sam   ]
         l_bed_rep = [ "%s/%s"        % (bed_rep_dir,brief_name )             for brief_name in l_brief ]
         
         sh_work += "sh %s  %s %s  %s %s  %s %s   %s \n" % ( sh_file, sh_sortbed,bed_mrg, bgzip_exe,tabix_exe,  samtools_exe,bedtools_exe,  " ".join(l_bed_rep) )
               
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="200m",maxjob=100 )


   def mrgbed_multi(self):
      home_dir     = os.path.abspath('./')
      
      bed_rep_dir= self['dir_name'].bed_rep
      bed_mrg_dir= self['dir_name'].bed_mrg

      if not os.path.isdir(  bed_mrg_dir ):
         os.mkdir( bed_mrg_dir )
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      sh_file      = "%s/s04.2.bedmrg.multi.sh"      % (script_dir)
      sh_work_file = "%s/s04.2.bedmrg.multi_work.sh" % (script_dir)
      
      sh_sortbed   = "%s/sort_bed.sh"        % (bin_dir)
      bgzip_exe = self['sftw_name'].bgzip
      tabix_exe = self['sftw_name'].tabix
      
      sh_info = """
sh_sortbed=$1
bed_mrg=$2
bgzip_exe=$3
tabix_exe=$4
shift
shift
shift
shift

for i in $@
do 
   cat $i.bed
done >$bed_mrg.bed

for i in $bed_mrg $@
do 
   sh $sh_sortbed $i.bed $i.sort.bed $bgzip_exe $tabix_exe
done

      """ 
      sh_work = ""
      for merge_sam in self['sam_info']['merge']:
         bed_mrg   = "%s/%s.multi"    % (bed_mrg_dir,merge_sam )
         
         l_sam     = self['sam_info']['merge_sam'][merge_sam]
         l_brief   = [ "%s"           % (self['sam_info']['samp_brief'][sam]) for sam        in l_sam   ]
         l_bed_rep = [ "%s/%s.multi"  % (bed_rep_dir,brief_name )             for brief_name in l_brief ]
         
         sh_work += "sh %s  %s %s  %s %s  %s \n" % ( sh_file, sh_sortbed,bed_mrg, bgzip_exe,tabix_exe,  " ".join(l_bed_rep) )
               
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_multi( cpu=8 )
#      my_job.running_SGE( vf="200m",maxjob=100 )
