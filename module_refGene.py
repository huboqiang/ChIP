from __future__ import division
import re,sys,os
import numpy as np
import subprocess
import time
from pybedtools import BedTool

def check_pos( chrom,beg,end,M_length ):
   if beg < 0:
      beg = 1
   if end > M_length[chrom]:
      end = M_length[chrom]-1
   return beg,end

def CpG_type(chrom,beg,end, in_fa,M_length):
   
   beg,end = check_pos(chrom,beg,end,M_length)
   seq = BedTool.seq( (chrom, beg, end),in_fa )
   
   seq = seq.upper()
   i       = 0;
   status  = "ICP";
   maxRcpg = 0;
   while  i+500 < len(seq) :
      tmp = seq[i:i+500]
      C  = tmp.count( "C"  )
      G  = tmp.count( "G"  )
      CpG= tmp.count( "CG" )
      
      Cgc  = (C+G)/500
      Rcpg = 0
      if C != 0 and G != 0:
         Rcpg = 500*CpG/(C*G)
      
      if Rcpg >= 0.75 and Cgc >= 0.55:
         status = "HCP"
         break
      else:
         if maxRcpg < Rcpg:
            maxRcpg = Rcpg
      i += 5
   
   if maxRcpg < 0.48:
      status = "LCP"
   return status

def sort_bed(  infile ):
   shell_info = """
bedtools sort -i %s >%s.tmp && mv %s.tmp %s
   """ % ( infile,infile,infile,infile )
   p = subprocess.Popen(shell_info,shell='True')
   while 1:
      run_cnt = 0
      if p.poll() is None:
         run_cnt += 1
         time.sleep(1)
      if run_cnt == 0:
         break
      
   

class refGeneLine(object):
   def __init__(self,line):
      f = line.split()
      self.bin	         = int(f[0])
      self.name	      =     f[1]
      self.chrom	      =     f[2]
      self.strand       =     f[3]
      self.txStart      = int(f[4])
      self.txEnd        = int(f[5])
      self.cdsStart     = int(f[6])
      self.cdsEnd       = int(f[7])
      self.exonCount    = int(f[8])
      self.exonStarts   = np.array( f[ 9].split(",")[:-1],dtype=int )
      self.exonEnds     = np.array( f[10].split(",")[:-1],dtype=int )
      self.score        = int(f[11])
      self.name2        =     f[12]
      self.cdsStartStat =     f[13]
      self.cdsEndStat   =     f[14]
      self.exonFrames   =     f[15]

class RefGeneTxt(object):
   def __init__(self,infile,genomeFa):
      self.infile   = infile
      self.genomeFa = genomeFa
      self.length   = {}
      self.__load_length()
      
      self.prefix = ".".join( self.infile.split(".")[:-1] )
      self.infile_sort = "%s.txt" % ( self.prefix )
      if not os.path.isfile( self.infile_sort ):
         self.__sort_refgene()
      
   def __sort_refgene(self):
      shell_info = """
sort -k1V -k2n -k3n %s >%s
      """ % ( self.infile, self.infile_sort )
      p = subprocess.Popen(shell_info,shell='True')
      while 1:
         run_cnt = 0
         if p.poll() is None:
            run_cnt += 1
            time.sleep(1)
         if run_cnt == 0:
            break
            
   def __load_length(self):
      f_inFai = open( "%s.fai" % (self.genomeFa),"r")
      for line in f_inFai:
         line = line.strip('\n')
         f    = line.split()
         self.length[ f[0] ] = int( f[1] )
      f_inFai.close()
      
   def refGene2bed(self,up_stream,down_stream,ext_type=""):
      self.refGene_ext_bed = "%s.up%d_down%d.%sBsorted.bed" % ( self.prefix, up_stream, down_stream, ext_type )
      f_infile_sort     =  open( self.infile_sort,"r" )
      f_refGene_ext_bed =  open( self.refGene_ext_bed,"w" )
      
      for line in f_infile_sort:
         line = line.strip('\n')
         line_info = refGeneLine( line )
         chrom = line_info.chrom
         
         beg   = line_info.txStart - up_stream
         end   = line_info.txEnd   + down_stream
         if line_info.strand == "-":
            beg   = line_info.txStart   - down_stream
            end   = line_info.txEnd     + up_stream
         
         if ext_type == "tss." or ext_type == "promoter.":
            beg   = line_info.txStart - up_stream
            end   = line_info.txStart + down_stream
            if line_info.strand == "-":
               beg   = line_info.txEnd - down_stream
               end   = line_info.txEnd + up_stream
            
         
         strand = line_info.strand
         tid    = line_info.name
         ltype  = "protein_coding"
         if tid[0:2] == "NR":
            ltype  = "noncoding"
         gid    = line_info.name2

         out    = "%s\t%d\t%d\t%s\t%s\t%s\t%s"  % ( chrom,beg,end, strand,tid,ltype,gid )
         print >>f_refGene_ext_bed, out
      f_infile_sort.close()
      f_refGene_ext_bed.close()
   
   def refGeneInfo(self,TSS_up_stream=1000,TSS_down_stream=500):
      self.tranInfo = {}
      self.geneInfo = {}
      
      f_infile_sort     =  open( self.infile_sort,"r" )
      
      for line in f_infile_sort:
         line = line.strip('\n')
         line_info = refGeneLine( line )
         chrom = line_info.chrom
         beg   = line_info.txStart
         end   = line_info.txEnd

         strand = line_info.strand
         tid    = line_info.name
         ltype  = "protein_coding"
         if tid[0:2] == "NR":
            ltype  = "noncoding"
         gid    = line_info.name2
         
         tss_beg     = line_info.txStart - TSS_up_stream
         tss_end     = line_info.txStart + TSS_down_stream
         if line_info.strand == "-":
            tss_beg  = line_info.txEnd   - TSS_down_stream
            tss_end  = line_info.txEnd   + TSS_up_stream
            
         tid_type = CpG_type( chrom,tss_beg,tss_end,  self.genomeFa, self.length )

         print chrom,beg,end,tid_type,gid,tid

         if gid not in self.geneInfo:
            self.geneInfo[ gid ] = { 'tran':[],  'line':{},  'max_len':0, 'max_len_tid':"None",'max_len_tid_type':""  }
         
         self.geneInfo[ gid ][ 'tran' ].append( tid )
         self.geneInfo[ gid ][ 'line' ][ tid ] = line
         if end-beg > self.geneInfo[ gid ][ 'max_len' ]:
            self.geneInfo[ gid ][ 'max_len'     ] = end-beg
            self.geneInfo[ gid ][ 'max_len_tid' ] = tid
            self.geneInfo[ gid ][ 'max_len_tid_type' ] = tid_type
            
         if tid not in self.tranInfo:
            self.tranInfo[ tid ] = { 'gene':gid, 'tmp_len':end-beg, 'exon_cnt':len(line_info.exonStarts), 'beg':beg, 'end':end, 'chr':chrom,\
                 'exon_pos':{ 'left':line_info.exonStarts,'right':line_info.exonEnds }, 'strand':strand,'tid_type':tid_type }

      f_infile_sort.close()
            
   def Only_LongestTid_Bed(self,up_stream,down_stream,ext_type=""):
      infile = "%s.up%d_down%d.%sBsorted.bed"            % ( self.prefix, up_stream, down_stream, ext_type )
      otfile = "%s.up%d_down%d.%sBsorted.longestTid.bed" % ( self.prefix, up_stream, down_stream, ext_type )
      f_infile = open( infile,"r" )
      f_otfile = open( otfile,"w" )
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         tid  = f[4]
         gid  = f[6]
         if tid == self.geneInfo[ gid ][ 'max_len_tid' ]:
            print >>f_otfile, "%s\t%s" % ( line,self.tranInfo[tid]['tid_type'] )
      f_infile.close()
      f_otfile.close()
      sort_bed( otfile )
      
