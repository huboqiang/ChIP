from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import scipy.stats
import time
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from scipy import stats

def generate_color(samples):
   colors = []
   Rfile = 'test.R'
   rfile = open(Rfile,"w")
   print >>rfile,'c<-rainbow(%d)\nwrite.table(c,"colors.xls",sep="\t",quote=FALSE)' % (len(samples))
   rfile.close()
   shell_info = 'Rscript %s' % (Rfile)
   p = subprocess.Popen(shell_info,shell='True')
   while 1:
      run_cnt = 0
      if p.poll() is None:
         run_cnt += 1
         time.sleep(1)
      if run_cnt == 0:
         break
   o_file = open('colors.xls',"r")
   colors = []
   in_h = o_file.readline()
   for line in o_file:
      line = line.strip('\n')
      f    = line.split()
      colors.append(f[1][0:7])
   shell_info = 'rm colors.xls test.R'
   p = subprocess.Popen(shell_info,shell='True')
   return colors

class PlotChIPRNA(object):
   def __init__( self, dir_name,TSS_promoter_up=5000,TSS_promoter_down=5000,width=500,bodybin=100 ):
      self.dir_name   = dir_name
      self.TSS_promoter_up   = TSS_promoter_up
      self.TSS_promoter_down = TSS_promoter_down
      self.width  =  width
      self.bodybin=  bodybin

      self.group_dens = {}

   def load_samp_class(self, merge_name):
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      
      infile =  "%s/%s.tss.up%d_down%d.Ranked_ChIP_Mean.xls" % ( Peak_mrg_TSS_Cor , merge_name, self.TSS_promoter_up, self.TSS_promoter_down )
      f_infile = open( infile,"r" )
      for line in f_infile:
         line = line.strip('\n')
         f    = line.split()
         group  = f[0]
         np_val = np.array(f[1:],dtype=float)
         self.group_dens[ group ] = np_val
      f_infile.close()
   
   def plot_class(self,merge_name, mapped_reads, read_len=101):
      Peak_mrg_TSS_Cor = self.dir_name.Peak_mrg_TSS_Cor
      
      len_up   = self.TSS_promoter_up
      len_down = self.TSS_promoter_down
      width    = self.width
      bodybin  = self.bodybin
      
      x = range(-int((len_up+1)/width),int((len_down+1)/width))
      x_lab = []
      x_pos = []
      for i in range(0,len(x)):
         num = x[i]*width
         if   i == 0:
            tag = '-5kb'
            x_lab.append(tag)
            x_pos.append(x[i])
         elif x[i] == 0:
            tag = 'TSS'
            x_lab.append(tag)
            x_pos.append(x[i])
         elif x[i] == x[-1]:
            tag = "+5kb"
            x_lab.append(tag)
            x_pos.append(x[i])
      
      
      x     = np.array(x)
      x_pos = np.array(x_pos)    
         
      fig = plt.figure(figsize=(5,5))
         
      ax = plt.subplot( 1, 1, 1 )      
      y_max = 0
      colors = generate_color( ['0','1','2'] )
      
      ax.set_title('%s' % (merge_name))
      
      xlab = 'Relative Distance(bp) with %d bp bins' % width
      
      print merge_name, read_len, mapped_reads
      
      y0 = self.group_dens[ "0" ] * width * 1000000 /  (read_len * mapped_reads)
      y1 = self.group_dens[ "1" ] * width * 1000000 /  (read_len * mapped_reads)
      y2 = self.group_dens[ "2" ] * width * 1000000 /  (read_len * mapped_reads)
###      y3 = self.samp_class[ samp ][ "3" ]              
###      ym = max( y0.max(),y1.max(),y2.max(),y3.max() )
      ym = max( y0.max(),y1.max(),y2.max() )
      if ym > y_max:
         y_max = ym

      ax.plot( x,y0, color=colors[0],linewidth=2,label="No   Exp" )
      ax.plot( x,y1, color=colors[1],linewidth=2,label="Low  Exp" )
      ax.plot( x,y2, color=colors[2],linewidth=2,label="High Exp" )

      ax.get_xaxis().set_ticks(x_pos)
      ax.get_xaxis().set_ticklabels(x_lab)
      
      ax.set_ylabel('Methylation Level')
      ax.set_xlabel(xlab)
         
###      plt.ylim( 0,1.0 )
      ax.legend(loc="upper right",fontsize=9)
      for tick in ax.xaxis.get_major_ticks():
         tick.tick1On = True
         tick.tick2On = False
      for tick in ax.yaxis.get_major_ticks():
         tick.tick1On = True
         tick.tick2On = False

      pdf_file = "%s/%s.tss.up%d_down%d.Ranked_ChIP_Mean.pdf" % ( Peak_mrg_TSS_Cor , merge_name, self.TSS_promoter_up, self.TSS_promoter_down )
      plt.savefig(pdf_file,format='pdf')

   def clean_class(self,merge_name):
      self.group_dens = {}
