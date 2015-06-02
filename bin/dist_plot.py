from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy   as np
import pandas  as pd

import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib        import pyplot as plt
import matplotlib as mpl

from optparse          import OptionParser
import pylab as P
from scipy.stats.kde   import gaussian_kde
from scipy.interpolate import spline
from sklearn           import linear_model

sns.set_style("white")

class PlotDens(object):

   def __init__(self,in_PeakFile, in_DensFile, out_prefix):
      self.in_PeakFile = in_PeakFile
      self.in_DensFile = in_DensFile
      self.outXls      = "%s.xls" % out_prefix
      self.outPdf      = "%s.pdf" % out_prefix
      self.outLgstPdf  = "%s.logestics.pdf" % out_prefix
      
      self.x_lim2 = 1
      self.x_lim1 = 1
      self.y_lim = 0.00005
      
   def load_dens( self ):

      pd_peak = pd.read_csv( self.in_PeakFile  ,sep="\t" )
      pd_dens = pd.read_csv( self.in_DensFile  ,sep="\t" )
      
      pd_dens_Peak = pd_dens[ pd_peak == 1]
      pd_dens_Norm = pd_dens[ pd_peak == 0]
      
      
      a  = pd.melt( pd_dens_Norm ,id_vars=['#chr','beg','end'])
      a1 = a.drop( a.columns[0:3],axis=1 )
      self.np_norm = a1[ pd.notnull(a1['value']) ]['value'].values
      
      a  = pd.melt( pd_dens_Peak ,id_vars=['#chr','beg','end'])
      a1 = a.drop( a.columns[0:3],axis=1 )
      self.np_peak = a1[ pd.notnull(a1['value']) ]['value'].values

   def calc_dens(self):      
      kde_np_peak  = gaussian_kde( self.np_peak )
      self.np_peak[ self.np_peak>self.x_lim2 ] = self.x_lim2
      
      kde_np_norm  = gaussian_kde( self.np_norm,bw_method=0.05 /  self.np_norm.std(ddof=1) )
      self.np_norm[ self.np_norm>self.x_lim2 ] = self.x_lim2
      
      dist_space = np.linspace( 0, self.x_lim2, 100 )
      x      = dist_space
      y_peak = kde_np_peak(dist_space)
      y_norm = kde_np_norm(dist_space)
      
      data = { 'x':x , 'y_peak':y_peak,  'y_norm':y_norm }
      pd_frame = pd.DataFrame( data )
      pd_frame.to_csv( self.outXls,sep="\t" )
      

   def plot_dens(self):
      pd_frame = pd.read_csv( self.outXls,sep="\t" )


      pd_peak = pd.read_csv( self.in_PeakFile  ,sep="\t",index_col=['#chr','beg','end'] )
      pd_dens = pd.read_csv( self.in_DensFile  ,sep="\t",index_col=['#chr','beg','end'] )
      
      x        = pd_frame[ 'x' ]
      y_peak   = pd_frame[ 'y_peak' ]
      y_norm   = pd_frame[ 'y_norm' ]

      fig = plt.figure(figsize=(10,8))
      ax = fig.add_subplot(1,1,1)
      
      self.np_peak[ self.np_peak>self.x_lim2 ] = self.x_lim2
      self.np_norm[ self.np_norm>self.x_lim2 ] = self.x_lim2
      
      n, bins, patches = ax.hist(self.np_peak, 25, normed=1, histtype='stepfilled')
      plt.setp(patches, 'facecolor', 'r', 'alpha', 0.25)    
      ax.plot( x,y_peak,"r",label="H3K27ac Peak    read-density" )
      
      n, bins, patches = ax.hist(self.np_norm, 25, normed=1, histtype='stepfilled')
      plt.setp(patches, 'facecolor', 'b', 'alpha', 0.25)
      ax.plot( x,y_norm,"b",label="H3K27ac no-Peak read-density" )
      


      pd_statDens = pd.DataFrame( {"stat":pd_peak[ pd_peak.columns[0] ], "dens":pd_dens[ pd_peak.columns[0] ]} )
      X = pd_statDens['dens'].values
      y = pd_statDens['stat'].values
      X[ X>self.x_lim2 ] = self.x_lim2
      X = np.array([X]).T
      X_test = np.linspace( 0,self.x_lim2,300 )

      clf = linear_model.LogisticRegression(solver='lbfgs')
      clf.fit( X,y )
            
      def model(x):
         return 1 / (1 + np.exp(-x))
         
      loss = model(X_test * clf.coef_ + clf.intercept_).ravel()
      plt.plot(X_test, loss, color='black', linewidth=3)
      
      idx_2 = X_test[ np.where( loss==loss[ loss>0.2 ][0] )[0] ]
      idx_3 = X_test[ np.where( loss==loss[ loss>0.3 ][0] )[0] ]
      idx_4 = X_test[ np.where( loss==loss[ loss>0.4 ][0] )[0] ]
      idx_5 = X_test[ np.where( loss==loss[ loss>0.5 ][0] )[0] ]
      idx_6 = X_test[ np.where( loss==loss[ loss>0.6 ][0] )[0] ]
      idx_7 = X_test[ np.where( loss==loss[ loss>0.7 ][0] )[0] ]
      idx_8 = X_test[ np.where( loss==loss[ loss>0.8 ][0] )[0] ]
      
      print idx_2, idx_3, idx_4, idx_5, idx_6, idx_7, idx_8 
      
      ax.set_title( "Reads density distribution", fontsize=16)
      ax.set_ylabel("Probability Density",size=16)

      x_cutoff = 0.0
      y_cutoff = 0.0
      for i,x_val in enumerate(x):
         if y_norm[i] > y_peak[i] and y_norm[i+1] <= y_peak[i+1]:
            x1  = x[i]
            x2  = x[i+1]
            y11 = y_norm[i]
            y12 = y_norm[i+1]
            y21 = y_peak[i]
            y22 = y_peak[i+1]
            x_cutoff = ( (y11-y21)*x2+(y22-y12)*x1 )/( (y11-y21)+(y22-y12) )
            y_cutoff = ( (x2-x_cutoff)*y11  +(x_cutoff-x1)*y12   )/(     x2   -    x1    )
            break


      ax.plot( [x_cutoff,x_cutoff],[y_cutoff,0],"k")
      ax.text( 0.5, 5,    "x = %1.4f" % ( x_cutoff ) )
      
      self.x_cutoff = x_cutoff
      
      ax.set_xlim(0,self.x_lim1)
      for tick in ax.xaxis.get_major_ticks():
         tick.tick1On = True
         tick.tick2On = False
         tick.label.set_fontsize(12) 
      
      
      for tick in ax.yaxis.get_major_ticks():
         tick.tick1On = True
         tick.tick2On = False
      
      
      plt.savefig( self.outPdf, format='pdf')



   def recal_Peak( self ):
      self.out_PeakFileRev = "%s.rev.xls" % ( ".".join( self.in_PeakFile.split(".")[:-1] ) )
      pd_peak  = pd.read_csv( self.in_PeakFile  ,sep="\t",index_col=['#chr','beg','end'] )
      pd_dens  = pd.read_csv( self.in_DensFile  ,sep="\t",index_col=['#chr','beg','end'] )
      pd_peak[ pd_dens >= self.x_cutoff ] = 1
      pd_peak[ pd_dens <  self.x_cutoff ] = 0
      
      pd_peak.to_csv( self.out_PeakFileRev, sep="\t" )
      
      



   def plot_logestic(self):
      pd_peak = pd.read_csv( self.in_PeakFile  ,sep="\t",index_col=['#chr','beg','end'] )
      pd_dens = pd.read_csv( self.in_DensFile  ,sep="\t",index_col=['#chr','beg','end'] )
      
      pd_statDens = pd.DataFrame( {"stat":pd_peak[ pd_peak.columns[0] ], "dens":pd_dens[ pd_peak.columns[0] ]} )
      
      X = pd_statDens['dens'].values
      y = pd_statDens['stat'].values

      
      X[ X>self.x_lim2 ] = self.x_lim2
         
      X = np.array([X]).T

      clf = linear_model.LogisticRegression(solver='lbfgs')
      clf.fit( X,y )
      
      
      fig = plt.figure(figsize=(10,8))
      ax = fig.add_subplot(1,1,1)
      
      ax.scatter( X,y,color='black' )
      X_test = np.linspace( 0,self.x_lim2,300 )
      
         
      def model(x):
         return 1 / (1 + np.exp(-x))
         
      print clf.coef_, clf.intercept_
         
      loss = model(X_test * clf.coef_ + clf.intercept_).ravel()
      plt.plot(X_test, loss, color='blue', linewidth=3)
         
      ax.set_xlim(0,self.x_lim1)
      for tick in ax.xaxis.get_major_ticks():
         tick.tick1On = True
         tick.tick2On = False
         tick.label.set_fontsize(12) 
      
      
      for tick in ax.yaxis.get_major_ticks():
         tick.tick1On = True
         tick.tick2On = False
      
      
      plt.savefig( self.outLgstPdf, format='pdf')
      
      


def prepare_optparser():
   usage ="""usage: %s [options]

   Using -h or --help for more information

Example:
   python %s  Merge_1kb_Peaks.enhancer.xls  Merge_1kb_Density.enhancer.xls  MarkerDensity

   """ % (sys.argv[0],sys.argv[0])

   description = " Select reads while converting bam to bed "

   optparser = OptionParser(version="%s v0.2 20150124" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   
   optparser.add_option("-r", "--region",   help="\nInput intervals.")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser
   
   
def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      in_PeakFile = args[0]
      in_DensFile = args[1]
      out_prefix  = args[2]

   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)

   m_dens = PlotDens( in_PeakFile, in_DensFile, out_prefix )
   m_dens.load_dens()
#   m_dens.calc_dens()
   m_dens.plot_dens()
#   m_dens.recal_Peak()
#   m_dens.plot_logestic()
   
if __name__ == '__main__':
   main()
