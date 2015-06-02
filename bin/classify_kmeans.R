args<-commandArgs(T)
if(length(args)<1)
{
  cat("Rscript classify_kmeans.R Merge_1kb_PeaksRev.enhancer.Rinput.xls\n")
  q()
}

join_str<-function(in_vec,join_str){
  out_str<-in_vec[1]
  if (length(in_vec) > 1){
    for (i in 2:length(in_vec))
    {
      out_str = paste(out_str,in_vec[i],sep=join_str  )  
    }
  }
  return(out_str)
}

z_score <- function(vec){
  ( vec-mean(vec) )/( sd(vec) )
}


get_args_infoV1<-function(args_in){
  split.arg<-strsplit(args_in[1],'/',perl=TRUE)
  prefix_strs<-split.arg[[1]][1:length(split.arg[[1]])]
  file_Info<-NULL
  file_Info$prefix<-join_str( prefix_strs[1:length(split.arg[[1]])-1],"/" )
  if ( is.na(file_Info$prefix) ){
    file_Info$prefix<-c("./")
  }
  file_Info$infile<-tail(prefix_strs,n=1)
  split.out<-strsplit(file_Info$infile,'\\.',perl=TRUE)
  prefix_out_strs<-split.out[[1]][1:length(split.out[[1]])]
  file_Info$outprefix<-join_str( prefix_out_strs[1:length(split.out[[1]])-1],"." )
  return( file_Info )
}

give.n <- function(x){
  return(c(y = median(x)+0.05, label = round(median(x),2)  ) ) 
}

callTag<-function( all,mat_tag, cor_min ){
  group_idx<-rep(0,dim(all)[1])
  for ( i in 1:dim(all)[1]){
    max_corr = -1
    max_j_idx= 0
    if ( sd(all[i,],na.rm=TRUE) != 0 & !is.na(sd(all[i,],na.rm=TRUE) ) ){
      for ( j in 1:dim(mat_tag)[1] ){
        corr = cor( as.numeric(all[i,]),as.numeric(mat_tag[j,]),use = "complete.obs",method="pearson")
        if ( !is.na(corr)  ){
          max_corr = corr
          max_j_idx = j
        }
      }
      if (max_corr > cor_min  ){
        group_idx[i]<-max_j_idx 
      }
    }
  }
  return( group_idx )
}

reorder_df<-function( all,group_idx,file_Info,tag ){
  all$idx<-group_idx
  o<-order( all$idx )
  all<-all[o,]
  index<-all$idx
  all.output<-all[ ,c(dim(all)[2],1:(dim(all)[2]-1)) ]
  all.output$gene<-row.names( all.output )
  all.output<-all.output[ ,c(dim(all.output)[2],1:(dim(all.output)[2]-1)) ]
  all$idx<-NULL
  
  out_xls<-paste(file_Info$outprefix, tag,"xls",sep=".")
  write.table( all.output, file=out_xls,sep="\t",quote = FALSE,row.names=F  )
  out_pdf<-paste(file_Info$outprefix, tag,"pdf",sep=".")
  pdf(out_pdf,width=20,height=20)
  pheatmap(  log10(1+all), cluster_rows=F,cluster_cols=F, color=color.palette,border_color=NA,show_rownames =F)
  dev.off()
  
  return(all)
}

library(plyr)
library(grid)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(gplots)

file_Info<-get_args_infoV1(args[1])
setwd(file_Info$prefix)
data<-read.table( file_Info$infile,header = TRUE )

row.names( data )<-data$pos
data$pos<-NULL

color.palette = colorRampPalette(c("green","black","red"))(100)

mat<-as.matrix(data)
km<-kmeans( mat,10,nstart=40 )
group_idx_stage<-as.numeric(km$cluster)
mat  <-reorder_df( as.data.frame(mat)  ,group_idx_stage,file_Info,"Plot_Kmeans"   )