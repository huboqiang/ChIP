sh /data/Analysis/huboqiang/project/human_sc/chIp_analysis/scripts/s11.PeakGeneRegion.sh   /data/Analysis/huboqiang/software/bedtools-2.17.0/bin/bedtools /data/Analysis/huboqiang/project/human_sc/chIp_analysis/03.2.Peak_mrg/mESC_NChIP_H3K27ac_1K/mESC_NChIP_H3K27ac_1K_VS_Input_treat_pileup.sort.bdg.gz /data/Analysis/huboqiang/project/human_sc/chIp_analysis/03.3.Peak_idr/mESC_NChIP_H3K27ac_1K/mESC_NChIP_H3K27ac_1K_VS_Input_peaks.conservative.regionPeak.gz /datc/huboqiang/test_chipPipe/Database/refGene.up5000_down5000.genebody.Bsorted.longestTid.bed /data/Analysis/huboqiang/project/human_sc/chIp_analysis/bin/get_bin/get_bin   5000 5000 5000 5000   50 100 1   /data/Analysis/huboqiang/project/human_sc/chIp_analysis/05.2.Peak_mrg.density.Genebody/mESC_NChIP_H3K27ac_1K/mESC_NChIP_H3K27ac_1K.genebody.up5000_down5000.xls /data/Analysis/huboqiang/project/human_sc/chIp_analysis/05.1.Peak_mrg.density.TSS/mESC_NChIP_H3K27ac_1K/mESC_NChIP_H3K27ac_1K.tss.up5000_down5000.xls 

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

zcat $peak_bed | cut -f 1-5 >$peak_bed.bed

zcat /data/Analysis/huboqiang/project/human_sc/chIp_analysis/03.2.Peak_mrg/mESC_NChIP_H3K27ac_1K/mESC_NChIP_H3K27ac_1K_VS_Input_treat_pileup.sort.bdg.gz | tail -n +2 | bedtools intersect -sorted -wo -a /dev/stdin -b  /datc/huboqiang/test_chipPipe/Database/refGene.up5000_down5000.genebody.Bsorted.longestTid.bed | /data/Analysis/huboqiang/project/human_embryo_sequencing/07.CHIP-seq/get_bin_ignore_peak/get_bin -U 5000 -D 5000 -T 5000 -b 50 -B 100 /dev/stdin stat1_genebody stat2_tss
      
