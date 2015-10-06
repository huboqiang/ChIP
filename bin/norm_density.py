# -*- coding:utf-8 -*-
# 这个模块用于和系统交互
import sys
import gzip
import MySQLdb
from optparse import OptionParser

__author__      = "Shuhui Bian"
__copyright__   = "Copyright 2015, Biopic, Peking University"


def prepare_optparser():
    usage ="""usage: %s [options] input.sort.bdg.gz


Using -h or --help for more information

Example:
    python %s -r hg19 /datd/huboqiang/ChIP_human/Week12/03.1.Peak_rep/Adult_brain_H3K4me3_rep1/Adult_brain_H3K4me3_rep1_VS_Input_control_lambda.sort.bdg.gz

    """ % (sys.argv[0], sys.argv[0])
    description = "Normalize bedGraph file."
    optparser = OptionParser(
        version="%s v0.1 2015.10" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )
    optparser.add_option("-r", "--ref", default="hg19",  help="\nreference of UCSC.")
    return optparser

def load_chrlen(ref):
    M_chr_len = {}
    db = MySQLdb.connect(host="genome-mysql.cse.ucsc.edu", db=ref, user="genome")
    c = db.cursor()
    c.execute("""SELECT chrom,size FROM chromInfo""")
    while 1:
        result = c.fetchone()
        if result != None:
            inf1 = result[0]
            inf2 = result[1]
            M_chr_len[inf1] = inf2
        else:
            break
        
    return M_chr_len


def parse_line(line):
    line = line.strip()
    f = line.split()
    return float(f[3])


def div_total(line, total_val):
    line = line.strip()
    f = line.split()
    normalized_value = float(f[3])/total_val
    return normalized_value



def main():
    prepare_optparser()
    (options, args) = prepare_optparser().parse_args()
    try:
        infile = args[0]
        ref = options.ref
    
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)

    M_chr_len = load_chrlen(ref)
    out_file1 = "%s.norm.bedGraph" % (".".join(infile.split(".")[:-2]))  # 定义输出的文件名称
    f_outFile1 = open(out_file1, "w")  # 打开输出文件

    total_val = 0.0
    with gzip.open(infile, "rb") as f_infile:
        for line in f_infile:
            line = line.strip()
            total_val += parse_line(line)
    
    with gzip.open(infile, "rb") as f_infile:
        for line in f_infile:
            line = line.strip()
            f = line.split()
            beg_site = int(f[1])
            end_site = int(f[2])
            chr_len = M_chr_len[f[0]]
           
            if end_site < chr_len:                 
                normalized_value = div_total(line, total_val)
                out = "%s\t%d\t%d\t%f" % (f[0], beg_site, end_site, normalized_value)
                print >>f_outFile1, out
            
    
            elif beg_site < chr_len:
                normalized_value = div_total(line, total_val)       
                out = "%s\t%d\t%d\t%f" % (f[0], beg_site, chr_len, normalized_value)
                print >>f_outFile1, out
            
            
           
    f_outFile1.close()



if __name__ == "__main__":
    main()
