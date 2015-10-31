#-*- coding:utf-8 -*-
import sys

class IntervalBed(object):
    "IntervalBed for a Queue structure of reading bed file."""
    def __init__(self, infile):
        self.chr = ""
        self.beg = 0
        self.end = 0
        self.val = 0.0
        self.infile = infile

    def read_bed(self):
        with open(self.infile, "r") as f_infile:
            for line in f_infile:
                self.__parse_line(line)
                if self.chr != self.line_chrom:
                    if len(self.chr) == 0:
                        # 初始化
                        self.__new_interval()
                    else:
                        # 换染色体
                        # 打印
                        self.__print_interval()
                        # 清空
                        self.__clean_interval()
                        # 初始化
                        self.__new_interval()
                    
                else:
                    if self.line_val == self.val:
                        # 延长
                        self.__elong_interval()
                    else:
                        # 打印
                        self.__print_interval()
                        # 清空
                        self.__clean_interval()
                        # 初始化
                        self.__new_interval()
        
            # 打印没被打印的最后一个interval
            self.__print_interval()
            # 清空
            self.__clean_interval()

    def __parse_line(self, line):
        line = line.strip()
        f = line.split()
        self.line_chrom = f[0]
        beg_1 = int(f[1])
        end_1 = int(f[2])
        beg_2 = int(f[5])
        end_2 = int(f[6])
    
        self.line_beg = max(beg_1, beg_2)
        self.line_end = min(end_1, end_2)
        self.line_val =  float(f[3]) - float(f[7])
        if self.line_val < 0:
            self.line_val = 0.0
    
    def __new_interval(self):
        self.chr = self.line_chrom
        self.beg = self.line_beg
        self.end = self.line_end
        self.val = self.line_val

    def __elong_interval(self):
        self.end = self.line_end

    def __print_interval(self):
        print "%s\t%d\t%d\t%1.6f" % (self.chr, self.beg, self.end, self.val)
    
    def __clean_interval(self):
        self.chr = ""
        self.beg = 0
        self.end = 0
        self.val = 0.0


        

#v3.2 简化类的使用 
infile = sys.argv[1]
Itv_bed = IntervalBed(infile)
Itv_bed.read_bed()
