A pipeline which could processing from raw fastq reads to peak analysis

First, before this pipeline in a server, make sure the required modules were installed. If not, running the following scripts for deploying.
```bash
mkdir install_packages

### install python anaconda 2.2.0
cd software/
wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda-2.2.0-Linux-x86_64.sh
bash Anaconda-2.2.0-Linux-x86_64.sh  # prefix=/path/for/anaconda
mv Anaconda-2.2.0-Linux-x86_64.sh install_packages

### install R 3.2.0
wget http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz
tar -zxvf R-3.2.0.tar.gz
cd R-3.2.0
./configure --prefix ~/software/R-3.2.0
make
make install
cd ..
mv R-3.2.0.tar.gz install_packages

### install samtools 0.1.18
### using old version because the latest one could have somewhat trouble with 
### other software like tophat.
wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
tar -jxvf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make
cd ..
mv samtools-0.1.18.tar.bz2 install_packages

### install bwa 0.7.5a
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.5a.tar.bz2
tar -jxvf bwa-0.7.5a.tar.bz2
cd bwa-0.7.5a
make
cd ..
mv bwa-0.7.5a.tar.bz2 install_packages

### install bedtools 2.24.0
git clone https://github.com/arq5x/bedtools2/
cd bedtools2
make

### install MACS2
wget https://pypi.python.org/packages/source/M/MACS2/MACS2-2.1.0.20150731.tar.gz
tar -zxvf MACS2-2.1.0.20150731.tar.gz
cd MACS2/
# sed -i 's/Ofast/O3/g' /data/Analysis/huboqiang/software/MACS/setup_w_cython.py if necessary
python setup_w_cython.py install

### install picard-tools
wget http://jaist.dl.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip
unzip picard-tools-1.119.zip

### install tabix and bgzip
wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
tar -jxvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
cd ..
mv tabix-0.2.6.tar.bz2 install_packages

### install igvtools
wget https://data.broadinstitute.org/igv/projects/downloads/igvtools_2.3.57.zip?
unzip igvtools_2.3.57.zip

### install idrcode and spp_package
wget http://www.broadinstitute.org/~anshul/softwareRepo/idrCode.tar.gz
wget http://www.broadinstitute.org/~anshul/softwareRepo/spp_package.tar.gz
tar -zxvf idrCode.tar.gz
tar -zxvf spp_package.tar.gz
cd spp_package
tar -zxvf spp_1.10.1.tar.gz

### install UCSC utilities
from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/


### install required python modules:
pip install ngslib
pip install Mysql
pip install svgwrite
pip install seaborn
pip install pysam==0.8.3
pip install pybedtools==0.6.9
pip install sklearn
```

After that, download this script:
```
cd $PYTHONPATH  # path for put the python packages. path/to/anaconda/lib/python2.7/site-packages/ for default
git clone https://github.com/hubqoaing/ChIP
```


Secondly, go to the ./setting file, and change the following values to your own path:
```python
self.Database            = "DIR/TO/DATABASE"                               #line 61
self.sftw_py             = "DIR/TO/SOFTWARE_EXE_FILE"  #line 77
self.sftw_pl             = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_java           = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_R              = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_MarkDup        = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bwa            = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_samtools       = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_macs2          = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_macs14         = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bedtools       = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_bgzip          = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_tabix          = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_igvtools       = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_batchIDR       = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_spp            = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_get_psudoCount = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_sort_bdg       = "DIR/TO/SOFTWARE_EXE_FILE"
self.sftw_ucsc_dir       = "DIR/TO/SOFTWARE_EXE_FILE"

```

Go to the analysis dictionary and copy the bin file here.
``` bash
cd PATH/FOR/ANALYSIS   # go to 
copy $PYTHONPATH/ChIP/run_chipseq.py ./
```

Next, make the input files. You can download these files in UCSC or so on and then using own-scripts to merge the ERCC information, and generate files in this format.###TAB split

``` bash
vim sample_input.xls
==> sample_input.xls <==
sample               stage   type        tissue      brief_name               merge_name             end_type   control
H3K27ac_mE105_brain  E105    H3K27ac     brain       H3K27ac_mE105_brain      H3K27ac_mE105_brain     SE        Input_mE105_brain
Input_mE105_brain1   E105    H3K27ac     brain       Input_mE105_brain1       Input_mE105_brain       SE        Input_mE105_brain
...
```
Notice that only NAME\_FOR\_RAW\_FQ were required that this NAME should be the same as 00.0.raw\_fq/NAME.
NAME\_FOR\_PROCESSING will be the name for the rest analysis's results.
NAME\_FOR\_READING    will be the name for files in statinfo.
stage and sample_group could be writen as anything. It was here only for make the downstream analysis easily.

Before running this pipeline, put the fastq reads in the ./00.0.raw_data dictionary.
```bash
mkdir 00.0.raw_data
for i in `tail -n +2 sample_input.xls | awk '{print $1}`
do
    mkdir 00.0.raw_data/$i && ln -s PATH/TO/RAW_DATA/$i/*gz 00.0.raw_data/$i
done
```

After that, running this pipeline:
```bash
python run_chipseq.py --ref YOUR_REF --TSS_genebody_up 5000 --TSS_genebody_down 5000 --TSS_promoter_up 5000  --TSS_promoter_down 5000 --Body_extbin_len 50 --Body_bincnt 100 --TSS_bin_len 1 --top_peak_idr 100000 sample_input.xls 
```

Wait for the results. 
Notice if you have to run it in a cluster, please do not running this scripts directly.
For example, if SGE system used, then:

Comments this command
```python
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
```
and using this command in modules in ./frame/*py
```
       my_job.running_SGE(vf="400m", maxjob=100, is_debug = self.is_debug)
```
Method for submit jobs in other system were still developing.
