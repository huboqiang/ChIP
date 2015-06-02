input=$1
out=$2
bgzip_exe=$3
tabix_exe=$4

#bgzip_exe=/usr/local/bin/bgzip
#tabix_exe=/usr/local/bin/tabix

sort  -k1V -k2n -k3n $input >$out

$bgzip_exe -f $out && $tabix_exe -f -p bed -s 1 -b 2 -e 3 -S 0 $out.gz
