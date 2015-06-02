fileName=$1

nlines=$( zcat ${fileName} | wc -l ) # Number of reads in the tagAlign file
nlines=$(( (nlines + 1) / 2 )) # half that number
zcat "${fileName}" | shuf | split -d -l ${nlines} - $fileName.div # This will shuffle the lines in the file and split it into two parts
bgzip -f $fileName.div00
bgzip -f $fileName.div01
mv $fileName.div00.gz $fileName.pr1.tagAlign.gz && tabix -f -p bed -s 1 -b 2 -e 3 -S 0 $fileName.pr1.tagAlign.gz
mv $fileName.div01.gz $fileName.pr2.tagAlign.gz && tabix -f -p bed -s 1 -b 2 -e 3 -S 0 $fileName.pr2.tagAlign.gz