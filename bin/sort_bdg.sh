inbed=$1
outbed=$2

for i in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X Y M
do 
	grep -w chr$i $inbed | sort -k 2n
done | bgzip -fc >$outbed.gz && tabix -f -p bed -s 1 -b 2 -e 3 -S 0 $outbed.gz