#!/bin/bash 

cd targets

for file in miRDB/SC2V-mir*
do
  var1=${file%.c*}
  var2=${var1#*miRDB/}
  file2=Diana/70/${var2}_gene.tab
  awk 'BEGIN{FS=";"};{if ($1 >= 70 && NR > 1){print $3}}' $file | sort > aux
  awk 'BEGIN{FS="\t"}{print $3}' $file2 | sort > aux2
  join aux aux2 > overlap_$var2.tab

done
rm aux aux2


