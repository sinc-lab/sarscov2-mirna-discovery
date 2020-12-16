#!/bin/bash 
cd targets
mkdir -p targetsDE

for file in overlap_SC2V*
do
  var1=${file%.t*}
  awk 'BEGIN{FS="\t"};FNR==NR{a[$1];next}($1 in a){print $1}' $file ../expressionData/DEgenes_SARS-CoV-2_Full_FC1.txt > aux2
  sort aux2 | uniq > targetsDE/overlapDE_$var1.tab
  rm aux2
  awk 'BEGIN{FS=","};{if ((($2 < 0 && $5=="TRUE") || ($3<0 && $6=="TRUE") || ($4 <0 && $7=="TRUE")) && NR > 1){print $1}}' ../expressionData/DE_SARS-CoV-2_Full_FC1.csv | sort > aux2
  sort $file >aux1
  join aux1 aux2 | uniq > targetsDE/overlapDOWN_$var1.tab
  rm aux1 aux2
done 

cat targetsDE/overlapDE* | sort | uniq > targetsDE/targetsDE.txt

cat targetsDE/overlapDOWN* | sort | uniq > targetsDE/targetsDEDOWN.txt

#generating targetsDE file

awk 'BEGIN{FS=",";OFS=","};FNR==NR{a[$1];next}($1 in a){print $0}' targetsDE/targetsDE.txt ../expressionData/DE_SARS-CoV-2_Full_FC1.csv > targetsDE/TargetsDEgenes_SARS-CoV-2_Full_FC1.csv

