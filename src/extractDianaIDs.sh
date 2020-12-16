#!/bin/bash

cd targets/Diana

mkdir -p 70
for file in SC2V*
do
  grep ">" $file | awk 'BEGIN{FS="|"}{if($3 >= 0.70)}{print $2, $3}' > 70/$file.tab
  
done
