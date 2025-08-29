#!/bin/bash
for i in $(seq 1 10)
do echo $i
cp ace_continuous_50t_BMOU_5fossil_2ace_ord_ID.R tmp

sed 's/<ID>/00'"$i"'/g' tmp > ace_continuous_50t_BMOU_5fossil_2ace_ord_${i}.R

done

