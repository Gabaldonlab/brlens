#!/bin/bash
for i in 1 2 3 4 5 6 7
do
	rm ../outputs/*
	start_time=$(date +%s.%3N)
	time python3 clade_sp_dist_v2.py -f ../data/0076_108.txt -g ../data/0076_norm_groups.csv -o ../outputs/ -c $i
	end_time=$(date +%s.%3N)
	elapsed=$(echo "scale=3; $end_time - $start_time" | bc)
	echo -e "$i\t$elapsed"
done
