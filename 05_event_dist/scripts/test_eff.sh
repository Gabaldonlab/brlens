#!/bin/bash
for i in 1 2 3 4 5 6 7
do
	rm ../outputs/*
	SECONDS=0
	time python3 clade_sp_dist_v2.py -f ../data/0076_108.txt -g ../data/0076_norm_groups.csv -o ../outputs/ -c $i
	duration=$SECONDS
	echo -e "$i\t$duration"
done
