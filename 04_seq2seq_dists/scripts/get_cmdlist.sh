#! /bin/bash

files='../data/*'

for file in $files
do
  ph_no=$(echo $file | rev | cut -d '/' -f1 | rev | cut -d '_' -f1)
  echo python3 seq2seq_dists.py -f $file -o ../outputs -p ../data/${ph_no}_all_protein_names.txt -c 4
done
