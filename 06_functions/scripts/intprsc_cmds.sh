#!/bin/bash

files="../data/splitted/*"

mkdir ../outputs

for file in $files
do
  filename=$(echo $file | rev | cut -d '/' -f1 | cut -d '.' -f2 | rev)
  echo interproscan.sh -i $file -f tsv -dp -goterms -o ../outputs/${filename}.tsv
done
