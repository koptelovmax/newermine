#!/bin/bash
# Mine frequent itemsets:
infile1='data//acids.txt'
infile2='data//acids_empty.txt' # empty text file
outfile='data//acids_freq.txt'
./testFrequentLinear 1 $infile1 $infile2 200 1 > $outfile
