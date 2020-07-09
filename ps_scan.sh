#!/bin/bash
# Mine Prosite motifs from amino-acid sequences
infile='data//acids.dat' # amino-acid sequences in .dat format (see README of Prosite)
outfile='data//acids.out'
perl ps_scan.pl -r $infile > $outfile

