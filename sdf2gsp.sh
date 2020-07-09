#!/bin/bash
# Convert .sdf to .gsp
# test data:
infile='data//drugs_sdf.sdf'
outfile='data//drugs_gsp.gsp'
logfile='data//log_sdf2gsp.txt'
perl sdf2gsp_nolabels_mod_v3.pl $infile $outfile > $logfile
