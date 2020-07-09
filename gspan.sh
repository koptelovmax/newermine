#!/bin/bash
# Mine frequent subgraphs from .gsp
infile='data//drugs_gsp.gsp'
./gSpan-64 -f $infile -s 0.3 -o -i
