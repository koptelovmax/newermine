#!/bin/bash
# Convert .smi to .sdf
infile='data//drugs_smiles.smi'
outfile='data//drugs_sdf.sdf'
obabel -ismi $infile -O $outfile

