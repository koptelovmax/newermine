# Link Prediction in Multi-layer Networks and Its Application to Drug Design
The code for the paper presented at IDA'18:
https://doi.org/10.1007/978-3-030-01768-2_15

## To reconstruct the IUPHAR network (the multi-layer graph with 6 layers):
1) The IUPHAR database needs to be downloaded from: https://www.guidetopharmacology.org/download.jsp and processed manualy with the use of tables editor such as MS Excel or OpenOffice Calc (see code for examples)
2) The DrugBank database needs to be downloaded from: https://www.drugbank.ca/releases, unzipped and placed into 'data' folder
3) The BioGrid network needs to be downloaded from: https://downloads.thebiogrid.org/BioGRID, unzipped and placed into 'data' folder
4) **datasets.py** -- parses DrugBank, BioGrid and construction of corresponding networks
5) **smi2sdf.sh** -- convert smiles codes representing drugs (taken from *ligands.csv* of IUPHAR) into .sdf format
6) **sdf2gsp.sh** -- convert .sdf file into .gsp format (with the use of *sdf2gsp_nolabels_mod_v3.pl*)
7) **gspan.sh** -- mine frequent subgraphs from .gsp file (requires *gSpan-64* to be downloaded from: https://www.cs.ucsb.edu/~xyan/software/gSpan.htm)
8) **fpsort.py** -- sort frequent subgraphs w.r.t. support
9) **parser.py** -- parse the NCBI database for amino-acids sequences
10) **freqlinear.sh** -- mine frequent subsequences from amino-acid sequences (requires *testFrequentLinear* to be downloaded from: https://www.bio.ifi.lmu.de/forschung/algorithmic_bioinformatics/succinct/)
11) **supsort.py** -- sort frequent subsequences w.r.t. support
12) **ps_scan.py** -- mine Prosite motifs from amino-acid sequences (requires *ps_scan_linux_x86* to be downloaded from: ftp://ftp.expasy.org/databases/prosite/ps_scan/)
13) **encoding.py** -- construct feature vectors for ligands and targets 
14) **similarities.py** -- construction of similarities networks

Already constructed network can be found here: https://zimmermanna.users.greyc.fr/supplementary/networks.zip

## To reconstruct the experiments:
1) Download, unzip and put into 'data' folder already constructed network or reconstruct the network yourself
2) **newermine.py** -- run the implementation of the NEVERMINE algorithm and leave-one-out cross-validation
3) **evaluation.py** -- run the framework for computing the quality measures

Environment requirements: Python 2.7, networkx 1.11, openbabel, scikit-learn, numpy
