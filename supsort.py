#%%
import numpy as np
import sys

infile1 = sys.argv[1] # "dataset//beta//acids_train_0.txt"
infile2= sys.argv[2] # "dataset//beta//acids_freq_0.txt"
outfile = sys.argv[3] # "dataset//beta//acids_freq_0.csv"

# Sort frequent subsequences w.r.t. support:
f = open(infile1, 'r')
seqs = ()
for line in f:
    try:
        seqs = seqs + (line.rstrip(),)
    except ValueError:
        print "Invalid input:", line
f.close()

f = open(infile2, 'r')
f_out = open(outfile, 'w')

items_list = []

for line in f:
    item = line.rstrip()
    num_item = np.sum([(item in seqs[i]) for i in range(len(seqs))])
    items_list.append((num_item,item))

for pair in sorted(items_list, reverse=True):
    f_out.write(pair[1]+','+str(pair[0])+'\n')

f.close()
f_out.close()
#%%