#%%
import numpy as np
import sys

infile = sys.argv[1]
# Sort frequent subgraphs w.r.t. support:
subgr_list = []

buff_str = ''
f = open(infile, 'r')
f_out = open(infile+".sorted", 'w')

for line in f:
    try:
        data = line.split(' ')
        if data[0] == 't':
            buff_str += line
            sup = np.int(data[4])
            stop = False
            while not stop:
                line = f.next()
                buff_str += line
                data = line.split(' ')
                if data[0] == 'x':
                    stop = True
            subgr_list.append((sup,buff_str))
            buff_str = ''
    except ValueError:
        print "Invalid input:", line

for pair in sorted(subgr_list, reverse=True):
    f_out.write(pair[1]+'\n')

f.close()
f_out.close()
#%%