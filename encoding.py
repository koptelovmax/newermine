import numpy as np

num_ligands = 6821
num_targets = 1882
d_lig = 200
d_tar = 200
#%%
X_lig = np.zeros((num_ligands,d_lig),bool)

f = open("data//ligands.gsp.fp.sorted", 'r')

i = 0
for line in f:
    try:
        data = line.split(' ')
        if (data[0] == 'x') and (i < d_lig):
            for item in data[1:-1]:
                X_lig[np.int(item),i] = True
            i+=1                        
    except ValueError:
        print "Invalid input:", line

f.close()
#%%
np.save("data//X_ligands_subgraphs.npy",X_lig)
#%%
X_tar = np.zeros((num_targets,d_tar),bool)

f = open("data//acids.txt", 'r')
acids_train = []

for line in f:
    try:
        acids_train.append(line.rstrip())
    except ValueError:
        print "Invalid input:", line

f.close()
#%%
f = open("data//acids_freq.csv", 'r')

i = 0
for line in f:
    try:
        data = line.split(',')
        if i < d_tar:
            for j in range(num_targets):
                if data[0] in acids_train[j]:
                    X_tar[j,i] = True
            i+=1
    except ValueError:
        print "Invalid input:", line

f.close()
#%%
np.save("data//X_targets_substrings.npy",X_tar)
#%%
X_tar2 = np.zeros((num_targets,d_tar),bool)

# Load biological descriptors:
patterns = []

f = open("data//patterns.csv", 'r')

line = f.next()
for line in f:
    try:
        data = line.split(',')
        patterns.append(data[0])
    except ValueError:
        print "Invalid input:", line
f.close()
#%%
# Encode targets with motifs:
f = open("data//acids.out", 'r')

i = 0
for line in f:
    try:
        if ':' in line:
            data = line.split(' ')
            if data[2] in patterns:
                X_tar2[np.int(data[0][3:])-1,patterns.index(data[2])] = True
            i+=1
    except ValueError:
        print "Invalid input:", line
f.close()
#%%
np.save("data//X_targets_motifs.npy",X_tar2)
#%%
