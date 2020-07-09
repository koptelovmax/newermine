import numpy as np
#%%
def vector_sim(vector1,vector2):
    # Jaccard/Tanimoto similarity:
    return float(np.sum(vector1 & vector2))/np.sum(vector1 | vector2)
#%%
# Load ligands ids:
f = open("data//ligands.txt", 'r')
ligands_id = []

for line in f:
    try:
        ligands_id.append(np.int(line))
    except ValueError:
        print "Invalid input:", line
    
f.close()
#%%
# Load targets ids:
f = open("data//targets.txt", 'r')
targets_id = []

for line in f:
    try:
        targets_id.append(np.int(line))
    except ValueError:
        print "Invalid input:", line
    
f.close()
#%%
# Remove not distinct targets:
#name = 'targets_substrings'
name = 'targets_motifs'

X_old = np.load('data//X_'+name+'.npy')
X = []

f = open("data//targets.txt", 'r')
targets_id = []

count = 0
for line in f:
    try:
        if np.int(line) not in targets_id:
            targets_id.append(np.int(line))
            X.append(X_old[count])
        count+=1
    except ValueError:
        print "Invalid input:", line
    
f.close()
#%%
#name = 'ligands_subgraphs'
#name = 'targets_substrings'
name = 'targets_motifs'

#X = np.load('data//X_'+name+'.npy')
f_out = open("data//"+name+"2.txt", "w")

num_vectors = np.size(X,0)
similarity = 0 #np.zeros((num_vectors,num_vectors),float)
#matrix = np.zeros((num_vectors,num_vectors),float)
sims = []

for i in range(num_vectors):
    print 'vector',i
    for j in range(i):
        if (i != j):
            similarity = np.round(vector_sim(X[i],X[j]),4)
            #matrix[i,j] = similarity
            #f_out.write('l'+str(ligands_id[i])+' l'+str(ligands_id[j])+' '+str(similarity)+'\n')
            f_out.write('t'+str(targets_id[i])+' t'+str(targets_id[j])+' '+str(similarity)+'\n')
            sims.append(similarity)
#np.save('data//'+name+'.npy',matrix)
f_out.close()
#%%
# some statistics:
np.mean(sims)
np.median(sims)
np.min(sims)
np.max(sims)
#%%
