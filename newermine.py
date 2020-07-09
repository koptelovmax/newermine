#%%
import numpy as np
import networkx as nx
import timeit
import sys

exp_id = sys.argv[1] # instance id
num_srv = 16 # number of instances to run in parallel
#%%
# matrix normalization by columns:
def norm_columns(M):
    M_res = np.zeros((np.size(M,0),np.size(M,1)),float)
    for i in range(np.size(M,1)):
        if np.sum(M[:,i]) != 0:
            M_res[:,i]=M[:,i]/float(np.sum(M[:,i]))
            
    return np.matrix(M_res)
#%%    
# column normalization:
def norm_column(M,i):
    if np.sum(M[:,i]) != 0:
        M[:,i]=M[:,i]/float(np.sum(M[:,i]))
#%%
def load_networks():

    G = nx.MultiGraph()
    
    lig = []
    tar = []
    networks = {}
    
    net = []
    count = 0
    
    # Load data from Drugbank (drug-drug interaction network):
    f = open("data//ligands_drugbank.txt", 'r')
    
    for line in f:
        try:
            data = line.split(' ')
            if data[0] != data[1].rstrip(): # check for loops
                G.add_edge(data[0],data[1].rstrip(),weight=1,key='drugbank',label='drugbank') # ligand-ligand interaction
                lig.append(data[0])
                lig.append(data[1].rstrip())
                net.append(data[0])
                net.append(data[1].rstrip())
            count+=1
        except ValueError:
            print "Invalid input:", line
            
    networks['drugbank'] = np.unique(net)
       
    f.close()
    
    print '1 out of 6 networks loaded'
    
    # Load data from BioGrid (protein-protein interaction network):
    f = open("data//targets_biogrid.txt", 'r')
    
    net = []
    count = 0
    
    for line in f:
        try:
            data = line.split(' ')
            if data[0] != data[1].rstrip(): # check for loops
                G.add_edge(data[0],data[1].rstrip(),weight=1,key='biogrid',label='biogrid') # target-target interaction        
                tar.append(data[0])
                tar.append(data[1].rstrip())
                net.append(data[0])
                net.append(data[1].rstrip())
            count+=1
        except ValueError:
            print "Invalid input:", line
            
    networks['biogrid'] = np.unique(net)
        
    f.close()
    
    print '2 out of 6 networks loaded'
    
    # Load data from IUPHAR (drug-protein interaction network):
    f = open("data//interactions_iuphar+.csv", 'r')
    
    line = f.next()
    net = []
    count = 0
    for line in f:
        try:
            data = line.split(',')
            if data[14].rstrip() == 'Positive':
                G.add_edge('l'+data[0],'t'+data[4],weight=1,key='iuphar',label='iuphar') # ligand-target positive interaction (interaction exists)
                lig.append('l'+data[0])
                tar.append('t'+data[4])
                net.append('l'+data[0])
                net.append('t'+data[4])
            elif data[14].rstrip() == 'Negative':
                G.add_edge('l'+data[0],'t'+data[4],weight=0,key='iuphar',label='iuphar') # ligand-target negative interaction (no interaction)
                lig.append('l'+data[0])
                tar.append('t'+data[4])
                net.append('l'+data[0])
                net.append('t'+data[4])
            count+=1
        except ValueError:
            print "Invalid input:", line
            
    networks['iuphar'] = np.unique(net)
        
    f.close()
    
    print '3 out of 6 networks loaded'
    
    # Load ligands similarities (based on subgraphs):
    f = open("data//ligands_subgraphs.txt", 'r')
    
    net = []
    count = 0
    for line in f:
        try:
            data = line.split(' ')
            G.add_edge(data[0],data[1],weight=float(data[2]),key='lig_sg',label='lig_sg') # ligand-ligand similarity
            lig.append(data[0])
            lig.append(data[1])
            net.append(data[0])
            net.append(data[1])
            #if count % 10000 == 0:
            #    print count
            count+=1
        except ValueError:
            print "Invalid input:", line
            
    networks['lig_sg'] = np.unique(net)
        
    f.close()
    
    print '4 out of 6 networks loaded'
    
    # Load targets similarities (based on substrings):
    f = open("data//targets_substrings.txt", 'r')
    
    net = []
    count = 0
    for line in f:
        try:
            data = line.split(' ')
            G.add_edge(data[0],data[1],weight=float(data[2]),key='tar_ss',label='tar_ss') # target-target similarity
            tar.append(data[0])
            tar.append(data[1])
            net.append(data[0])
            net.append(data[1])
            count+=1
        except ValueError:
            print "Invalid input:", line
            
    networks['tar_ss'] = np.unique(net)
        
    f.close()
    
    print '5 out of 6 networks loaded'
    
    # Load targets similarities (based on motifs):
    f = open("data//targets_motifs.txt", 'r')
    
    net = []
    count = 0
    for line in f:
        try:
            data = line.split(' ')
            G.add_edge(data[0],data[1],weight=float(data[2]),key='tar_mot',label='tar_mot') # target-target similarity
            tar.append(data[0])
            tar.append(data[1])
            net.append(data[0])
            net.append(data[1])
            count+=1
        except ValueError:
            print "Invalid input:", line
            
    networks['tar_mot'] = np.unique(net)
        
    f.close()
    
    print '6 out of 6 networks loaded\n'
    
    print 'Graph has been loaded!'
    print 'Nodes: ',G.number_of_nodes()
    print 'Edges: ',G.number_of_edges()
    
    return G,list(np.unique(lig)),list(np.unique(tar)),networks
#%%
def construct_M_N(G,ligands,targets,layers,networks,num_tt,num_dd):
    
    n = len(targets) # number of targets
    m = len(ligands) # number of ligands
    num_dt = len(networks)-(num_tt+num_dd) # number of dt networks
    
    # adjecency matrix initialization:
    M = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),num_tt*n+num_dd*m+num_dt*(n+m)),float)
    
    # network transitions matrix initialization:
    N = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),num_tt*n+num_dd*m+num_dt*(n+m)),bool)
    
    # parse graph G:
    for e in G.edges_iter(data=True):

        if e[2]['label'] in networks:
            if networks.index(e[2]['label']) < num_tt:
                # parse tt network:
                node1_index = networks.index(e[2]['label'])*n+targets.index(e[0])
                node2_index = networks.index(e[2]['label'])*n+targets.index(e[1])
            elif networks.index(e[2]['label']) < (num_tt+num_dd):
                # parse dd network:
                node1_index = num_tt*n+(networks.index(e[2]['label'])-num_tt)*m+ligands.index(e[0])
                node2_index = num_tt*n+(networks.index(e[2]['label'])-num_tt)*m+ligands.index(e[1])
            else:
                # parse dt network:
                if ('t' in e[0]) and ('l' in e[1]):
                    node1_index = num_tt*n+num_dd*m+(networks.index(e[2]['label'])-num_tt-num_dd)*(n+m)+targets.index(e[0])
                    node2_index = num_tt*n+num_dd*m+(networks.index(e[2]['label'])-num_tt-num_dd)*(n+m)+n+ligands.index(e[1])
                else:
                    node1_index = num_tt*n+num_dd*m+(networks.index(e[2]['label'])-num_tt-num_dd)*(n+m)+n+ligands.index(e[0])                    
                    node2_index = num_tt*n+num_dd*m+(networks.index(e[2]['label'])-num_tt-num_dd)*(n+m)+targets.index(e[1])

            # fill in matrix M:
            M[node1_index,node2_index] = e[2]['weight']
            M[node2_index,node1_index] = e[2]['weight']
            
            # fill in matrix N:           
            N[node1_index,node1_index] = True
            N[node2_index,node2_index] = True       
            
            # parse all networks:
            for i in range(len(networks)):

                # fill in node1: 
                if (e[0] in layers[networks[i]]) and networks[i] != e[2]['label']:
                    if i < num_tt:
                        # parse tt network:
                        node_index = i*n+targets.index(e[0])
                    elif i < (num_tt+num_dd):
                        # parse dd network:
                        node_index = num_tt*n+(i-num_tt)*m+ligands.index(e[0])
                    else:
                        # parse dt network:
                        if 't' in e[0]:
                            node_index = num_tt*n+num_dd*m+(i-num_tt-num_dd)*(n+m)+targets.index(e[0])
                        else:
                            node_index = num_tt*n+num_dd*m+(i-num_tt-num_dd)*(n+m)+n+ligands.index(e[0])                            
                    N[node_index,node1_index] = True
                
                # fill in node2:
                if (e[1] in layers[networks[i]]) and networks[i] != e[2]['label']:
                    if i < num_tt:
                        # parse tt network:
                        node_index = i*n+targets.index(e[1])
                    elif i < (num_tt+num_dd):
                        # parse dd network:
                        node_index = num_tt*n+(i-num_tt)*m+ligands.index(e[1])
                    else:
                        # parse dt network:
                        if 't' in e[1]:
                            node_index = num_tt*n+num_dd*m+(i-num_tt-num_dd)*(n+m)+targets.index(e[1])
                        else:
                            node_index = num_tt*n+num_dd*m+(i-num_tt-num_dd)*(n+m)+n+ligands.index(e[1])                            
                    N[node_index,node2_index] = True
    
    return norm_columns(M),N
#%%
def init_node(V,start_node,G,ligands,targets,networks,num_tt,num_dd):
    
    n = len(targets) # number of targets
    m = len(ligands) # number of ligands
    num_dt = len(networks)-(num_tt+num_dd) # number of dt networks
    
    if start_node in targets:
        index_start_node = targets.index(start_node)
        for i in range(num_tt):
            if filter(lambda x:networks[i] in x, G.edges(start_node,data='label')):
                V[i*n+index_start_node] = 1
        for i in range(num_dt):
            if filter(lambda x:networks[num_tt+num_dd+i] in x, G.edges(start_node,data='label')):
                V[num_tt*n+num_dd*m+i*(n+m)+index_start_node] = 1
    elif start_node in ligands:
        index_start_node = ligands.index(start_node)
        for i in range(num_dd):
            if filter(lambda x:networks[num_tt+i] in x, G.edges(start_node,data='label')):
                V[num_tt*n+i*m+index_start_node] = 1
        #for i in range(num_dt):
        #    if filter(lambda x:networks[num_tt+num_dd+i] in x, G.edges(start_node,data='label')):
        #        V[num_tt*n+num_dd*m+i*(n+m)+n+index_start_node] = 1        
#%%
def init_vector(start_ligand,G,V_tar,ligands,targets,networks,num_tt,num_dd):
    
    eta = 0 # importance of targets networks
    n = len(targets) # number of targets
    m = len(ligands) # number of ligands
    num_dt = len(networks)-(num_tt+num_dd) # number of dt networks
    
    # temporary vectors initialization:
    V1 = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),1),float)
    #V2 = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),1),float)
    
    # initialize given ligand:    
    init_node(V1,start_ligand,G,ligands,targets,networks,num_tt,num_dd)
    #V1 = norm_columns(V1)
    #V1 = V1*mju
        
    # initialize targets:
    #for t in targets:
    #    init_node(V2,t,G,ligands,targets,networks,num_tt,num_dd)
    #V2 = norm_columns(V2)
    #V2 = V2*(1-mju)
    
    return np.dot(norm_columns(V1),(1-eta))+np.dot(V_tar,eta)
#%%
def newermine(M,V):

    V_prev = np.matrix(np.zeros((np.size(V,0),np.size(V,1)),float))
    V_cur = np.matrix(np.zeros((np.size(V,0),np.size(V,1)),float))
    
    V_cur[:] = V[:]
    
    betta = 0.7 # taxation  
    count = 0
    count_stop = 40
    stop = False
    
    # convergence:
    while not stop:
        V_prev[:] = V_cur[:]
        V_cur = np.dot(betta,np.dot(M,V_prev))+np.dot((1-betta),V)
        count += 1
        if (np.sum(np.abs(V_cur - V_prev)) < 1.0e-03) or (count >= count_stop):
            stop = True
            
    return V_cur,count
#%%
def merge_V(V_avg,ligands,targets,num_tt,num_dd):
    
    n = len(targets) # number of targets
    m = len(ligands) # number of ligands
    V_res = np.zeros((n+m,1),float) # resulting vector
      
    # merge everything back:
    for i in range(len(V_avg)):
        if i < (num_tt*n):
            V_res[i % n,0] += V_avg[i,0]
        elif i < (num_tt*n+num_dd*m):
            V_res[n + (i-num_tt*n) % m,0] += V_avg[i,0]
        else:
            V_res[(i-num_tt*n-num_dd*m) % (n+m),0] += V_avg[i,0]
    
    return V_res[0:n]
#%%
# select networks:
networks = ['tar_ss','tar_mot','biogrid','lig_sg','drugbank','iuphar']

num_tt = 3 # number of target-target networks
num_dd = 2 # number of drug-drug networks
num_dt = len(networks)-(num_tt+num_dd) # number of dt networks
#%%
# load data:
G,ligands,targets,layers = load_networks()

n = len(targets) # number of targets
m = len(ligands) # number of ligands

# build matrices M_norm, N:
M_norm,N = construct_M_N(G,ligands,targets,layers,networks,num_tt,num_dd)

# save for the experiments:
#np.save('matrix_M_norm.npy',M_norm)
#np.save('matrix_N.npy',N)

# load matrices:
#M_norm = np.load('matrix_M_norm.npy')
#N = np.load('matrix_N.npy')

# compute final matrix M:
M = norm_columns(np.array(np.dot(M_norm,N)))

## save final M:
#M_prev = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),num_tt*n+num_dd*m+num_dt*(n+m)),float)
#M_prev[:,:] = M[:,:]
#%%     
# leave-one-out cross-validation:
dt_network = 'iuphar' # drug-target network name for experiments
dt_nodes_index = num_tt*n+num_dd*m+(networks.index(dt_network)-num_tt-num_dd)*(n+m)

top_k = 20 # precision at top 20

# list of existing links (from dt_network):
existing_links = filter(lambda x:dt_network in x, G.edges(data='label'))

# data classes initialization:
Y_test = np.zeros(len(existing_links),bool)
Y_pred = np.zeros(len(existing_links),bool)

# initialize targets:
V_tar = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),1),float)
for t in targets:
    init_node(V_tar,t,G,ligands,targets,networks,num_tt,num_dd)
V_tar = norm_columns(V_tar)
#%%
# define job list for the given machine:
num_experiments = (len(existing_links) // num_srv + 1)
exp_start = num_experiments*(int(exp_id)-1)
if int(exp_id) != num_srv:
    exp_end = num_experiments*int(exp_id)
else:
    exp_end = len(existing_links)
#%%
f_out = open('experiments//y_pred_'+str(exp_id)+'.txt','w')
# iterate over the links:
for i in range(exp_start,exp_end):
    print 'Prediction of a link ',i-exp_start,' out of ',exp_end-exp_start,' exp #',exp_id
    f_out.write('Prediction of a link '+str(i-exp_start)+' out of '+str(exp_end-exp_start)+' exp #'+str(exp_id)+'\n')
    
    start = timeit.default_timer() # Initialize the timer to compute the time
    
    # memorize an edge:
    if ('l' in existing_links[i][0]) and ('t' in existing_links[i][1]):
        node_start = existing_links[i][0]
        node_end = existing_links[i][1]
    elif ('l' in existing_links[i][1]) and ('t' in existing_links[i][0]):
        node_start = existing_links[i][1]
        node_end = existing_links[i][0]
    
    # load link class:    
    Y_test[i] = G.get_edge_data(node_start,node_end)[dt_network]['weight']
        
    # save edge data before removal:
    edge_rem_weight = G.get_edge_data(node_start,node_end)[dt_network]['weight']

    # remove edge from G:
    G.remove_edge(node_start,node_end)

    # remove edge from M_norm:
    node_start_index = dt_nodes_index+n+ligands.index(node_start)
    node_end_index = dt_nodes_index+targets.index(node_end)
    M_norm[node_start_index,node_end_index] = 0
    M_norm[node_end_index,node_start_index] = 0
    norm_column(M_norm,node_start_index)
    norm_column(M_norm,node_end_index)
   
    #start = timeit.default_timer() # Initialize the timer to compute the time
   
    # list of ligand-target indices of the removed edge:
    nodes_indices = []
    
    # add ligand indices to the list:
    for l in range(len(networks)):
        if node_start in layers[networks[l]]:
            if l < (num_tt+num_dd):
                # parse dd network:
                nodes_indices.append(num_tt*n+(l-num_tt)*m+ligands.index(node_start))
            else:
                # parse dt network:
                nodes_indices.append(num_tt*n+num_dd*m+(l-num_tt-num_dd)*(n+m)+n+ligands.index(node_start))
    
    # add target indices to the list:
    for l in range(len(networks)):
        if node_end in layers[networks[l]]:
            if l < num_tt:
                # parse tt network:
                nodes_indices.append(l*n+targets.index(node_end))
            else:
                # parse dt network:
                nodes_indices.append(num_tt*n+num_dd*m+(l-num_tt-num_dd)*(n+m)+targets.index(node_end))
    
    # memorize columns of M and update it by recomputing the columns:
    M_columns = []
    for j in range(len(nodes_indices)):
        temp_column = np.zeros((num_tt*n+num_dd*m+num_dt*(n+m),1),float)
        temp_column[:] = M[:,nodes_indices[j]]
        M_columns.append(temp_column)
        M[:,nodes_indices[j]] = np.dot(M_norm,N[:,nodes_indices[j]]).reshape(-1,1)
        norm_column(M,nodes_indices[j])
        
    # try to predict removed edge by NEWERMINE:
    V = init_vector(node_start,G,V_tar,ligands,targets,networks,num_tt,num_dd)
    V_cur,count = newermine(M,V)
    V_res = merge_V(V_cur,ligands,targets,num_tt,num_dd)
    
    print 'Steps:',count
    f_out.write('Steps:'+str(count)+'\n')
   
    # list of edges starting from node_start:
    verified_links = filter(lambda x:dt_network in x, G.edges(node_start,data='label'))
    ver_targs_indexes = []
    for j in range(len(verified_links)):
        if 't' in verified_links[j][0]:
            ver_targs_indexes.append(targets.index(verified_links[j][0]))
        elif 't' in verified_links[j][1]:
            ver_targs_indexes.append(targets.index(verified_links[j][1]))

    # dictionary of targets without targets connected to node_start:
    targets_probs = {}
    for j in range(len(targets)):
        if j not in ver_targs_indexes:
            targets_probs[targets[j]]=float(V_res[j]) 

    # take top "top_k" (top 20) targets:
    top_targets_probs = sorted(targets_probs.items(),key=lambda x: x[1],reverse=True)[0:top_k]
    top_targets = [t[0] for t in top_targets_probs]

    # perform prediction:
    if node_end in top_targets:
        Y_pred[i] = True
        print 'Prediction: True'
    else:
        print 'Prediction: False'
    
    # put removed link back to graph G:
    G.add_edge(node_start,node_end,weight=edge_rem_weight,key=dt_network,label=dt_network)

    # put removed link back to matrix M_norm:
    if np.max(M_norm[:,node_end_index]) != 0:
        M_norm[node_start_index,node_end_index] = np.max(M_norm[:,node_end_index])
        norm_column(M_norm,node_end_index)
    else:
        M_norm[node_start_index,node_end_index] = 1
        
    if np.max(M_norm[:,node_start_index]) != 0:
        M_norm[node_end_index,node_start_index] = np.max(M_norm[:,node_start_index])
        norm_column(M_norm,node_start_index)
    else:
        M_norm[node_end_index,node_start_index] = 1
        
    # put memorised columns back to M:
    for j in range(len(nodes_indices)):
        M[:,nodes_indices[j]] = M_columns[j]
        
    stop = timeit.default_timer() # stop the time counting
    print 'Time prediction, sec: '+str(stop-start)+'\n'
    f_out.write('Time spent, sec: '+str(stop-start)+'\n\n')
#%%            
# save results:            
np.save('experiments//y_pred_'+str(exp_id)+'.npy',Y_pred)

f_out.close()

print 'Experiment is done!'
#%%
