#%%
import numpy as np
#%%
# Load ligands INN:
ligands_names = {}

f = open("data//ligands_INN+.csv", 'r')

line = f.next()
for line in f:
    try:
        data = line.split(',')
        ligands_names[np.int(data[0])] = data[1].rstrip()
    except ValueError:
        print "Invalid input:", line
    
f.close()
#%%
# Load ligands groups and some targets properties:
ligands_groups = {}
targets_types = {}
targets_families = {}
targets_swissprot = {}
targets_entrez = {}

f1 = open("data//lig_types.csv", 'r')
f2 = open("data//tar_types2.csv", 'r')

line=f1.next()
for line in f1:
    try:
        data = line.split(',')
        ligands_groups[np.int(data[0])] = data[1].rstrip()
    except ValueError:
        print "Invalid input:", line

line=f2.next()        
for line in f2:
    try:
        data = line.split(',')
        targets_types[np.int(data[2])]=data[0].replace('"', '')
        targets_families[np.int(data[2])]=data[1].replace('"', '')
        if data[3] != '':
            targets_swissprot[np.int(data[2])]=data[3].replace('"', '')
        if data[4] != '':
            targets_entrez[np.int(data[2])]=data[4].replace('"', '').rstrip()
    except ValueError:
        print "Invalid input:", line
    
f1.close()
f2.close()
#%%
# Parse DrubBank database:
f = open("data//full_database.xml", 'r')
f_out = open("data//drugbank.txt", "w")

count = 0
for line in f:
    try:
        if 'drug type' in line:
            while '<name>' not in line:
                line = f.next()
                count+=1
            drug_name = line.strip()[6:-7]
            while ('<drug-interactions>' not in line) and ('</drug>' not in line):
                line = f.next()
                count+=1
            while ('</drug-interactions>' not in line) and ('</drug>' not in line):
                if '<name>' in line:
                    f_out.write(drug_name+'\t'+line.strip()[6:-7]+'\n')
                line = f.next()
                count+=1
            #f_out.write(drug_name+'\n')
        else:
            count+=1                
    except ValueError:
        print "Invalid input:", line

f.close()
f_out.close()
#%%
# Load DrugBank drugs and create list of DrugBank drugs:
f = open("data//drugbank.txt", 'r')

drugbank_list = []
count = 0
for line in f:
    try:
        data = line.split('\t')
        #drugbank_list.append(data[0].strip())
        #drugbank_list.append(data[1].strip())
        drugbank_list.append(data[0].strip().lower())
        drugbank_list.append(data[1].strip().lower())
        count+=1
    except ValueError:
        print "Invalid input:", line

f.close()

drugbank_drugs = np.unique(drugbank_list)
#%%
# Parse BioGrid database and create list of BioGrid proteins:
f = open("data//BIOGRID-ALL-3.4.154.tab2.txt", 'r')

line = f.next()
biogrid_list = []
count = 0
for line in f:
    try:
        data = line.split('\t')
        biogrid_list.append(data[1].strip())
        biogrid_list.append(data[2].strip())
        count+=1
    except ValueError:
        print "Invalid input:", line

f.close()

biogrid_proteins = np.unique(biogrid_list)
#%%
# Construct dataset:
f = open("data//all_data.csv", 'r')
f_out = open("data//iuphar_interactions.csv", "w")
#f = open("dataset//all_int.csv", 'r')
#f_out = open("dataset//networks//all_int_data.csv", "w")

f_out.write('ligand id,ligand group,smiles code,INN,target id,target type,target family,acids seq,swissprot id,entrez id,drugbank,biogrid,affinity,units,class\n')

line = f.next()
for line in f:
    try:
        data = line.split(',')
        
        if np.int(data[0]) in ligands_groups:
            lig_group = ligands_groups[np.int(data[0])]
        else:
            lig_group = ''

        if np.int(data[0]) in ligands_names:
            lig_INN = ligands_names[np.int(data[0])]
            
            #if ligands_names[np.int(data[0])].lower() in drugbank_drugs:
            if any(ligands_names[np.int(data[0])].lower() in item for item in drugbank_drugs):            
                lig_db = '+'
            else:
                lig_db = ''
        else:
            lig_INN = ''
            lig_db = ''
        
        if np.int(data[2]) in targets_types:
            tar_type = targets_types[np.int(data[2])]
        else:
            tar_type = ''
        
        if np.int(data[2]) in targets_families:
            tar_fam = targets_families[np.int(data[2])]
        else:
            tar_fam = ''
            
        if np.int(data[2]) in targets_swissprot:
            tar_prot = targets_swissprot[np.int(data[2])]
        else:
            tar_prot = ''
            
        if np.int(data[2]) in targets_entrez:
            tar_entr = targets_entrez[np.int(data[2])]
            
            if targets_entrez[np.int(data[2])] in biogrid_proteins:
                tar_biogr = '+'
            else:
                tar_biogr = ''
        else:
            tar_entr = ''
            tar_biogr = ''
            
        f_out.write(data[0]+','+lig_group+','+data[1]+','+lig_INN+','+data[2]+','+tar_type+','+tar_fam+','+data[3]+','+tar_prot+','+tar_entr+','+lig_db+','+tar_biogr+','+data[4]+','+data[5]+','+data[6])
    except ValueError:
        print "Invalid input:", line
    
f.close()
f_out.close()
#%%
f = open("data//drugbank.txt", 'r')
f_out = open("data//ligands_drugbank.txt", "w")

for line in f:
    try:
        data = line.split('\t')
        if (data[0].strip().lower() in ligands_names.values()) and (data[1].strip().lower() in ligands_names.values()):
            f_out.write("l"+str(ligands_names.keys()[ligands_names.values().index(data[0].strip().lower())])+" l"+str(ligands_names.keys()[ligands_names.values().index(data[1].strip().lower())])+"\n")
    except ValueError:
        print "Invalid input:", line

f.close()
f_out.close()
#%%
f = open("data//BIOGRID-ALL-3.4.154.tab2.txt", 'r')
f_out = open("data//targets_biogrid.txt", "w")

for line in f:
    try:
        data = line.split('\t')
        if (data[1] in targets_entrez.values()) and (data[2] in targets_entrez.values()):
            f_out.write("t"+str(targets_entrez.keys()[targets_entrez.values().index(data[1])])+",t"+str(targets_entrez.keys()[targets_entrez.values().index(data[2])])+","+data[17]+","+data[18]+"\n")
    except ValueError:
        print "Invalid input:", line

f.close()
f_out.close()
#%%
