#%%
import numpy as np
from lxml import html
import requests
#%%
def get_protein_id(protein_name):
    
    request_str = 'https://www.ncbi.nlm.nih.gov/protein/'+protein_name+'?report=fasta&format=text'
    page = requests.get(request_str)
    tree = html.fromstring(page.content)
    
    return tree.xpath('//div[@class="seq gbff"]')[0].get('val')
#%%
f = open("data//proteins.txt", 'r')
f_out = open("data//proteins_id.txt", "w")

for line in f:
    try:
        data = line.split('|')
        buff_str = get_protein_id(data[0].rstrip('\n'))
        print data
        if len(data) > 1:
            for i in range(1,len(data)):
                buff_str += '|'
                buff_str += get_protein_id(data[i].rstrip('\n'))
        f_out.write(buff_str+'\n')
    except ValueError:
        print "Invalid input:", line
    
f.close()
f_out.close()
#%%    
def get_acid_sequence(protein_id):
    
    request_str = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id='+protein_id+'&db=protein&report=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000'
    page = requests.get(request_str)
    seq = page.text.splitlines()
    
    return  ''.join(seq[1:])
#%%
f = open("data//proteins_ids.txt", 'r')
f_out = open("data//acids.txt", "w")

for line in f:
    try:
        data = line.split('|')
        print data
        buff_str = ''
        for i in range(len(data)):
            buff_str += get_acid_sequence(data[i].rstrip('\n'))
        f_out.write(buff_str+'\n')
    except ValueError:
        print "Invalid input:", line
    
f.close()
f_out.close()
#%%
ligands_dict = {}
ligands = []

f1 = open("data//ligands.txt", 'r')
f2 = open("data//smiles.txt", 'r')

for line in f1:
    try:
        ligands.append(line.rstrip('\n'))
    except ValueError:
        print "Invalid input:", line
        
i = 0
for line in f2:
    try:
        ligands_dict[ligands[i]]=line.rstrip('\n')
        i+=1
    except ValueError:
        print "Invalid input:", line
    
f1.close()
f2.close()
#%%
targets_dict = {}
targets = []

f1 = open("data//targets.txt", 'r')
f2 = open("data//acids.txt", 'r')

for line in f1:
    try:
        targets.append(line.rstrip('\n'))
    except ValueError:
        print "Invalid input:", line
        
i = 0
for line in f2:
    try:
        targets_dict[targets[i]]=line.rstrip('\n')
        i+=1
    except ValueError:
        print "Invalid input:", line
    
f1.close()
f2.close()
#%%
f = open("data//int_edited.csv", 'r')
f_out = open("data//int_selected.csv", "w")

f.next()
f_out.write('ligand id,smiles code,target id,acids seq,affinity,units\n')
for line in f:
    try:
        data = line.split(',')
        if data[1] in ligands_dict.keys():
            smile_code = ligands_dict[data[1]]
        else:
            smile_code = ''
        if data[0] in targets_dict.keys():
            acids_seq = targets_dict[data[0]]
        else:
            acids_seq = ''
        f_out.write(data[1]+','+smile_code+','+data[0]+','+acids_seq+','+data[4]+','+data[2]+'\n')
    except ValueError:
        print "Invalid input:", line
    
f.close()
f_out.close()
#%%
f = open("data//int_selected_edited.csv", "r")

pairs = []

f.next()
for line in f:
    try:
        data = line.split(',')
        pairs.append(data[0]+'-'+data[2])
    except ValueError:
        print "Invalid input:", line
    
f.close()

#from collections import Counter
#print Counter(pairs).most_common(100)
#%%
f = open("data//int_edited_new.csv", 'r')
f_out = open("data//new_dataset.csv", "w")

pairs = []

f.next()
f_out.write('ligand id,smiles code,target id,acids seq,affinity,units\n')
for line in f:
    try:
        data = line.split(',')
        if data[1] in ligands_dict.keys():
            smile_code = ligands_dict[data[1]]
        else:
            smile_code = ''
        if data[0] in targets_dict.keys():
            acids_seq = targets_dict[data[0]]
        else:
            acids_seq = ''
        
        pair = data[1]+'-'+data[0]
        
        if pair not in pairs:
            pairs.append(pair)
        
            f_out.write(data[1]+','+smile_code+','+data[0]+','+acids_seq+','+data[4]+','+data[2]+'\n')
    except ValueError:
        print "Invalid input:", line
    
f.close()
f_out.close()
#%%
