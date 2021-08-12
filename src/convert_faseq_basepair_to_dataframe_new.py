import glob
import sys,os
import numpy as np
import pandas as pd


os.chdir(sys.argv[1])
fa_file=sys.argv[2]


RBP=fa_file.split('/')[-1].split('.')[0]
edges_files=glob.glob("*_edges")
    
LIST_NAME=[]
LIST_SEQ=[]
LIST_EDGE=[]
with open(fa_file,'r') as fa:
    lines=fa.readlines() 
    list_=[]
    for i in range(0,len(lines),2):
        if '(' in lines[i]:
            name=(lines[i].strip().split('(')[0].replace(':','_'))[1:]
        else:
            name=(lines[i].strip())[1:]
        #name=(lines[i].strip().replace("(","_").replace(")","").replace(":","_"))[1:]
        #name=(lines[i].split(':')[0]+"_"+lines[i].split(':')[1][:-4])[1:]
        LIST_NAME.append(name)
        LIST_SEQ.append(lines[i+1])
        ind=edges_files.index(name+"_edges")
        LIST_edge_temp=[]
        with open(edges_files[ind]) as ed:
            for line in ed:
                values = line.split()[0:2]  
                new_values = ['A' + str(x) for x in values]
                LIST_edge_temp.append(new_values)
        LIST_EDGE.append(LIST_edge_temp)
df = pd.DataFrame(list(zip(LIST_NAME, LIST_SEQ,LIST_EDGE)), columns =['chr_start_end', 'sequence','base_pair_0.5cutoff']) 
df.to_pickle(RBP+'_seq_basepair.pkl')           

