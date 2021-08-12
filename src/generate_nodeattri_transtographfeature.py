#2021-01-12 by Hongli Ma
#This python script is crated as a new version of generating the nodeattri by removing the onehot encoding; 3mer,4mer, as well as 5mer 

import numpy as np
import sys,os
import pandas as pd


kmer=int(sys.argv[1])

# change working directory
os.chdir(sys.argv[2])
# input file
fafile=sys.argv[3]

RBP=fafile.split('.')[0]
strufile="unpaired_probability_"+RBP+"_u1.txt"
name_list=[]
SEQ_list=[]
with open(fafile,'r') as fa:
    lines=fa.readlines()
    for line in lines:
        if '>' in line:
            name_list.append(line.strip())
        else:
            node_seq_list=[]
            line=line.replace('a','A').replace('t','U').replace('T','U').replace('c','C').replace('g','G')
            na_list=list(line.strip())
            for i in range(0,100-kmer+1):
                node_seq=[]
                for ii in range(0,kmer):
                    node_seq.append(na_list[i+ii])
                node_seq_list.append(node_seq)
            SEQ_list.append(node_seq_list)

SEC_list=[]
with open(strufile,'r') as stru:
    lines2=stru.readlines()
    for i in range(0,len(lines2),6):
        node_sec_list=[]
        p1_list=lines2[i+1].strip().split('\t')
        p2_list=lines2[i+2].strip().split('\t')
        p3_list=lines2[i+3].strip().split('\t')
        p4_list=lines2[i+4].strip().split('\t')
        p5_list=lines2[i+5].strip().split('\t')

        for j in range(0,100-kmer+1):
            node_sec=[]
            for k in range(0,kmer):
                node_sec.extend((p1_list[j+k],p2_list[j+k],p3_list[j+k],p4_list[j+k],p5_list[j+k]))
            node_sec_list.append(node_sec)
        SEC_list.append(node_sec_list)

feature_list=[]
for i in range(0,len(SEQ_list)):
    A=np.array(SEQ_list[i])
    B=np.array(SEC_list[i])
    C=np.concatenate((A,B),axis=1)
    feature_list.append(C)

df = pd.DataFrame(zip(name_list, feature_list), columns = ['name', 'nodeattri'])


file=pd.read_pickle('graphk_updated_'+RBP+'_'+str(kmer)+'mer.pkl')
name_list=file['name']
edge_list=file['result']

nafile=df

name_list2=nafile['name']
na_list=nafile['nodeattri']

name_list_new=[]
for i in list(name_list2):
    i=i.rsplit("(",1)[0][1:].replace(":","_")
    name_list_new.append(i)
    
Final_list=[]

#print(list(name_list))
#print(name_list_new)

for i in range(0,len(name_list)):
    
    nn=name_list[i]
    

    try:
        index=name_list_new.index(nn)
        na=na_list[index]
    
        
        temp=edge_list[i]      
        n1=int(temp[0])
        n2=int(temp[1])
        temp.extend(na[n1-1])
        temp.extend(na[n2-1])
        #print(temp)

    except:
           pass

    Final_list.append(temp)

df2=pd.DataFrame(zip(name_list,Final_list),columns=['name','feature'])
df2.to_pickle(RBP+'_'+str(kmer)+'mer_graphedgewithfeature.pkl')

