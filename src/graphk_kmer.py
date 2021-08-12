import sys, os
import os.path as osp
import glob
import numpy as np
from diameter import Diameter as D
import pandas as pd
import multiprocessing
from multiprocessing import Pool


# change working directory
os.chdir(sys.argv[1])

RBP=sys.argv[2]
k=int(sys.argv[3])

file_list = glob.glob('*_edges') 

def Graphk(file):
    d = {}
    with open(file,'r') as f:
        for line in f:
            words = str.split(line)
            d[(int(words[0]),int(words[1]))] = float(words[2])
    
    result={}
    
    for key in d.keys(): 
    
            i=key[0]
            j=key[1]    
            temp_list = []
            V=[]
            for addnum in range(0,k):
                V.append(i+addnum)
                V.append(j-addnum)
                if addnum <k-1:
                    temp_list.append([i+addnum,i+addnum+1,1])
                    temp_list.append([j-addnum-1,j-addnum,1])
                      
            for ii in range(i,i+k):
                for jj in range(j,j-k,-1):
                    try:
                        temp_list.append([ii,jj,d[ii,jj]])
                    except:
                        continue
            dia=D(len(temp_list),V,temp_list).floyd()
            weight=1/dia
            result['\t'.join([str(i),str(j-k+1),str(weight)])] = [i,j-k+1,weight,-k,-k] 
            
            
            
            if key[0]>k-1 and key[1]<100-k+2:
                i=key[0]
                j=key[1]    
                temp_list = []
                V=[]
                for addnum in range(0,k):
                    V.append(i-addnum)
                    V.append(j+addnum)
                    if addnum <k-1:
                        temp_list.append([i-addnum-1,i-addnum,1])
                        temp_list.append([j+addnum,j+addnum+1,1])
                      
                for ii in range(i,i-k,-1):
                    for jj in range(j,j+k):
                        try:
                            temp_list.append([ii,jj,d[ii,jj]])
                        except:
                            continue
                dia=D(len(temp_list),V,temp_list).floyd()
                weight=1/dia
                result['\t'.join([str(i-k+1),str(j),str(weight)])] = [i-k+1,j,weight,-1,-1]           
                
                
                
    df = pd.DataFrame(list(result.values()), columns=('node1','node2','weight','typeA','typeB'))
    df_final=df[df.weight < 1/(2*(k-1)+0.5)]
    final_list=df_final.values.tolist()         
    return final_list 


p=Pool(multiprocessing.cpu_count())
result_list = p.map(Graphk, file_list)

def from_iterable_indexadded(iterable):
    index=0
    for it in iterable:
        name_of_edges_file=file_list[index].split("_")
        if len(name_of_edges_file)>1:
            ele_chr="_".join(file_list[index].split("_")[0:-1])
        else:
            ele_chr=name_of_edges_file[0]
        for element in it:
            yield element, ele_chr
        index+=1
        
def flat(iterable):
    return map(list, zip(*from_iterable_indexadded(iterable)))
    
flat_result, corresponding_file = flat(result_list)

df = pd.DataFrame(zip(corresponding_file, flat_result), columns = ['name', 'result'])
df.to_pickle('graphk_updated_'+str(RBP)+'_'+str(k)+'mer.pkl')
