import pandas as pd
import numpy as np
import heapq
from multiprocessing import Pool
import multiprocessing
import sys,os
from scipy.stats import rankdata
from heapq import nlargest
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist, jaccard
import networkx as nx
from vdm3 import ValueDifferenceMetric
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import logomaker
import itertools
import collections
import operator
from operator import itemgetter


os.chdir(sys.argv[1])
file_name = str(sys.argv[2])
kmer=str(sys.argv[3])
segment = int(sys.argv[4])
total_seg = int(sys.argv[5])
par1=float(sys.argv[6])
par2=int(sys.argv[7])
fafile=str(sys.argv[8])

k=int(kmer)

dfp=pd.read_pickle(file_name+"_"+kmer+"mer_graphedgewithfeature.pkl")
dfpf=dfp['feature']
dfpn=dfp['name']

listpf=dfpf.tolist()
listpn=dfpn.tolist()

X2=[]
y2=[]
namename=[]

##### centered
for i in range(0,len(listpf)):
    if (listpf[i][0]>=30 and listpf[i][1]+k-1<=70 and int(listpf[i][3])==-k) or (listpf[i][0]-k+1>=30 and listpf[i][1]<=70 and int(listpf[i][3]==-1)):
        X2.append(listpf[i])
        namename.append(listpn[i])

centered_feature=X2

        
pre1="#".join('P{}#E{}#H{}#I{}#M{}'.format(x+1,x+1,x+1,x+1,x+1) for x in range(0,k))
final1=pre1.split("#") 
pre2="#".join('P{}#E{}#H{}#I{}#M{}'.format(x+1,x+1,x+1,x+1,x+1) for x in range(k,2*k))
final2=pre2.split("#") 
df=pd.DataFrame(X2,columns=['pos1','pos2','weight','type1','type2']+['seq{}'.format(x+1) for x in range(0,k)]+final1+['seq{}'.format(x+1) for x in range(k,2*k)]+final2)        
    
pos1_cat=pd.cut(df["pos1"], np.arange(30.0, 71.0, 10.0),labels=False)
pos2_cat=pd.cut(df["pos2"], np.arange(30.0, 71.0, 10.0),labels=False)
pos1_str=[str(pos1_label)for pos1_label in pos1_cat]
pos2_str=[str(pos2_label)for pos2_label in pos2_cat]
pos_str=[pos1_str[i]+pos2_str[i] for i in range(0,len(pos1_str))]
y2=np.array(pos_str)
#y2=np.array(list(zip(np.array(pos1_cat), np.array(pos2_cat))))//cannot be 2-d list

X2=np.array(df.iloc[:,5:])

case = ValueDifferenceMetric(X=X2, y=y2, continuous=list(range(k,6*k))+list(range(7*k,12*k)))
case.fit()

def process(indices):
    i,j = indices[0],indices[1]
    return case.get_distance(ins_1=X2[i],ins_2=X2[j])
 
index_list = []
 
start = int((segment)*len(X2)/total_seg)
end = int((segment+1)*len(X2)/total_seg)
 
for i in range(start,end):
    for j in range(i+1,len(X2)):
        index_list.append((i,j))

p = Pool(multiprocessing.cpu_count())    
result = p.map(process, index_list)

Distance_Matrix=np.zeros((len(X2), len(X2)))
index = 0

for i in range(start,end):
    for j in range(i+1,len(X2)):
        Distance_Matrix[i][j]=result[index]
        Distance_Matrix[j][i]=result[index]
        index +=1

np.save("result_VDM3_"+file_name+"_"+kmer+"_"+str(segment)+".npy",Distance_Matrix)
#print(Distance_Matrix)


p_file=np.load("result_VDM3_"+file_name+"_"+kmer+"_0.npy")
row, col = np.diag_indices_from(p_file)
p_file[row,col] = 1000
df=pd.DataFrame(p_file)

df_newnew=pd.DataFrame(np.where(df.rank(axis=1,method='min')>par1*0.005*len(df), 0, df),columns=df.columns)
overlap=(df_newnew.T==df_newnew)*df_newnew
indices = np.argwhere(overlap.values != 0)

res = pdist((df_newnew != 0).astype(int), 'jaccard')
WW=1-squareform(res)
UU=np.transpose(indices)
weight_new=WW[UU[0],UU[1]]


#indices_new=np.c_[indices,weight_new]  
#indices_new=np.concatenate((indices,np.reshape(weight_new,(-1, 1))),axis=1)
indices_new=np.column_stack([indices,np.reshape(weight_new,(-1, 1))])
df_df=pd.DataFrame(indices_new,columns=('node1','node2','weight'))
df_df_df=df_df[df_df.weight > np.percentile(df_df.weight,100-par2*10)]

G = nx.Graph([e for e in df_df_df.values[:,0:2]])
cliques=list(nx.find_cliques(G))

sorted_c=sorted(cliques,key=len)
final_c = [i for i in sorted_c if len(i)>=0.5*max(map(len,sorted_c))]

np.save(file_name+"_positive_"+str(par1)+"_"+str(par2)+"_"+kmer+"mer_nopaircontrol.npy",final_c)


final_A11=final_c
#final_A11=np.load(file_name+"_positive_"+str(par1)+"_"+str(par2)+"_"+kmer+"mer_nopaircontrol.npy",allow_pickle=True)

D=dict(collections.Counter([x for sublist in final_A11 for x in sublist]))
top50D=list(dict(sorted(D.items(), key=lambda item: item[1])).keys())[-50:]
#print(D)
#print(top50D)
f = operator.itemgetter(*list(map(int,top50D)))


LIST_NAME=f(namename)
LIST_F=f(centered_feature)

seq_bp=pd.read_pickle(file_name+"_seq_basepair.pkl")

LIST_N=[]
LIST_S=[]
LIST_PO=[]
EKS_label=[]
EKS_type=[]
SEG1_SEQ=[]
SEG2_SEQ=[]

SEG1_Structure_Prob=[]
SEG2_Structure_Prob=[]


for i in range(len(LIST_NAME)):
    name_pre=LIST_NAME[i]
    LIST_N.append(name_pre)
    ind=seq_bp.index[seq_bp['chr_start_end'] == name_pre].tolist()[0]
    LIST_S.append(seq_bp['sequence'][ind])
    LIST_PO.append((int(LIST_F[i][0]),int(LIST_F[i][1])))
    EKS_label.append(LIST_F[i][2])
    if LIST_F[i][3]==-1.0:
        EKS_type.append("left-opened")
    else:
        EKS_type.append("right-opened")
    SEG1_seq=LIST_F[i][5]
    for nn in range(1,k):
        SEG1_seq=SEG1_seq+LIST_F[i][5+nn]
    SEG1_SEQ.append(SEG1_seq)
    SEG1_Structure_prob=[]
    for mm in range(0,5*k):
        SEG1_Structure_prob.append(LIST_F[i][5+k+mm])
    SEG1_Structure_Prob.append(SEG1_Structure_prob)
    
    SEG2_seq=LIST_F[i][5+6*k]
    for nn in range(1,k):
        SEG2_seq=SEG2_seq+LIST_F[i][5+6*k+nn]
    SEG2_SEQ.append(SEG2_seq)
    SEG2_Structure_prob=[]
    for mm in range(0,5*k):
        SEG2_Structure_prob.append(LIST_F[i][5+7*k+mm])
    SEG2_Structure_Prob.append(SEG2_Structure_prob)

dfdf=pd.DataFrame(list(zip(LIST_N,LIST_S,LIST_PO,EKS_label,EKS_type,SEG1_SEQ,SEG1_Structure_Prob,SEG2_SEQ,SEG2_Structure_Prob)),columns=['chr_start_end','sequence','(first'+str(k)+'merstartpos,second'+str(k)+'merstartpos)','EKS_label','EKS_type','segment1_sequence','segment1_Structure_Prob(P1,E1,H1,I1,M1,P2,E2,H2,I2,M2,...)','segment2_sequence','segment2_Structure_Prob(P1,E1,H1,I1,M1,P2,E2,H2,I2,M2,...)'])
#dfdf.to_pickle(RBP+"_Motif_"+str(k)+".pkl")
dfdf.to_csv(file_name+"_EKS_Motif_k="+str(k)+".csv",index=False)
