import numpy as np
import sys
import glob
rbp=sys.argv[1]
kmer=sys.argv[2]

pfile_list=glob.glob("result_VDM3_"+rbp+"_positive_"+kmer+"_*.npy")

pfile1=np.load(pfile_list[0])
psha=np.shape(pfile1)
pmatrix=np.zeros(psha)
for pfile in pfile_list:
    file=np.load(pfile)
   # file=np.fromfile(pfile,dtype=np.float32)
    pmatrix+=file
np.save("positive_"+rbp+"_vdm3_nopaircontrol_distance_matrix_"+kmer+"mer.npy",pmatrix)

