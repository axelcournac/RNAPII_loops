# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 19:03:08 2018
@author: axel KournaK 
do an agglomerated plot as we did  
cool: adapted to cool data 
"""
import matplotlib.pyplot as plt
import pandas as pd
import random as rando
import sys
import os
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts"))
import hicstuff as hcs
import numpy as np
import cooler

# Contact data:
#cool_file = "/home/axel/Bureau/all_cool_files_renamed/aragon_2019_march_seqs_out__AT327_tmp_valid_idx_pcrfree.pairs.cool"
cool_file = sys.argv[1]   # name of the cool file 
name_bank = sys.argv[2]   # name of the bank

# Input of pairs: 
#file_bg2='/home/axel/Bureau/YEAST/pairs_peaks_condensins3.txt.bg2.10kb.50kb.2'
file_bg2 = sys.argv[3]    # bg2 file containing the positions of pairs of peaks
name_peak = sys.argv[4]    # name we want to give to the protein 


df=pd.read_table(file_bg2,header=None, delimiter="\t")  

# Parameters
bin_matrice = 2000 # size of bin of matrice (in bp)
n_random = 10
area = 10  #  number of bins around the sites

# Loading of contact data:
c = cooler.Cooler(cool_file)
#cooler.balance_cooler(c, store=True, mad_max=10)   # Normalisation 

# Loading and normalisation of contact matrices:
list_all_chrms=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9'
                ,'chr10','chr11','chr12','chr13','chr14','chr15','chr16')
matscn = {}
matraw = {}
th_sum = {}
indices_matrice = {}
for chr1 in list_all_chrms : 
    print(chr1)
    matraw_t = c.matrix(balance=False).fetch(chr1, chr1)
#    matscn[c] = hcs.normalize_sparse(matraw,norm="SCN")
#    matscn[chr1] = matraw
    matscn[chr1] = hcs.normalize_dense(matraw_t, norm="SCN",order=2)
#    matscn[c] = matscn[c].tolil()
    matraw_t = c.matrix(balance=False).fetch(chr1, chr1)
    matraw[chr1] = matraw_t
    th_sum[chr1] = np.median( np.array(matraw_t.sum(axis=0))) - 1.* np.std(matraw_t.sum(axis=0))
    indices_matrice[chr1] = np.where( matraw_t.sum(axis=0) > th_sum[chr1] )

    
#  Computation of agglo 
print(len(df) )
df1=df
n_pos_set = 0

MAT_SUM1 = np.zeros(   (area*2+1,area*2+1)  )
MAT_SUM11 = np.zeros(   (area*2+1,area*2+1)  )
MAT_OCC1 = np.zeros(   (area*2+1,area*2+1)  )

for chr1 in list_all_chrms :
    print(chr1)
    n1 = matscn[chr1].shape[0]
    print("Number of bins and after filtering poor interacting bins:")
    print(n1)
    print(len(indices_matrice[chr1][0]))
    pos_set = df1.loc[(df1[0] == chr1)]
    pos_set = np.array(pos_set)
    n_pos_set = n_pos_set + pos_set.shape[0]
    ns =0
    for pos_i in range(pos_set.shape[0] ) :
        site1 = int(pos_set[pos_i,1]/bin_matrice)
        site2 = int(pos_set[pos_i,4]/bin_matrice)
        ns +=1
        pi =0
        for i in range(site1-area,site1+area+1) :
            pj =0
            for j in range(site2-area,site2+area+1) :
                if i >=0 and j>=i and i<n1 and j<n1 and (i in 
                indices_matrice[chr1][0]) and (j in indices_matrice[chr1][0]):
                    MAT_SUM1[pi,pj] = MAT_SUM1[pi,pj] + matscn[chr1][i,j]
                    MAT_SUM11[pi,pj] = MAT_SUM11[pi,pj] + matraw[chr1][i,j]
                    MAT_OCC1[pi,pj] = MAT_OCC1[pi,pj] + 1.0
                pj +=1
            pi +=1     
MAT_SUM1 = MAT_SUM1 / MAT_OCC1
MAT_SUM11 = MAT_SUM11 / MAT_OCC1
      
# Random part: 
MAT_SUM2 = np.zeros(   (area*2+1,area*2+1)  )
MAT_SUM22 = np.zeros(   (area*2+1,area*2+1)  )
for r in range(n_random) :
    print(r)
    MAT_SUM_scn = np.zeros(   (area*2+1,area*2+1)  )
    MAT_SUM_raw = np.zeros(   (area*2+1,area*2+1)  )
    MAT_OCC = np.zeros(   (area*2+1,area*2+1)  )
    for chr1 in list_all_chrms :
        n1 = matscn[chr1].shape[0]
        pos_set = df1.loc[(df1[0] == chr1)]
        pos_set = np.array(pos_set)
        ns =0
        for i in range(pos_set.shape[0] ) :
            site1 =  rando.sample(range(0,n1),1)[0] # we picked a bin in the chr
            site2 = site1 + int((pos_set[i,4] - pos_set[i,1] )/bin_matrice) 
            ns +=1
            pi =0
            for i in range(site1-area,site1+area+1) :
                pj =0
                for j in range(site2-area,site2+area+1) :
                    if i >=0 and j>=i and i<n1 and j<n1 and (i in
                    indices_matrice[chr1][0]) and (j in indices_matrice[chr1][0]) :
                        MAT_SUM_scn[pi,pj] = MAT_SUM_scn[pi,pj] + matscn[chr1][i,j]
                        MAT_SUM_raw[pi,pj] = MAT_SUM_raw[pi,pj] + matraw[chr1][i,j]
                        MAT_OCC[pi,pj] = MAT_OCC[pi,pj] + 1.0
                    pj +=1
                pi +=1    
    MAT_SUM_scn = MAT_SUM_scn / MAT_OCC            
    MAT_SUM_raw = MAT_SUM_raw / MAT_OCC
    
    MAT_SUM2 = MAT_SUM2 + MAT_SUM_scn
    MAT_SUM22 = MAT_SUM22 + MAT_SUM_raw
    
MAT_SUM2 = MAT_SUM2 / n_random
MAT_SUM22 = MAT_SUM22 / n_random
        
print("n_pos_set")
print(n_pos_set)

# --------------PLots:--------------------------------------------------------
plt.figure(1)
filename = name_bank +"_"+ name_peak + "_scn" 
print(filename)

plt.imshow(np.log(MAT_SUM1/MAT_SUM2),interpolation="none",vmin=-0.2,vmax=0.2,cmap="bwr")
#    plt.imshow(MAT_SUM1**0.1,interpolation="none",cmap="afmhot_r",vmin=0.335,vmax=0.370)
plt.colorbar()
tick_locs=(0,area*2)
tick_lbls=('-'+str( int(area*bin_matrice/1000) )+' kb',
           '+ '+ str( int(area*bin_matrice/1000) )+'kb')

plt.xticks(tick_locs, tick_lbls,fontsize=20)
plt.yticks(tick_locs, tick_lbls,fontsize=20)
plt.xlabel("Nb= "+str(n_pos_set),fontsize=10  )
plt.title(name_bank +" "+ name_peak+"\n",fontsize=10)

plt.savefig(filename+".pdf")
np.savetxt(filename+"_"+str( int(area*bin_matrice/1000)) +".txt", MAT_SUM1/MAT_SUM2)

plt.figure(2)
filename = name_bank +"_"+ name_peak + "_raw" 
print(filename)

plt.imshow(np.log(MAT_SUM11/MAT_SUM22),interpolation="none",vmin=-0.2,vmax=0.2,cmap="bwr")
#    plt.imshow(MAT_SUM1**0.1,interpolation="none",cmap="afmhot_r",vmin=0.335,vmax=0.370)
plt.colorbar()
tick_locs=(0,area*2)
tick_lbls=('-'+str( int(area*bin_matrice/1000) )+' kb',
           '+ '+str( int(area*bin_matrice/1000) )+'kb')

plt.xticks(tick_locs, tick_lbls,fontsize=10)
plt.yticks(tick_locs, tick_lbls,fontsize=10)
plt.xlabel("Nb= "+str(n_pos_set)  )
plt.title(name_bank +" "+ name_peak+"\n",fontsize=10)

plt.savefig(filename+".pdf")
#np.savetxt(filename+"_"+str( int(area*bin_matrice/1000)) +".txt", MAT_SUM11/MAT_SUM22)
#    np.savetxt(filename+"_"+str( int(area*bin_matrice/1000) )+".txt", MAT_SUM1)    
plt.close('all')

print("ALL the home-made agglomerated plots are finished!")
