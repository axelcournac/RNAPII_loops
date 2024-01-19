# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 18:15:46 2018
@author: axel KournaK 
To create all possible pairs from genomic positions.
"""
import numpy as np
import pandas as pd
import itertools
import sys

# Input: list of genomic positions :
#file_prot='/home/axel/Bureau/figure_all_chrms_sup1/redone_norm3/th_08/HSC_plasmids_in_Micro-C_WT_log_SC288_genome.txt.sort.formated3.sorted'
file_prot=sys.argv[1]
name_prot=sys.argv[2] 

df=pd.read_table(file_prot,header=None, delimiter="\t")

print("Number of peaks in the initial file:")
print(len(df))
#name_prot="plasmid_2mu"

BIN=2000
df2 = df[[0, 1]]
df2[1] = df2[1].apply(lambda x: int(x/BIN))
df3=df2.drop_duplicates()
print("Number of peaks in the initial file after removing of duplicates:")
print(len(df3))
df3[1] = df3[1].apply(lambda x: int(x*BIN))
 
#------------------------------------------------------------------------------
# Intra and inter pairs :

f_out1 = open("pairs_intra_"+name_prot+".txt","w+")
f_out2 = open("pairs_inter_"+name_prot+".txt","w+")

len(df)
list_pairs = list(itertools.combinations( range(len(df3)) ,2 ))
len(list_pairs)    

distance_min = 10000    # minimal distance to consider for the pair
distance_max = 50000    # maximal distance to consider for the pair

for e in list_pairs :
    e1, e2 = e
    chr1, pos1 = df3.iloc[e1][0],df3.iloc[e1][1]
    chr2, pos2 = df3.iloc[e2][0],df3.iloc[e2][1]
    if chr1 == chr2 and np.abs(pos2-pos1)< distance_max and np.abs(pos2-pos1)>distance_min: 
        f_out1.write(chr1 + '\t'+ str( int(pos1)) + '\t' + str( int(pos1+BIN))+ '\t' + 
                     chr1 + '\t'+ str( int(pos2)) + '\t' + str( int(pos2+BIN))+ '\t' + "1" +'\n')
    if chr1 != chr2:
        f_out2.write(chr1 + '\t'+ str( int(pos1) )+ '\t' + str( int(pos1+BIN))+ '\t' + 
                     chr2 + '\t'+ str( int(pos2)) + '\t' + str( int(pos2+BIN))+ '\t' + "1" +'\n')

f_out1.close()
f_out2.close()


