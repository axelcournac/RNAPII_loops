# -*- coding: utf-8 -*-
"""
@author: Axel KournaK
Select Chip-seq peaks to do after pileup plots with contact data like Hi-C or microC
"""
import numpy as np
import matplotlib
import pandas as pd
#import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

import os

os.chdir('/home/axel/Bureau/z_python_scripts_copy')
import itertools
import matplotlib.backends.backend_pdf

# Various 1D analog signals from ChIP-seq or other already tried:
# ---------------------------------------------------------------------------------------------------------------------- 
df4=pd.read_table('/home/axel/Bureau/polII_peaks_2024/SRR2065097.fastq.sam.MQ30.FigS5_Scc1PK9_IP_G1_releasing_60min',header=None, delimiter=" ")
df8=pd.read_table('/home/axel/Bureau/polII_peaks_2024/SRR2065092.fastq.sam.MQ30.FigS5_Scc1PK9_WCE_G1_releasing_60min',header=None, delimiter=" ")
name_bank="cohesin_classic"

df44=pd.read_table('/home/axel/Bureau/polII_peaks_2024/SRR7175393.fastq.sam.MQ30.Rpb3_in_Log_All_Replicates_Merged_IP',header=None, delimiter=" ")
df88=pd.read_table('/home/axel/Bureau/polII_peaks_2024/SRR7175394.fastq.sam.MQ30.Rpb3_in_Log_All_Replicates_Merged_Input',header=None, delimiter=" ")
name_bank="PolII_classic"

#------------
len(df4)
list_all_chrms= ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                 "chr10","chr11","chr12","chr13","chr14","chr15","chr16")

# we take the positions of the signal for each chromosome 
ka={}
kb={}
kaa={}
kbb={}

for c in  list_all_chrms:
    print(c)
    ka[c] = df4.loc[(df4[0] == c)]
    kb[c] = df8.loc[(df8[0] == c)]
    kaa[c] = df44.loc[(df44[0] == c)]
    kbb[c] = df88.loc[(df88[0] == c)]

# Histogram of the signals 
BIN_histo = 2000 #  bin for 1D histogram
hist1={}
hist2={}
bins1={}
bins2={}

hist11={}
hist22={}
bins11={}
bins22={}

for c in list_all_chrms:
    print(c)
    va= np.array(ka[c][1])
    va = [ int(x) for x in va ]
    vb= np.array(kb[c][1])
    vb = [ int(x) for x in vb ]
    max_reached = int( max(max(vb), max(va) ) )
    
    hist1[c], bins1[c] = np.histogram(va,bins= range(0,max_reached,BIN_histo),density="True")
    hist2[c], bins2[c] = np.histogram(vb,bins= range(0,max_reached,BIN_histo),density="True")
    
    va= np.array(kaa[c][1])
    va = [ int(x) for x in va ]
    vb= np.array(kbb[c][1])
    vb = [ int(x) for x in vb ]
    hist11[c], bins11[c] = np.histogram(va,bins= range(0,max_reached,BIN_histo),density="True")
    hist22[c], bins22[c] = np.histogram(vb,bins= range(0,max_reached,BIN_histo),density="True")

# Filtering & plots: 
c = np.loadtxt("/home/axel/Bureau/YEAST/centro1.dat3")      
i=0 
threshold = 1.5 # above which we consider to have a peak 
threshold2 = 1.5

area = 5   # area around centro 
f_out = open("positions_peaks_"+str(threshold)+"_"+name_bank+".txt","w+")
f_out2 = open("pairs_peaks_short_"+str(threshold)+"_"+name_bank+".txt","w+")
f_out3 = open("pairs_peaks_long_"+str(threshold)+"_"+name_bank+".txt","w+")
BIN= 2000

for chr in list_all_chrms:
    i+=1
    plt.figure(i)
    print(chr)  
    b=   bins1[chr][:-1]
    b2=  bins11[chr][:-1]    
    h =  hist1[chr] / hist2[chr]
    h2 = hist11[chr] / hist22[chr]
    
    plt.title(chr)
    
    plt.plot(b,h,label="cohesin")  #    plot(b,h,label="ChIP-IP-183-2")
    plt.axhline(y= threshold, color='black', linestyle='dashed')
    plt.xlabel("Positions along the chromosome")
    plt.ylabel("ChIP / input")
    plt.legend()
#    plt.ylim(0,10)
    
    # basic Selection of peaks:
    centro = i
    centro_binned = int(c[centro-1][1]/ BIN)      
    bs2 = b[ h > threshold ]   #  selection here !!
    hs2 = h[ h > threshold ]
    plt.plot(bs2,hs2,'o',label="near centro")
    
    indices = np.where(bs2>=0)
    
#   removing of bins clothed to centro :
    indices = np.where( np.logical_or( (bs2 < (centro_binned -area)* BIN ), (bs2 > (centro_binned + area) *BIN ) ) ) 
    if centro ==12:
        rdna = 458991
        rdna_binned = int( rdna /BIN)
        indices1 = np.where( np.logical_or( (bs2 < (centro_binned -area)* BIN ), (bs2 > (centro_binned + area) *BIN ) ) ) 
        indices2 = np.where( np.logical_or( (bs2 < (rdna_binned -area * 3 )* BIN ), (bs2 > (rdna_binned + area * 3 ) *BIN ) ) ) 
        indices = np.intersect1d(indices1, indices2 )
        
    bs2 =  bs2[indices]
    hs2 =  hs2[indices]
    plt.plot(bs2,hs2,'o',color="yellow",label="selected")
    
    plt.xlabel("Positions along the chromosome")
    plt.ylabel("ChIP / input")
    plt.legend()
#    plt.ylim(0,10)
    
    # Generation of pairs positions files:
    distance_min1 = 10000    # minimal distance to consider for the pair
    distance_max1 = 50000    # maximal distance to consider for the pair

    distance_min2 = 50000    # minimal distance to consider for the pair
    distance_max2 = 450000    # maximal distance to consider for the pair
    
    for p in bs2 :  #  to write positions of peaks
        f_out.write(chr +'\t'+ str( int([p][0]) ) + '\t' +str( int([p][0])+BIN )  +'\t'+
                    chr +'\t'+ str( int([p][0]) ) + '\t' +str( int([p][0])+BIN )  +'\n')
    combi_pos = list(itertools.combinations(bs2,2) )
    for p in range(0, len(combi_pos) ):    #  to write pairs of positions of peaks and potential loops
        pos1=int(combi_pos[p][0])
        pos2=int(combi_pos[p][1]) 
        if np.abs(pos2-pos1)< distance_max1 and np.abs(pos2-pos1)>distance_min1: 
            f_out2.write(chr +'\t'+ str(pos1) + '\t' +str(pos1+BIN)  +'\t'+
                        chr +'\t'+ str(pos2) + '\t' +str(pos2+BIN)  +'\n')
        
        if np.abs(pos2-pos1)< distance_max2 and np.abs(pos2-pos1)>distance_min2: 
            f_out3.write(chr +'\t'+ str(pos1) + '\t' +str(pos1+BIN)  +'\t'+
                        chr +'\t'+ str(pos2) + '\t' +str(pos2+BIN)  +'\n')
        
        
pdf = matplotlib.backends.backend_pdf.PdfPages(name_bank+"_"+str(threshold)+".pdf")
for fig in range(1, plt.figure().number): ## will open an empty extra figure :(
    pdf.savefig( fig )
pdf.close()
plt.close('all')
    
f_out.close() 
f_out2.close()
f_out3.close()   

