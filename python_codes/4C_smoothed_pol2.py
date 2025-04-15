#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 15:58:19 2024
@author: axel
Plot 4C like signals from Hi-C maps with some smart smoothing (lowess)
Plots for Cohesion paper
_2: with cohesin Chip-seq plot 
"""

import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import cooler
import matplotlib.gridspec as gridspec
import pyBigWig
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess as  sm_lowess

# Files

# # #fig1
# cool_file1="/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/fastq_Hi-C/out_NB4/tmp/valid_idx_pcrfree.pairs.2000.cool"
# name1="HiC pGALScc1V137K, Scc1WT"
# file="/media/axel/RSG51/diverse_yeast_data_copy/cohesion_paper/results_pGALScc1V137K-HA/tracks/NB22_nxq/NB22_nxq^mapped_W303_2micron_glabrata^HKU3QY.vs-NB26_nxq.bw"
# name11="ChIP-seq Scc1-pk-aid, pGALScc1V137K-HA, expressed inG1, IP anti-HA"

# cool_file1="/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/fastq_Hi-C/out_NB5/tmp/valid_idx_pcrfree.pairs.2000.cool"
# name1="HiC pGALScc1V137K, Scc1-pk-aid"
# file="/media/axel/RSG51/diverse_yeast_data_copy/cohesion_paper/results_pGALScc1WT-HA_calib/tracks/NB42_nxq/NB42_nxq_unmapped_glabrata_mapped_W303_2micron_AVU13Q.vs-NB46_nxq.bw"
# name11="ChIP-seq Scc1-pk-aid, pGALScc1WT-HA, expressed inG1, IP anti-HA"

# #fig2
# cool_file1="/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/fastq_Hi-C/out_NB31_39/tmp/valid_idx_pcrfree.pairs.2000.cool"
# name1="HiC pGALScc1V137K, Scc1-pk-aid, noco + oestradiol"
# file="/media/axel/RSG51/diverse_yeast_data_copy/cohesion_paper/last_version_tiny/results_with_Scc1_calib/tracks/NB45_nxq/NB45_nxq_unmapped_glabrata_mapped_W303_2micron_HL0VL8.CPM.calibrated.bw"
# name11="ChIP-seq pGALScc1V137K, Scc1-pk-aid, expressed inG2, IP anti-HA"

# cool_file1="/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/fastq_Hi-C/out_NB33_40/tmp/valid_idx_pcrfree.pairs.2000.cool"
# name1="HiC pGALScc1V137K, Scc1WT, noco + oestradiol"
# file="/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/last_version_tiny/results_without_Scc1_calib/tracks/NB44_nxq/NB44_nxq_unmapped_glabrata_mapped_W303_2micron_O5DJXA.CPM.calibrated.bw"
# name11="ChIP-seq SCC1WT, pGALScc1V137K, expressed inG2, IP anti-HA"




# sup fig 1  for GAL7 gene 
cool_file1="/media/axel/RSG5/disk/copy_diverse_yeast/gal_libraries/out_CH403_on_W303/tmp/valid_idx_pcrfree.pairs.2000.cool"
name1="CH403 FB220-1b nocodazole +IAA + oestradiol from G1 HiC"

file1="/media/axel/RSG5/disk/copy_diverse_yeast/gal_libraries/results_Pol2_b_from_G1/tracks/CH397_nxq/CH397_nxq^mapped_W303_2micron_glabrata^BWUIZ7.vs-CH401_nxq.bw"
name11="CH397 FB218-8d nocodazole +oestradiol from G1 IP pol2 ChIP"

cool_file2="/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/fastq_Hi-C/out_NB2/tmp/valid_idx_pcrfree.pairs.2000.cool"
name2="NB2 W303-1A nocodazole +IAA + oestradiol from G1  HiC"

file2="/media/axel/RSG5/disk/copy_diverse_yeast/gal_libraries/results_Pol2_from_G1/tracks/CH396_nxq/CH396_nxq^mapped_W303_2micron_glabrata^806WR1.vs-CH400_nxq.bw"
name22="CH396 FB218-4a nocodazole +oestradiol from G1 IP pol2 ChIP"

file3="//media/axel/RSG5/disk/copy_diverse_yeast/gal_libraries/results_W303_Scc1_b_from_G1/tracks/CH395_nxq/CH395_nxq^mapped_W303_2micron_glabrata^NHNSGP.vs-CH399_nxq.bw"
name33="Scc1 Galactose"


# # # fig 1  for centro 4C 
# cool_file1= "/media/axel/RSG5/disk/copy_diverse_yeast/data_LChaptel/out_LCH11/tmp/valid_idx_pcrfree.pairs.2000.W303.cool"   # W303
# name1="LCH11"

# file1="/media/axel/RSG5/disk/copy_diverse_yeast/data_LChaptel/chIP-seq/results_G2_arrested_cells_25C_1h_no_rapa_25C_1h_shift_35.5C/tracks/LCH25_nxq/LCH25_nxq^unmapped_glabrata^mapped_W303_2micron^1S68AG.vs-LCH24_nxq.bw"
# name11="Pol II"

# cool_file2="/media/axel/RSG5/disk/copy_diverse_yeast/data_LChaptel/out_LCH12/tmp/valid_idx_pcrfree.pairs.2000.W303.cool"
# name2="LCH12"

# file2="/media/axel/RSG5/disk/copy_diverse_yeast/data_LChaptel/chIP-seq/results_G2_arrested_cells_25C_1h_no_rapa_25C_1h_shift_35.5C/tracks/LCH25_nxq/LCH25_nxq^unmapped_glabrata^mapped_W303_2micron^1S68AG.vs-LCH24_nxq.bw"
# name22="Pol II"

# file3="/home/axel/Bureau/pds5_chip_seq/results_on_W303/tracks/NB63/NB63^mapped_W303_2micron^1TNB6J.vs-NB73.bw"    # cohesin in Pds5 mutant 
# name33="Scc1"


# centro file with their positions 
centro = pd.read_table("/media/axel/RSG5/disk/copy_diverse_yeast/cohesion_paper/genomes/centromeres.txt2",
                       header=None, delimiter="\t")
# Create a zip object from two lists
zipbObj = zip(centro[0], centro[1])
# Create a dictionary from zip object
centro_dict = dict(zipbObj)
centro_dict['chrM'] = 0
centro_dict['plasmid_p2-micron'] = 0


# Contact data
BIN=2000
c1 = cooler.Cooler(cool_file1)
c2 = cooler.Cooler(cool_file2)

cooler.balance_cooler(c1,store=True, mad_max=10)   # Normalisation
cooler.balance_cooler(c2,store=True, mad_max=10)   # Normalisation

# chr1="chrXI"
# chr2="chrXI"
# position=446464

# function for nan sum
def nan_sum(row):
    total = np.nanmean(row)*len(row[~np.isnan(row)])
    return total

# GAL7
chr1="chrII"
chr2="chrII"
position=274427

# chr1="chrX"
# chr2="chrX"
# position=centro_dict[chr1]

# chr1="chrII"
# chr2="chrII"
# position=centro_dict[chr1]

# binning and simple 4C with the matrice
pos1=int(position/BIN)
m1 = c1.matrix(balance=True).fetch(chr1, chr2)
s1=m1[pos1,:]

m2 = c2.matrix(balance=True).fetch(chr1, chr2)
s2=m2[pos1,:]

# region of focus and use of fetch function
area=2000
chr3 = (chr1,position-area,position+area) # region of focus  

m12 = c1.matrix(balance=True).fetch(chr1, chr3)
reads_chr12= m12.sum()
coverage_plasmid = np.apply_along_axis(nan_sum, 1, m12)
# coverage_plasmid [ np.isnan(coverage_plasmid1)] = np.nan

m12 = c2.matrix(balance=True).fetch(chr1, chr3)
reads_chr12= m12.sum()
coverage_plasmid2 = np.apply_along_axis(nan_sum, 1, m12)

s1=coverage_plasmid
s2=coverage_plasmid2

# ChIP-seq data
bw = pyBigWig.open(file1)
bw2 = pyBigWig.open(file2)
bw3 = pyBigWig.open(file3)

i=np.array(bw.intervals(chr1) )
i2=np.array(bw2.intervals(chr1) )
i3=np.array(bw3.intervals(chr1) )
 
x1=i[:,0]
x2=i2[:,0]
x3=i3[:,0]

v1=i[:,2]
v2=i2[:,2]
v3=i3[:,2]

# # Smooth 
# v1s=np.convolve(v1, np.ones(50)/50, mode='same')
# v2s=np.convolve(v2, np.ones(50)/50, mode='same')
# v3s=np.convolve(v3, np.ones(50)/50, mode='same')

#------------------------------------------------------------------------------
# Plot
plt.figure(figsize=(24,8))

l1=range(len(s1))
l2=range(len(s2))

gs = gridspec.GridSpec(3, 1,height_ratios=[1,1,1]);
ax1 = plt.subplot(gs[0])

# ax1.bar(l1,s1,width=1, alpha=.99,color='purple')
# ax1.bar(l1,s2,width=1, alpha=.6,color='orchid')

# ax1.bar(l1,s1,width=1, alpha=.99,color='darkorange')
# ax1.bar(l1,s2,width=1, alpha=.6,color='bisque')

ax1.bar(l1,s1,width=1, alpha=.99,color='orange')
ax1.bar(l1,s2,width=1, alpha=.6,color='royalblue')

ax1.spines[['right', 'top']].set_visible(False)
plt.xlabel("Genomic coordinates (bin 2kb)")
plt.ylabel("4C like score")
plt.title(chr1+" from "+name1+"\n"+"4Clike from centro")
plt.xlim(0,np.max(l1)*1.)
# plt.ylim(0,0.05)
plt.ylim(0,0.007)
plt.legend(loc='upper right')



ax2 = plt.subplot(gs[1], sharex=ax1)
ax2.spines[['right', 'top']].set_visible(False)

plt.plot(x1/BIN, v1, color="tomato", alpha=.99)
# plt.fill_between(x1/BIN,0,v1, color="tomato",  alpha=.99, label="Pol II ChIP-seq")
# plt.bar(x1/BIN,v1, color="tomato",  alpha=.99, label="Pol II ChIP-seq")


plt.xlabel("Genomic coordinates (bin 2kb)")
plt.ylabel("ChIP-seq IP/input")
plt.title(chr1+" from "+name11)

plt.xlim(0,np.max(x1/BIN)*1.) 
plt.ylim(0,np.max(v1)*1.1)
plt.tight_layout()
plt.legend(loc='upper right')

ax2 = plt.subplot(gs[2], sharex=ax1)
ax2.spines[['right', 'top']].set_visible(False)

plt.plot(x3/BIN, v3, color="tomato", alpha=.99)
# plt.fill_between(x3/BIN,0,v3s, color="grey",  alpha=.99, label="Scc1 ChIP-seq")
# plt.bar(x3/BIN,v3, color="grey",  alpha=.99, label="Scc1 ChIP-seq")
   
plt.xlabel("Genomic coordinates (bin 2kb)")
plt.ylabel("ChIP-seq IP/input")
plt.title(chr1+" from "+name33)

plt.xlim(0,np.max(x3/BIN)*1.) 
plt.ylim(0,np.max(v3)*1.1)
plt.tight_layout()
plt.legend(loc='upper right')

plt.tight_layout()

plt.savefig(name1+".pdf", dpi=400, format='pdf')
plt.savefig(name1+".svg", dpi=400, format='svg')






