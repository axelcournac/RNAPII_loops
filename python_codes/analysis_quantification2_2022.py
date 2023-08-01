#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:24:15 2019
@author: axel KournaK
To have some plots after quantifications with ChromoSight
je reprends cette analyse en 2022 pour la répliquer avec Pol II
"""
import numpy as np
import scipy
from scipy.sparse import coo_matrix 
from scipy.sparse import csr_matrix 
import time
import itertools
import matplotlib.pylab as plt
import json
import pandas as pd
from scipy import stats
import matplotlib.gridspec as gridspec
import seaborn as sns
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from skmisc.loess import loess

# before: 
#chromosight quantify --pattern=loops_small /home/axel/Bureau/YEAST/pairs_peaks_condensins3.txt.bg2.all_sizes.2 /home/axel/Bureau/loops_quantification_2022/cool_files/SRR7706227_SRR7706226_hic_scer_mitotic_2kb.cool SRR7706227_SRR7706226_hic_scer_mitotic_2kb
#chromosight quantify --pattern=loops /home/axel/Bureau/YEAST/pairs_peaks_1.5_PolII_Log_1.5.txt.bg2 /media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_Micro-C_WT_log_redone/tmp/valid_idx_pcrfree.pairs.cool quiescence_paper_log_2

# input contact data: 
df=pd.read_table('/home/axel/Bureau/loops_quantification_2022/SRR7706227_SRR7706226_hic_scer_mitotic_2kb.tsv',header=0, delimiter="\t") 
df=pd.read_table('/home/axel/Bureau/loops_quantification_2022/SRR7706227_SRR7706226_hic_scer_mitotic_2kb_polII_2.tsv',header=0, delimiter="\t") 
df=pd.read_table('/home/axel/Bureau/loops_quantification_2022/quiescence_paper_log_2.tsv',header=0, delimiter="\t") 
df=pd.read_table('/home/axel/Bureau/loops_quantification_2022/SRR8718853_Hi-C_Wildtype_small.tsv',header=0, delimiter="\t")
df=pd.read_table('/home/axel/Bureau/loops_quantification_2022_spectre/all_centro_paper.tsv',header=0, delimiter="\t")

dist = df["start2"] - df["start1"]
print( min(dist) )
print( max(dist) )
print( max(df['score']) )

bin_matrice=10000

N=100
mean_tab = np.zeros(N)
std_tab = np.zeros(N)
size_tab = np.zeros(N)
no_tab = np.zeros(N)

for i in range(N) :
    kb = df.loc[ ( (df['chrom1'] == df['chrom2']) & 
                  ((df["start2"] - df["start1"] >= i*bin_matrice ) & (df["start2"] - df["start1"] < (i+1)*bin_matrice ) )) ]
    print(i, len(kb) )
    size_tab[i] = i*bin_matrice
    mean_tab[i] = np.nanmedian( kb['score'] )
    std_tab[i] = np.std( kb['score'] )
    no_tab[i] = len( kb['score'] )

mean_tab[np.isnan(mean_tab)] = 0
    
# plot with lowess
x = size_tab
y = mean_tab

plt.plot(x,y,'o')
plt.plot(x,y)

l = loess(x,y, span=0.15)   # do not take a big percentage
l.fit()
pred = l.predict(x, stderror=True)  #  time consuming step!
conf = pred.confidence()

lowess = pred.values
ll = conf.lower
ul = conf.upper

#------------- plot: 
color1="tomato"
T="all_centro_paper"
plt.title("Loop Spectrum")
plt.grid()
plt.plot(x, lowess, color = color1, label=T, linewidth=3.)
plt.fill_between(x,ll,ul,alpha=.33, color = color1)

plt.xlabel("Size between two PolII sites")
plt.ylabel("Median Loop score")

plt.legend()
plt.show()


