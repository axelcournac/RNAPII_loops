#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:24:15 2019
@author: axel KournaK
Spectral plots after quantifications with ChromoSight
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
# chromosight quantify --perc-undetected=100  --perc-zero=100  /home/axel/Bureau/YEAST/pairs_peaks_1.5_PolII_Log_1.5.txt.bg2 /media/axel/RSG51/diverse_yeast_data_copy/copy_cool_files_diverse_yeast_data/tsuki2_paper_out_repos_out_SRR13736654_55_tmp_valid_idx_pcrfree.pairs.2000.cool SRR13736654_55

df=pd.read_table('/home/axel/Bureau/Pol2_Project/loops_quantification_2022_spectre/SRR13736654_55.tsv',header=0, delimiter="\t")

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

min_number=15 # minimal number of events to do the average or median 

for i in range(N) :
    kb = df.loc[ ( (df['chrom1'] == df['chrom2']) & 
                  ((df["start2"] - df["start1"] >= i*bin_matrice ) & (df["start2"] - df["start1"] < (i+1)*bin_matrice ) )) ]
#    print(i, len(kb) )
    if len(kb) > min_number:
        size_tab[i] = i*bin_matrice
        mean_tab[i] = np.nanmean( kb['score'] )
        std_tab[i] = np.std( kb['score'] )
        no_tab[i] = len( kb['score'] )
    else:
        size_tab[i] = i*bin_matrice
        mean_tab[i] = np.nan
        std_tab[i] = np.nan
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


