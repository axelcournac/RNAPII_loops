#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

To have spectrum plot after quantifications with ChromoSight
"""
import numpy as np
import scipy
import matplotlib.pylab as plt
import pandas as pd
from skmisc.loess import loess

# before quantify the cohesin loops with chromosight : 
# chromosight quantify --pattern=loops_small --perc-zero=100 --perc-undetected=20  /home/axel/Bureau/YEAST/pairs_peaks_cohesins3.txt.bg2.all_sizes.2  /media/axel/RSG5/disk/copy_diverse_yeast/data_LChaptel/out_LCH10/tmp/valid_idx_pcrfree.pairs.2000.S288C.cool  LCH10

 
bin_matrice = 2000
df=pd.read_table('/home/axel/Bureau/python_codes/LCH10.tsv',header=0, delimiter="\t") 
T="LCH12"
           

dist = df["start2"] - df["start1"]
print( min(dist) )
print( max(dist) )
print( np.nanmax(df['score'])  )


N=80
mean_tab = np.zeros(N)
std_tab = np.zeros(N)
size_tab = np.zeros(N)
no_tab = np.zeros(N)

perc=0.10  # pourcentage of points that have a value

for i in range(N) :
    kb = df.loc[ ( (df['chrom1'] == df['chrom2']) & 
                  ((df["start2"] - df["start1"] >= i*bin_matrice ) & (df["start2"] - df["start1"] < (i+1)*bin_matrice ) )) ]
    print(i, len(kb) )
    size_tab[i] = i*bin_matrice
    mean_tab[i] = np.nanmedian( kb['score'] )
    std_tab[i] = np.std( kb['score'] )
    no_tab[i] = len( kb['score'] )
    s=kb['score']
    if len(s[np.isfinite(s)]) < perc * len( kb['score'] ):
        mean_tab[i] = np.nan
    
x = size_tab[np.isfinite(mean_tab)]
y = mean_tab[np.isfinite(mean_tab)]

#
#plt.plot(x,y,linewidth=3., color="red",label="N431")
#plt.plot(x,y, 'o')

#plt.xlabel("Size between two cohesin peaks")
#plt.ylabel("Median Loop score")
#
#plt.legend()
#plt.show()
#plt.grid()


# lowess from scikit
l = loess(x,y, span=0.18)   # do not take a big percentage
l.fit()
pred = l.predict(x, stderror=True)  #  time consuming step!
conf = pred.confidence()

lowess = pred.values
ll = conf.lower
ul = conf.upper

# 
#------------- plot: 
color1="red"

plt.title("Loop Spectrum")
plt.grid()
plt.plot(x, lowess, label=T, linewidth=3.0)
plt.fill_between(x,ll,ul,alpha=.33)

plt.xlabel("Size between two cohesin peaks")
plt.ylabel("Median Loop score")

plt.legend()
plt.savefig("spectrum.pdf")
plt.show()


