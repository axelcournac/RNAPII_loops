#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:08:57 2020
@author: axel
To extract plasmid p2 micron vs chrms from cool files. 
2: adding of ChIPseq signals. 
"""
import time
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random
import sys
import os
from scipy.stats.stats import spearmanr
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))
import scn 
import ice_mirny3
import scipy
import scipy.ndimage
import scipy.io as sio
import distance_law_human
import hicstuff as hcs
import numpy as np
import json
import sys
import chromosight.utils.plotting as cup
import pandas as pd
import cooler
import function_geno

# contact data:

#cool_file1= "/media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_Micro-C_WT_log_redone/tmp/valid_idx_pcrfree.pairs.cool"
#name_bank1= "log tsuky1"
#
##cool_file1= "/media/axel/RSG4/diverse_yeast_data/AC/fastq_files/seqs/out2_WT_AC_2micron/tmp/valid_idx_pcrfree.pairs.cool"
##name_bank1= "AC"
#
#cool_file1= "/media/axel/RSG4/diverse_yeast_data/quiescence_2019/fastq/out_SMC4-off/tmp/valid_idx_pcrfree.pairs.cool"
#name_bank1= "SMC_OFF"

# Tsuky data 
cool_file1= "/media/axel/RSG4/diverse_yeast_data/tsuki2_paper/out_repos/out_SRR13736655/tmp/valid_idx_pcrfree.pairs.2000.cool"
name_bank1= "SRR13736655-log-tsuky2"

cool_file2= "/media/axel/RSG4/diverse_yeast_data/tsuki2_paper/out_repos/out_SRR13736659/tmp/valid_idx_pcrfree.pairs.2000.cool"
name_bank2= "SRR13736659_GSM5090900_Q"

#cool_file2= "/home/axel/Bureau/COOL_files/all_cool_files_koshland/res_2000/cool_files/out3_SRR11893114_tmp_valid_idx_pcrfree.pairs.2000.cool"
#name_bank2= "SRR11893114"

c1 = cooler.Cooler(cool_file1)
cooler.balance_cooler(c1, store=True)   # Normalisation 
d=c1.info
total_reads = d['sum']

c2 = cooler.Cooler(cool_file2)
cooler.balance_cooler(c2, store=True)   # Normalisation 
d=c2.info
total_reads = d['sum']

BIN=2000
bin_histo=2000
list_chr =  ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
             "chr16"]

centro = pd.read_table("/home/axel/Bureau/YEAST/centro1.dat4",
                       header=None, delimiter=" ")
# Create a zip object from two lists
zipbObj = zip(centro[0], centro[1])
# Create a dictionary from zip object
centro_dist = dict(zipbObj)
centro_dist['chrM'] = 0
centro_dist['plasmid_p2-micron'] = 0

sizes = pd.read_table("/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/sacCer3.chr_sizes.txt",
                      header=None, delimiter="\t") 
zipbObj = zip(sizes[0], sizes[1])
sizes_dist = dict(zipbObj)

#
chr2 = "plasmid_p2-micron"   # plasmid micron
chrM = "chrM"   # mito 

sizes_dist[chrM] = 100000
sizes_dist['plasmid_p2-micron']=6300
sizes_dist['pRS413'] = 4970  # the one of Fabien 
sizes_dist['pSH47'] = 6979
#
genes = pd.read_csv('/home/axel/Bureau/YEAST/GENES_SC288/orf_genomic_all.txt2', 
                    sep=" ", header= None)
long_genes = genes[ abs(genes[3]-genes[2]) > 7000 ]
len(long_genes)

# Genomic data
ori = pd.read_csv('/home/axel/Bureau/YEAST/origins/origins_Alvino.txt2', 
                    sep=" ", header= None)
ars = pd.read_csv('/home/axel/Bureau/YEAST/ARS/ars.txt2', 
                    sep=" ", header= None)
acetyl = pd.read_csv('/home/axel/Bureau/YEAST/H4_acetylation/GSM3445778_YSB787-1-H4ac_peaks.narrowPeak3', 
                    sep=" ", header= None)
motif = pd.read_csv('/home/axel/Bureau/2micron_plasmid_PROJECT/TGCATTTTT.sam2', 
                    sep=" ", header= None)
g4 = pd.read_csv('/home/axel/Bureau/YEAST/G4/pcbi.1000861.s001.txt2.lifted.renamed2', 
                    sep="\t", header= None)

# ChIPseq data:
# 1rst condition 

chip_ip1_a = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175393.fastq.sam.MQ30.Rpb3_in_Log_All_Replicates_Merged_IP', 
                    sep=" ", header= None)
chip_input1_a  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175394.fastq.sam.MQ30.Rpb3_in_Log_All_Replicates_Merged_Input', 
                    sep=" ", header= None)

chip_ip2_a  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175373_Brn1_in_Log_All_Replicates_Merged_IP.fastq.sam.MQ30', 
                    sep=" ", header= None)
chip_input2_a  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175374.fastq.sam.MQ30.Brn1_in_Log_All_Replicates_Merged_Input', 
                    sep=" ", header= None)

chip_ip3_a  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR2065097.fastq.sam.MQ30.FigS5_Scc1PK9_IP_G1_releasing_60min', 
                    sep=" ", header= None)     # cohesin
chip_input3_a  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR2065092.fastq.sam.MQ30.FigS5_Scc1PK9_WCE_G1_releasing_60min', 
                    sep=" ", header= None)

chip_ip4_a  = pd.read_csv('/home/axel/Bureau/smc6_chip/SRR1555047_IP.fastq.sam.MQ0', 
                    sep=" ", header= None)     # smc6
chip_input4_a  = pd.read_csv('/home/axel/Bureau/smc6_chip/SRR1555048_input.fastq.sam.MQ0', 
                    sep=" ", header= None)

# 2d condition 
chip_ip1_b = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175399.fastq.sam.MQ30.Rpb3_in_Q_All_Replicates_Merged_IP', 
                    sep=" ", header= None)
chip_input1_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175382.fastq.sam.MQ30.Brn1_in_Q_All_Replicates_Merged_Input', 
                    sep=" ", header= None)

chip_ip2_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175381.fastq.sam.MQ30.Brn1_in_Q_All_Replicates_Merged_IP', 
                    sep=" ", header= None)
chip_input2_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR7175382.fastq.sam.MQ30.Brn1_in_Q_All_Replicates_Merged_Input', 
                    sep=" ", header= None)

chip_ip3_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR11466784.1_1.fastq.sam.MQ30', 
                    sep=" ", header= None)     # cohesin in Q phase .. 
chip_input3_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR13866613_1_18504_NoTag_BY4741_ChIP-exo.fastq.sam.MQ30', 
                    sep=" ", header= None)

chip_ip4_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR11466784.1_1.fastq.sam.MQ30', 
                    sep=" ", header= None)     # H3
chip_input4_b  = pd.read_csv('/media/axel/RSG4/diverse_yeast_data/CHip_seq_2018/SRR13866613_1_18504_NoTag_BY4741_ChIP-exo.fastq.sam.MQ30', 
                    sep=" ", header= None)


# Multi plot: 
name_prot="Rpb3"
name_prot2="Condensin"
name_prot3="Cohesin"
name_prot4="Smc6"
#
limit_y =  8*10**-7

if not os.path.exists(name_bank1+"_files") :
    os.makedirs(name_bank1+"_files")

list_all_contact = []
list_all_chip = []

list_all_contact_lg = []
list_all_chip_lg = []
i_fig=0

list_chr = ["chr15"]
for chr1 in list_chr :
    i_fig=i_fig+1
    genes_chr = genes[genes[0] ==chr1]
    long_genes_chr = long_genes[long_genes[0] ==chr1]
    ori_chr = ori[ori[0] ==chr1]
    ars_chr = ars[ars[0] ==chr1]
    acetyl_chr = acetyl[acetyl[0] ==chr1]
    g4_chr = g4[g4[0] ==chr1]
    
    chip_ip_chr1_a = chip_ip1_a[chip_ip1_a[0] ==chr1]
    chip_input_chr1_a = chip_input1_a[chip_input1_a[0] ==chr1]

    chip_ip_chr2_a = chip_ip2_a[chip_ip2_a[0] ==chr1]
    chip_input_chr2_a = chip_input2_a[chip_input2_a[0] ==chr1]
       
    chip_ip_chr3_a = chip_ip3_a[chip_ip3_a[0] ==chr1]
    chip_input_chr3_a = chip_input3_a[chip_input3_a[0] ==chr1]
    
    chip_ip_chr4_a = chip_ip4_a[chip_ip4_a[0] ==chr1]
    chip_input_chr4_a = chip_input4_a[chip_input4_a[0] ==chr1]
    
    chip_ip_chr1_b = chip_ip1_b[chip_ip1_b[0] ==chr1]
    chip_input_chr1_b = chip_input1_b[chip_input1_b[0] ==chr1]

    chip_ip_chr2_b = chip_ip2_b[chip_ip2_b[0] ==chr1]
    chip_input_chr2_b = chip_input2_b[chip_input2_b[0] ==chr1]
       
    chip_ip_chr3_b = chip_ip3_b[chip_ip3_b[0] ==chr1]
    chip_input_chr3_b = chip_input3_b[chip_input3_b[0] ==chr1]
    
    chip_ip_chr4_b = chip_ip4_b[chip_ip4_b[0] ==chr1]
    chip_input_chr4_b = chip_input4_b[chip_input4_b[0] ==chr1]
    
    # matrice 
    matscn1 = c1.matrix(balance=True).fetch(chr1, chr1) 
    matscn1[np.isnan(matscn1)] = 0 
    coverage = matscn1.sum(axis=0)

    mat = c1.matrix(balance=False).fetch(chr1, chr2)
    mat[np.isnan(mat)] = 0
    coverage_plasmid = mat.sum(axis=1)
#    coverage_plasmid[coverage_plasmid==0] = np.nan
    
    mat2 = c1.matrix(balance=True).fetch(chr1, chrM)
    mat2[np.isnan(mat2)] = 0   
    coverage_chrM = mat2.sum(axis=1)
#    coverage_chrM[ coverage_chrM==0 ] = np.nan
    
    m1 = c1.matrix(balance=False).fetch(chr1, chr1)
    reads_chr1= m1.sum()
    m2 = c1.matrix(balance=False).fetch(chr2, chr2)
    reads_chr2= m2.sum()
    m12 = c1.matrix(balance=False).fetch(chr1, chr2)
    reads_chr12= m12.sum()
    
    matscn2 = c2.matrix(balance=True).fetch(chr1, chr1) 
    matscn2[np.isnan(matscn2)] = 0 
    coverage = matscn2.sum(axis=0)

    mat = c2.matrix(balance=False).fetch(chr1, chr2)
    mat[np.isnan(mat)] = 0
    coverage_plasmid = mat.sum(axis=1)
#    coverage_plasmid[coverage_plasmid==0] = np.nan
    
    mat2 = c1.matrix(balance=False).fetch(chr1, chrM)
    mat2[np.isnan(mat2)] = 0   
    coverage_chrM = mat2.sum(axis=1)
#    coverage_chrM[ coverage_chrM==0 ] = np.nan
    
    m1 = c2.matrix(balance=False).fetch(chr1, chr1)
    reads_chr1= m1.sum()
    m2 = c2.matrix(balance=False).fetch(chr2, chr2)
    reads_chr2= m2.sum()
    m12 = c2.matrix(balance=False).fetch(chr1, chr2)
    reads_chr12= m12.sum()
    
    # computation of correlation:
    maxi= matscn2.shape[0]*bin_histo
    
    values_chip1_a = function_geno.chip_compute(chip_ip_chr1_a, chip_input_chr1_a, bin_histo, maxi)
    values_chip2_a = function_geno.chip_compute(chip_ip_chr2_a, chip_input_chr2_a, bin_histo, maxi)
    values_chip3_a = function_geno.chip_compute(chip_ip_chr3_a, chip_input_chr3_a, bin_histo, maxi)
    values_chip4_a = function_geno.chip_compute(chip_ip_chr4_a, chip_input_chr4_a, bin_histo, maxi)
    
    values_chip1_b = function_geno.chip_compute(chip_ip_chr1_b, chip_input_chr1_b, bin_histo, maxi)
    values_chip2_b = function_geno.chip_compute(chip_ip_chr2_b, chip_input_chr2_b, bin_histo, maxi)
    values_chip3_b = function_geno.chip_compute(chip_ip_chr3_b, chip_input_chr3_b, bin_histo, maxi)
    values_chip4_b = function_geno.chip_compute(chip_ip_chr4_b, chip_input_chr4_b, bin_histo, maxi)

    corr = spearmanr(coverage_plasmid, values_chip1_a,
                     nan_policy="omit")
    corr_coeff = corr.correlation
    #------------------------------------------------------------------------- 
    # MULTIPLOT:  
    plt.figure()
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8, 11.8)
    fig.tight_layout(pad=3.0)
    gs = gridspec.GridSpec(5, 2,height_ratios=[9,1,1,1,1])

    ax1 = plt.subplot(gs[0])
    ax1.imshow(matscn1**0.15,interpolation="none", cmap="afmhot_r", 
             aspect="auto")
    plt.title(chr1+" "+str(total_reads)+" "
             +str(reads_chr1)+" "+str(reads_chr2)+" "+str(reads_chr12)+"\n"+
             name_bank1)

    ax11 = plt.subplot(gs[1], sharex=ax1, sharey=ax1)
    ax11.imshow(matscn2**0.15,interpolation="none", cmap="afmhot_r", 
             aspect="auto")
    plt.title(chr1+" "+str(total_reads)+" "
             +str(reads_chr1)+" "+str(reads_chr2)+" "+str(reads_chr12)+"\n"+
             name_bank2)
    
    ax2 = plt.subplot(gs[2], sharex=ax1)
    ax2.plot(values_chip1_a, color="tomato",linewidth=3.0)
    ax2.fill_between(range(len(values_chip1_a)),values_chip1_a, facecolor='tomato')
    plt.title("PolII ChIPseq "+name_prot)
    plt.xticks([], [])
    plt.ylim(0,8)
    
    coverage_plasmid = coverage_plasmid /(sizes_dist[chr2])
    median_plasmid = np.nanmedian(coverage_plasmid) 
    print(median_plasmid)
    
    ax3 = plt.subplot(gs[4], sharex=ax1)
    ax3.plot(values_chip2_a, color="royalblue")
    ax3.fill_between(range(len(values_chip2_a)),values_chip2_a, facecolor='royalblue')
    plt.plot(centro_dist[chr1]/BIN,0.0,'o', color="orange",label="Centro")
    plt.ylim(0,2.5)
    plt.title("ChIP-seq "+name_prot2)
    plt.xticks([], [])
    
    coverage_chrM = coverage_chrM/(sizes_dist[chrM]) 
    median_chrM = np.nanmedian(coverage_chrM) 
    
    ax4 = plt.subplot(gs[6], sharex=ax1)
    ax4.plot(values_chip3_a, color="orange" )
    ax4.fill_between(range(len(values_chip3_a)),values_chip3_a, facecolor='orange')
    plt.title("Cohesin occupancy")
    plt.ylim(0,3)
    list_all_contact=np.concatenate((list_all_contact, coverage_plasmid), axis=0)
    list_all_chip=np.concatenate((list_all_chip, values_chip1_a), axis=0)
    
    ax4 = plt.subplot(gs[8], sharex=ax1)
    ax4.plot(values_chip4_a, color="green" )
    ax4.fill_between(range(len(values_chip4_a)),values_chip4_a, facecolor='green')
    plt.title("Smc5 occupancy")
    plt.ylim(0,3)

    # only for long genes: 
    l1=long_genes_chr[3]/BIN
    l1=list(l1)
    l1 = [int(x) for x in l1]
    list_all_contact_lg=np.concatenate((list_all_contact_lg, 
                                        coverage_plasmid[l1]), axis=0)
    list_all_chip_lg=np.concatenate((list_all_chip_lg, 
                                     values_chip1_a[l1]), axis=0)
    

    ax2 = plt.subplot(gs[3], sharex=ax1)
    ax2.plot(values_chip1_b, color="tomato",linewidth=3.0)
    ax2.fill_between(range(len(values_chip1_b)),values_chip1_b, facecolor='tomato')
    plt.title("PolII ChIPseq "+name_prot)
    plt.xticks([], [])
    plt.ylim(0,5.0)
    
    coverage_plasmid = coverage_plasmid /(sizes_dist[chr2])
    median_plasmid = np.nanmedian(coverage_plasmid) 
    print(median_plasmid)
    
    ax3 = plt.subplot(gs[5], sharex=ax1)
    ax3.plot(values_chip2_b, color="royalblue" )
    ax3.fill_between(range(len(values_chip2_b)),values_chip2_b, facecolor='royalblue')
    plt.plot(centro_dist[chr1]/BIN,0.0,'o', color="orange",label="Centro")

    plt.ylim(0,2.5)
    plt.title("ChIP-exo "+name_prot2)
    plt.xticks([], [])
  
    coverage_chrM = coverage_chrM/(sizes_dist[chrM]) 
    median_chrM = np.nanmedian(coverage_chrM) 
    
    ax4 = plt.subplot(gs[7], sharex=ax1)
#    ax4.plot(coverage_plasmid/median_plasmid, color="orange", linewidth=3.0)
    ax4.plot(values_chip3_b, color="orange" )
    ax4.fill_between(range(len(values_chip3_b)),values_chip3_b, facecolor='orange')
#    plt.title(str(chrM)+" Median="+str(round(median_chrM*10**8,2))+"x 10^8")
    plt.title("Cohesin occupancy")
    plt.ylim(0,3)
    
    ax4 = plt.subplot(gs[9], sharex=ax1)
#    ax4.plot(coverage_plasmid/median_plasmid, color="orange", linewidth=3.0)
    ax4.plot(values_chip3_b, color="pink" )
    ax4.fill_between(range(len(values_chip3_b)),values_chip3_b, facecolor='pink')
#    plt.title(str(chrM)+" Median="+str(round(median_chrM*10**8,2))+"x 10^8")
    plt.title("Cohesin occupancy")
    plt.ylim(0,3)
    
    
    list_all_contact=np.concatenate((list_all_contact, coverage_plasmid), axis=0)
    list_all_chip=np.concatenate((list_all_chip, values_chip1_a), axis=0)
    
    plt.savefig(name_bank1+"_files"+"/MAT_SCN_"+chr1+"_"+
            name_bank1+"_"+str(BIN/1000)+"bp"+".pdf")
    

plt.close("all")
# Correlation plot 
corr = spearmanr(list_all_contact, list_all_chip,
                     nan_policy="omit")
corr_coeff = corr.correlation  
 
plt.plot(list_all_chip, list_all_contact, 'o', 
     color="royalblue", alpha=0.2) 

plt.plot(list_all_chip_lg, list_all_contact_lg, 'o', 
     color="red", alpha=0.8) 

plt.title("Spearman= "+
              str(round(corr.correlation,2))+
              " pvalue="+str(corr.pvalue))

plt.xlabel("ChIP / Input "+name_prot)
plt.ylabel("Contact of the plasmid with host genome")    
    
name_set="HSC_log_2000"
plt.savefig("1D_ENRICHMENT_HSC_"+name_prot+"_"+
            name_set+".pdf",  dpi=600, format='pdf')
plt.close("all")   
    
np.savetxt("1D_ENRICHMENT_HSC_"+name_prot+"_"+
            name_set+".txt",X=s)    
    
    
