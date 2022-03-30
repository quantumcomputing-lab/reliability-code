#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 17:08:23 2022
@author: samudra
"""

import matplotlib.pyplot as plt
import calendar
import pandas as pd
import numpy as np
import scipy
import scipy.stats

A4ratio = 297/210

#### POINT IN TIME
df = pd.read_csv('/home/samudra/Desktop/sjv_stability_feb_2022/data/data_wash_dummy.csv')
# change below to select a specific date
df = df.sample(frac=1,axis=1).reset_index(drop=True)

sns_df = pd.DataFrame()
sns_df["qubit"]=np.arange(0,127)
sns_df["FI"]=np.nanmean(df,axis=0)

# delta_f
delta_f_mat = np.zeros((127,127))

for i in range(127):
    delta_f_mat[i,i]=np.nan
    for j in np.arange(i+1, 127):
        delta_f_mat[j,i]=np.abs( sns_df["FI"][j]-sns_df["FI"][i] )
        delta_f_mat[i,j] = np.nan

delta_f_mat[delta_f_mat==0]=np.nan
        
row_col_min = np.argwhere( delta_f_mat == np.nanmin(delta_f_mat))
row_col_min = row_col_min[0]
row_col_max = np.argwhere( delta_f_mat == np.nanmax(delta_f_mat))
row_col_max = row_col_max[0]
y = ( delta_f_mat[row_col_min[0], row_col_min[1]], delta_f_mat[row_col_max[0], row_col_max[1]])

sns.set()
sns_plot = sns.heatmap(delta_f_mat, robust=True, fmt="f", cmap='Oranges',
                       cbar_kws={'label': 'Difference in Fidelity (%)'})
sns_plot.figure.set_size_inches(16, int(16*A4ratio))
sns_plot.set_facecolor('white')
labels = [i for i in range(0,127,6)]
sns_plot.axes.set_xticks(labels)
#sns_plot.axes.set_yticks(labels)
#sns_plot.axes.set_xticklabels(labels, rotation=0, position=(0,2.8), fontsize="10", va="center")
sns_plot.axes.set_xticklabels(labels, rotation=0, fontsize="10", va="center")
#sns_plot.axes.set_yticklabels(labels, rotation=0, position=(0,0.28), fontsize="10", va="center")
sns_plot.axes.set_ylabel('Qubit #', fontsize=12)
sns_plot.axes.set_xlabel('Qubit #', fontsize=12)
figure = sns_plot.get_figure()

left, bottom, width, height = [0.54, 0.6, 0.2, 0.24]
ax = figure.add_axes([left, bottom, width, height])
x_pos = (1,2)
ax.bar(x_pos, y, align='center', alpha = 0.8, color='darkorange',
       error_kw=dict(ecolor='k', lw=2, capsize=5, capthick=2))
ax.set_xticks(x_pos)
ax.set_ylim(-0.01,30)
min_pair = 'Min\n(' + str(row_col_min[0]) + ', ' + str(row_col_min[1]) +')'
max_pair = 'Max\n(' + str(row_col_max[0]) + ', ' + str(row_col_max[1]) +')'
ax.set_xticklabels((min_pair, max_pair), fontsize=10)
ax.set_xlabel('Qubit Pair',fontsize=12)
ax.set_ylabel('Difference in Fidelity (%)',fontsize=12)
ax.yaxis.grid(True)
ax.set_facecolor('lightgray')
ax.set_frame_on(True)
figure.savefig('FI_raw_lattice.png', bbox_inches = "tight", dpi=300)

#### ACROSS TIME
# Style 1: time-series
sns_df = pd.read_csv('/home/samudra/Desktop/sjv_stability_feb_2022/data/data_wash_dummy.csv')
x  = np.arange(0,df.shape[0])
q=0
for k in range(8):
    print(k)
    fig, ax = plt.subplots(nrows=16, sharex=True, subplot_kw=dict(frameon=False)) # frameon=False removes frames
    fig.set_size_inches(21, 16)
    plt.subplots_adjust(hspace=.08)

    for i in range(16):
        if q==127:
            ax[i].get_yaxis().set_ticks([])
            break

        ax[i].xaxis.label.set_visible(False)
        ax[i].get_yaxis().set_ticks([])
        ax[i].set_ylabel('q'+str(q))
        ax[i].plot(x, sns_df[str(q)], color='k')

        if i==15 or q==126:
            ax[i].xaxis.label.set_visible(True)

        q=q+1

    fig.savefig('FI_raw_time_'+str(k)+'.png', bbox_inches = "tight", dpi=300)

# Style 2: scatter

# Prepare data format
df = pd.read_csv('/home/samudra/Desktop/sjv_stability_feb_2022/data/data_wash_dummy.csv')

x =[]
y = []
sns_df = pd.DataFrame()
for q  in range(df.shape[1]):
    x = np.append(x, np.zeros(df.shape[0])+q)
    y = np.append(y, df[str(q)])
sns_df["qubit"] = x.astype('int')
sns_df["FI"] = y

sns.set_style(style="whitegrid")
#tips = sns.load_dataset("tips")
#ax = sns.stripplot(x=sns_df["qubit"])
#ax = sns.stripplot(x="qubit", y="FI", data=sns_df)
#ax = sns.stripplot(x="qubit", y="FI", data=sns_df, linewidth=1)

for k in range(8):
    print(k)
    q0=k*16
    if k==7:
        q1=126
    else:
        q1=k*16+15

    row_start = q0*df.shape[0]
    row_end = (q1+1)*df.shape[0]-1
    
    dataset = sns_df.loc[row_start:row_end,:]    
    ax = sns.stripplot(x="qubit", y="FI", data=dataset, jitter=0.01, linewidth=0.8, alpha=0.8)
    ax.figure.set_size_inches(21, 16)
    ax.set_xlabel("Qubit #", fontsize=12)
    ax.set_ylabel("Initialization Fidelity (%)", fontsize=12)
    ax.grid(True, which='both', alpha=0.4)
    fig=ax.get_figure()
    fig.savefig('FI_raw_scatter_'+str(k)+'.png', bbox_inches = "tight", dpi=300)
