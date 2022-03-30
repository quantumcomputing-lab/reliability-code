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
import seaborn as sns

A4ratio = 297/210

filename = 'washington_readout_error_mar10_2022.csv'
foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'

#### POINT IN TIME
df = pd.read_csv(foldername +filename)
# change below to select a specific date
# df = df.sample(frac=1,axis=1).reset_index(drop=True)

columns = []
for i in range(127):
    columns = columns+['q'+str(i)+'_readout_err']

df0 = df[columns]
df1=df0.loc[3:,:].reset_index(drop=True)

sns_df = pd.DataFrame()
sns_df["qubit"]=np.arange(0,127)
sns_df["FI"]=np.array(df1.loc[df1.shape[0]-1,:])*100 #np.nanmean(df1,axis=0)

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
ax.set_ylim(-0.01,40)
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
x  = np.arange(0,df1.shape[0])

df1_min = np.floor(np.nanmin(df1)*100)
df1_max = np.ceil(np.nanmax(df1)*100)

q=0
for k in range(8):
    print(k)
    fig, ax = plt.subplots(nrows=16, sharex=True, subplot_kw=dict(frameon=False)) # frameon=False removes frames
    fig.set_size_inches(21, 16)
    plt.subplots_adjust(hspace=.08)

    for i in range(16):

        ax[i].xaxis.label.set_visible(False)
        ax[i].get_yaxis().set_ticklabels([])
        ax[i].get_xaxis().set_ticklabels([])
        ax[i].grid(True, color='lightgrey', alpha=0.8)

        if q<127:
            ax[i].set_ylabel('q'+str(q))    
            ax[i].plot(x, df1['q'+str(q)+'_readout_err'], color='k')

        if i==15 or q==126:
            ax[i].xaxis.label.set_visible(True)
            ax[i].set_xticks([0,df1.shape[0]])
            ax[i].set_xticklabels(['Dec-21','Mar-22'])

        q=q+1

    fig.savefig('FI_raw_time_'+str(k)+'.png', bbox_inches = "tight", dpi=300)