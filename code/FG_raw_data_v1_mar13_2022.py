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
import math

A4ratio = 297/210

filename = 'washington_rawdata_mar12_2022.csv'
foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'
nqubits = 127

#### POINT IN TIME
df = pd.read_csv(foldername +filename)
# change below to select a specific date
# df = df.sample(frac=1,axis=1).reset_index(drop=True)

header_list = np.array(df.columns)
columns = ['query_date', 'last_update_date']
no_connection = []
for i in range(nqubits):
    for j in np.arange(i+1, nqubits):
        colname1 = 'cx'+str(i)+'_'+str(j)+'_gate_error_value'
        colname2 = 'cx'+str(j)+'_'+str(i)+'_gate_error_value'
        if colname1 in header_list:
            columns = columns+[colname1]
        elif  colname2 in header_list:
            columns = columns+[colname2]
        else:
            no_connection.append(colname1)
df=df[columns]
last_df = df.iloc[-1][2:]
keys = np.array( last_df.index.astype('str')  )
vals = np.array( last_df.values.astype('float') )

N_rows = [not math.isnan(item) for item in vals ]
N_rows = [a and b for a, b in zip( N_rows  , vals !=0.0) ]
N_rows = [a and b for a, b in zip( N_rows  , vals !=1.0) ]
last_df = last_df[N_rows]
npairs = len(last_df)

sns_df = pd.DataFrame()
sns_df["qubit_pair"]=np.array(last_df.index)
sns_df["FG"]=np.array(last_df.values)*100 #np.nanmean(df1,axis=0)

# delta_f
delta_f_mat = np.zeros((npairs,npairs))
for i in range(npairs):
    delta_f_mat[i,i]=np.nan
    for j in np.arange(i+1, npairs):
        delta_f_mat[j,i]=np.abs( sns_df["FG"][j]-sns_df["FG"][i] )
        #delta_f_mat[i,j] = delta_f_mat[j,i]
        delta_f_mat[i,j] = np.nan
       
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
ticks = [i for i in range(0, npairs, 6)]
labels = [item.replace('_gate_error_value', '') for item in np.array(sns_df["qubit_pair"][ticks])]
sns_plot.axes.set_xticks(ticks)
sns_plot.axes.set_xticklabels(labels, rotation=90, fontsize="12", va="center", position=(0,-.028))

sns_plot.axes.set_yticks(ticks)
sns_plot.axes.set_yticklabels(labels, rotation=0, fontsize="12", va="center")

#sns_plot.axes.set_yticks(labels)
#sns_plot.axes.set_xticklabels(labels, rotation=0, position=(0,2.8), fontsize="10", va="center")
#sns_plot.axes.set_yticklabels(labels, rotation=0, position=(0,0.28), fontsize="10", va="center")
sns_plot.axes.set_ylabel('Pair #', fontsize=12)
sns_plot.axes.set_xlabel('Pair #', fontsize=12)
figure = sns_plot.get_figure()
figure.savefig('FG_raw_lattice.png', bbox_inches = "tight", dpi=300)

# comments: only a subest of the 144 pairs is labelled. 
# Gates which had nan for their calibration are not shown. 
# There were 9 such gates.
# The data for 135 is shown. 
# The ticks are obviously not ordered as order has no meaning.
# Only 23 are shown
# Only lower triangular shown

#### ACROSS TIME
df1 = df.iloc[:,2:]
x  = np.arange(0,df1.shape[0])
#vals = df.iloc[:,2:].values.tolist()
#v=[item for item in vals if (item!=0 and item!=1.0) ]
#df1_min = np.floor( np.nanmin(v)*100 )
#df1_max = np.ceil( np.nanmax(v)*100 )

q=0
for k in np.arange(0,8):
    print(k)
    fig, ax = plt.subplots(nrows=18, sharex=True, subplot_kw=dict(frameon=False)) # frameon=False removes frames
    fig.set_size_inches(21, 16)
    plt.subplots_adjust(hspace=.08)

    for i in range(18):

        ax[i].xaxis.label.set_visible(False)
        ax[i].get_yaxis().set_ticklabels([])
        ax[i].get_xaxis().set_ticklabels([])
        ax[i].grid(True, color='lightgrey', alpha=0.8)

        if q<144:
            y = np.array( df1.iloc[:,q] )
            y[y==0.0]=np.nan
            y[y==1.0]=np.nan
            ax[i].set_ylabel(df1.columns[q], rotation=0, fontsize="12",)
            ax[i].yaxis.set_label_coords(-0.07,0)
            ax[i].plot(x, y, color='k')

        if i==18-1 or q==144-1:
            ax[i].xaxis.label.set_visible(True)
            ax[i].set_xticks([0,df1.shape[0]])
            ax[i].set_xticklabels(['Dec-21','Mar-22'])

        q=q+1

    fig.savefig('FG_raw_time_'+str(k)+'.png', bbox_inches = "tight", dpi=150)

'''
notice that there are gaps in the plot. 
these are times when the calibration data came back erroneous as denoted by 
nan/ 0/ 1
'''