#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 23, 2022
@author: samudra
Purpose: creates historical database for any online IBM device

TODOS:
data_update():
    - from device_data.csv, select the ones which are online and do the following for each
    - check which dates exist (call it last_available_date). 
      collect the static data from online_since date or last_available_date to today
    - collect the latest eta data; 
      it is fine if eta is not latest and greatest; 
      but trigger run for eta

#this needs to be combined with plot_heatmap:
    plot_nmi_heatmap(nmi_df, nqubits, colorbar_label, xylabel, figname )

WRITEUP:
- we have dropped dates pre-dec 9 because that's when it was commissioned
- we have dropped dates nan, 1, 0 as they are cleary erroneous entries
- for all entries which have multiple entries on the same date, we have taken the latest hourly value
translate all timestamps to date only
- FG TS: 
    notice that there are gaps in the plot. 
    these are times when the calibration data came back erroneous as denoted by nan/ 0/ 1

- for FG, in the spatial graph:
    only a subest of the 144 pairs is labelled. 
    gates which had nan for their calibration are not shown (there were 9 such gates).
    data for 135 is shown. 
    ticks are obviously not ordered as order has no meaning.
    only 23 are shown
    only lower triangular shown
"""

import matplotlib.pyplot as plt
import calendar
import pandas as pd
import numpy as np
import scipy
import scipy.stats
import seaborn as sns
import sys
import math

device_list = {
'washington': 127,
'brooklyn': 65,
'toronto': 27,
'mumbai': 27,
'cairo': 27,
'montreal': 27,
'auckland': 27,
'hanoi': 27,
'guadalupe': 16,
'perth': 7,
'lagos': 7,
'jakarta': 7,
'manila': 5,
'bogota': 5,
'santiago': 5,
'quito': 5,
'lima': 5,
}

devicename = 'brooklyn'
nqubits = device_list[devicename]
A4ratio = 297/210
filename = 'device_data.csv'
foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'+devicename+'/'
df_all = pd.read_csv(foldername +filename)
df_all = df_all.loc[:, ~df_all.columns.str.contains('^Unnamed')]
tolerance = np.float64(0.01)

# Initialization Fidelity
df = FI_df(df_all)
# CI ACROSS TIME
F = time_array_FI_FG(df)
plot_CI(F, xlabel="Initialization fidelity", 
        ylabel="Normalized histogram", figname = 'FI_across-time.png')
# CI POINT-IN-TIME
F = latest_time_array_FI_FG(df)
plot_CI(F, xlabel="Inter-qubit difference in initialization fidelity", 
        ylabel="Normalized histogram", figname = 'FI_point-in-time.png')
# RAW DATA POINT IN TIME
last_df = latest_time_array_FI_FG(df)
ticks = [i for i in range(0, nqubits, 6)]
labels = [i for i in range(0,nqubits,6)]
plot_heatmap(last_df, ticks, labels, xylabel = 'Qubit #',
                 figname = 'FI_raw_lattice.png',
                 colorbar_label = dict({'label': 'Difference in Fidelity (%)'}))
# RAW DATA ACROSS TIME
df1 = df.iloc[:,2:]
plot_raw_ts(df1, rows_per_sublpot=int(df1.shape[1]/8), total_ts=df1.shape[1], figname = 'FI_raw_time_')

# Gate Fidelity
df = FG_df(df_all)
# CI ACROSS TIME
F = time_array_FI_FG(df)
plot_CI(F, xlabel="CNOT fidelity", 
        ylabel="Normalized histogram", figname = 'FG_across-time.png')
# CI POINT-IN-TIME
F = latest_time_array_FI_FG(df)
plot_CI(F, xlabel="Inter-gate difference in CNOT fidelity", 
        ylabel="Normalized histogram", figname = 'FG_point-in-time.png')
# RAW DATA POINT IN TIME
last_df = latest_time_array_FI_FG(df)
npairs = len(last_df)
ticks = [i for i in range(0, npairs, 6)]
labels = [item.replace('_gate_error_value', '') for item in np.array( last_df.index )][::6]
plot_heatmap(last_df, ticks, labels, xylabel = 'Connection pair #',
                 figname = 'FG_raw_lattice.png',
                 colorbar_label = dict({'label': 'Difference in Fidelity (%)'}))
# RAW DATA ACROSS TIME
df1 = df.iloc[:,2:]
plot_raw_ts(df1, rows_per_sublpot=int(df1.shape[1]/8), total_ts=df1.shape[1], figname = 'FG_raw_time_')

# Duty Cycle
df = tau_df(df_all, nqubits)
# CI ACROSS TIME
Flist=list()
for col in df.columns[2:]:
    Flist = Flist + [item for item in df[col] if not math.isnan(item)]
F = np.array(Flist)
plot_CI(F, xlabel="CNOT duty cycle", 
        ylabel="Normalized histogram", figname = 'tau_across-time.png')
# CI POINT-IN-TIME
F = [item for item in df.iloc[-1,2:] if not math.isnan(item)]
plot_CI(F, xlabel="Inter-gate difference in duty cycle", 
        ylabel="Normalized histogram", figname = 'tau_point-in-time.png')
# RAW DATA POINT IN TIME
last_df = pd.DataFrame(df.iloc[-1,2:])
Nrows = [not math.isnan(item) for item in last_df.values]
last_df = last_df[Nrows]
last_df = last_df.sort_values(by=last_df.columns[0])
last_df = last_df/1000
npairs = len(last_df)
ticks = [i for i in range(0, npairs, int(npairs/5))]
labels = [item.replace('_tau_value', '') for item in np.array( last_df.index )][::35]
xlabel = 'Gate Pair #'
ylabel = 'Point-in-time variation\n of CNOT duty cycle\n (x1000)'
figname = 'tau_raw_lattice.png'
plot_sorted_ts(last_df, ticks, labels, xlabel, ylabel, figname)

# RAW DATA ACROSS TIME
df1 = df.iloc[:,2:]
plot_raw_ts(df1, rows_per_sublpot=int(df1.shape[1]/8), total_ts=df1.shape[1], figname = 'tau_raw_time_')

# Addressability
csvfile = 'nmi_data.csv'
nmi_df  = get_nmi_df(foldername+csvfile, devicename, data_refresh=0) # if you want error bars, we will need to run multiple times
plot_nmi_heatmap(nmi_df, 
                 nqubits, 
                 colorbar_label = 'Mutual Information\n (normalized, %)', 
                 xylabel = "Qubit #", 
                 figname = 'nmi.png') #this needs to be combined with plot_heatmap

#### GENERAL ROUTINES
def remove_nan_1_0(series_s):
    N_rows  = [not math.isnan(item) for item in series_s]
    N_rows = [a and b for a, b in zip( N_rows  , series_s!=0.0)]
    N_rows = [a and b for a, b in zip(N_rows, series_s!=1.0)]
    series_s = series_s[N_rows]
    return series_s

def tau_df(df, nqubits):
    #filename = 'washington_rawdata_mar12_2022.csv'
    #foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'
    #df = pd.read_csv(foldername +filename)
    #nqubits = 127
    #df.columns[df.columns.str.contains('cx0_1')]
    #df.columns[df.columns.str.contains('T2')]


    ind = np.logical_and(df.columns.str.contains('_gate_length'), df.columns.str.contains('cx'))
    ind = np.logical_or(ind, df.columns.str.contains('T2'))
    cols = ['query_date', 'last_update_date']
    cols.extend( np.array(df.columns[ind]) )
    df1 = df.loc[:, cols]

    tau_colnames=['query_date', 'last_update_date']
    for q1 in range(nqubits):
        for q2 in range(nqubits):
            q1str = str(q1)
            q2str = str(q2)

            gate_length_colname = 'cx' + q1str + '_' + q2str + '_gate_length_value'
            q1_colname = 'q_' + q1str + '__T2_value'
            q2_colname = 'q_' + q2str + '__T2_value'

            gate_length_unitname = 'cx' + q1str + '_' + q2str + '_gate_length_unit'
            q1_unitname = 'q_' + q1str + '__T2_unit'
            q2_unitname = 'q_' + q2str + '__T2_unit'


            if ((gate_length_colname in df1.columns) and 
                (q1_colname in df1.columns) and
                (q2_colname in df1.columns)):

                df1[gate_length_unitname]=df1[gate_length_unitname].astype('str')
                df1[q1_unitname]=df1[q1_unitname].astype('str')
                df1[q2_unitname]=df1[q2_unitname].astype('str')

                tau_colname = 'cx_' + q1str + '_' + q2str + '_tau_value'
                tau_colnames.append(tau_colname)

                if (set(df1[gate_length_unitname].dropna().unique())!= set(['ns']) or 
                    set(df1[q1_unitname].dropna().unique())!= set(['us']) or
                    set(df1[q2_unitname].dropna().unique())!= set(['us'])):
                    L=df1.shape[0]
                    tau = []
                    for i in range(L):
                        if df1.loc[i,gate_length_unitname]=='us':
                            gate_length = df1.loc[i, gate_length_colname]*1000
                        elif df1.loc[i,gate_length_unitname]=='ns':
                            gate_length = df1.loc[i, gate_length_colname]                            
                        else:
                            gate_length = np.nan
                            
                        if df1.loc[i,q1_unitname]=='ns':
                            q1_T2=df1.loc[i, q1_colname]/1000
                        elif df1.loc[i,q1_unitname]=='us':
                            q1_T2=df1.loc[i, q1_colname]                            
                        else:
                            q1_T2 = np.nan
                            
                        if df1.loc[i,q2_unitname]=='ns':
                            q2_T2=df1.loc[i, q2_colname]/1000    
                        elif df1.loc[i,q2_unitname]=='us':
                            q2_T2=df1.loc[i, q2_colname]                            
                        else:
                            q2_T2 = np.nan

                        avg_T2 = 1 / ( 1/q1_T2 + 1/q2_T2 )                                           
                        tau.append( gate_length / avg_T2 *1000)
                    #sys.exit('unit error')
                else:
                    gate_length = df1.loc[:, gate_length_colname]
                    avg_T2 = 1 / ( 1/df1.loc[:, q1_colname] + 1/df1.loc[:, q2_colname] )               
                    tau = gate_length / avg_T2 *1000

                df1[tau_colname] = tau

    tau_df = df1.loc[:, tau_colnames]

    for q1 in range(nqubits):
        for q2 in np.arange(q1+1, nqubits):
            q1str = str(q1)
            q2str = str(q2)
            tau_colname1 = 'cx_' + q1str + '_' + q2str + '_tau_value'
            tau_colname2 = 'cx_' + q2str + '_' + q1str + '_tau_value'
            if (tau_colname1 in tau_df.columns) and (tau_colname2 in tau_df.columns):
                if np.sum(np.abs(tau_df[tau_colname1]-tau_df[tau_colname2]))==0:
                    tau_df = tau_df.drop(tau_colname2,axis=1)
                else:
                    tau_df[tau_colname1] = tau_df[[tau_colname1, tau_colname2]].mean(axis=1)
                    tau_df = tau_df.drop(tau_colname2,axis=1)

    return tau_df

def FI_df(df_all):
    columns = ['query_date', 'last_update_date']
    for i in range(nqubits):
        columns = columns+['q_'+str(i)+'__readout_error_value']
    df = df_all[columns]
    times = pd.to_datetime(df['last_update_date'].str.strip('+'), format='%Y-%m-%dT%H:%M:%S')
    good=list()
    bad=list()
    for i in range(len(times)-1):    
        if (times[i+1] == times[i]): #or df.loc[i,data_label] <=0 or df.loc[i,data_label]>=1:
            bad.append(i)
        else:
            good.append(i)
    
    df = df.loc[good,:].reset_index(drop=True)
    return df

def FG_df(df):
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
    
    df = df[columns]
    times = pd.to_datetime(df['last_update_date'].str.strip('+'), format='%Y-%m-%dT%H:%M:%S')
    good=list()
    bad=list()
    for i in range(len(times)-1):    
        if (times[i+1] == times[i]): #or df.loc[i,data_label] <=0 or df.loc[i,data_label]>=1:
            bad.append(i)
        else:
            good.append(i)
    
    df = df.loc[good,:].reset_index(drop=True)
    return df

def find_region_prob(x,x1, x2):
    x=np.sort(x)
    N = len(x)
    (I1, val) = min(enumerate(x), key=lambda x: abs(x[1]-x1))
    (I2, val) = min(enumerate(x), key=lambda x: abs(x[1]-x2))
    p = (I2-I1)/N
    return p

def empirical_CI(x):
    N = len(x)
    sorted_estimates = np.sort(np.array(x))
    conf_interval = [sorted_estimates[int(0.025 * N)], sorted_estimates[int(0.975 * N)]]
    return conf_interval

def bootstrapped_CI(x, N = 10000):
    mean_estimates = []
    for _ in range(N):
        print(_)
        re_sample_idx = np.random.randint(0, len(x), 1)
        mean_estimates.append(np.mean(x[re_sample_idx]))
    
    sorted_estimates = np.sort(np.array(mean_estimates))
    conf_interval = [sorted_estimates[int(0.025 * N)], sorted_estimates[int(0.975 * N)]]
    return conf_interval

def beta_fitted_pdf(F): # this is not used?
    x=np.arange(start=np.min(F)*.99, stop=np.min((1.0, np.max(F)*1.1)), step=0.01)
    N = len(x)
    dist = getattr(scipy.stats, 'beta')
    params = dist.fit(F) # beningn warning
    a=params[0]
    b=params[1]
    loc = params[2]
    scale = params[3]
    pdf_fitted = dist.pdf(x, a=a, b=b, loc=loc, scale=scale) * N
    return beta_fitted_pdf

def get_hist(F):
    nbins = int( (max(F) - min(F))/( 2*(np.quantile(F,0.75)-np.quantile(F,0.25))/len(F)**(1/3)) )
    (hist_y, hist_bins, patches) = plt.hist(F, density=True, bins=nbins)
    hist_x = [np.mean( [hist_bins[i], hist_bins[i+1]]) for i in range(len(hist_y))]
    hist_x.append(1.0)
    hist_y = np.append(hist_y, np.array([0]))
    hist_y_normalized  = hist_y / np.max(hist_y)
    return (hist_x, hist_y_normalized )

def time_array_FI_FG(df):
    Flist=list()
    for col in df.columns[2:]:
        #print(col)
        y=remove_nan_1_0(df[col])
        y = list( [1.0-float(item) for item in y] )
        Flist = Flist + y
    F = np.array(Flist)
    return F

def latest_time_array_FI_FG(df):
    F = remove_nan_1_0(df.iloc[-1,2:])  
    return F

def plot_CI(F, xlabel, ylabel, figname):
    # three methods
    #confidence_interval3 = scipy.stats.beta.interval(0.95, a, b, loc=loc, scale=scale) # this is dangerous
    #confidence_interval1 = bootstrapped_CI(F)
    confidence_interval2 = empirical_CI(F)
    #my_mean, my_var = scipy.stats.beta.stats(a, b, loc=loc, scale=scale, moments='mv')
    my_mean=np.mean(F)
    confidence_interval = confidence_interval2
    
    (hist_x, hist_y) = get_hist(F)
    orangeline_max = hist_x[np.nanargmin( np.abs( hist_x - confidence_interval[1]) )]
    orangeline_min = hist_x[np.nanargmin( np.abs( hist_x - confidence_interval[0]) )]
    blueline_min = hist_x[np.nanargmin( np.abs( hist_x - (my_mean*(1-tolerance)) ) )]
    blueline_max = hist_x[np.nanargmin( np.abs( hist_x - (my_mean*(1+tolerance))) )]
    
    plt.clf()
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_facecolor('white')
    ax.grid(True, color='lightgrey', alpha=0.8)
    
    plt.plot(hist_x, hist_y, marker = '', 
             linestyle='-', alpha=0.8, 
             color='k', label='Empirical PDF')
    plt.fill_between(hist_x, hist_y, 
                     where=(orangeline_min <= np.array(hist_x)) & (np.array(hist_x) <= orangeline_max ), 
                     color='darkorange', alpha=0.4, label='95% C.I.')
    plt.fill_between(hist_x, hist_y, 
                     where=( blueline_min <= np.array(hist_x)) & (np.array(hist_x) <= blueline_max ), 
                     color='darkblue', alpha=0.4, label='Tolerance')
    plt.grid(True, which='both', alpha=0.4)
    plt.axvline(x=my_mean, color='k', linestyle = '--', label='Mean')
    orangeline_max = hist_x[np.nanargmin( np.abs( hist_x - confidence_interval[1]) )]
    orangeline_min = hist_x[np.nanargmin( np.abs( hist_x - confidence_interval[0]) )]
    plt.axvline(x=orangeline_min, color='darkorange', linestyle = '--')
    plt.axvline(x=orangeline_max, color='darkorange', linestyle = '--')
    plt.legend(loc="top left")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(figname, bbox_inches = "tight", dpi=300)
    
    print("(Mean, Tol+, Tol-) = ", my_mean, my_mean*(1+tolerance), my_mean*(1-tolerance))
    print("Confidence interval = ", confidence_interval)
    print("Tolerance region probability = ", find_region_prob(F,my_mean*(1-tolerance),my_mean*(1+tolerance))*100)
    
def plot_heatmap(last_df, ticks, labels, xylabel, figname, colorbar_label):
    sns_df = pd.DataFrame()
    #sns_df["qubit_pair"]=np.array(last_df.index)
    sns_df["FG"]=np.array(last_df.values)*100 #np.nanmean(df1,axis=0)
    n_elements = len(last_df)
    
    delta_f_mat = np.zeros((n_elements,n_elements))
    for i in range(n_elements):
        delta_f_mat[i,i]=np.nan
        for j in np.arange(i+1, n_elements):
            delta_f_mat[j,i]=np.abs( sns_df["FG"][j]-sns_df["FG"][i] )
            delta_f_mat[i,j] = np.nan
               
    sns.set()
    sns_plot = sns.heatmap(delta_f_mat, robust=True, fmt="f", cmap='Oranges',
                           cbar_kws=colorbar_label)
    sns_plot.figure.set_size_inches(16, int(16*A4ratio))
    sns_plot.set_facecolor('white')
    sns_plot.axes.set_xticks(ticks)
    sns_plot.axes.set_xticklabels(labels, rotation=90, fontsize="12", va="center", position=(0,-.028))    
    sns_plot.axes.set_yticks(ticks)
    sns_plot.axes.set_yticklabels(labels, rotation=0, fontsize="12", va="center")    
    sns_plot.axes.set_ylabel(xylabel, fontsize=12)
    sns_plot.axes.set_xlabel(xylabel, fontsize=12)
    figure = sns_plot.get_figure()
    figure.savefig(figname, bbox_inches = "tight", dpi=300)
    
def plot_raw_ts(df1, rows_per_sublpot, total_ts, figname):
    x  = np.arange(0,df1.shape[0])
    q=0
    for k in np.arange(0,8):
        print(k)
        fig, ax = plt.subplots(nrows=rows_per_sublpot, 
                               sharex=True, 
                               subplot_kw=dict(frameon=False)) # frameon=False removes frames
        fig.set_size_inches(21, 16)
        plt.subplots_adjust(hspace=.08)
    
        for i in range(rows_per_sublpot):
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
    
            if i==rows_per_sublpot-1 or q==total_ts-1:
                ax[i].xaxis.label.set_visible(True)
                ax[i].set_xticks([0,df1.shape[0]])
                ax[i].set_xticklabels(['Dec-21','Mar-22'])
    
            q=q+1
    
        fig.savefig(figname + str(k) + '.png', bbox_inches = "tight", dpi=150)

def plot_sorted_ts(last_df, ticks, labels, xlabel, ylabel, figname):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.set_facecolor('white')
    ax.grid(True, color='lightgrey', alpha=0.8)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation = 90, ha="center")
    plt.plot(last_df.values)
    plt.legend(loc="top left")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(figname, bbox_inches = "tight", dpi=300)

def retrive_data_from_json(filename, foldername):
    import json    
    scale = 16 ## equals to hexadecimal
    f = open(foldername+filename)
    data = json.load(f) 
    f.close()
    
    #counts=data['results'][0]['data']['counts']
    memory=data['results'][0]['data']['memory']

    str_array = []    
    for i in range(len(memory)):
        #print( memory[i][2:] )
        my_hexdata = memory[i][2:]
        x=bin(int(my_hexdata, scale))[2:].zfill(nqubits)
        #print( x )
        str_array.append(x)
    return str_array

def get_entropy(x):
    alphabets = set(x)
    N = len(x)
    h = 0
    for a in alphabets:
        p = np.sum([1 for i in range(len(x)) if (x[i]==a)])/N
        if p != 0:
            h = h + p*np.log(p)/np.log(2)

    return -h

def get_nmi(x,y):
    if len(x)!=len(y):
        sys.exit('Error: Unequal x, y')
    xy = [(x[i],y[i]) for i in range(len(x))]
    nmi = 0 
    hx = get_entropy(x)
    hy = get_entropy(y)
    hxy = get_entropy(xy)
      
    num = hx+hy-hxy
    den = 0.5*(hx+hy)
    if den == 0:
        nmi = np.nan
    else:
        nmi = num/den
    
    return nmi

def get_nmi_df(csvfile, devicename, data_refresh=0):
    if data_refresh:
        filename = '622d13ca25b05bf21476fa4e-output.json'
        foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'\
                     + devicename + '/job-622d13ca25b05bf21476fa4e/'
        results_str = retrive_data_from_json(filename, foldername)
        
        qubit_data_matrix = np.zeros((8192,nqubits))
        for i in range(len(results_str)):
            for k in range(nqubits):
                qubit_data_matrix[i,k]=int(results_str[i][k])
        
        nmi = np.zeros((nqubits, nqubits))+np.nan
        for i in range(nqubits):
            for j in np.arange(i+1, nqubits):
                nmi[i,j] = get_nmi( qubit_data_matrix[:,i], qubit_data_matrix[:,j] )
                print(i,j)
        
        import csv
        with open(csvfile, 'w') as csv_file:  
            writer = csv.writer(csv_file)
            writer.writerow(np.arange(0,nqubits))
            for i in range(nqubits):
               writer.writerow(nmi[i,:])
           
    nmi_df = pd.read_csv(csvfile)
    return nmi_df

def plot_nmi_heatmap(nmi_df, nqubits, colorbar_label, xylabel, figname ):
    sns.set()
    sns_plot = sns.heatmap(nmi_df*100, robust=True, fmt="f", cmap='copper_r',
                           cbar_kws=dict({'label': colorbar_label }))
    #sns_plot.figure.set_size_inches(16, int(16*A4ratio))
    sns_plot.set_facecolor('white')
    ticks = [i for i in range(0, nqubits, max(1,int(nqubits/10)))]
    labels = [i for i in range(0, nqubits, int(nqubits/10))]
    sns_plot.axes.set_xticks(labels)
    sns_plot.axes.set_xticklabels(labels, rotation=90, position=(0,0), fontsize="12", va="center")
    sns_plot.axes.set_yticks(labels)
    sns_plot.axes.set_yticklabels(labels, rotation=0, position=(0,0.28), fontsize="12", va="center")
    sns_plot.axes.set_ylabel(xylabel, fontsize=12)
    sns_plot.axes.set_xlabel(xylabel, fontsize=12)
    figure = sns_plot.get_figure()
    figure.savefig(figname, bbox_inches = "tight", dpi=300)
