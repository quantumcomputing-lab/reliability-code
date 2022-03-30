import matplotlib.pyplot as plt
import calendar
import pandas as pd
import numpy as np
import scipy
import scipy.stats
import seaborn as sns
import sys
import math

'''
TODOS:
run tau using this script
run eta using this script
adjust dataframe: all entries which have the same date stamp, take their last value
run other devices
slot into html nicely

def tau_df(df): # duty cycle
def eta_df(df): # addressability
def data_update():

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
'''

A4ratio = 297/210
filename = 'washington_readout_error_mar10_2022.csv'
foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'
df_all = pd.read_csv(foldername +filename)
nqubits = 127
tolerance = 0.01

# CI ACROSS TIME
# Initialization Fidelity
df = FI_df(df_all)
plot_CI(df, xlabel="Initialization fidelity", ylabel="Normalized histogram", 
        figname = 'FI_across-time.png')

# Gate Fidelity
df = FG_df(df_all)
plot_CI(F, xlabel="CNOT fidelity", ylabel="Normalized histogram", 
        figname = 'FG_across-time.png')

# CI POINT-IN-TIME
df = FG_df(df_all)
error_rate = remove_nan_1_0(df.iloc[-1,2:])  
N = len(error_rate)
delta_f = np.zeros(int(N*(N-1)/2))*np.nan
s=0
for i in range(N):
    for j in np.arange(i+1,N):
        delta_f[s] = np.abs(error_rate[j]-error_rate[i])
        s=s+1

tolerance = np.float64(0.02)
F=delta_f
(hist_x, hist_y) = beta_fitted_pdf(F)
plot_CI(F, xlabel="Inter-gate difference in CNOT fidelity", 
        ylabel="Normalized histogram", figname = 'FG_point-in-time.png')

#### RAW DATA POINT IN TIME
df = FG_df(df_all)
last_df = remove_nan_1_0(df.iloc[-1][2:])
npairs = len(last_df)
ticks = [i for i in range(0, npairs, 6)]
labels = [item.replace('_gate_error_value', '') for item in np.array(sns_df["qubit_pair"][ticks])]
xylabel = 'Pair #'
figname = 'FG_raw_lattice.png'
colorbar_label = 'label': 'Difference in Fidelity (%)'

plot_heatmap(last_df, ticks, labels, xylabel, figname, colorbar_label):
        
#### RAW DATA ACROSS TIME
df = FG_df(df_all)
df1 = df.iloc[:,2:]
rows_per_sublpot=18
total_ts=144
figname = 'FG_raw_time_'
plot_raw_ts(df1, rows_per_sublpot, total_ts, figname)

#### GENERAL ROUTINES
def remove_nan_1_0(series_s):
    N_rows  = [not math.isnan(item) for item in series_s]
    N_rows = [a and b for a, b in zip( N_rows  , series_s!=0.0)]
    N_rows = [a and b for a, b in zip(N_rows, series_s!=1.0)]
    series_s = series_s[N_rows]
    return series_s

def FI_df(df_all):
    columns = ['input_date', 'last_update_date']
    for i in range(nqubits):
        columns = columns+['q'+str(i)+'_readout_err']
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
    nbins = int( (max(F) - min(F))/( 2*(np.quantile(F,0.75)-np.quantile(F,0.25))/len(F)**(1/3)) )
    dist = getattr(scipy.stats, 'beta')
    params = dist.fit(F) # beningn warning
    a=params[0]
    b=params[1]
    loc = params[2]
    scale = params[3]
    pdf_fitted = dist.pdf(x, a=a, b=b, loc=loc, scale=scale) * N
    return beta_fitted_pdf

def get_hist(F):
    (hist_y, hist_bins, patches) = plt.hist(F, density=True, bins=nbins)
    hist_x = [np.mean( [hist_bins[i], hist_bins[i+1]]) for i in range(len(hist_y))]
    hist_x.append(1.0)
    hist_y = np.append(hist_y, np.array([0]))
    hist_y_normalized  = hist_y / np.max(hist_y)
    return (hist_x, hist_y_normalized )

def plot_CI(F, xlabel="CNOT fidelity", ylabel="Normalized histogram", figname = 'FG_across-time.png'):
    Flist=list()
    for col in df.columns[2:]:
        print(col)
        y=remove_nan_1_0(df[col])
        y = list( [1.0-float(item) for item in y] )
        Flist = Flist + y
    F = np.array(Flist)
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
    blueline_min = hist_x[np.nanargmin( np.abs( hist_x - (my_mean-tolerance) ) )]
    blueline_max = hist_x[np.nanargmin( np.abs( hist_x - (my_mean+tolerance)) )]
    
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
    
    print("(Mean, Tol+, Tol-) = ", my_mean, my_mean+tolerance, my_mean-tolerance)
    print("Confidence interval = ", confidence_interval)
    print("Tolerance region probability = ", find_region_prob(F,my_mean-tolerance,my_mean+tolerance )*100)
    
def plot_heatmap(last_df, ticks, labels, xylabel = 'Pair #',
                 figname = 'FG_raw_lattice.png',
                 colorbar_label = {'label': 'Difference in Fidelity (%)'}):
    sns_df = pd.DataFrame()
    sns_df["qubit_pair"]=np.array(last_df.index)
    sns_df["FG"]=np.array(last_df.values)*100 #np.nanmean(df1,axis=0)
    
    delta_f_mat = np.zeros((npairs,npairs))
    for i in range(npairs):
        delta_f_mat[i,i]=np.nan
        for j in np.arange(i+1, npairs):
            delta_f_mat[j,i]=np.abs( sns_df["FG"][j]-sns_df["FG"][i] )
            delta_f_mat[i,j] = np.nan
               
    sns.set()
    sns_plot = sns.heatmap(delta_f_mat, robust=True, fmt="f", cmap='Oranges',
                           cbar_kws={colorbar_label})
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
