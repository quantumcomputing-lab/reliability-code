import matplotlib.pyplot as plt
import calendar
import pandas as pd
import numpy as np
import scipy
import scipy.stats

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

# POINT-IN-TIME
filename = 'washington_readout_error_mar10_2022.csv'
foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'
df = pd.read_csv(foldername +filename)
columns = []
for i in range(127):
    columns = columns+['q'+str(i)+'_readout_err']
df0 = df[columns]
df1=df0.iloc[-1,:].reset_index(drop=True)
error_rate = np.array(df1)
N = len(error_rate)
delta_f = np.zeros(int(N*(N-1)/2))*np.nan
s=0
for i in range(N):
    for j in np.arange(i+1,N):
        delta_f[s] = np.abs(error_rate[j]-error_rate[i])
        s=s+1

user_defined_tol = np.float64(0.02)

F=delta_f
x=np.arange(start=np.min(F)*.99, stop=np.min((1.0, np.max(F)*1.1)), step=0.01)
N = len(x)

dist = getattr(scipy.stats, 'beta')
params = dist.fit(F) # beningn warning
a=params[0]
b=params[1]
loc = params[2]
scale = params[3]
pdf_fitted = dist.pdf(x, a=a, b=b, loc=loc, scale=scale) * N

# Freedman-Diaconis rule (an imporvemeny over Sturge's rule)
nbins = (max(F) - min(F))/( 2*(np.quantile(F,0.75)-np.quantile(F,0.25))/len(F)**(1/3))
(hist_y, hist_bins, patches) = plt.hist(F, density=True, bins=16)
hist_y = list(hist_y)
hist_x = [np.mean( [hist_bins[i], hist_bins[i+1]]) for i in range(len(hist_y))]
hist_x.insert(0, 0.0)
hist_y.insert(0, 0)
hist_y_normalized  = hist_y / np.max(hist_y)

# three methods
confidence_interval3 = scipy.stats.beta.interval(0.95, a, b, loc=loc, scale=scale) # this is dangerous
confidence_interval1 = bootstrapped_CI(F)
confidence_interval2 = empirical_CI(F)
#my_mean, my_var = scipy.stats.beta.stats(a, b, loc=loc, scale=scale, moments='mv')
my_mean = np.nanmean(F)
confidence_interval = confidence_interval2

plt.clf()
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.set_facecolor('white')
ax.grid(True, color='lightgrey', alpha=0.8)
#plt.plot(x, pdf_fitted/max(pdf_fitted), label='pdf', linewidth = 1.1, alpha=0.8, color='k')
plt.plot(hist_x, hist_y_normalized, marker = '', 
         linestyle='-', alpha=0.8, 
         color='k', label='Empirical PDF')
plt.fill_between(hist_x, hist_y_normalized, 
                 where=(confidence_interval[0] <= hist_x) & (hist_x <= confidence_interval[1] ), 
                 color='darkorange', alpha=0.3, label='95% C.I.')
plt.fill_between(hist_x, hist_y_normalized, 
                 where= (hist_x < user_defined_tol), 
                 color='darkblue', alpha=0.4, label='Tolerance')
#plt.fill_between(x, pdf_fitted/max(pdf_fitted), 
#                 where=(confidence_interval[0] < x) & (x < confidence_interval[1] ), 
#                 color='grey', alpha=0.4, label='CI')
#plt.xlim(0.8, 1.01)
plt.grid(True, which='both', alpha=0.4)
plt.axvline(x=my_mean, color='k', linestyle = '--', label='Mean')
plt.axvline(x=confidence_interval[1], color='darkorange', linestyle = '--')
plt.axvline(x=confidence_interval[0], color='darkorange', linestyle = '--')
plt.legend(loc="top right")
plt.xlabel("Inter-qubit difference in initialization fidelity")
plt.ylabel("Normalized histogram")
#plt.plot(hist_x, hist_y/np.max(hist_y), marker = '', linestyle='-', alpha=0.8, color='darkorange')
#plt.scatter(hist_x, hist_y/max(hist_y), color='darkorange', alpha=0.8)
plt.savefig('FI_point-in-time.png', facecolor=fig.get_facecolor(), bbox_inches = "tight", dpi=300)

#print mean
print(my_mean)
# user-defined tolerance region
print(0, user_defined_tol)
# empirical CI
print(confidence_interval) 
# empirical probability that val is within tolerance
print(find_region_prob(F,0,user_defined_tol )*100)

# ACROSS TIME
filename = 'washington_readout_error_mar10_2022.csv'
foldername = '/home/samudra/Desktop/sjv_stability_feb_2022/data/'
df = pd.read_csv(foldername +filename)

Flist=list()
nqubits = 27
user_defined_tol = 0.01

for q in range(nqubits):
    q=str(q)
    time_label = 'q'+str(q)+'_readout_err_cal_time'
    data_label = 'q'+str(q)+'_readout_err'
    
    times = pd.to_datetime(df[time_label].str.strip('+'), format='%Y-%m-%dT%H:%M:%S')
    good=list()
    bad=list()
    for i in range(len(times)-1):    
        if (times[i+1] == times[i]): #or df.loc[i,data_label] <=0 or df.loc[i,data_label]>=1:
            bad.append(i)
        else:
            good.append(i)
    
    df_temp = df.loc[good,:].reset_index(drop=True)
    df_new = df_temp.loc[:,(data_label, time_label)].reset_index(drop=True)
    Flist.append(1-df_new[data_label])

F = np.array([item for sublist in Flist for item in sublist])

N = len(x)
#x=np.arange(start=np.min(F)*.99, stop=np.min((1.0, np.max(F)*1.1)), step=user_defined_tol)
nbins = int( (max(F) - min(F))/( 2*(np.quantile(F,0.75)-np.quantile(F,0.25))/len(F)**(1/3)) )
#x = np.linspace(start=np.min(F)*.99, stop=np.min((1.0, np.max(F)*1.1)), bin=nbins)

dist = getattr(scipy.stats, 'beta')
params = dist.fit(F) # beningn warning
a=params[0]
b=params[1]
loc = params[2]
scale = params[3]
pdf_fitted = dist.pdf(x, a=a, b=b, loc=loc, scale=scale) * N

(hist_y, hist_bins, patches) = plt.hist(F, density=True, bins=nbins)
hist_x = [np.mean( [hist_bins[i], hist_bins[i+1]]) for i in range(len(hist_y))]
hist_x.append(1.0)
hist_y = np.append(hist_y, np.array([0]))
hist_y_normalized  = hist_y / np.max(hist_y)

# three methods
confidence_interval3 = scipy.stats.beta.interval(0.95, a, b, loc=loc, scale=scale) # this is dangerous
confidence_interval1 = bootstrapped_CI(F)
confidence_interval2 = empirical_CI(F)
#my_mean, my_var = scipy.stats.beta.stats(a, b, loc=loc, scale=scale, moments='mv')
my_mean=np.mean(F)
confidence_interval = confidence_interval2

plt.clf()
#plt.plot(x, pdf_fitted/max(pdf_fitted), label='pdf', linewidth = 1.1, alpha=0.8, color='k')
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.set_facecolor('white')
ax.grid(True, color='lightgrey', alpha=0.8)

plt.plot(hist_x, hist_y_normalized, marker = '', 
         linestyle='-', alpha=0.8, 
         color='k', label='Empirical PDF')
plt.fill_between(hist_x, hist_y_normalized, 
                 where=(confidence_interval[0] < hist_x) & (hist_x < confidence_interval[1] ), 
                 color='darkorange', alpha=0.4, label='95% C.I.')
plt.fill_between(hist_x, hist_y_normalized, 
                 where=(my_mean-user_defined_tol < hist_x) & (hist_x < my_mean+user_defined_tol ), 
                 color='darkblue', alpha=0.4, label='Tolerance')
#plt.fill_between(x, pdf_fitted/max(pdf_fitted), 
#                 where=(confidence_interval[0] < x) & (x < confidence_interval[1] ), 
#                 color='grey', alpha=0.4, label='CI')
#plt.xlim(0.8, 1.01)
plt.grid(True, which='both', alpha=0.4)
plt.axvline(x=my_mean, color='k', linestyle = '--', label='Mean')
plt.axvline(x=confidence_interval[1], color='darkorange', linestyle = '--')
plt.axvline(x=confidence_interval[0], color='darkorange', linestyle = '--')
plt.legend(loc="top left")
plt.xlabel("Initialization fidelity")
plt.ylabel("Normalized histogram")
#plt.plot(hist_x, hist_y/np.max(hist_y), marker = '', linestyle='-', alpha=0.8, color='darkorange')
#plt.scatter(hist_x, hist_y/max(hist_y), color='darkorange', alpha=0.8)
plt.savefig('FI_across-time.png', bbox_inches = "tight", dpi=300)

print(my_mean, my_mean+user_defined_tol, my_mean-user_defined_tol)
print(confidence_interval)
print("Probability = ", find_region_prob(F,my_mean-user_defined_tol,my_mean+user_defined_tol )*100)