import matplotlib.pyplot as plt
import calendar
import pandas as pd
import numpy as np
import scipy
import scipy.stats


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


q='0'
time_label = 'q'+str(q)+'_readout_err_cal_time'
data_label = 'q'+str(q)+'_readout_err'

df = pd.read_csv('/home/samudra/Desktop/sjv_dissertation/stability-project/data/yorktown_daily_readout_error.csv')
times = pd.to_datetime(df[time_label].str.strip('+'), format='%Y-%m-%dT%H:%M:%S')
good=list()
bad=list()
for i in range(len(times)-1):    
    if (times[i+1] == times[i]): #or df.loc[i,data_label] <=0 or df.loc[i,data_label]>=1:
        bad.append(i)
    else:
        good.append(i)

df = df.loc[good,:].reset_index(drop=True)
df_new = df.loc[:,(data_label, time_label)].reset_index(drop=True)
F = 1-df_new[data_label]

x=np.arange(start=np.min(F)*.99, stop=np.min((1.0, np.max(F)*1.1)), step=0.01)
N = len(x)
#x = np.linspace(start=np.min(F)*.99, stop=np.min((1.0, np.max(F)*1.1)), num=N)

dist = getattr(scipy.stats, 'beta')
for i in (0,1):
    params = dist.fit(F) # beningn warning
a=params[0]
b=params[1]
loc = params[2]
scale = params[3]
pdf_fitted = dist.pdf(x, a=a, b=b, loc=loc, scale=scale) * N

(hist_y, hist_bins, patches) = plt.hist(F, density=True, bins=x)
hist_x = [np.mean( [hist_bins[i], hist_bins[i+1]]) for i in range(len(hist_y))]

# three methods
confidence_interval3 = scipy.stats.beta.interval(0.95, a, b, loc=loc, scale=scale)
confidence_interval1 = bootstrapped_CI(F)
confidence_interval2 = empirical_CI(F)
my_mean, my_var = scipy.stats.beta.stats(a, b, loc=loc, scale=scale, moments='mv')

confidence_interval = confidence_interval2

plt.clf()
plt.plot(x, pdf_fitted/max(pdf_fitted), label='pdf', linewidth = 1.1, alpha=0.8, color='k')
plt.fill_between(x, pdf_fitted/max(pdf_fitted), 
                 where=(confidence_interval[0] < x) & (x < confidence_interval[1] ), 
                 color='grey', alpha=0.4, label='CI')
plt.xlim(0.8, x[-1])
plt.grid(True, which='both', alpha=0.4)
plt.axvline(x=my_mean, color='r', linestyle = '--', label='mean')
plt.legend(loc="top left")
plt.xlabel("xlabel")
plt.ylabel("ylabel")
#plt.plot(hist_x, hist_y/np.max(hist_y), marker = '', linestyle='--', alpha=0.8, color='darkorange')
#plt.scatter(hist_x, hist_y/max(hist_y), color='darkorange', alpha=0.8)
plt.savefig('FI_temporal_q'+str(q)+'.png', bbox_inches = "tight", dpi=300)
