import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set_style(style='darkgrid')

time_df = pd.read_csv('time_measures_35000.csv', index_col = 0)
time_df.loc[time_df['Null Image']]
N = np.unique(time_df['N'].values)
time_df['Total time'] = time_df.loc[:,'Voronoi':'Saving'].sum(axis=1)
means_std = time_df.loc[:,'Voronoi':'Total time']
means_std['N'] = time_df['N']
means_std = means_std.groupby('N').agg(['mean', 'std'])
means_std

#%%-----------------------------------------------------------------------------
#PLOT of TIMES w/ dispersion - GOOD

plt.figure(figsize=(12,8))
plt.title('Section Timing', fontsize = 24)
plt.xlabel('N points in ROI', fontsize = 16)
plt.ylabel('Time (s)', fontsize = 16)

cols = list(means_std.columns.levels[0])
cols
colors = ['b','y','g','r','c']

for col, c in zip(cols, colors):
    plt.fill_between(N,
                     means_std[col]['mean'] - means_std[col]['std'],
                     means_std[col]['mean'] + means_std[col]['std'],
                     alpha = 0.6,
                     color = c)
    plt.plot(N, means_std[col]['mean'], 'o-', c = c, markeredgecolor = 'k', label = col)
plt.gca().legend(loc = 'upper left', fontsize = 15)
plt.savefig('time_plot_new.png')


#%%-----------------------------------------------------------------------------
#BOXPLOT of TIMES - BAD
plt.figure(figsize = (12,8))
g = sns.boxplot('N','Total time', data =time_df, hue ='N')
g.legend_.remove()


#%%-----------------------------------------------------------------------------
#BARCHART of TIMES - BAD
cols = list(means_std.columns.levels[0])
barWidth = 1000

plt.figure(figsize = (12,8))
plt.bar(N, means_std[cols[0]]['mean'], width = barWidth)
plt.yscale('log')
for n in range(1,len(cols)):
    plt.bar(N, means_std[cols[n]]['mean'], bottom = means_std[cols[n-1]]['mean'], width = barWidth)
