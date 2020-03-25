import pandas as pd

time_df = pd.read_csv('time_measures.csv', index_col=False)
time_df['Total time'] = time_df.iloc[:,3:7].sum(axis=1)
time_df.iloc[:,[1,3,4,5,6,8]].groupby('N').agg(['mean', 'std'])
