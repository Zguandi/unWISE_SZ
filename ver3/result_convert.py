import numpy as np
import pandas as pd

result0_path = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/namaster_comparison0.txt'
result1_path = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/namaster_comparison1.txt'
result2_path = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/namaster_comparison2.txt'

OUTPATH = '/mnt/c/Users/gdzhao/projects/unwise_sz/ver3/result/nobeam/'

data0 = pd.read_csv(result0_path,sep=' ')
header0 = data0.columns
data1 = pd.read_csv(result1_path,sep=' ')
header1 = data1.columns
data2 = pd.read_csv(result2_path,sep=' ')
header2 = data2.columns

data_gy = []
header_gy = []
data_yy = []
header_yy = []

values_ell = data0.values[:,0]

################################
### Add first column
################################
header_gy.append('ell')
data_gy.append(values_ell)
header_yy.append('ell')
data_yy.append(values_ell)

for i in range(len(header0)):
    if header0[i].endswith('_gy'):
        header_gy.append(header0[i])
        data_gy.append(data0.values[:,i])
    if header1[i].endswith('_gy'):
        header_gy.append(header1[i])
        data_gy.append(data1.values[:,i])
    if header2[i].endswith('_gy'):
        header_gy.append(header2[i])
        data_gy.append(data2.values[:,i])
    if header0[i].endswith('_yy'):
        header_yy.append(header0[i])
        data_yy.append(data0.values[:,i])
    if header1[i].endswith('_yy'):
        header_yy.append(header1[i])
        data_yy.append(data1.values[:,i])
    if header2[i].endswith('_yy'):
        header_yy.append(header2[i])
        data_yy.append(data2.values[:,i])

df_gy = pd.DataFrame(data_gy).T
df_gy.columns = header_gy
df_yy = pd.DataFrame(data_yy).T
df_yy.columns = header_yy

df_gy.to_csv(OUTPATH+'namaster_comparison_gy.csv',sep=',',index=False)
df_yy.to_csv(OUTPATH+'namaster_comparison_yy.csv',sep=',',index=False)