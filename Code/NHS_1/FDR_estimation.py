#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os
import scipy.stats as st
import matplotlib.pyplot as plt

get_ipython().magic(u'matplotlib inline')


# In[2]:


fdr = pd.read_csv('./permutations.csv')
del fdr['Unnamed: 0']

pval_path = './EWAS_results.csv'
out_df = pd.read_csv(pval_path)

out_df.index = out_df.varname.tolist()
del out_df['varname']
out_df = out_df[['p_val']]

fdr = fdr.transpose()
out_df = pd.concat([out_df, fdr], axis=1)
out_df.columns = ['p_val']+np.arange(0,len(fdr.columns),1).tolist()


# In[3]:


out_df3 = out_df.copy()
out_df3 = out_df3.sort_values(by = ['p_val'])

raw = []
thre = 0.05
n_perm = len(out_df3.columns)-1
print 'NUMBER OF PERMUTATIONS: ', n_perm


for i in range(len(out_df3)):
    num = 0
    for j in out_df3.columns[1:]:
        num += np.sum(np.where(out_df3[j] <= out_df3.iloc[i].p_val,1,0))
        
    denom = np.sum(np.where(out_df3.p_val <= out_df3.iloc[i].p_val,1,0))
    raw.append((np.float(num)/n_perm)/denom)


# In[4]:


raw2 = []
for i in range(len(raw)):
    raw2.append(np.min(raw[i:]))
out_df3['fdr'] = raw2
out_df3 = out_df3[['p_val','fdr']]

# an FDR of 0.05 is specified to determine significant assoications
thre = 0.05
print 'NUMBER OF SIGNIFICANT TESTS : ', thre , len(out_df3[out_df3.fdr < thre])


# In[5]:


out_df3.to_csv('./fdr.csv')

