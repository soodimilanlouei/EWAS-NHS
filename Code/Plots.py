#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import os
from matplotlib import font_manager as fm, rcParams

get_ipython().magic(u'matplotlib inline')


# In[2]:


fpath = "./HelveticaNeue-Light.otf" 
prop = fm.FontProperties(fname=fpath)


# # Figure 1

# In[3]:


df = pd.read_excel('./KG-data.xlsx', sheet_name='ALL')
g = nx.MultiGraph()

for i in range(len(df)):
    g.add_node(df.lab[i], mode = 'exposure', cat = df['CAT'][i])
    g.add_node(df.phenotype[i], mode = 'phenotype')
    
for i in range(len(df)):
    g.add_edges_from([(df.lab[i], df.phenotype[i], dict(effect = df.association[i]))])

nx.write_gml(g, './full-KG.gml')    


# # Figure 3-a

# In[4]:


df = pd.read_csv('./ewas.csv')
df = df.sort_values(by = 'exp_coef')
df['Y'] = -np.log(df.p_val)
desc = pd.read_excel('./variable_description.xlsx')
df = pd.merge(df, desc, on='varname', how='outer')

fpath = "./HelveticaNeue-Light.otf" 
prop = fm.FontProperties(fname=fpath)
fname = os.path.split(fpath)[1]
plt.rcParams['xtick.major.pad']='8'
plt.rcParams['ytick.major.pad']='8'
plt.rcParams['axes.labelpad']='20'

plt.figure(figsize = (20,6))
c = 'gray'

thre_df = pd.read_csv('./fdr.csv')
thre_df = thre_df[thre_df.fdr < 0.05]
thre = thre_df.iloc[-1].p_val

for i in range(len(df)):
    if (df.p_val[i] <= thre ) and (df.exp_coef[i] < 1):
        c = 'green'#'#0099bf'
        a = 0.7
        
    elif (df.p_val[i] <= thre ) and (df.exp_coef[i] > 1):
        c = 'red'
        a = 0.7
    else: 
        c = 'gray'
        a = 0.6
        
    if df['Type'][i] == 'Food':
        plt.scatter(df.exp_coef[i], df.Y[i], color = c, label = None, alpha = a, marker = 'D')#'-',)
    else:
        
        plt.scatter(df.exp_coef[i], df.Y[i], color = c, label = None, alpha = a)#'-',)
        
    plt.bar(df.exp_coef[i], df.Y[i], width=0.0004, color = c, alpha = a)
       
    
plt.plot( [np.min(df.exp_coef)-0.01, np.max(df.exp_coef)+0.01], [-np.log(0.05), -np.log(0.05)], 
         linestyle = '--', color = '#fa8072', linewidth = 0.6)

plt.plot( [np.min(df.exp_coef)-0.01, np.max(df.exp_coef)+0.01], [-np.log(thre), -np.log(thre)], 
         linestyle = '-', color = '#f3553c',  linewidth = 0.6)

plt.yticks(fontsize = 15, fontproperties=prop, fontname = 'HelveticaNeue-Light' )
plt.xticks(fontsize = 15, fontproperties=prop, fontname = 'HelveticaNeue-Light' )
plt.ylabel('-log(p-value)', fontsize = 20, fontproperties=prop, fontname = 'HelveticaNeue-Light' )
plt.xlabel('hazard ratio', fontsize = 20, fontproperties=prop, fontname = 'HelveticaNeue-Light' )
plt.xlim(np.min(df.exp_coef)-0.015, np.max(df.exp_coef) + 0.015)
plt.ylim(0,25)

plt.show()


# # Figure 3-b

# In[5]:


df = pd.read_excel('./Food_composition.xlsx')

g = nx.Graph()
w = []

for i in range(4, len(df)):
    for j in range(6, len(df.columns)):
        d = 'grey'
        if df[df.columns[j]][i] != 0:
        
            if df['direction'][i] == 0:
                g.add_node(df.Label[i], lab =  'beneficial', mod = 'nutrient', HR = np.abs(np.log(df['HR'][i])))
            else:
                g.add_node(df.Label[i], lab =  'harmful', mod = 'nutrient', HR =   np.abs(np.log(df['HR'][i])))
            
            if df.iloc[2][j] == 0:
                

                g.add_node(df.iloc[1][j], lab = 'beneficial', mod = 'food', HR =  np.abs(np.log(df.iloc[3][j])))
            else:
                g.add_node(df.iloc[1][j], lab = 'harmful', mod = 'food', HR =  np.abs(np.log(df.iloc[3][j])))
                
                
            if df['HR'][i] > 1 and df.iloc[3][j] > 1:
                d = 'red'
            if df['HR'][i] < 1 and df.iloc[3][j] < 1:
                d = 'green'

                
            g.add_edge(df.Label[i], df.iloc[1][j], weight = df[df.columns[j]][i], direction = d)
            w.append(df[df.columns[j]][i])
            
nx.write_gml(g, 'food-nutrient-network.gml')


# # Figure 4-b

# In[6]:


weighted_HR = {}
for i in range(4, len(df)):
    try:
        aa = df.iloc[i][6:]
        bb = np.sum(df.iloc[i][6:])
        w = [np.float(x)/bb for x in aa]
        food_hr =  df.iloc[3][6:].tolist()
        
        power = []
        for j in range(len(w)):
            power.append(np.power(food_hr[j], w[j]))
        
        weighted_HR[df.Label[i]] = np.prod(power)
    except: 
        ZeroDivisionError
        print i,df.Label[i]
        
HR = {}
for i in range(4, len(df)):
    HR[df.Label[i]] = df.HR[i]


# In[7]:


plt.rcParams['xtick.major.pad']='8'
plt.rcParams['ytick.major.pad']='8'
plt.rcParams['axes.labelpad']='20'

plt.figure(figsize = (8,8))
for i in weighted_HR:
    if HR[i] > 1 and weighted_HR[i]>1:
        plt.scatter(HR[i],weighted_HR[i], alpha = .8, color = 'red')
        
    elif HR[i] <1 and weighted_HR[i] < 1:
        plt.scatter(HR[i],weighted_HR[i], alpha = .8, color = 'green')
    
    else:
        plt.scatter(HR[i],weighted_HR[i], alpha = .8, color = 'orange')
        

plt.plot([0.85,1.21],[1,1],'--',color = 'grey', linewidth = 1)
plt.plot([1,1],[0.85,1.21], '--', color = 'grey', linewidth = 1)
plt.xlabel('Actual HR', fontsize = 15, fontproperties=prop, fontname = 'HelveticaNeue-Light')
plt.ylabel('Werighted food HR', fontsize = 15, fontproperties=prop, fontname = 'HelveticaNeue-Light')
plt.xticks(fontproperties=prop, fontname = 'HelveticaNeue-Light', fontsize = 12)
plt.yticks(fontproperties=prop, fontname = 'HelveticaNeue-Light', fontsize = 12)
plt.xlim(0.85,1.21)
plt.ylim(0.85,1.21)
plt.show()

