#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 21:16:04 2025

@author: joannafernandez
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os 

#import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes

siperc = Find("/path/to/siH2AX_percentages.csv")
siperc= siperc.dropna() #remove empty rows
sipval = Find("/path/to/siH2AX_pvals.csv")
sipval= sipval.dropna() #remove empty rows

#where p>0.001 (non-significant indel) remove from consideration
siperc.loc[:, ~siperc.columns.isin(['replicate', 'indels'])] = siperc.loc[:, ~siperc.columns.isin(['replicate', 'indels'])].mask(sipval > 0.001, 0)

#formating for plotting
singles = siperc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['siSCR_6hr', 'siSCR_12hr', 'siSCR_24hr',
                                      'siSCR_48hr', 'siH2AX_6hr', 'siH2AX_12hr', 'siH2AX_24hr',
                                      'siH2AX_48hr'],
                               value_name='percentage')
singles['percentage'] = singles['percentage'].astype(int)

#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    else:
        return None


singles['class'] = singles['indels'].apply(classify_indel)
singles = singles[singles['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    singles
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df[summed_df['class'].isin(['NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


custom_palette = {
    'NHEJ': '#575757',
    'uncut': '#EFB54F'
}

time_order = ['6hr', '12hr', '24hr', '48hr']
condition_order = ['siSCR', 'siH2AX']
class_order = ['NHEJ', 'uncut']

df = subset.copy()
df['time'] = pd.Categorical(df['time'], categories=time_order, ordered=True)
df['condition'] = pd.Categorical(df['condition'], categories=condition_order, ordered=True)
df['class'] = pd.Categorical(df['class'], categories=class_order, ordered=True)

# palette as dict {class: color}
if not isinstance(custom_palette, dict):
    pal_list = sns.color_palette(n_colors=len(class_order))
    custom_palette = dict(zip(class_order, pal_list))

means = (
    df.groupby(['time','condition','class'], observed=True)['relative_percentage']
      .mean()
      .unstack('class')
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=['time','condition']),
               columns=class_order)
      .fillna(0.0)
)


ncols = len(time_order)
fig, axes = plt.subplots(1, ncols, figsize=(4.8*ncols, 4.5), sharey=True)
axes = np.atleast_1d(axes)

x_pos = {'siSCR': -0.35, 'siH2AX': 0.35}
bar_width = 0.55

for j, t in enumerate(time_order):
    ax = axes[j]

    for cond in condition_order:
        hrow = means.loc[(t, cond)]
        bottom = 0.0
        for cls in class_order:
            h = float(hrow.get(cls, 0.0))
            ax.bar(x_pos[cond], h, width=bar_width, bottom=bottom,
                   color=custom_palette.get(cls, 'C0'),
                   edgecolor='black', linewidth=1.5, alpha=0.8, zorder=2)
            bottom += h

        dft = df[(df['time'] == t) & (df['condition'] == cond) & (df['class'] == 'NHEJ')]
        vals = dft['relative_percentage'].to_numpy()
        if vals.size:
            jitter = (np.random.rand(vals.size) - 0.5) * (bar_width * 0.25)
            ax.scatter(np.full(vals.size, x_pos[cond]) + jitter, vals, s=180,
                       facecolors='white', edgecolors='black',
                       linewidths=1.2, zorder=4, clip_on=False)

            if vals.size > 1:
                mean_val = float(vals.mean())
                sem = float(vals.std(ddof=1) / np.sqrt(vals.size))

    ax.set_title(t, fontsize=12)
    ax.set_xlim(-0.9, 0.9)
    ax.set_xticks([x_pos['siSCR'], x_pos['siH2AX']])
    ax.set_xticklabels(['siSCR', 'siH2AX'], fontsize=11)
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.axhline(50, linestyle='--', color='grey', linewidth=1)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(axis='both', width=1.6, length=4)
for ax in np.ravel(axes):
    ax.tick_params(axis='y', labelsize=22)  
    

sns.despine()
plt.tight_layout()
plt.show()
    
plt.savefig("sih2axfig1.svg", dpi=300, bbox_inches='tight')
#%%

#now for ATMI
atmiperc = Find("/path/to/ATMi_percentages.csv")
atmiperc= atmiperc.dropna() #remove empty rows
atmipval = Find("/path/to/ATMi_pvals.csv")
atmipval= atmipval.dropna() #remove empty rows

#where p>0.001 (non-significant indel) remove from consideration
atmiperc.loc[:, ~atmiperc.columns.isin(['replicate', 'indels'])] = atmiperc.loc[:, ~atmiperc.columns.isin(['replicate', 'indels'])].mask(atmipval > 0.001, 0)

#formating for plotting
inhibitor = atmiperc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['DMSO_6hr', 'DMSO_12hr', 'DMSO_24hr',
                                      'DMSO_48hr', 'ATMi_6hr', 'ATMi_12hr', 'ATMi_24hr', 'ATMi_48hr'],
                               value_name='percentage')
inhibitor['percentage'] = inhibitor['percentage'].astype(int)

#assign indel type 
inhibitor['class'] = inhibitor['indels'].apply(classify_indel)
inhibitor = inhibitor[inhibitor['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    inhibitor
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)

subset.to_csv('atmi260925.csv')


#%%
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests



#%%
subsetnhej = subset[(subset['class'] == 'NHEJ')]
df = subsetnhej.copy()

df['condition'] = pd.Categorical(df['condition'], ['DMSO','ATMi'])
df['time'] = pd.Categorical(df['time'], ['6hr','12hr','24hr','48hr'], ordered=True)

# 1) TWO-WAY ANOVA (overall)
model = smf.ols(
    'relative_percentage ~ C(condition, Treatment(reference="DMSO")) * '
    'C(time, Treatment(reference="6hr"))',
    data=df
).fit()
anova_tbl = sm.stats.anova_lm(model, typ=2)
print(anova_tbl)

# 2) POST-HOC: ATMi vs DMSO within each time, Šidák-corrected
term_cond = 'C(condition, Treatment(reference="DMSO"))[T.ATMi]'

def contrast_str(t):
    if t == '6hr':
        return f'{term_cond} = 0'
    else:
        return f'{term_cond} + {term_cond}:C(time, Treatment(reference="6hr"))[T.{t}] = 0'

rows = []
for t in df['time'].cat.categories:
    tt = model.t_test(contrast_str(t))
    ci_low, ci_high = tt.conf_int().ravel()
    rows.append({
        'time': t,
        'diff_ATMi_minus_DMSO': float(tt.effect),
        'SE': float(tt.sd),
        't': float(tt.tvalue),
        'p_raw': float(tt.pvalue),
        'CI95_low': ci_low,
        'CI95_high': ci_high
    })

posthoc = pd.DataFrame(rows)
posthoc['p_sidak'] = multipletests(posthoc['p_raw'], method='sidak')[1]
print(posthoc[['time','diff_ATMi_minus_DMSO','p_raw','p_sidak','CI95_low','CI95_high']])

#%%


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch


time_order = ['6hr', '12hr', '24hr', '48hr']
condition_order = ['DMSO', 'ATMi']
class_order = ['NHEJ', 'uncut']  

df = subset.copy()
df['time'] = pd.Categorical(df['time'], categories=time_order, ordered=True)
df['condition'] = pd.Categorical(df['condition'], categories=condition_order, ordered=True)
df['class'] = pd.Categorical(df['class'], categories=class_order, ordered=True)


if not isinstance(custom_palette, dict):
    pal_list = sns.color_palette(n_colors=len(class_order))
    custom_palette = dict(zip(class_order, pal_list))

means = (
    df.groupby(['time','condition','class'], observed=True)['relative_percentage']
      .mean()
      .unstack('class')
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=['time','condition']),
               columns=class_order)
      .fillna(0.0)
)

ncols = len(time_order)
fig, axes = plt.subplots(1, ncols, figsize=(4.8*ncols, 4.5), sharey=True)
axes = np.atleast_1d(axes)


x_pos = {'DMSO': -0.35, 'ATMi': 0.35}
bar_width = 0.55

for j, t in enumerate(time_order):
    ax = axes[j]


    for cond in condition_order:
        hrow = means.loc[(t, cond)]
        bottom = 0.0
        for cls in class_order:
            h = float(hrow.get(cls, 0.0))
            ax.bar(x_pos[cond], h, width=bar_width, bottom=bottom,
                   color=custom_palette.get(cls, 'C0'),
                   edgecolor='black', linewidth=1.5, alpha=0.8, zorder=2)
            bottom += h

        # NHEJ: replicate points + SEM
        dft = df[(df['time'] == t) & (df['condition'] == cond) & (df['class'] == 'NHEJ')]
        vals = dft['relative_percentage'].to_numpy()
        if vals.size:
            jitter = (np.random.rand(vals.size) - 0.5) * (bar_width * 0.25)
            ax.scatter(np.full(vals.size, x_pos[cond]) + jitter, vals, s=180,
                       facecolors='white', edgecolors='black',
                       linewidths=1.2, zorder=4, clip_on=False)

            if vals.size > 1:
                mean_val = float(vals.mean())
                sem = float(vals.std(ddof=1) / np.sqrt(vals.size))
                ax.errorbar(x_pos[cond], mean_val, yerr=sem, fmt='none',
                            ecolor='black', elinewidth=1.8, capsize=6, capthick=1.8, zorder=5)


    ax.set_title(t, fontsize=12)
    ax.set_xlim(-0.9, 0.9)
    ax.set_xticks([x_pos['DMSO'], x_pos['ATMi']])
    ax.set_xticklabels(['DMSO', 'ATMi'], fontsize=11)
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.axhline(50, linestyle='--', color='grey', linewidth=1)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(axis='both', width=1.6, length=4)

for ax in np.ravel(axes):
    ax.tick_params(axis='y', labelsize=22) 
   


sns.despine()
plt.tight_layout()
plt.show()
#plt.savefig("atmifig1.svg", dpi=300, bbox_inches='tight')  

#%%

