#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 16:42:02 2025

@author: joannafernandez
"""

""

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


R0perc = Find("/path/to/tideR0_percentages.csv")
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find("/path/to/tideR0_pvals.csv")
R0pval= R0pval.dropna() #remove empty rows

#where p>0.001 (non-significant indel) remove from consideration
R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)

#formating for plotting
melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['sgLAD_6hr',
                                      'sgLAD_24hr', 'sgLAD_48hr'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)

#%%
#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df['class'] = melted_df['indels'].apply(classify_indel)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)



#%%

hue_order = ['6hr', '24hr', '48hr']

custom_palette = {
    'MMEJ': '#D33873',
    'NHEJ': '#575757',
    'uncut': '#EFB54F'
}

time_order = ['6hr','24hr', '48hr']
condition_order = ['sgLAD']
class_order = ['MMEJ','NHEJ', 'uncut']  


#%%

g = sns.catplot(
    data=subset,
    x="class",
    y="relative_percentage",
    hue="class",
    col="time",
    row="condition",
    kind="bar",
    estimator=np.mean,                 
    errorbar=("se", 1),              
    capsize=0.25,                     
    err_kws=dict(color="black", lw=1.8),  
    edgecolor="black",
    linewidth=1.5,
    width=0.7,
    palette=custom_palette,
    dodge=False,                       
    legend=False,
    alpha=0.9,
    height=4,
    aspect=1.8,
)


for (row_val, col_val), ax in g.axes_dict.items():
    subsety = subset[(subset['condition'] == row_val) & (subset['time'] == col_val)]

    sns.swarmplot(
        data=subsety,
        x="class",
        y="relative_percentage",
        color="white",
        edgecolor="black",
        linewidth=1.5,
        size=10,
        dodge=False,
        ax=ax,
        zorder=4,
        clip_on=False
    )

    ax.axhline(y=50, color='black', linestyle='--', linewidth=1)


for ax in g.axes.flat:
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', which='both', labelleft=True)
    ax.set_yticks([0, 20, 40, 60, 80, 100])

    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='x', which='both', labelbottom=True)

    ax.axhline(50, linestyle='--', color='grey', linewidth=1)

    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(axis='both', width=1.6, length=4)


for ax in g.axes.flat:
    for p in ax.patches:
        p.set_edgecolor("black")
        p.set_linewidth(1.5)
        p.set_joinstyle("miter")

g.set_xlabels(fontsize=10)
g.set_ylabels(fontsize=10)
g.set_titles(size=10)

for ax in g.axes.flat:
    ax.tick_params(axis='y', labelsize=22)

sns.despine()
plt.tight_layout()
plt.show()


#%%
'''
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


x_pos = {'sgLAD': 0}
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
  #  ax.set_xticks([x_pos['DMSO'], x_pos['ATMi']])
    ax.set_xticklabels(['sgLAD'], fontsize=11)
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

'''
#%%


import itertools as it
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats


subsetnhej = subset[(subset['class'] == 'NHEJ')]

# ---- 1) Filter to MMEJ and set time order (adjust the order if you like) ----

# If you want a specific order, set it here; otherwise it will infer from appearance.
time_order = ['6hr','24hr','48hr']  # <- uncomment & edit if needed
subsetnhej['time'] = pd.Categorical(subsetnhej['time'], categories=time_order, ordered=True)


# ---- 2) Ordinary one-way ANOVA: relative_percentage ~ time ----
model = smf.ols('relative_percentage ~ C(time)', data=subsetnhej).fit()
anova_tbl = sm.stats.anova_lm(model, typ=2)
print("\n=== One-way ANOVA (MMEJ % across time) ===")
print(anova_tbl)

# ---- 3) Bonferroni multiple comparisons between all timepoints (GraphPad-like) ----
# Use the ANOVA residual MSE as pooled variance, like GraphPad's "Bonferroni's test"
levels = list(subsetnhej['time'].cat.categories) if pd.api.types.is_categorical_dtype(subsetnhej['time']) \
         else sorted(subsetnhej['time'].unique())

# group means and sizes
group_stats = (
    subsetnhej.groupby('time')['relative_percentage']
    .agg(['mean','count'])
    .reindex(levels)
)

MSE = model.mse_resid
df_resid = int(model.df_resid)
pairs = list(it.combinations(levels, 2))
m = len(pairs)  # number of pairwise tests for Bonferroni

rows = []
for a, b in pairs:
    mean_a, n_a = group_stats.loc[a, 'mean'], group_stats.loc[a, 'count']
    mean_b, n_b = group_stats.loc[b, 'mean'], group_stats.loc[b, 'count']
    diff = mean_b - mean_a  # (b - a), change order if you prefer

    # SE based on pooled (ANOVA) variance:
    SE = np.sqrt(MSE * (1.0/n_b + 1.0/n_a))

    # t statistic, 2-sided p using residual df
    t_stat = diff / SE
    p_raw = 2 * stats.t.sf(np.abs(t_stat), df_resid)

    # Bonferroni-adjusted p-value
    p_bonf = min(p_raw * m, 1.0)

    # Bonferroni-adjusted CI: use alpha/m for each comparison
    alpha = 0.05
    tcrit = stats.t.ppf(1 - (alpha/(2*m)), df_resid)
    ci_low = diff - tcrit * SE
    ci_high = diff + tcrit * SE

    rows.append({
        'comparison': f'{b} vs {a}',
        'mean_b_minus_a': diff,
        'SE_pooled': SE,
        't': t_stat,
        'df_resid': df_resid,
        'p_raw': p_raw,
        'p_bonf': p_bonf,
        'CI95_low_Bonf': ci_low,
        'CI95_high_Bonf': ci_high,
        'n_a': int(n_a),
        'n_b': int(n_b)
    })

posthoc = pd.DataFrame(rows)
print("\n=== Bonferroni multiple comparisons (GraphPad-style; pooled variance) ===")
print(posthoc[['comparison','p_raw','p_bonf','CI95_low_Bonf','CI95_high_Bonf','n_a','n_b']])

