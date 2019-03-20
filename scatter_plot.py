from __future__ import division
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import seaborn as sns
sns.set(color_codes=True)


# IC50 data of each cell line (Table S2)
public_pharma = 'public_pharmacology_IC50.csv'

# name of file containing output data from segmentation (Table S3)
cell_lines = 'Exp10_Wrapper_Mapk_125_ExportedSubpopDataALL copy.csv'

# unique drug ID
drug_ids = 'drug_ids.csv'

# minimum and maximum IC50 values of each drug
min_max = '000_min_max_conc.csv'

# additional drug info including target and market name
drug_info = 'drug_info.csv'

# read csv files of the above
df1 = pd.read_csv(cell_lines, index_col = 0, error_bad_lines=False)
df2 = pd.read_csv(public_pharma, index_col = 0, error_bad_lines=False)
df3 = pd.read_csv(drug_ids, index_col = 0, error_bad_lines=False)
df4 = pd.read_csv(min_max, index_col = 0, error_bad_lines=False)
df5 = pd.read_csv(drug_info, index_col = 0, error_bad_lines=False)

drug1 = 'SB590885'
drug2 = 'CI-1040'
unit = 'log10(micromolar)'
sensitive = 'sensitive'
resistant = 'resistant'

# target of each drug
drug1_target = df5['TARGET'].loc[92]
drug2_target = df5['TARGET'].loc[88]

# biomarker file after significance testing
biomarker_file = drug1 + '_' + drug2 + '.csv'
names = ["Subpop.", "Number of cells in Subpop.", "Feature", "Percent in Subpop.", "Number of cells with Feature", "P-value", "Adjusted P-value"]
df6 = pd.read_csv(biomarker_file, error_bad_lines=False)

# biomarkers needed to be annotated based on significance
anon_feat1 = df6['Feature'].loc[0]
anon_feat2 = df6['Feature'].loc[1]

# IC50 values of all cell lines
all_cells1 = df2[drug1].dropna()
all_cells2 = df2[drug2].dropna()

# 20th percentile cutoff lines
drug1_20th_percentile = np.percentile(all_cells1, 20)
drug2_20th_percentile = np.percentile(all_cells2, 20)

# minimum and maximum IC50 values for each drug
drug1_min = str(df4['min_conc_log10_microMol'].loc[246])
drug1_max = str(df4['max_conc_log10_microMol'].loc[246])
drug2_min = str(df4['min_conc_log10_microMol'].loc[196])
drug2_max = str(df4['max_conc_log10_microMol'].loc[196])


final_ind1_ic50 = []
final_ind2_ic50 = []
pltaverage1_ic50 = []
pltaverage2_ic50 = []
pltspop_num = []
cbtick_label = []

fig = plt.figure(figsize=(8, 8))
plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
title = 'POET IC50 Data Visualisation'

for column in df1:
    cells = df1[column].dropna()
    ind1_ic50 = df2.loc[cells, drug1].dropna()
    ind2_ic50 = df2.loc[cells, drug2].dropna()

    #
    query1 = (df2.loc[cells, drug1].dropna().sum())
    query2 = (df2.loc[cells, drug2].dropna().sum())

    # number of cell in subpopulation
    spop_num = (float(df1[column].loc[2]))
    pltspop_num.append(spop_num)

    # label of each subpopulation
    label = df1[column].loc[1]
    cbtick_label.append(label)

    average1_ic50 = query1/spop_num
    pltaverage1_ic50.append(average1_ic50)

    average2_ic50 = query2/spop_num
    pltaverage2_ic50.append(average2_ic50)

    final_ind1_ic50.append(ind1_ic50)
    final_ind2_ic50.append(ind2_ic50)

    # annotate biomarker on scatter plot
    plt.annotate(anon_feat1, xy=(0.6, -0.3), xytext=(3, 1.5),
        fontsize=8,
        xycoords='data',
        textcoords='offset points')

    plt.annotate(anon_feat2, xy=(1.6, 0.1), xytext=(3, 1.5),
        fontsize=8,
        xycoords='data',
        textcoords='offset points')

    # annotate subpopulation number
    # plt.annotate(label, (average1_ic50 + 0.05, average2_ic50 - 0.05),
    #     fontsize=8,
    #     xycoords='data',
    #     textcoords='offset points')

# final individual IC50 values of cell lines in each subpopulation
final_ind1_ic50 = pd.concat(final_ind1_ic50)
final_ind2_ic50 = pd.concat(final_ind2_ic50)

def plot_hist_scatter():

    # set up colors
    label = df1.loc[1]
    cb_labels = map(int, label)

    cmap = plt.cm.get_cmap('jet',max(cb_labels)-min(cb_labels)+1)
    bounds = range(min(cb_labels),max(cb_labels)+2)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # obtain and plot each individual cell lines in scatter plot
    # cell1 = df1.cell1.dropna()
    # cell2 = df1.cell2.dropna()
    # cell3 = df1.cell3.dropna()
    # cell4 = df1.cell4.dropna()
    # cell5 = df1.cell5.dropna()
    # cell6 = df1.cell6.dropna()
    # cell7 = df1.cell7.dropna()
    #
    # query1_cell1 = df2.loc[cell1, drug1].dropna()
    # query2_cell1 = df2.loc[cell1, drug2].dropna()
    #
    # query1_cell2 = df2.loc[cell2, drug1].dropna()
    # query2_cell2 = df2.loc[cell2, drug2].dropna()
    #
    # query1_cell3 = df2.loc[cell3, drug1].dropna()
    # query2_cell3 = df2.loc[cell3, drug2].dropna()
    #
    # query1_cell4 = df2.loc[cell4, drug1].dropna()
    # query2_cell4 = df2.loc[cell4, drug2].dropna()
    #
    # query1_cell5 = df2.loc[cell5, drug1].dropna()
    # query2_cell5 = df2.loc[cell5, drug2].dropna()
    #
    # query1_cell6 = df2.loc[cell6, drug1].dropna()
    # query2_cell6 = df2.loc[cell6, drug2].dropna()
    #
    # query1_cell7 = df2.loc[cell7, drug1].dropna()
    # query2_cell7 = df2.loc[cell7, drug2].dropna()

    # plt.scatter(query1_cell1, query2_cell1, c='#f58200', cmap=cmap, norm=norm, alpha=0.5, edgecolor='none')
    # plt.scatter(query1_cell2, query2_cell2, c='#939393', cmap=cmap, norm=norm, alpha=0.5, edgecolor='none')
    # plt.scatter(query1_cell3, query2_cell3, c='#339d00', cmap=cmap, norm=norm, alpha=0.5, edgecolor='none')
    # plt.scatter(query1_cell4, query2_cell4, c='#0039c7', cmap=cmap, norm=norm, alpha=0.5, edgecolor='none')
    # plt.scatter(query1_cell5, query2_cell5, c='#0039c7', cmap=cmap, norm=norm, alpha=0.6, edgecolor='none')
    # plt.scatter(query1_cell6, query2_cell6, c='#339d00', cmap=cmap, norm=norm, alpha=0.5, edgecolor='none')
    # plt.scatter(query1_cell7, query2_cell7, c='#339d00', cmap=cmap, norm=norm, alpha=0.5, edgecolor='none')

    # plot average IC50 value of subpopulation
    plt.scatter(pltaverage1_ic50, pltaverage2_ic50, c=cb_labels, s=pltspop_num, alpha=0.7, cmap=cmap, norm=norm)
    plt.xlim(-2, 3)
    plt.ylim(-2, 3)
    plt.yticks(rotation=90)
    plt.axvline(drug1_20th_percentile, color='#212121', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.axhline(drug2_20th_percentile, color='#212121', linestyle='--', linewidth=0.5, alpha=0.5)

    plt.xlabel(drug1 + ' ' + '[' + drug1_target + ']' + ' ' + '-' + ' ' + unit)
    plt.ylabel(drug2 + ' ' + '[' + drug2_target + ']' + ' ' + '-' + ' ' + unit)
    plt.annotate(sensitive, xy=(-1.5, -2.8),
            color='red',
            annotation_clip=False)
    plt.annotate(resistant, xy=(1.5, -2.8),
            color='blue',
            annotation_clip=False)
    plt.annotate(sensitive, xy=(-2.9, -0.5),
            rotation=90,
            color='red',
            annotation_clip=False)
    plt.annotate(resistant, xy=(-2.9, 1.5),
            rotation=90,
            va='center',
            color='blue',
            annotation_clip=False)
    plt.tick_params(
        axis='both',
        which='both',
        left='on',
        right='off',
        top='off',
        bottom='on',
        length=0,
        size=0)

    # plot histograms of individual cell lines
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=1)
    plt.title(title)
    plt.yticks(rotation=90)
    plt.xlim(-2, 3)
    plt.axvline(z, color='#212121', linestyle='--', linewidth=0.5)
    plt.xlabel('')
    plt.ylabel('Number of cell lines')
    plt.annotate('Min conc.',
        xy=(drug1_min, 0),
        xytext=(drug1_min, 30),
        fontsize=8, textcoords='data',
        ha='center',
        arrowprops=dict(arrowstyle="-"),
        )
    plt.annotate('Max conc.',
        xy=(drug1_max, 0),
        xytext=(drug1_max, 30),
        fontsize=8, textcoords='data',
        ha='center',
        arrowprops=dict(arrowstyle="-"),
        )
    ax1.tick_params(
        axis='both',
        which='both',
        left='on',
        right='off',
        top='off',
        bottom='on',
        length=0,
        size=0)

    plt.hist(all_cells1, color='b', alpha=0.5, bins=20)
    plt.hist(final_ind1_ic50, color='g', alpha=0.5, bins=20)

    ax2 = plt.subplot2grid((3, 3), (1, 2), colspan=1, rowspan=2)
    plt.ylim(-2, 3)
    plt.yticks(rotation=90)
    plt.axhline(a, color='#212121', linestyle='--', linewidth=0.5)
    plt.ylabel('')
    plt.xlabel('Number of cell lines')
    plt.annotate('Min conc.',
        xy=(0, drug2_min),
        xytext=(30, drug2_min),
        fontsize=8, textcoords='data',
        va='center',
        arrowprops=dict(arrowstyle="-"),
        )
    plt.annotate('Max conc.',
        xy=(0, drug2_max),
        xytext=(30, drug2_max),
        fontsize=8, textcoords='data',
        va='center',
        arrowprops=dict(arrowstyle="-"),
        )
    ax2.tick_params(
        axis='both',
        which='both',
        left='on',
        right='off',
        top='off',
        bottom='on',
        length=0,
        size=0)

    plt.hist(all_cells2, alpha=0.5, color='b', orientation="horizontal", bins=20)
    plt.hist(final_ind2_ic50, alpha=0.5, color='g', orientation="horizontal", bins=20)

    # plot legend
    axes = plt.subplot2grid((3, 3), (0, 2), axisbg='white')
    axes.tick_params(axis='x', colors='white')
    axes.tick_params(axis='y', colors='white')
    divider = make_axes_locatable(axes)

    size = [10, 50, 100]
    x1 = [-1, -1, -1]
    y1 = [0.9, 1.1, 1.3]
    cells10 = '10 cells'
    cells50 = '50 cells'
    cells100 = '100 cells'

    plt.scatter(x1, y1, s=size, alpha=0.5)

    plt.annotate(cells10, xy=(-1, 0.9), xytext=(4, 1.5),
        fontsize=8,
        xycoords='data',
        textcoords='offset points')

    plt.annotate(cells50, xy=(-1, 1.1), xytext=(4, 1.5),
        fontsize=8,
        xycoords='data',
        textcoords='offset points')

    plt.annotate(cells100, xy=(-1, 1.3), xytext=(4, 1.5),
        fontsize=8,
        xycoords='data',
        textcoords='offset points')

    leg_title = 'Legend: Number of cell lines in subpopulation'
    plt.title(leg_title, fontsize=8, ha='center')

    plt.show()

plot_hist_scatter()
