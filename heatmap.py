from __future__ import division
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd

# import list of drug pairs and matrix to visualize after classification of pair-wise drug responses
matrix = "MAPK_AKT_drug-drug_matrix.csv"

df1 = pd.read_csv(matrix, index_col = 0, error_bad_lines=False)

# convert to data frame
updated_matrix = pd.DataFrame(df1)

# count the types of categories
val_counts = df1.apply(pd.value_counts)
val_counts['Count of cats'] = val_counts.sum(axis=1)

# convert and print to json for supplemental website S1
json_matrix = updated_matrix.to_json(orient='split')
print json_matrix

# define colour map
cmap = colors.ListedColormap(['#444444', '#f96881', '#fdff26','#7fbf7f','#0066ff'])

# plot bar chart
bar_count = pd.DataFrame(val_counts['Count of cats'])
print bar_count
bar_count.plot.barh(color=['#444444', '#f96881', '#fdff26','#7fbf7f','#0066ff'], alpha=0.8, width=0.5)

# initialize heatmap
fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.4,left=0.35)
ax.set_axis_bgcolor("#444444")
heatmap = ax.imshow(updated_matrix, cmap=cmap, interpolation='nearest')
x_labels = pd.Series(updated_matrix.columns)
y_labels = pd.Series(updated_matrix.index)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(updated_matrix.shape[1]) - 0.05, minor=False)
ax.set_yticks(np.arange(updated_matrix.shape[0]) - 0.05, minor=False)
ax.tick_params(
    axis='both',
    which='both',
    left='off',
    right='off',
    top='off',
    bottom='off')

# set x and y tick labels
ax.set_xticklabels(x_labels[0:], rotation=90, fontsize=8)
ax.set_yticklabels(y_labels[0:], fontsize=8)

# set title and label parameters
plt.title('Different responses between MAPK and AKT pathway inhibitors')
ttl = ax.title
ttl.set_position([.5, 1.05])
plt.ylabel('MAPK pathway signalling')
plt.xlabel('AKT/PI3K pathway signalling')
plt.show()
