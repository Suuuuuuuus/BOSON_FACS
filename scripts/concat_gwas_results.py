import sys
import numpy as np
import FACSus as fs
import pandas as pd

indir = '/well/ansari/users/gjx698/boson/boson_vcf/results/linear_allCovariates/'
head = 20

def get_file(marker, indir = indir):
    df = pd.read_csv(indir + marker + '.baseline.gwas.txt', sep = ' ').dropna()
    df['marker'] = marker
    return df

markers = pd.read_csv('markers.txt', header = None, names = ['markers'])
markers = list(markers['markers'].values)

for marker in markers:
    if 'df_all' in globals():
        tmp = get_file(marker)
        df_all = pd.concat([df_all, tmp])
        top_signals = pd.concat([top_signals, tmp.head(head)])
    else:
        df_all = get_file(marker)
        top_signals = df_all.head(head)

df_all.to_csv(indir + 'all_markers.csv', sep = ',', index = None)
top_signals = top_signals.sort_values(by = 'p', ascending = True)
top_signals.to_csv(indir + 'all_markers_top_signals.csv', sep = ',', index = None)

