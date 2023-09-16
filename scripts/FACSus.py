import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import scipy as sp
import csv
import io 
import os
from scipy.stats import poisson
import itertools
import collections
import seaborn as sns
from sklearn import preprocessing

### filters ###

# Exclude single panel sum outliers
def sanity_check_filter(df, cols, newname, lower_bound = 95, upper_bound = 105, plot = False, save = False, save_path = None):
    df[newname] = df[cols].sum(axis=1)
    n = df.shape[0]
    if plot:
        plot_hist(df, newname, 'Percentage', save = save, save_path = save_path)
    # Filtering
    df = df[(df[newname] > lower_bound) & (df[newname] < upper_bound)]
    df = df.drop(columns = [newname])
#     print('Number of samples removed for filter', newname, ':', n-df.shape[0])
    return df
# Exclude cross panel std dev outliers 
def sanity_check_filter2(df, cols, newname, plot = False, save = False, save_path = None):
    n = df.shape[0]
    df[newname] = df[cols].std(axis=1)
    df['normalised']=(df[newname]-df[newname].mean())/df[newname].std()
    if plot:
        plot_hist(df, 'normalised', 'Std Dev (Normalised)', save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df[newname] = df[cols].mean(axis=1)
#     df = df.drop(columns = cols)
    df = df.drop(columns = ['normalised'])
#     print('Number of samples removed for filter', newname, ':', n-df.shape[0])
    return df
# Exclude cross panel abs diff outliers
def sanity_check_filter3(df, colname, newname, plot = False, save = False, save_path = None):
    n = df.shape[0]
    df[newname] = (df[colname + '_panel1'] - df[colname + '_panel2']).abs()
    df['normalised']=(df[newname]-df[newname].mean())/df[newname].std()
    if plot:
        plot_hist(df, 'normalised', 'Abs Diff (Normalised)', save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df[newname] = (df[colname + '_panel1'] + df[colname + '_panel2'])/2
    df = df.drop(columns = ['normalised'])
    df = df.drop(columns = [colname + '_panel1', colname + '_panel2'])
#     print('Number of samples removed for filter', newname, ':', n-df.shape[0])
    return df
# Exclude single panel single marker outliers
def sanity_check_filter4(df, col, lower_bound = 3, upper_bound = -3, plot = False, save = False, save_path = None):
    n = df.shape[0]
    df['normalised']=(df[col]-df[col].mean())/df[col].std()
    if plot:
        plot_strip(df, col, save = save, save_path = save_path)
    # Filtering
    df = df[(df['normalised'] > -3) & (df['normalised'] < 3)]
    df = df.drop(columns = ['normalised'])
#     print('Number of samples removed for filter', col, ':', n-df.shape[0])
    return df

### end of filters ###
### plot ###

def plot_strip(df, column, xlabel = None, title = None, upper = 0.75, lower = 0.25, save = False, save_path = None):
    data = df[column]
    upper_quartile = data.quantile(upper)
    lower_quartile = data.quantile(lower)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    sns.stripplot(data=data, jitter = 0.01, ax = ax, alpha = 0.6, zorder = 1)
    ax.hlines(upper_quartile, xmin=-0.05, xmax=0.05, colors="black", linestyles="-", linewidth=2, alpha = 0)
    ax.hlines(lower_quartile, xmin=-0.05, xmax=0.05, colors="black", linestyles="-", linewidth=2, alpha = 0)
    ax.hlines(upper_quartile, xmin=-0.005, xmax=0.005, colors="black", linestyles="-", linewidth=2)
    ax.hlines(lower_quartile, xmin=-0.005, xmax=0.005, colors="black", linestyles="-", linewidth=2)
    ax.vlines(x = 0, ymin = upper_quartile, ymax = lower_quartile, colors = "black", linestyles = "-", linewidth = 2)
    ax.scatter(0, data.mean(), marker = 'x', alpha = 1, s = 40, c = 'black')
    
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_ylabel('Percentage (%)')
    if save:
        plt.savefig(save_path, bbox_inches = "tight", dpi=300)
def plot_hist(df, column, xlabel = None, bins = 20, label = True, save = False, save_path = None):
    counts, edges, bars = plt.hist(df[column], alpha=0.75, bins = bins, histtype='bar', ec='black')
    if label:
        plt.bar_label(bars)
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title(column)
    if save:
        plt.savefig(save_path, bbox_inches = "tight", dpi=300)
    plt.show()
def plot_box(df, column):
    counts, edges, bars = plt.hist(df[column], alpha=0.75, histtype='bar', ec='black')
    plt.bar_label(bars)
    plt.xlabel('Percentage')
    plt.ylabel('Counts')
    plt.title(column)

### end of plots ###
### auxiliaries ###

def format_float(f, dec = 4):
    return f'{f:.{dec}f}'
def quantile_normalisation(df, col):
    return preprocessing.quantile_transform(df[col].values.reshape(-1, 1), output_distribution = 'normal', n_quantiles = df.shape[0])
    
### end of auxiliaries ###

### phe ###

def get_phenotype(df, col, quantile_transformation = True, convert_decimal = False, trailing_str = 'BOSON', save = False, save_path = None):
    df['Patient ID'] = trailing_str + df['Patient ID']
    df['IID'] = df['Patient ID']
    if quantile_transformation:
        df[col] = quantile_normalisation(df, col)
    if convert_decimal:
        df[col] = df[col]/100
    phe = df.loc[:, ['Patient ID', 'IID', col]]
    phe[col] = phe[col].apply(format_float)
    if save:
        phe.to_csv(save_path, header = None, index = None, sep = ' ')
    return None

### end of phe ###
### parse plink GWASs results ###

def parse_plink(file, save = False, save_path = None):
    gwas = []
    for line in open(file): 
        lst = line.split('\n')[0].strip().split(' ')
        lstnew = [i for i in lst if i != '']
        gwas.append(lstnew)
    gwas.pop(0)
    gwas = pd.DataFrame(gwas, columns=['chr', 'ID', 'pos', 'alt', 'add', 'sample_size', 'beta', 'SE', 'p'])
    if save:
        gwas.to_csv(save_path, index = None, sep = ' ')
    return gwas