import io
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
# import statsmodels.api as sm
import csv
import pycircos
from scipy.stats import poisson
import itertools
import collections
import seaborn as sns

import FACSus as fs

sns.set(style="whitegrid")

def encode_sex(s):
    if s == 'Female':
        return 2
    elif s == 'Male':
        return 1
    else:
        return 0
pca = pd.read_csv('../boson_vcf/archived/covar.txt', sep = ' ', 
                  header = None, names = ['Patient ID','ID']+['pca'+str(i) for i in range(1,11)])
df = pd.merge(pca, tmp, on = 'Patient ID')
covar = pd.read_csv('../boson_metadata/Boson_data.csv', sep = ',', 
                    usecols = ['patient_id', 'gender', 'age'], dtype = {'patient_id': 'str'})
covar['patient_id'] = 'BOSON' + covar['patient_id']
covar = covar.rename(columns = {'patient_id': 'Patient ID'})
df = pd.merge(covar, df, on = 'Patient ID').dropna()
df['ID'] = df['Patient ID']
df['gender'] = df['gender'].apply(encode_sex)
covariates = df[['Patient ID', 'ID'] + ['pca'+str(i) for i in range(1,11)]].drop_duplicates()
covariates.to_csv('../boson_vcf/covar.factor.txt', sep = ' ', index = False, header = False)