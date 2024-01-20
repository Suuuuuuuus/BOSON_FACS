configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import FACSus as fs

panel_number = [i for i in range(1,6)]
panel_time = ['baseline', 'pt']
subjects = ['host'] #subject = ['host', 'virus']

rule merge_facs_dfs:
    input:
        facs_csvs = expand("data/facs_in_use/{time}/Panel{panel}_Baseline.csv", panel = panel_number)
    output:
        facs_merged = "data/facs_in_use/{time}/merged.csv"
    run: 
        df_ary = []
        for file in input.facs_csvs:
            df = pd.read_csv(file, sep = ',')
            df['Patient_ID'] = df['Patient_ID'].astype(str)
            df_ary.append(df)
        df[0].iloc[:,1:-9] = df[0].iloc[:,1:-9].apply(pd.to_numeric)
        df[0] = df[0].drop(columns = ['Azim_data.SUBJID'])
        for i in range(1, 5):
            df_ary[i].iloc[:, 1:] = df_ary[i].iloc[:, 1:].apply(pd.to_numeric)

        new_dfs = [i.copy() for i in df_ary]
        for i in range(5):
            new_columns = [df_ary[i].columns[0]] + [col + '_panel' + str(i+1) for col in df_ary[i].columns[1:]]
            new_dfs[i].columns = new_columns
        """
           write a command here to remove the postfix '_panel' or modify the previous one if the marker measurement is not repeated
           modify original csv marker names
        """
        tmp = pd.merge(new_dfs[0], new_dfs[1], on = 'Patient_ID')
        tmp = pd.merge(tmp, new_dfs[2], on = 'Patient_ID')
        tmp = pd.merge(tmp, new_dfs[3], on = 'Patient_ID')
        df = pd.merge(tmp, new_dfs[4], on = 'Patient_ID')
        df.to_csv(output.facs_merged, sep = ',', index = False)
        # df['Patient_ID'] = df['Patient_ID'].astype(str)
        # df = df.dropna() 
        # df['Patient ID'] = 'BOSON' + df['Patient ID']

rule get_phenotypes:
    input:
        merged_csv = rules.merge_facs_dfs.output.facs_merged
    output:
        plink_phenotype = "results/{subject}/{time}/phenotypes/{marker}.txt", 
        marker_distribution = "results/{subject}/{time}/{marker}/strip.{marker}.png"
    run: 
        df = pd.read_csv(input.merged_csv, sep = ',') # Need to check 1) dtypes 2) missing values
        if wildcards.marker == 'CD3+':
            tmp = df[['Patient_ID', 'CD3-_panel1'] + [col for col in df.columns if 'CD3+_panel' in col]].dropna()
            CD3 = fs.sanity_check_filter(tmp, ['CD3+_panel1', 'CD3-_panel1'], 'CD3')
            CD3 = fs.sanity_check_filter2(CD3, ['CD3+_panel' + str(i+1) for i in range(5)], 'CD3+')
            fs.get_phe('CD3+', CD3, 'CD3+', plot = True, save_plot = True, save_phe = True, save_path_plot = output.marker_distribution, save_path_phe = output.plink_phenotype)
        elif wildcards.marker == 'CD3+CD4+':
            tmp = df[['Patient_ID'] + ['CD3+CD4+_panel' + str(i+1) for i in range(5)]].dropna()
            CD4 = fs.sanity_check_filter2(tmp, ['CD3+CD4+_panel' + str(i+1) for i in range(5)], 'CD3+CD4+')
            fs.get_phe('CD3+CD4+', CD4, 'CD3+CD4+', plot = True, save_plot = True, save_phe = True, save_path_plot = output.marker_distribution, save_path_phe = output.plink_phenotype)
        elif wildcards.marker == 'CD3+CD8+':
            tmp = df[['Patient_ID'] + ['CD3+CD8+_panel' + str(i+1) for i in range(5)]].dropna()
            CD8 = fs.sanity_check_filter2(tmp, ['CD3+CD8+_panel' + str(i+1) for i in range(5)], 'CD3+CD8+')
            fs.get_phe('CD3+CD8+', CD8, 'CD3+CD8+', plot = True, save_plot = True, save_phe = True, save_path_plot = output.marker_distribution, save_path_phe = output.plink_phenotype)
        elif wildcards.marker == 'CD3+CD8+CD161++':
            tmp = df[['Patient_ID'] + ['CD3+CD8+CD161++_panel' + str(i+1) for i in range(2)]].dropna()
            CD161 = fs.sanity_check_filter3(tmp, 'CD3+CD8+CD161++', 'CD3+CD8+CD161++')
            fs.get_phe('CD3+CD8+CD161++', CD161, 'CD3+CD8+CD161++', plot = True, save_plot = True, save_phe = True, save_path_plot = output.marker_distribution, save_path_phe = output.plink_phenotype)
        elif wildcards.marker == 'CD3+CD8+CD161+Va7.2+':
            tmp = df[['Patient_ID'] + ['CD3+CD8+CD161+Va7.2+_panel' + str(i+1) for i in range(2)]].dropna()
            Va72 = fs.sanity_check_filter3(tmp, 'CD3+CD8+CD161+Va7.2+', 'CD3+CD8+CD161+Va7.2+')
            fs.get_phe('CD3+CD8+CD161+Va7.2+', Va72, 'CD3+CD8+CD161+Va7.2+', plot = True, save_plot = True, save_phe = True, save_path_plot = output.marker_distribution, save_path_phe = output.plink_phenotype)
        else:
            tmp = df[['Patient_ID'] + [wildcards.marker]].dropna()
            fs.get_phe(wildcards.marker, df, wildcards.marker, plot = True, save_plot = True, save_phe = True, save_path_plot = output.marker_distribution, save_path_phe = output.plink_phenotype)
        