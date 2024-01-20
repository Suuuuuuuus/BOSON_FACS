configfile: "pipelines/config.json"
include: "auxiliary.smk"

include: "preprocess.smk"
# include: "gwas.smk"
# include: "post_gwas.smk"

# include: "test.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
import collections
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")
sys.path.append("scripts")
import FACSus as fs

markers = read_tsv_as_lst(config["markers"])
panel_number = [i for i in range(1,6)]
panel_time = ['baseline', 'pt']
model = ["additive", "recessive", "dominant"]
subjects = ['host'] # subject = ['host', 'virus']

rule preprocess_all:
    input:
        facs_merged = expand("data/facs_in_use/{time}/merged.csv", time = panel_time),
        plink_phenotype = expand("results/{subject}/{time}/phenotypes/{marker}.txt", time = panel_time, marker = markers, subject = subjects), 
        marker_distribution = expand("results/{subject}/{time}/{marker}/strip.{marker}.png", time = panel_time, marker = markers, subject = subjects)

rule gwas_all:
    input:
        plink_result = expand("results/{subject}/{time}/{marker}/plink_reports/{marker}.assoc.linear", subject = subjects, time = panel_time, marker = markers),
        plink_parsed = expand("results/{subject}/{time}/{marker}/plink_reports/{marker}.{time}.gwas.txt", subject = subjects, time = panel_time, marker = markers)

rule post_gwas_all:
    input:
        top_signals = expand("results/{subject}/{time}/top_signals.csv", subject = subjects, time = panel_time)

rule test_all:
    input:
        # vcf = "results/imputation/tmp/res.txt"

# rule all:
#     input:
#         lcwgs_wrap_up = "results/lcwgs_results.csv"

