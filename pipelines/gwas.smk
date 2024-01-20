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
subjects = ['host'] # subject = ['host', 'virus']
model = ["additive", "recessive", "dominant"] # Modify rule gwas such that plink can run different models

rule gwas:
    input:
        plink_phenotype = "results/{subject}/{time}/phenotypes/{marker}.txt",
        bed = "data/vcf/{subject}.bed", # Modify plink file name (on cluster) 
        bim = "data/vcf/{subject}.bim",
        fam = "data/vcf/{subject}.fam",
        covar = "data/auxiliaries/covar.all.txt" # Need to modify this later such that different gwas runs can use different covar files
    output:
        plink_result = "results/{subject}/{time}/{marker}/plink_reports/{marker}.assoc.linear"
    shell: """
        plink # Go to cluster see gwas.sh
    """

rule parse_gwas_result:
    input:
        plink_original = rules.gwas.output.plink_result
    output:
        plink_parsed = "results/{subject}/{time}/{marker}/plink_reports/{marker}.{time}.gwas.txt"
    run: 
        fs.parse_plink(input.plink_original, True, output.plink_parsed)

"""
1. See how to run R in smk
2. Write gwas.R and dag.R into pipeline
"""

rule plot_QQ:

rule plot_Manhattan:

rule calculate_gif: