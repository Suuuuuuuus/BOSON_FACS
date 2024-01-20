configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import FACSus as fs

markers = read_tsv_as_lst(config["markers"])
panel_number = [i for i in range(1,6)]
panel_time = ['baseline', 'pt']
subjects = ['host'] # subject = ['host', 'virus']

# Using numbers as a cut-off is not ideal, at least need to add an extra condition to iteratively check if the following signal (like 21st or 22nd) is associated or not
rule extract_top_signals:
    input:
        plink_parsed = "results/{subject}/{time}/{marker}/plink_reports/{marker}.{time}.gwas.txt"
    output:
        top_signal = temp("results/{subject}/{time}/{marker}/plink_reports/tmp/{marker}.{time}.top_signal.txt")
    params:
        num_top = config["num_top"],
        check_range = config["check_range"]
    run: 
        df = pd.read_csv(input.plink_parsed, sep = ' ').dropna()
        df['marker'] = wildcards.marker
        top_signals = df.head(params.num_top).sort_values(by = 'p', ascending = True)
        top_signals.to_csv(output.top_signal, sep = ',', index = None)

rule aggregate_top_signals:
    input:
        top_signal = expand("results/{subject}/{time}/{marker}/plink_reports/tmp/{marker}.{time}.top_signal.txt", marker = markers)
    output:
        top_signals = "results/{subject}/{time}/top_signals.csv"
    shell: """
        cat {input.top_signal} >> {output.top_signals}
    """

rule locuszoom:
    input:
        plink_original = rules.gwas.output.plink_result
    output:
        plink_parsed = "results/{subject}/{time}/{marker}/plink_reports/{marker}.{time}.gwas.txt"
    run: 
        """
            Need to modify this as it:
                1. takes different wildcards
                2. didn't save fig at this time
                3. code is not optimised
                4. may need some pre-locuszoom step to filter out the original gwas results
        """
        df = pd.read_csv('boson_vcf/results/all/CD3+iNKT+.baseline.gwas.txt', sep = ' ').dropna()
        df['p'] = df['p'].astype(float)
        df['REF'] = df['alt']
        df = df.rename(columns = {'alt': 'ALT', 'p': 'P', 'pos': 'POS', 'chr': 'CHROM'})
        df = df.sort_values(by = 'P', ascending = True)

        chr = str(df.iloc[0, 0])
        pos = int(df.iloc[0, 2])
        df.iloc[0,3] = 'G'
        window = 150000
        LocusZoom(
            df[['CHROM', 'POS', 'REF', 'ALT', 'P']],
            chrom=chr,   
            start=pos - window,
            end=pos + window,
            build='GRCh38'
            # width=400
        )

        chr = str(df.iloc[1, 0])
        pos = int(df.iloc[1, 2])
        df.iloc[1,3] = 'T'
        window = 100000
        LocusZoom(
            df[['CHROM', 'POS', 'REF', 'ALT', 'P']],
            chrom=chr,   
            start=pos - window,
            end=pos + window,
            build='GRCh38'
        )