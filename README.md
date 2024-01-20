# BOSON_FACS
Sus' pipeline to process and analyse FACS data in a GWAS approach. The pipeline for now includes:
* Preparing and formatting FACS data
* GWAS
* Post-GWAS analyses

To Do:
* Please note that this project is still in progress and currently held, although updates and analyses may be added time by time.

Inputs:
* A `config.json` file to specify sample names, etc. You should modify this file which is under directory `pipelines`
* A `markers.tsv` file to save all markers
* `host` and/or `virus` genomic files in plink format
* `FACS` inputs in `.csv` format (see pipeline) - For my case I have FACS data collected at different timepoints, so analyses for multiple FACS profiles are made possible. However, one can simply hack this pipeline as a normal GWAS pipeline

Explanation of Entries in the Config File:
* `markers`: A `.tsv` file to store all marker names. Note that in all pipelines I have developped, `.tsv` files will always be a single column file without index and header, and they are always used as snakemake wildcards
* `num_top`: Number of top signals to extract from each GWAS run
* `check_range`: Check if the unincluded signals are within this range, if so include it

Outputs:
* `preprocess.smk`:
    * Processed FACS files and basic description of the data
* `gwas.smk`:
    * GWAS runs, QQ/Manhattan plots, and genomic inflation factor calculation
* `post_gwas.smk`:
    * Extracting top signals for further look-up
    * Post_GWAS analyses
* `dump.smk`:
    * Currently unused rules. These rules should NEVER be included in `master.smk` and only used to store potentially useful codes, as dependencies and names can be awful or even detrimental to the main pipeline.
* `test.smk`:
    * For the purpose of short testing
* `auxiliary.smk`:
    * Auxiliary `Python` snakemake functions

Run the Pipeline:
* For now, the whole pipeline is separated into different snakemake files that groups a bunch of jobs together. To run a specific file, use, for example, `snakemake -s pipelines/master.smk -c 1 gwas_all` (needless to say, don't forget to dry-run first by `-nr`).
* Alternatively, a `submit_snakemake.sh` submission script is provided. This file takes multiple parameters, with the first to be number of cores required, and the second to be the name of the rules one'd like to run. For example, you can submit by `./submit_snakemake.sh 8 gwas_all`. Additional flags can be passed by the third parameter, like `-nr`.

Stay Tuned :)


