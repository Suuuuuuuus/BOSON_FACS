#!/bin/bash
#SBATCH -J boson_gwas
#SBATCH -A ansari.prj
#SBATCH -p short
#SBATCH -o output.out
#SBATCH -D /well/ansari/users/gjx698/boson/boson_vcf/
#SBATCH -c 4
#SBATCH --array=1-6

module purge
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate sus

marker=$(cat PCs.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
covar="covar.all.txt"
indir="results/linear_allCovariates_PCs/"
plotdir="../graphs/facs-host/baseline/linear_allCovariates_PCs/"
plot="TRUE"

mkdir -p $indir
mkdir -p $plotdir

touch genomic_inflation_factor_all.txt

plink \
--file boson.renamed \
--pheno phenotypes/"$marker".txt \
--allow-no-sex \
--linear hide-covar \
--covar $covar \
--out "$indir"/"$marker"

python ../scripts/parse_plink.py $marker $indir

Rscript ../scripts/gwas.R $marker $plot $indir $plotdir

