import sys
import pandas as pd
import session_info
from pyhere import here
from pathlib import Path
import statsmodels.stats.multitest

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

data_directory = "/users/jsundstr/hicks_home/eQTL/"

run_mode = sys.argv[1]
if run_mode not in ['nominal', 'cis', 'independent', 'interaction']:
    print("'run_mode', the first command-line argument, must be one of 'nominal', 'cis', 'interaction', or 'independent'.")
    sys.exit()
if run_mode == "interaction":
    if len(sys.argv) != 3:
        print("Since 'interaction' mode was selected, exactly one covariate name was expected to be specified at the command line.")
        sys.exit()
    interaction_cov = sys.argv[2]
    out_dir_suffix = f'{run_mode}_{interaction_cov}'
else:
    interaction_cov = "none"
    out_dir_suffix = run_mode

#Specify paths and variables
in_dir = Path(here(data_directory, "tensorQTL_input"))
out_dir = Path(here(data_directory, "tensorQTL_output", out_dir_suffix))
plink_prefix_path = "/dcs05/lieber/liebercentral/libdGenotype_LIBD001/BrainGenotyping/subsets/habenula_new/plink/habenula_maf05"
covariates_file = str(in_dir / "covariates.txt")
expression_bed = str(in_dir / "logcounts.bed.gz")

prefix = "habenula"
add_chr = True

out_dir.mkdir(exist_ok = True)

print(f"Using 'run_mode'={run_mode} and 'interaction_cov'={interaction_cov}.")

# load phenotypes and covariates
print("Reading Expression file: " + expression_bed)
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
print("Phenotype dimensions:")
print(phenotype_df.shape)

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

print("." *20 )
# PLINK reader for genotypes
print("Reading Plink files: " + plink_prefix_path)
pr = genotypeio.PlinkReader(plink_prefix_path)
print("Loading Genotypes...", end='')
genotype_df = pr.load_genotypes().rename({'Br0983': 'Br983'}, axis = 1)
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
print("Genotype dimensions:", end='')
print(genotype_df.shape)