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
plink_prefix_path = "/dcs04/hicks/data/data_proc/bulkDNAgenotype_habenula_human_pilot/habenula_maf05"
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

genotype_ids = set(genotype_df.columns)
covariates_ids = set(covariates_df.index)

#   Exclude any samples that have genotypes but not covariates; halt in the
#   reverse situation
if (genotype_ids > covariates_ids):
    extra_donors = list(genotype_ids - covariates_ids)
    print("Excluding extra genotyped donors:", ','.join(extra_donors))
    genotype_df.drop(extra_donors, axis = 1, inplace = True)
elif (genotype_ids != covariates_ids):
    print("Not all covariate samples in genotyped samples.")
    sys.exit()

## Fix chr names (maybe fix in plink?)
if add_chr:
    print("Adding 'chr' to genotype positions")
    variant_df.chrom = [s.split(':')[0] for s in list(variant_df.index)]

variant_chrom = set(variant_df.chrom)
express_chrom = set(phenotype_pos_df.chr)

#   Exclude any expressed chromosomes that aren't in genotyping data; halt in
#   the reverse situation
if express_chrom > variant_chrom:
    print("Excluding phenotypes from these chromosomes:")
    print(express_chrom - variant_chrom)
    
    chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
    print("Phenotypes with chr in variants: " )
    print(chrom_filter.value_counts())
    
    phenotype_df = phenotype_df[chrom_filter]
    phenotype_pos_df = phenotype_pos_df[chrom_filter]
elif express_chrom != variant_chrom:
    print("Expression data must contain the same or a superset of chromosomes as the genotyping data.")
    sys.exit()

################################################################################
#   Run tensorQTL
################################################################################

if run_mode == "nominal":
    cis.map_nominal(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        prefix = prefix, covariates_df = covariates_df, maf_threshold = 0.05,
        interaction_df = None, maf_threshold_interaction = 0, group_s = None,
        window = 500000, run_eigenmt = True, output_dir = out_dir,
        write_top = False, verbose = False
    )
elif run_mode == "cis":
    cis_out = cis.map_cis(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        covariates_df = covariates_df, group_s = None, maf_threshold = 0.05,
        beta_approx = True, nperm = 10000, window = 500000,
        random_tiebreak = False, logger = None, seed = 118, verbose = True
    )
    cis_out.to_csv(out_dir / "cis_out.csv")
elif run_mode == "independent":
    cis_out = pd.read_csv(
        out_dir.parent / "cis" / "cis_out.csv", index_col = 0
    )

    #   Compute q value, an expected input column
    cis_out['qval'] = statsmodels.stats.multitest.fdrcorrection(cis_out['pval_beta'])[1]

    ind_out = cis.map_independent(
        genotype_df = genotype_df, variant_df = variant_df, cis_df = cis_out,
        phenotype_df = phenotype_df, phenotype_pos_df = phenotype_pos_df,
        covariates_df = covariates_df, group_s = None, maf_threshold = 0.05,
        nperm = 10000, window = 500000, random_tiebreak = False, logger = None,
        seed = 119, verbose = True
    )

    ind_out.to_csv(out_dir / "independent_out.csv")
else:
#'run_mode' must be 'interaction' based on check at the top of script
    col_data = pd.read_csv(in_dir / "colData.csv", index_col = 'BrNum')
    if interaction_cov not in col_data.columns:
        print(f'Expected {interaction_cov} to be a valid colData column.')
        sys.exit()

    nominal_out = cis.map_nominal(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        prefix = prefix, covariates_df = covariates_df, output_dir = out_dir, 
        interaction_df = col_data[[interaction_cov]],
        maf_threshold_interaction = 0.05, group_s = None, window = 500000,
        run_eigenmt = False, write_stats = True
    )

session_info.show()