library(here)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(snpStats)
library(bigsnpr)
library(sessioninfo)
library(cowplot)
library(data.table)
library(jaffelab)
library(getopt)
library(MRutils)
library(Polychrome)
data(palette36)

#   Read in which tensorQTL run mode is being used
spec <- matrix(
    c("mode", "m", 1, "character", "tensorQTL run mode"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

accepted_modes = c('nominal', 'cis', 'independent')
if (!(opt$mode %in% accepted_modes)) {
    stop(
        sprintf(
            "'opt$mode' must be in '%s'.",
            paste(accepted_modes, collapse = "', '")
        )
    )
}

eqtl_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', opt$mode, 'FDR05.csv'
)
deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_narrow_path = here(
    'processed-data', '17_eQTL', 'trubetskoy_gwas_supp_tab1.xls'
)
gwas_narrow_filt_path = here(
    'processed-data', '17_eQTL', 'gwas_narrow_filtered.csv'
)
gwas_wide_path = here(
    "processed-data", "13_MAGMA","GWAS", "scz2022",
    "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
)
gwas_wide_filt_path = here(
    'processed-data', '17_eQTL', 'gwas_wide_filtered.csv'
)
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
gene_pcs_path = here('processed-data', '03_bulk_pca', 'PCs.rds')
snp_pcs_path = here(
    'processed-data', '08_bulk_snpPC', 'v3',
    'habenula_R.9_MAF.05.RSann_filt.snpPCs.tab'
)
plink_path = here(
   'processed-data', '08_bulk_snpPC', 'v3', 'habenula_R.9_MAF.05.RSann.bed'
)
raw_geno_path = here(
    'processed-data', '08_bulk_snpPC', 'v3',
    'habenula_R.9_MAF.05.RSann_filt.traw'
)
paired_variants_path = here(
    "processed-data", "17_eQTL", "DEA_paired_variants.txt"
)
rs_path = here(
    "processed-data", "17_eQTL", "rsID_independent_deg_or_gwas_wide.csv"
)
plot_dir = here('plots', '17_eQTL', opt$mode)

sig_cutoff_deg_explore = c(0.1, 0.05)