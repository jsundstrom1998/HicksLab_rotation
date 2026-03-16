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

#For testing
#opt <- list(mode = "cis")

accepted_modes = c('nominal', 'cis', 'independent')
if (!(opt$mode %in% accepted_modes)) {
    stop(
        sprintf(
            "'opt$mode' must be in '%s'.",
            paste(accepted_modes, collapse = "', '")
        )
    )
}

eqtl_data_directory = "/users/jsundstr/hicks_home/"
lieber_data_directory = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/"
genotype_directory = "/dcs04/hicks/data/data_proc/bulkDNAgenotype_habenula_human_pilot/v3"


eqtl_path = here(
    eqtl_data_directory, 'eQTL',  'tensorQTL_output', opt$mode, 'FDR05.csv'
)
deg_path = here(lieber_data_directory,
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_narrow_path = here(lieber_data_directory,
    'processed-data', '17_eQTL', 'trubetskoy_gwas_supp_tab1.xls'
)
gwas_narrow_filt_path = here(lieber_data_directory,
    'processed-data', '17_eQTL', 'gwas_narrow_filtered.csv'
)
gwas_wide_path = here(lieber_data_directory,
    "processed-data", "13_MAGMA","GWAS", "scz2022",
    "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
)
gwas_wide_filt_path = here(lieber_data_directory,
    'processed-data', '17_eQTL', 'gwas_wide_filtered.csv'
)
rse_path = here(lieber_data_directory,
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
gene_pcs_path = here(lieber_data_directory, 'processed-data', '03_bulk_pca', 'PCs.rds')
snp_pcs_path = here(genotype_directory,
    'habenula_R.9_MAF.05.RSann_filt.snpPCs.tab'
)
plink_path = here(genotype_directory,
 'habenula_R.9_MAF.05.RSann.bed'
)
raw_geno_path = here(genotype_directory,
    'habenula_R.9_MAF.05.RSann_filt.traw'
)
paired_variants_path = here(lieber_data_directory,
    "processed-data", "17_eQTL", "DEA_paired_variants.txt"
)
rs_path = here(eqtl_data_directory, 'eQTL', 
     "rsID_independent_deg_or_gwas_wide.csv"
)
plot_dir = here(eqtl_data_directory, 'eQTL',
    'plots', opt$mode)

sig_cutoff_deg_explore = c(0.1, 0.05)

#   We can use a stricter cutoff for nominal, where we have more results
if (opt$mode == 'nominal') {
    sig_cutoff_deg_plot = 0.05
} else {
    sig_cutoff_deg_plot = 0.1
}

sig_cutoff_gwas = 5e-8

#   Covariates used in the model for tensorQTL and for DEA, respectively
eqtl_covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5',
    paste0('PC', 1:13)
)
deg_covariates = c(
    'PrimaryDx', 'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN',
    'totalAssignedGene', 'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4',
    'qSV5', 'qSV6', 'qSV7', 'qSV8', 'tot.Hb', 'tot.Thal'
)

geno_colors = c(palette36[6], "#DB9813", palette36[8])
names(geno_colors) = c("0", "1", "2")

dx_colors = c("#2a9d8f", "#f77f00")
names(dx_colors) = c("Control", "SCZD")

lift_over_path = system('which liftOver', intern = TRUE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   Return a merged tibble of genotyping data, residualized expression, and
#   colData
merge_exp_df = function(
        rse_gene, mod_deg, mod_eqtl, eqtl, plink, mismatched_snps
    ) {
    #   Grab the expression for all genes paired with significant eQTLs;
    #   convert to long format. Residualize using both the tensorQTL and DEA
    #   models
    express_eqtl = assays(rse_gene[unique(eqtl$phenotype_id),])$logcounts |>
        cleaningY(mod = mod_eqtl, P = 1) |>
        as.data.frame() |>
        rownames_to_column("gene_id") |>
        pivot_longer(
            cols = -gene_id, names_to = "sample_id",
            values_to = "resid_logcount_eqtl"
        )
    express_deg = assays(rse_gene[unique(eqtl$phenotype_id),])$logcounts |>
        cleaningY(mod = mod_deg, P = 2) |>
        as.data.frame() |>
        rownames_to_column("gene_id") |>
        pivot_longer(
            cols = -gene_id, names_to = "sample_id",
            values_to = "resid_logcount_deg"
        )
    express_both = left_join(
        express_eqtl, express_deg, by = c('gene_id', 'sample_id')
    )
    
    geno_df = plink$genotypes[, unique(eqtl$variant_id)] |>
        as.data.frame() |>
        rownames_to_column("sample_id") |>
        pivot_longer(
            cols = -sample_id, names_to = "snp_id", values_to = "genotype"
        ) |>
        mutate(
            sample_id = ifelse(sample_id == "Br0983", "Br983", sample_id),
            genotype = factor(
                #   Ensure 0 is reference, 1 is heterozygous, and 2 is
                #   homozygous for the ALT
                ifelse(
                    snp_id %in% mismatched_snps,
                    as.integer(genotype) - 1,
                    3 - as.integer(genotype)
                )
            )
        )

    #   Join genotyping data, expression, and colData
    exp_df = eqtl |>
        dplyr::select(phenotype_id, variant_id) |>
        dplyr::rename(gene_id = phenotype_id, snp_id = variant_id) |>
        #   Add expression data residualized from both models
        full_join(express_both, by = 'gene_id') |>
        #   Add genotypes for each SNP and sample
        left_join(geno_df, by = c('snp_id', 'sample_id')) |>
        #   Add colData
        left_join(
            colData(rse_gene) |> as_tibble(), by = join_by(sample_id == BrNum)
        ) |>
        #   Add gene symbol
        mutate(
            gene_symbol = rowData(rse_gene)$Symbol[
                match(gene_id, rownames(rse_gene))
            ]
        )

    return(exp_df)
}

#   Return a copy of 'exp_df' with an additional 'rs_id' column
add_rs_id = function(exp_df) {
    a = tibble(snp_id = unique(exp_df$snp_id)) |>
        separate(
            snp_id, into = c("chr", "pos", "A1", "A2"), sep = ":",
            remove = FALSE
        ) |>
        mutate(chr = str_extract(chr, '^chr(.*)', group = 1))
    
    #   Individually search for rs IDs
    a$rs_id = sapply(
        seq(nrow(a)), function(i) {
            get_rsid_from_position(
                chrom = a$chr[i], pos = a$pos[i], ref = a$A1[i], alt = a$A2[i],
                assembly = "hg38"
            )
        }
    )

    return(exp_df |> left_join(a |> select(snp_id, rs_id), by = "snp_id"))
}

plot_triad_exploratory = function(eqtl, exp_df, plot_dir, plot_prefix) {
    #   Plot expression by genotype of each variant with one gene per page and
    #   potentially several variants per gene
    plot_list_geno = list()
    plot_list_dx = list()
    plot_list_fraction = list()
    plot_list_fraction_no_labels = list()
    for (this_gene in unique(exp_df$gene_id)) {
        this_symbol = exp_df$gene_symbol[
            match(this_gene, exp_df$gene_id)
        ]

        label_df = eqtl |>
            filter(
                phenotype_id == this_gene,
                variant_id %in% exp_df$snp_id
            ) |>
            mutate(sig_label = sprintf(" p = %s", signif(pval_nominal, 3))) |>
            dplyr::rename(snp_id = variant_id)

        plot_list_geno[[this_gene]] = exp_df |>
            filter(gene_id == this_gene) |>
            ggplot() +
                geom_boxplot(
                    mapping = aes(
                        x = genotype, y = resid_logcount_eqtl, color = genotype
                    ),
                    outlier.shape = NA) +
                geom_jitter(
                    mapping = aes(
                        x = genotype, y = resid_logcount_eqtl, color = genotype
                    )
                ) +
                geom_text(
                    data = label_df,
                    mapping = aes(label = sig_label, x = -Inf, y = Inf),
                    hjust = 0,
                    vjust = 1
                ) +
                facet_wrap(~ snp_id) +
                scale_color_manual(values = geno_colors) +
                labs(
                    x = "Genotype", y = "Residualized eQTL Expression",
                    title = this_symbol
                ) +
                theme_bw() +
                theme(legend.position = "none")
        
        plot_list_dx[[this_gene]] = exp_df |>
            filter(gene_id == this_gene) |>
            #   Don't take duplicate rows if there are multiple SNPs per gene
            filter(snp_id == snp_id[1]) |>
            ggplot(
                    mapping = aes(
                        x = PrimaryDx, y = resid_logcount_deg, color = PrimaryDx
                    )
                ) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter() +
                scale_color_manual(values = dx_colors) +
                labs(y = "Residualized DEG Expression", title = this_symbol) +
                theme_bw(base_size = 20) +
                theme(legend.position = "none")
        
        #   Grab all SNPs associated with this gene
        these_snp_ids = exp_df |>
            filter(gene_id == this_gene) |>
            pull(snp_id) |>
            unique()

        #   Each page will consist of one SNP-gene pair and two plots: one for
        #   habenula and one for thalamus fraction
        for (this_snp_id in these_snp_ids) {
            temp = list()
            temp_no_labels = list()
            this_title = sprintf('%s: %s', this_symbol, this_snp_id)
            for (x_var_name in c("tot.Hb", "tot.Thal")) {
                temp[[x_var_name]] = exp_df |>
                    filter(gene_id == this_gene, snp_id == this_snp_id) |>
                    ggplot(
                        mapping = aes(
                            x = get({{ x_var_name }}), y = resid_logcount_eqtl,
                            color = genotype
                        )
                    ) +
                    geom_point(size = 3) +
                    geom_smooth(method = lm) +
                    scale_color_manual(values = geno_colors) +
                    coord_cartesian(xlim = c(0, 1)) +
                    theme_bw(base_size = 20) +
                    theme(plot.title = element_text(size = 20)) +
                    labs(y = "Residualized eQTL Expression", color = "Genotype")
                
                #   Only want one title and legend, not 2
                if (x_var_name == "tot.Hb") {
                    temp[[x_var_name]] = temp[[x_var_name]] +
                        labs(x = "Habenula Fraction", title = this_title) +
                        theme(legend.position = "none")
                } else {
                    temp[[x_var_name]] = temp[[x_var_name]] +
                        labs(x = "Thalamus Fraction", title = " ")
                }

                #   Create a version easier to plot as a large grid in
                #   illustrator for some manuscript figures
                temp_no_labels[[x_var_name]] = temp[[x_var_name]] +
                    theme_bw(base_size = 30) +
                    theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = "none",
                        plot.title = element_blank()
                    )
            }
            plot_list_fraction[[this_title]] = plot_grid(
                plotlist = temp, ncol = 2, rel_widths = 4:5
            )
            plot_list_fraction_no_labels[[this_title]] = plot_grid(
                plotlist = temp_no_labels, ncol = 2
            )
        }
    }

    pdf(file.path(plot_dir, sprintf('expr_by_geno_%s.pdf', plot_prefix)))
    print(plot_list_geno)
    dev.off()

    pdf(file.path(plot_dir, sprintf('expr_by_dx_%s.pdf', plot_prefix)))
    print(plot_list_dx)
    dev.off()

    pdf(
        file.path(
            plot_dir, sprintf('expr_by_geno_fraction_%s.pdf', plot_prefix)
        ),
        width = 14, height = 7
    )
    print(plot_list_fraction)
    dev.off()

    pdf(
        file.path(
            plot_dir,
            sprintf('expr_by_geno_fraction_%s_no_labels.pdf', plot_prefix)
        ),
        width = 14, height = 7
    )
    print(plot_list_fraction_no_labels)
    dev.off()
}

#   Residualized expression vs. genotype boxplots faceted by SNP ID
exp_vs_geno_manuscript_plot = function(
        eqtl, exp_df, plot_dir, plot_suffix, note_risk_allele, color_by_dx,
        facet_nrow = 1, pdf_width = 8, pdf_height = 6
    ) {
    a = exp_df |>
        #   Add 'pval_nominal', 'slope' columns from 'eqtl'
        left_join(
            eqtl |>
                dplyr::rename(snp_id = variant_id, gene_id = phenotype_id) |>
                select(snp_id, gene_id, pval_nominal, slope) |>
                mutate(
                    sig_label = sprintf(
                        "beta = %s \np = %s \n",
                        signif(slope, 3),
                        signif(pval_nominal, 3)
                    )
                ),
            by = c("snp_id", "gene_id")
        )

    #   Add informative facet title: include 2 forms of SNP IDs, gene
    #   symbol, and if specified, the risk allele
    if (note_risk_allele) {
        a = a |> mutate(
            anno_label = sprintf(
                "SNP: %s\n%s\n(risk allele: %s)\nGene: %s",
                rs_id, snp_id, risk_allele, gene_symbol
            )
        )
    } else {
        a = a |> mutate(
            anno_label = sprintf(
                "SNP: %s\n%s\nGene: %s",
                rs_id, snp_id, gene_symbol
            )
        )
    }
    
    #   Avoid duplicate labels
    label_df = a |>
        group_by(anno_label) |>
        slice_head(n = 1) |>
        ungroup()

    p = ggplot(a) +
        facet_wrap(~ anno_label, nrow = facet_nrow) +
        labs(x = "Genotype", y = "Residualized eQTL Expression") +
        theme_bw(base_size = 20) +
        theme(
            legend.position = "none", strip.text.x = element_text(size = 13)
        )
    
    #   Either color by genotype or diagnosis
    if (color_by_dx) {
        p = p +
            geom_boxplot(
                mapping = aes(x = genotype, y = resid_logcount_eqtl),
                outlier.shape = NA
            ) +
            geom_jitter(
                mapping = aes(
                    x = genotype, y = resid_logcount_eqtl, color = PrimaryDx
                )
            ) +
            scale_color_manual(values = dx_colors)
    } else {
        p = p +
            geom_boxplot(
                mapping = aes(
                    x = genotype, y = resid_logcount_eqtl, color = genotype
                ),
                outlier.shape = NA) +
            geom_jitter(
                mapping = aes(
                    x = genotype, y = resid_logcount_eqtl, color = genotype
                )
            ) +
            scale_color_manual(values = geno_colors)
    }

    #   Add text last to place it in front of points/ boxes in a couple
    #   unavoidable cases of overlap
    p = p +
        geom_text(
            data = label_df,
            mapping = aes(label = sig_label, x = Inf, y = -Inf),
            hjust = 1,
            vjust = 0,
            size = 6
        )
    
    pdf(
        file.path(plot_dir, paste0('expr_by_geno_', plot_suffix)),
        width = pdf_width, height = pdf_height
    )
    print(p)
    dev.off()

    #   Also plot a version where the y scale is free in the facet
    pdf(
        file.path(plot_dir, paste0('expr_by_geno_free_y_', plot_suffix)),
        width = pdf_width, height = pdf_height
    )
    print(p + facet_wrap(~ snp_id, scales = "free_y"))
    dev.off()
}

################################################################################
#   Read in eQTL, DEA, and GWAS results, and the RSE
################################################################################

eqtl = read_csv(eqtl_path, show_col_types = FALSE)

deg_full = read_tsv(deg_path, show_col_types = FALSE)

rse_gene = get(load(rse_path))

colData(rse_gene) = colData(rse_gene) |>
    as_tibble() |>
    #   Overwrite SNP PCs with the most recent values computed from the properly
    #   filtered genotyping data
    select(!matches('^snpPC')) |>
    left_join(
        read_tsv(snp_pcs_path, show_col_types = FALSE) |>
            #   Fix the name of one donor
            mutate(BrNum = ifelse(BrNum == "Br0983", "Br983", BrNum)),
        by = 'BrNum'
    ) |>
    #   Add in gene PCs
    left_join(readRDS(gene_pcs_path), by = "RNum") |>
    DataFrame()

colnames(rse_gene) = rse_gene$BrNum
rse_gene$PrimaryDx[rse_gene$PrimaryDx == "Schizo"] = "SCZD"

#   Due to repeated issues with the FTP servers performing the liftover,
#   save filtered copies and load from that rather than lifting over
#   in each interactive session

# gwas_narrow = read_excel(gwas_narrow_path) |>
#     filter(P < sig_cutoff_gwas) |>
#     dplyr::rename(chr = CHR, pos = BP) |>
#     #   Map from hg19 to hg38 and drop anything that fails
#     snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
#     filter(!is.na(pos)) |>
#     #   Construct variant_id from SNP info
    # mutate(
    #     variant_id = sprintf(
    #         'chr%s:%s:%s:%s',
    #         chr,
    #         pos,
    #         str_split_i(A1A2, '/', 1),
    #         str_split_i(A1A2, '/', 2)
    #     )
    # )
#
# write_csv(gwas_narrow, gwas_narrow_filt_path)
gwas_narrow = read_csv(gwas_narrow_filt_path, show_col_types = FALSE)

# gwas_wide = fread(gwas_wide_path) |>
#     as_tibble() |>
#     filter(PVAL < sig_cutoff_gwas) |>
#     dplyr::rename(chr = CHROM, pos = POS) |>
#     #   Map from hg19 to hg38 and drop anything that fails
#     snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
#     filter(!is.na(pos)) |>
#     #   Construct variant_id from SNP info
#     mutate(variant_id = sprintf('chr%s:%s:%s:%s', chr, pos, A1, A2))

# write_csv(gwas_wide, gwas_wide_filt_path)
gwas_wide = read_csv(gwas_wide_filt_path, show_col_types = FALSE)

################################################################################
#   Compare each type of result
################################################################################

#-------------------------------------------------------------------------------
#   Info about eQTLS alone
#-------------------------------------------------------------------------------

eqtl_gene = eqtl |>
    group_by(phenotype_id) |>
    summarize(n = n())

#   Significance hardcoded here because it was pre-filtered in 03_adj_p_vals.R
message(
    sprintf("%s significant SNP-gene pairs at FDR<0.05", nrow(eqtl))
)
message(
    sprintf("%s unique SNPs from these eQTLs", length(unique(eqtl$variant_id)))
)
message(
    sprintf(
        "%s Genes total and median %s SNPs per gene",
        nrow(eqtl_gene),
        median(eqtl_gene$n)
    )
)

#-------------------------------------------------------------------------------
#   Compare eQTLs with GWAS SNPs
#-------------------------------------------------------------------------------

gwas_paired_genes = list()
for (breadth in c('narrow', 'wide')) {
    if (breadth == 'narrow') {
        gwas = gwas_narrow
    } else {
        gwas = gwas_wide
    }

    message(
        sprintf(
            "Number of %s Trubetskoy et al. GWAS SNPs matching sig eQTLs:",
            breadth
        )
    )
    print(table(gwas$variant_id %in% eqtl$variant_id))

    gwas_paired_genes[[breadth]] = eqtl |>
        filter(variant_id %in% gwas$variant_id) |>
        pull(phenotype_id) |>
        unique()

    message(
        sprintf("Genes paired with %s-GWAS-significant eQTL SNPs:", breadth)
    )
    rowData(rse_gene)$Symbol[
        match(gwas_paired_genes[[breadth]], rownames(rse_gene))
    ] |>
        paste(collapse = ', ') |>
        message()
}

#-------------------------------------------------------------------------------
#   Compare eQTLS with DE genes
#-------------------------------------------------------------------------------

for (sig_cutoff in sig_cutoff_deg_explore) {
    deg =  deg_full |> filter(adj.P.Val < sig_cutoff)
    message(
        sprintf(
            "Number and name of DE genes at FDR<%s implicated in any eQTLs:",
            sig_cutoff
        )
    )
    print(table(deg$gencodeID %in% eqtl_gene$phenotype_id))
    deg$Symbol[deg$gencodeID %in% eqtl_gene$phenotype_id] |>
        paste(collapse = ', ') |>
        message()

    dea_paired_genes = deg$gencodeID[deg$gencodeID %in% eqtl_gene$phenotype_id]

    for (breadth in c('narrow', 'wide')) {
        message(
            sprintf(
                "Genes implicated in the %s GWAS, DEA at FDR<%s, and eQTLs:",
                breadth,
                sig_cutoff
            )
        )

        dea_paired_genes[dea_paired_genes %in% gwas_paired_genes[[breadth]]] |>
            paste(collapse = ', ') |>
            message()
    }
}

################################################################################
#   Read in and preprocess genotype info
################################################################################

#   Read in genotypes
plink = read.plink(plink_path)

#   Read in genotype metadata and line up with plink genotypes
geno_raw = read_delim(raw_geno_path, delim = '\t')
stopifnot(setequal(geno_raw$SNP, plink$map$snp.name))
geno_raw = geno_raw[match(plink$map$snp.name, geno_raw$SNP),]

#   Keep track of which genotypes should be flipped later, based on code from
#   https://github.com/LieberInstitute/brainseq_phase2/blob/d6779afe6f1509165b4c4eecdfdb7ff7a5d16a19/misc/pull_genotype_data.R#L76-L81
mismatched_snps = plink$map$snp.name[plink$map$allele.2 == geno_raw$COUNTED]

################################################################################
#   Plots exploring how genotype affects expression at select eQTLs
################################################################################

#   Model covariates to later call jaffelab::cleaningY() to residualize
#   expression for downstream plots
pd = as.data.frame(colData(rse_gene))
mod_eqtl <- model.matrix(
        as.formula(paste('~', paste(eqtl_covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(eqtl_covariates))]
mod_deg <- model.matrix(
        as.formula(paste('~', paste(deg_covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(deg_covariates))]

#   After exploring multiple significance cutoffs, choose one for plotting
deg = deg_full |> filter(adj.P.Val < sig_cutoff_deg_plot)
dea_paired_variants = eqtl |>
    filter(phenotype_id %in% deg$gencodeID) |>
    pull(variant_id)

#   Write paired variants to a text file. This will be read used to subset the
#   big VCF to ensure the below method for reading in genotypes (reading in the
#   plink bed file) works as VariantAnnotation::readVcf() does
if (opt$mode == "nominal") {
    writeLines(dea_paired_variants, paired_variants_path)
}

#   For all modes, plot any SNPs in an eQTL pair paired with a DEG
exp_df = merge_exp_df(
    rse_gene,
    mod_deg,
    mod_eqtl,
    eqtl |> filter(phenotype_id %in% deg$gencodeID),
    plink,
    mismatched_snps
)
plot_triad_exploratory(
    eqtl,
    exp_df,
    plot_dir,
    plot_prefix = sprintf(
        "eqtls_paired_with_dea_genes_FDR%s",
        as.character(sig_cutoff_deg_plot) |> str_split_i('\\.', 2)
    )
)

if (opt$mode == "independent") {
    #---------------------------------------------------------------------------
    #   Manuscript-ready plots for 3 DEG-paired eQTLs we want to showcase
    #---------------------------------------------------------------------------

    showcase_degs = c('TSPAN11', 'RNASEL', 'MYRFL')

    plot_suffix = sprintf(
        'eqtls_3_dea_genes_manuscript_FDR%s.pdf',
        as.character(sig_cutoff_deg_plot) |> str_split_i('\\.', 2)
    )
    
    this_exp_df = exp_df |>
        filter(gene_symbol %in% showcase_degs) |>
        #   Re-order genes to match geno plot below
        mutate(gene_symbol = factor(gene_symbol, levels = showcase_degs)) |>
        add_rs_id()
    
    #   Expression vs diagnosis for DEGs
    p = ggplot(
            this_exp_df,
            mapping = aes(
                x = PrimaryDx, y = resid_logcount_deg, color = PrimaryDx
            )
        ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        facet_wrap(~gene_symbol) +
        scale_color_manual(values = dx_colors) +
        labs(y = "Residualized DEG Expression") +
        theme_bw(base_size = 20) +
        theme(legend.position = "none")
    pdf(
        file.path(plot_dir, paste0('expr_by_dx_', plot_suffix)),
        width = 8, height = 6
    )
    print(p)
    dev.off()

    #   Expression by genotype for SNPs paired with DEGs
    exp_vs_geno_manuscript_plot(
        eqtl, this_exp_df, plot_dir, plot_suffix, note_risk_allele = FALSE,
        color_by_dx = TRUE
    )

    #---------------------------------------------------------------------------
    #   Manuscript-ready supplementary plots for remaining 4 DEG-paired eQTLs
    #---------------------------------------------------------------------------

    plot_suffix = sprintf(
        'eqtls_other_dea_genes_manuscript_FDR%s.pdf',
        as.character(sig_cutoff_deg_plot) |> str_split_i('\\.', 2)
    )

    #   The names of the remaining 4 of 7 DEGs paired in eQTLs
    showcase_degs = c('DPY19L2', 'RP11-624M8.1', 'ACCS', 'SMAGP')
    
    this_exp_df = exp_df |>
        filter(gene_symbol %in% showcase_degs) |>
        #   Re-order genes to match geno plot below
        mutate(gene_symbol = factor(gene_symbol, levels = showcase_degs)) |>
        add_rs_id()
    
    #   Expression vs diagnosis for DEGs
    p = ggplot(
            this_exp_df,
            mapping = aes(
                x = PrimaryDx, y = resid_logcount_deg, color = PrimaryDx
            )
        ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        facet_wrap(~gene_symbol, nrow = 1) +
        scale_color_manual(values = dx_colors) +
        labs(y = "Residualized DEG Expression") +
        theme_bw(base_size = 20) +
        theme(legend.position = "none")
    pdf(
        file.path(plot_dir, paste0('expr_by_dx_', plot_suffix)),
        width = 10, height = 6
    )
    print(p)
    dev.off()

    #   Expression by genotype for SNPs paired with DEGs
    exp_vs_geno_manuscript_plot(
        eqtl, this_exp_df, plot_dir, plot_suffix, note_risk_allele = FALSE,
        color_by_dx = TRUE, pdf_width = 10
    )

    #---------------------------------------------------------------------------
    #   For independent, plot (15) SNPs overlapping wider GWAS and their paired
    #   (16) genes
    #---------------------------------------------------------------------------

    gwas_variants = gwas_wide |>
        filter(variant_id %in% eqtl$variant_id) |>
        pull(variant_id)
    
    exp_df = merge_exp_df(
        rse_gene = rse_gene,
        mod_deg = mod_deg,
        mod_eqtl = mod_eqtl,
        eqtl = eqtl |> filter(variant_id %in% gwas_wide$variant_id),
        plink = plink,
        mismatched_snps = mismatched_snps
    )
    plot_triad_exploratory(
        eqtl,
        exp_df,
        plot_dir = plot_dir,
        plot_prefix = "wide_gwas_eqtls"
    )

    #---------------------------------------------------------------------------
    #   For a manuscript plot, we'll also want to sample 3 of these eQTLs and
    #   produce an expression-by-genotype plot faceted by eQTL. The remaining
    #   13 also become a supplementary figure
    #---------------------------------------------------------------------------

    gwas_3_snps = eqtl |>
        mutate(
            gene_symbol = rowData(rse_gene)$Symbol[
                match(phenotype_id, rowData(rse_gene)$gencodeID)
            ]
        ) |>
        filter(
            variant_id %in% gwas_wide$variant_id,
            gene_symbol %in% c('DND1P1', 'CRHR1-IT1', 'LRRC37A4P')
        ) |>
        pull(variant_id)
    
    exp_df_gwas = exp_df |>
        add_rs_id() |>
        #   Grab risk allele from the GWAS data
        left_join(
            gwas_wide |>
                mutate(risk_allele = ifelse(BETA > 0, A1, A2)) |>
                dplyr::rename(snp_id = variant_id) |>
                dplyr::select(snp_id, risk_allele),
            by = 'snp_id'
        )
    
    exp_vs_geno_manuscript_plot(
        eqtl,
        exp_df_gwas |> filter(snp_id %in% gwas_3_snps),
        plot_dir,
        plot_suffix = '3_gwas_wide_manuscript.pdf',
        note_risk_allele = TRUE,
        color_by_dx = FALSE
    )

    exp_vs_geno_manuscript_plot(
        eqtl,
        exp_df_gwas |> filter(!(snp_id %in% gwas_3_snps)),
        plot_dir,
        plot_suffix = 'other_gwas_wide_manuscript.pdf',
        note_risk_allele = TRUE,
        color_by_dx = FALSE,
        facet_nrow = 3,
        pdf_width = 12,
        pdf_height = 12
    )

    #---------------------------------------------------------------------------
    #   Also plot the top 10 (by significance) eQTLs for independent
    #---------------------------------------------------------------------------
    
    exp_df = merge_exp_df(
        rse_gene = rse_gene,
        mod_deg = mod_deg,
        mod_eqtl = mod_eqtl,
        eqtl = eqtl |> arrange(FDR) |> slice_head(n = 10),
        plink = plink,
        mismatched_snps
    )
    plot_triad_exploratory(
        eqtl,
        exp_df,
        plot_dir = plot_dir,
        plot_prefix = "top_10_eqtls"
    )

    #---------------------------------------------------------------------------
    #   For eQTLs containing a DEG or a wide GWAS SNP, compare expression
    #   residualized by the DEG vs. eQTL models 
    #---------------------------------------------------------------------------

    exp_df = merge_exp_df(
            rse_gene,
            mod_deg,
            mod_eqtl,
            eqtl |>
                dplyr::filter(
                    (variant_id %in% gwas_wide$variant_id) |
                    (phenotype_id %in% deg$gencodeID)
                ),
            plink,
            mismatched_snps
        ) |>
        mutate(
            symbol = rowData(rse_gene)$Symbol[
                match(gene_id, rownames(rse_gene))
            ]
        )

    p = ggplot(exp_df, aes(x = resid_logcount_deg, y = resid_logcount_eqtl)) +
        geom_point() +
        theme_bw(base_size = 15) +
        labs(
            x = "Residualized DEG Expression",
            y = "Residualized eQTL Expression"
        )
    pdf(file.path(plot_dir, 'resid_comparison.pdf'))
    print(p)
    dev.off()
    
    #---------------------------------------------------------------------------
    #   Next, find SNP ID ("rs ID") for SNPs overlapping the wide GWAS or
    #   paired with a DEG
    #---------------------------------------------------------------------------

    deg = deg_full |>
        filter(adj.P.Val < sig_cutoff_deg_plot)

    a = eqtl |>
        separate(
            variant_id, into = c("chr", "pos", "A1", "A2"), sep = ":",
            remove = FALSE
        ) |>
        mutate(chr = str_extract(chr, '^chr(.*)', group = 1)) |>
        dplyr::filter(
            (variant_id %in% gwas_wide$variant_id) |
            (phenotype_id %in% deg$gencodeID)
        )

    #   Individually search for rs IDs
    a$rs_id = sapply(
        seq(nrow(a)), function(i) {
            get_rsid_from_position(
                chrom = a$chr[i], pos = a$pos[i], ref = a$A1[i], alt = a$A2[i],
                assembly = "hg38"
            )
        }
    )

    #   Write to a CSV of results for manually querying BSP1+BSP2 eQTL browser
    a |>
        mutate(
            source = ifelse(
                variant_id %in% gwas_wide$variant_id,
                "gwas_wide",
                "independent_deg_snp"
            ),
            symbol = rowData(rse_gene)$Symbol[
                match(phenotype_id, rownames(rse_gene))
            ]
        ) |>
        dplyr::select(phenotype_id, symbol, variant_id, rs_id, source) |>
        write_csv(rs_path)
}

################################################################################
#   Add dummy plot with legends for genotype and diagnosis colors
################################################################################

if(opt$mode == 'independent') {
    p = ggplot(
            tibble(x = 1:3, y = 1:3, geno = as.character(0:2)),
            mapping = aes(x = x, y = y, color = geno)
        ) +
        geom_point() +
        scale_color_manual(values = geno_colors)
    
    pdf(file.path(plot_dir, 'geno_legend.pdf'))
    print(p)
    dev.off()

    p = ggplot(
            tibble(x = 1:2, y = 1:2, dx = c("Control", "SCZD")),
            mapping = aes(x = x, y = y, color = dx)
        ) +
        geom_point() +
        scale_color_manual(values = dx_colors)
    
    pdf(file.path(plot_dir, 'dx_legend.pdf'))
    print(p)
    dev.off()
}

session_info()
