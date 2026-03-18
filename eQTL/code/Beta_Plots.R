library(here)
library(tidyverse)
library(SummarizedExperiment)
library(sessioninfo)
library(data.table)
library(jaffelab)
library(cowplot)
library(ggrepel)
library(bigsnpr)
library(MRutils)

eqtl_independent_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'independent', 'FDR05.csv'
)
eqtl_int_path = here(
    'processed-data', '17_eQTL', 'tensorQTL_output', 'combined_interaction_subset.csv'
)
deg_path = here(
    'processed-data', '10_DEA', '04_DEA',
    'DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
)
gwas_wide_filt_path = here(
    'processed-data', '17_eQTL', 'gwas_wide_filtered.csv'
)
gwas_wide_path = here(
    "processed-data", "13_MAGMA","GWAS", "scz2022",
    "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
)
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
bsp2_path = here('processed-data', '17_eQTL', 'BSP2_cleaned_eqtl.csv')
supp_tab_path = here(
    'processed-data', '17_eQTL', 'independent_eqtl_supp_tab.csv'
)

plot_dir = here('plots', '17_eQTL')

sig_cutoff_deg = 0.1

source_colors = c("#14599D", "#78574C")
names(source_colors) = c("GWAS SNP", "DEG")

larger_bsp2_colors = c("#808080", "#000000")

lift_over_path = system('which liftOver', intern = TRUE)

################################################################################
#   Functions
################################################################################

# MRutils::get_rsid_from_position() sometimes fails due to temporary network
# issues, and returns NULL instead of NA when this happens. Write a wrapper
# which re-runs the function several times upon errors and returns a
# character value or NA
get_rsid_from_position_robust = function(variant_id, num_tries = 5) {
    #   If the ref or alt is more than one basepair, it isn't a SNP and
    #   therefore doesn't have an RS ID
    ref = str_split_i(variant_id, ':', 3)
    alt = str_split_i(variant_id, ':', 4)
    if ((nchar(ref) > 1) || (nchar(alt) > 1)) {
        return(NA)
    }

    rs_id = NULL
    attempt_num = 0
    while(is.null(rs_id) && attempt_num < num_tries) {
        rs_id = get_rsid_from_position(
            chrom = str_split_i(variant_id, ':', 1),
            pos = as.integer(str_split_i(variant_id, ':', 2)),
            ref = str_split_i(variant_id, ':', 3),
            alt = str_split_i(variant_id, ':', 4),
            assembly = "hg38"
        )
        attempt_num = attempt_num + 1
    }
    if (is.null(rs_id)) return(NA)
    return(rs_id)
}

#   Since MRutils::get_rsid_from_position() is so slow, we can't afford to
#   run it on the same SNP twice. Run on the unique set and then duplicate
#   the results where necessary
get_rsid_from_position_robust_fast = function(variant_id) {
    rs_df = tibble(
        variant_id = unique(variant_id),
        rs_id = map_chr(unique(variant_id), get_rsid_from_position_robust)
    )
    return(rs_df$rs_id[match(variant_id, rs_df$variant_id)])
}

################################################################################
#   Read in and preprocess data
################################################################################

#   For adding gene symbol
rse_gene = get(load(rse_path))

#-------------------------------------------------------------------------------
#   DEGs and ~20k GWAS SNPs
#-------------------------------------------------------------------------------

#   FDR < 0.1 DEGs
deg = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < sig_cutoff_deg)

#   ~20k "wide" GWAS SNPs
gwas_wide = read_csv(gwas_wide_filt_path, show_col_types = FALSE)

#-------------------------------------------------------------------------------
#   Independent and interaction habenula eQTLs
#-------------------------------------------------------------------------------

#   Read in all significant independent eQTLs
eqtl_independent = read_csv(eqtl_independent_path, show_col_types = FALSE) |>
    mutate(pair_id = paste(phenotype_id, variant_id, sep = '_'))

#   The same eQTLs but taking stats from the interaction model
eqtl_int = read_csv(eqtl_int_path, show_col_types = FALSE)

#   Warn that not all eQTLs were present in interaction models
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in habenula interaction model",
        nrow(eqtl_int |> filter(interaction_var == "hb")),
        nrow(eqtl_independent)
    )
)
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in thalamus interaction model",
        nrow(eqtl_int |> filter(interaction_var == "thal")),
        nrow(eqtl_independent)
    )
)

#-------------------------------------------------------------------------------
#   BSP2 eQTLs with genotypes
#-------------------------------------------------------------------------------

bsp2 = read_csv(bsp2_path, show_col_types = FALSE)

#   Warn that not all eQTLs were present in BSP2
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in BSP2 DLPFC",
        sum(bsp2$region == "DLPFC"),
        nrow(eqtl_independent)
    )
)
message(
    sprintf(
        "Only %s of %s independent eQTLs observed in BSP2 hippo",
        sum(bsp2$region == "Hippocampus"),
        nrow(eqtl_independent)
    )
)

################################################################################
#   Thalamus vs. habenula interaction plots for eQTLs overlapping DEGs or
#   wide GWAS SNPs
################################################################################

#   Keep a filtered copy to those overlapping either the wide GWAS or FDR < 0.1
#   DEGs
filt_eqtl_independent = eqtl_independent |>
    filter(
        phenotype_id %in% deg$gencodeID | 
        variant_id %in% gwas_wide$variant_id
    )

#   Note that there are 11 SNPs overlapping the GWAS data, but these SNPs appear
#   in 13 eQTLs total
message(
    "Num SNPs overlapping wide GWAS data: ",
    eqtl_independent |>
        filter(variant_id %in% gwas_wide$variant_id) |>
        pull(variant_id) |>
        unique() |>
        length()
)
message(
    "Num eQTLs including those SNPs: ",
    eqtl_independent |>
        filter(variant_id %in% gwas_wide$variant_id) |>
        nrow()
)

hb_pairs = eqtl_int |> filter(interaction_var == "hb") |> pull(pair_id)
thal_pairs = eqtl_int |> filter(interaction_var == "thal") |> pull(pair_id)
eqtl_int_both = eqtl_int |>
    dplyr::filter(
        pair_id %in% hb_pairs,
        pair_id %in% thal_pairs,
        pair_id %in% filt_eqtl_independent$pair_id
    )

#   Make note of filtered independent pairs not measured in each interaction
#   model
message("Pairs present in independent but absent from at least one interaction model:")
filt_eqtl_independent |>
    mutate(
        symbol = rowData(rse_gene)$Symbol[
            match(phenotype_id, rowData(rse_gene)$gencodeID)
        ]
    ) |>
    filter(!pair_id%in% eqtl_int_both$pair_id) |>
    select(symbol, variant_id) |>
    print()

p = eqtl_int_both |>
    dplyr::select(phenotype_id, variant_id, b_gi, pval_gi, interaction_var) |>
    mutate(
        gene_symbol = rowData(rse_gene)$Symbol[
            match(phenotype_id, rowData(rse_gene)$gencodeID)
        ],
        source = factor(
            ifelse(variant_id %in% gwas_wide$variant_id, "GWAS SNP", "DEG")
        )
    ) |>
    pivot_wider(
        values_from = c("b_gi", "pval_gi"), names_from = "interaction_var"
    ) |>
    mutate(
        font_face = ifelse(
            (b_gi_hb > 0) & (b_gi_thal < 0), "bold.italic", "plain"
        )
    ) |>
    ggplot(mapping = aes(x = b_gi_hb, y = b_gi_thal, color = source)) +
        geom_point(size = 3) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_label_repel(
            aes(label = gene_symbol, fontface = font_face), max.overlaps = 15,
            show.legend = FALSE
        ) +
        scale_color_manual(values = source_colors) +
        theme_bw(base_size = 23) +
        labs(
            x = "Beta: Habenula Interaction", y = "Beta: Thalamus Interaction",
            color = "eQTL Overlap"
        )

pdf(file.path(plot_dir, 'interaction_beta.pdf'), width = 10, height = 7)
print(p)
dev.off()

################################################################################
#   BSP2 vs habenula beta value plots at the same eQTLs
################################################################################

p_list = list()
for (this_region in c("DLPFC", "Hippocampus")) {
    this_bsp2 = bsp2 |>
        filter(region == this_region) |>
        select(pair_id, beta) |>
        dplyr::rename(beta_bsp2 = beta) |>
        left_join(
            eqtl_independent |>
                select(pair_id, slope) |>
                dplyr::rename(beta_habenula = slope),
            by = "pair_id"
        ) |>
        mutate(
            larger_bsp2 = (beta_bsp2 < abs(beta_habenula)) &
                (sign(beta_bsp2) == sign(beta_habenula))
        )
    this_bsp2_lm = lm(beta_bsp2 ~ beta_habenula, this_bsp2)
    message(
        sprintf(
            'BSP2 %s: y = %sx + %s',
            this_region,
            signif(this_bsp2_lm$coefficients[2], 3),
            signif(this_bsp2_lm$coefficients[1], 3)
        )
    )
    p_list[[this_region]] = ggplot(
            this_bsp2,
            mapping = aes(x = beta_habenula, y = beta_bsp2)
        ) +
        geom_point(mapping = aes(color = larger_bsp2)) +
        geom_abline(slope = 1, color = "#a98743") +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        geom_smooth(method = "lm", formula = y ~ x) +
        scale_color_manual(values = larger_bsp2_colors) +
        guides(color = "none") +
        theme_bw(base_size = 20) +
        labs(x = "Beta: Habenula", y = paste("Beta:", this_region))
}

pdf(file.path(plot_dir, 'BSP2_vs_habenula_beta.pdf'), width = 12, height = 6)
print(plot_grid(plotlist = p_list))
dev.off()

pdf(
    file.path(plot_dir, 'BSP2_vs_habenula_beta_coord_fixed.pdf'),
    width = 12, height = 6
)
print(
    plot_grid(
        plotlist = list(
            p_list[["DLPFC"]] + coord_fixed(),
            p_list[["Hippocampus"]] + coord_fixed()
        )
    )
)
dev.off()

################################################################################
#   Supplementary table including significant independent eQTLs and interaction
#   stats at those eQTLs as well
################################################################################

eqtl_int_man = eqtl_int |>
    #   Most of these columns will be identical with those found in the
    #   independent SNPs, and need not be duplicated
    dplyr::select(-c(pair_id, start_distance, ma_samples, ma_count, af)) |>
    dplyr::rename(FDR_gi = FDR) |>
    mutate(
        interaction_var = ifelse(
            interaction_var == "hb", "habenula", "thalamus"
        )
    ) |>
    #   Pivot wider so thalamus and habenula interaction stats have their own
    #   columns
    pivot_wider(
        values_from = !matches('^(phenotype_id|variant_id|interaction_var)$'),
        names_from = "interaction_var",
        names_glue = "{interaction_var}_{.value}"
    ) |>
    #   Clarify these stats come from interaction runs
    rename_with(
        ~ paste("inter", .x, sep = "_"),
        !matches('^(phenotype_id|variant_id)$')
    )

#   Get variant ID and p value for all wide GWAS SNPs
gwas_wide = fread(gwas_wide_path) |>
    as_tibble() |>
    dplyr::rename(chr = CHROM, pos = POS, GWAS_p_val = PVAL, snp_rs_id = ID) |>
    #   Map from hg19 to hg38 and drop anything that fails
    snp_modifyBuild(lift_over_path, from = 'hg19', to = 'hg38') |>
    filter(!is.na(pos)) |>
    #   Construct variant_id from SNP info
    mutate(variant_id = sprintf('chr%s:%s:%s:%s', chr, pos, A1, A2)) |>
    dplyr::select(GWAS_p_val, variant_id, snp_rs_id)

eqtl_independent = eqtl_independent |>
    #   Remove unnecessary columns and 'end_distance', which is just a duplicate
    #   of 'start_distance'
    dplyr::select(-c(`...1`, pair_id, end_distance)) |>
    dplyr::rename(FDR_beta = FDR) |>
    mutate(
        gene_symbol = rowData(rse_gene)$Symbol[
            match(phenotype_id, rownames(rse_gene))
        ]
    ) |>
    #   Clarify which stats come from independent run
    rename_with(
        ~ paste("indep", .x, sep = "_"),
        !matches('^(phenotype_id|gene_symbol|variant_id|snp_rs_id|start_distance|af|ma_samples|ma_count)$')
    ) |>
    #   Add interaction stats, where they exist, at independent SNPs
    left_join(eqtl_int_man, by = c('phenotype_id', 'variant_id')) |>
    #   Add FDR from differential expression for each gene
    left_join(
        read_tsv(deg_path, show_col_types = FALSE) |>
            dplyr::rename(phenotype_id = gencodeID, DEG_FDR = adj.P.Val) |>
            dplyr::select(phenotype_id, DEG_FDR),
        by = 'phenotype_id'
    ) |>
    #   Add p value from GWAS (and RS ID) for each SNP
    left_join(gwas_wide, by = 'variant_id')

#   Query RS ID where it's missing
eqtl_independent$snp_rs_id[is.na(eqtl_independent$snp_rs_id)] =
    get_rsid_from_position_robust_fast(
        eqtl_independent$variant_id[is.na(eqtl_independent$snp_rs_id)]
    )

message(
    sprintf(
        '%s of %s RS IDs still missing.',
        length(which(is.na(eqtl_independent$snp_rs_id))),
        nrow(eqtl_independent)
    )
)

eqtl_independent |>
    #   Gene and SNP first
    relocate(phenotype_id, gene_symbol, variant_id, snp_rs_id) |>
    write_csv(supp_tab_path)

session_info()