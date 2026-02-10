#Load libraries 
library("SummarizedExperiment")
library("sessioninfo")
library("tidyverse")
library("jaffelab")
library("here")
library("data.table")

#Define data locations relative to project root
data_directory <- "/dcl02/lieber/ajaffe/Roche_Habenula"

rse_path = here(
    data_directory, "processed-data", "rse_objects", "rse_gene_Habenula_Pilot.rda"
)

pca_path = here(
    data_directory, "processed-data", "03_bulk_pca", "PCs.rds"
)

snp_pcs_path = here(
    data_directory, "processed-data", "08_bulk_snpPC", "v3", 
    "habenula_R.9_MAF.05.RSann_filt.snpPCs.tab"
)

covariates = c(
    "PrimaryDx", 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5',
    paste0('PC', 1:13)
)

#Create output directory
output_dir = here("eQTL", "tensorQTL_input")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Convert SummarizedExperiment to bed format

rse_to_bed <- function(rse, assay_name = "logcounts") {
    stopifnot(assay_name %in% assayNames(rse))

    rr_df <- as.data.frame(rowRanges(rse)) |>
        tibble::rownames_to_column("ID") |>
        dplyr::mutate(
            start = ifelse(strand == "+", start, end),
            end = start + 1
        ) |>
        dplyr::arrange(seqnames, start) |>
        dplyr::select(`#Chr` = seqnames, start, end, ID)
  
    message(
        Sys.time(), " - Converting assay '", assay_name, "' to bed format: (",
        nrow(rse), ", ", ncol(rse),")"
    )
    counts <- SummarizedExperiment::assays(rse)[[assay_name]]

    counts <- rr_df |> 
        dplyr::left_join(
            as.data.frame(counts) |> 
                tibble::rownames_to_column("ID"), 
                by = "ID"
            )
    return(counts)
}

#Write covariate file

rse = get(load(rse_path, verbose = TRUE))

colData(rse) = colData(rse) |>
    as_tibble() |>
    #Overwrite SNP PCs with the most recent values computed from the properly
    #filtered genotyping data
    select(!matches('^snpPC')) |>
    left_join(
        read_tsv(snp_pcs_path, show_col_types = FALSE) |>
            #Fix the name of one donor
            mutate(BrNum = ifelse(BrNum == "Br0983", "Br983", BrNum)),
        by = 'BrNum'
    ) |>
    #Add in gene PCs
    left_join(readRDS(pca_path), by = "RNum") |>
    DataFrame()
colnames(rse) = rse$BrNum

#Also write colData to CSV, to have easy access to potential interaction covariates
col_data = as_tibble(colData(rse))
colnames(col_data) = gsub('\\.', '_', colnames(col_data)) 
write_csv(col_data, file.path(output_dir, 'colData.csv'))

message(Sys.time(), " - Format covariates")

#Model covariates and save
pd = as.data.frame(colData(rse)[, covariates])
covars <- model.matrix(
        as.formula(paste('~', paste(covariates, collapse = " + "))),
        data = pd
    )[, 2:(1 + length(covariates))] |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id")
corner(covars)

write_tsv(covars, file = file.path(output_dir, "covariates.txt"))

#Write logcounts to a ".bed.gz" file

message(Sys.time(), " - Format logcounts as bed in memory")
expression_bed <- rse_to_bed(rse)
corner(expression_bed)

message(Sys.time(), " - Write bed table to .bed.gz file")
data.table::fwrite(
    expression_bed, here(output_dir, "logcounts.bed.gz"), sep = "\t",
    quote = FALSE, row.names = FALSE
)

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()