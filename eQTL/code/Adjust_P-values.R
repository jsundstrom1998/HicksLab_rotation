library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")
library("qs2")
library("getopt")

data_directory = "/users/jsundstr/hicks_home/eQTL/"
#pd_nominal = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/17_eQTL/tensorQTL_output/nominal"
#pd_inter_thal = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/17_eQTL/tensorQTL_output/interaction_tot_Thal"
#pd_inter_Hb = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/17_eQTL/tensorQTL_output/interaction_tot_Hb"

#Read in which tensorQTL run mode is being used
spec <- matrix(
    c(
        "mode", "m", 1, "character", "tensorQTL run mode",
        "covariate", "c", 2, "character", "covariate applicable for interaction model"
    ),
    byrow = TRUE, ncol = 5
)

#Parse the command line arguments based on the specification
opt <- getopt(spec)

#Check if the provided mode is valid
accepted_modes = c('nominal', 'cis', 'independent', 'interaction')
if (!(opt$mode %in% accepted_modes)) {
    stop(
        sprintf(
            "'opt$mode' must be in '%s'.",
            paste(accepted_modes, collapse = "', '")
        )
    )
}

#Set output folder name suffix based on the mode and covariate (if interaction)
if (opt$mode == "interaction") {
    out_dir_suffix = sprintf('%s_%s', opt$mode, opt$covariate)
} else {
    out_dir_suffix = opt$mode
}

#Set output directory path
out_dir = here(data_directory, "tensorQTL_output", out_dir_suffix)
p_val_cutoff = 0.05

#Retrieve the relevant parquet files based on the mode and read them in, adjusting the p-values accordingly
if (opt$mode == 'nominal') {
    parquet_files <- list.files(
        out_dir,
        pattern = "\\.parquet$",
        full.names = TRUE
    )

    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) |>
        mutate(FDR = p.adjust(pval_nominal, "fdr"))
} else if (opt$mode == "interaction") {
    parquet_files <- list.files(
        out_dir,
        pattern = "\\.parquet$",
        full.names = TRUE
    )

    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) |>
        mutate(FDR = p.adjust(pval_gi, "fdr"))
    
#
    write_csv(eqtl_out, file.path(out_dir, 'all.csv'))
} else if (opt$mode %in% c('cis', 'independent')) {
    eqtl_out = read_csv(
            file.path(out_dir, paste0(opt$mode, '_out.csv')),
            show_col_types = FALSE
        ) |>
        mutate(FDR =  p.adjust(pval_beta, "fdr"))
} else { 
    stop(sprintf("Unsupported 'opt$mode' = '%s'.", opt$mode))
}

message("n pairs: ", nrow(eqtl_out))

#Filter
eqtl_out_filtered <- eqtl_out |>
    filter(FDR < p_val_cutoff)
message("n pairs FDR<", p_val_cutoff, ": ", nrow(eqtl_out_filtered))

fn <- here(
    out_dir,
    paste0(
        "FDR", as.character(p_val_cutoff) |> str_extract('\\.(.*)$', group = 1)
    )
)

#Save as CSV and the faster qs
write_csv(eqtl_out_filtered, file = paste0(fn, ".csv"))
qs2::qs_save(eqtl_out_filtered, file = paste0(fn, ".qs"))

#Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()