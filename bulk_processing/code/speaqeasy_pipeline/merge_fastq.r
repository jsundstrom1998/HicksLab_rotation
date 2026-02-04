#Load libraries
library('getopt')
library('BiocParallel')
library('devtools')

#Store file locations as variables in a list
optloc <- list(
    samples = '/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest',
    output = '/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/merged_fastq'
)

#Disable testing
testing <- FALSE

#Read in manifest 
manifest <- read.table(optloc$samples, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE)

#Check if paired end
paired <- ncol(manifest) > 3

#Find file extensions
files <- manifest[, 1]
extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
ext_found <- sapply(files, function(file) {
    extensions[names(unlist(sapply(patterns, grep, file))) == patterns]
})

if(any(is.na(ext_found))) {
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}

extensions <- split(ext_found, manifest[, ncol(manifest)])

#Check if extensions are the same per group
if(any(sapply(extensions, function(x) { length(unique(x)) }) != 1)) {
    stop("For each sample name, the extensions of the fastq files to be merged have to be the same")
}

#Create output directory
dir.create(optloc$output, showWarnings = FALSE, recursive = TRUE)

#Split according to the sample names
file_groups <- split(manifest, manifest[, ncol(manifest)])
extensions <- sapply(extensions, '[', 1)

merge_files <- function(file_names, new_file) {
    message(paste(Sys.time(), 'creating', new_file))
    call <- paste('cat', paste(file_names, collapse = ' '), '>', new_file)
    print(call)
    if(!testing) system(call)
}

res <- bpmapply(function(common, new_name, extension) {
    merge_files(common[, 1],
        file.path(optloc$output, paste0(new_name, '.', extension)))
    if(paired) {
        merge_files(common[, 3],
            file.path(optloc$output, paste0(new_name, '_read2.', extension)))
    }
}, file_groups, names(file_groups), extensions, 
    BPPARAM = MulticoreParam(1))

#Reproducibility info
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()