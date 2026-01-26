#Load libraries
library('getopt')
library('devtools')

#Store file locations as variables in a list
optloc <- list(
    samples = '/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest',
    output = '/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/'
)

#optloc <- getopt(loc)

#Read in manifest 
manifest <- read.table(optloc$samples, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE)

#Check if paired end
paired <- ncol(manifest) > 3
if(paired) system('touch .paired_end')

## Find file extensions & create error file if unknown extension found
files <- manifest[, 1]
extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
ext_found <- sapply(files, function(file) {
    extensions[names(unlist(sapply(patterns, grep, file))) == patterns]
})

if(any(is.na(ext_found))) {
    system('touch find_sample_error')
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}

extensions <- split(ext_found, manifest[, ncol(manifest)])

## Check if merging is required
merged <- length(unique(manifest[, ncol(manifest)])) == nrow(manifest)
if(!merged) {
    system('touch .requires_merging')
    message(paste0(Sys.time(), ' creating .samples_unmerged.manifest'))
    system(paste('mv', optloc$sample, file.path(dirname(optloc$sample),
        '.samples_unmerged.manifest')))

    message(paste(Sys.time(), 'creating the new samples.manifest file with the merged samples'))

    ## Split according to the sample names
    file_groups <- split(manifest, manifest[, ncol(manifest)])
    extensions <- sapply(extensions, '[', 1)
    
    if(paired) {
        new_manifest <- data.frame(
            file.path(optloc$output, paste0(names(file_groups), '.', extensions)),
            rep(0, length(file_groups)),
            file.path(optloc$output, paste0(names(file_groups), '_read2.',
                extensions)),
            rep(0, length(file_groups)),
            names(file_groups), stringsAsFactors = FALSE
        )
    } else {
        new_manifest <- data.frame(
            file.path(optloc$output, paste0(names(file_groups), '.', extensions)),
            rep(0, length(file_groups)),
            names(file_groups), stringsAsFactors = FALSE
        )
    }
    ## Make names short, in case you want to interactively check the new manifest
    colnames(new_manifest) <- paste0('V', seq_len(ncol(new_manifest)))

    write.table(new_manifest, file = optloc$sample, row.names = FALSE,
        col.names = FALSE, quote = FALSE, sep = '\t')
}

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()