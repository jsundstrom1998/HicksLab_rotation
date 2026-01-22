#Load libraries
library('getopt')
library('devtools')

#Store file locations as variables in a list
loc <- list(
    samples = '/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest',
    output = '/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/'
)

#Read in manifest 
manifest <- read.table(loc$samples, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE)
G
#Check if paired end
paired <- ncol(manifest) > 3
if(paired) system('touch .paired_end')

## Find file extensions
files <- manifest[, 1]
extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
ext_found <- sapply(files, function(file) {
    extensions[names(unlist(sapply(patterns, grep, file))) == patterns]
})