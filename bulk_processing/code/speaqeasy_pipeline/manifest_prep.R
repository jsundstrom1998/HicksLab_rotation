###
#Load packages
library(jaffelab)

#Define FASTQ directory and list all files
fqPath = "/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/"
fqFiles = list.files(fqPath, recur=TRUE)

#Filter out undefined reads and split remaining into left + right
fqFiles = fqFiles[!grepl("^Und", fqFiles)]

leftReads = fqFiles[grep("_R1_", fqFiles)]
rightReads = fqFiles[grep("_R2_", fqFiles)]

#Create manifest table
man = data.frame(leftReads = paste0(fqPath, leftReads),
	leftMd5 = 0,rightReads = paste0(fqPath, rightReads), 
	rightMd5 = 0, SampleID = ss(leftReads, "_"),stringsAsFactors=FALSE) 
	
write.table(man, file="/dcs04/hicks/data/jsundstr/bulk_processing/fastq_test/samples.manifest", 
	sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)
print("Done")

