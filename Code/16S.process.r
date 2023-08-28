
lapply(c('ggfortify', 'psych','msa', 'ggplot2','RColorBrewer', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn'), require,character.only=TRUE) #add as necessary

fully.process<-function(forward_read_names=".fastq.gz", seed=1, displayx=10, filetype=".jpg", maxeef=2, trainingset="silva_nr_v138_train_set.fa.gz", truncLenF=240, tree=FALSE,multithread=TRUE, loopnr=1){

path= getwd()
fnFs= sort(list.files(path, pattern= forward_read_names))

sample.names = sapply(strsplit(fnFs, "_"), `[`, 1)  #extract sample names
fnFs = file.path(path, fnFs) #specify global path to the sequence files

filt_path= file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs=file.path(filt_path, paste0(sample.names, forward_read_names))#Create a subdirectory and file names for the filtered files

out=filterAndTrim(fnFs, filtFs, truncLen=c(truncLenF), maxN=0, maxEE=c(maxeef), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

errF = learnErrors(filtFs, multithread=TRUE)

derepFs = derepFastq(filtFs, verbose=TRUE)

names(derepFs) = sample.names

dadaFs = dada(derepFs, err=errF, multithread=TRUE)

seqtab = makeSequenceTable(dadaFs)

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) = c("Original Reads", "Filtered Reads", "Denoised Reads", "tabled", "Non-chimeric Reads")
rownames(track) = sample.names
write.table(track, paste0(loopnr,"tracking_table.txt"), sep="\t")


taxa = assignTaxonomy(seqtab.nochim, trainingset, multithread=TRUE)

merged = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  tax_table(taxa))
merged@sam_data=sample_data(as.data.frame(sample_names(merged)))
sample_names(sample_data(merged))=sample_names(merged)
sample_data(merged)[ , 2] <- sample_data(merged)[ ,1]
colnames(merged@sam_data)=c("names", "dummy")
saveRDS(merged, paste0(loopnr,"merged.rds"))

write.table(sort(sample_sums(merged)), paste0(loopnr,"reads.per.sample.txt"), sep="\t")}


tests=seq(from = 50, to = 200, by = 10)

for (i in 1:length(tests)){ 
	fully.process(forward_read_names=".fastq.gz", maxeef=2, truncLenF=tests[i], loopnr=tests[i])}

files <- list.files(pattern="tracking_table.txt", recursive=TRUE)

tracking_all= NULL

#merge all tracking tables
for (f in 1:length(files)) {
   dat <- read.table(files[f], header=T, sep="\t", na.strings="", row.names=NULL)
   dat$forward.trim=paste0(gsub(".tracking_table.txt", "", files[f]),"0")
   tracking_all <- rbind(tracking_all, dat)
}


write.table(tracking_all, "tracking_all.txt")

