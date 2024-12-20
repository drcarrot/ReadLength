
lapply(c('ggfortify', 'psych','msa', 'ggplot2','RColorBrewer','phangorn', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn', 'DECIPHER', 'hillR'), require,character.only=TRUE) #add as necessary

print("A log of processing is saved as log.txt")
sink(file="log.txt")


fully.process<-function(forward_read_names="_1.fastq.gz", seed=1, trimleft=25, maxeef=2, trainingset="silva_nr_v138_train_set.fa.gz", sptrainingset="silva_species_assignment_v138.1.fa.gz", truncLenF=220, tree=TRUE,multithread=FALSE, loopnr=1){
   
  path= getwd()
  fnFs= sort(list.files(path, pattern= forward_read_names))
  
  sample.names = sapply(strsplit(fnFs, "_"), `[`, 1)  #extract sample names
  fnFs = file.path(path, fnFs) #specify global path to the sequence files
  
  filt_path= file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
  filtFs=file.path(filt_path, paste0(sample.names, forward_read_names))#Create a subdirectory and file names for the filtered files
  
  Forwardpqp=plotQualityProfile(fnFs[1])
  ggsave("Qualityprof.jpg",Forwardpqp)
  
  out=filterAndTrim(fnFs, filtFs,  truncLen=truncLenF,trimLeft=10, maxN=0, maxEE=maxeef, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
  
  errF = learnErrors(filtFs, multithread=FALSE)

  derepFs = derepFastq(filtFs, verbose=TRUE)

  names(derepFs) = sample.names

  dadaFs = dada(derepFs, err=errF, multithread=FALSE)

  seqtab = makeSequenceTable(dadaFs)
  
  #seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
  
  taxa = assignTaxonomy(seqtab, trainingset, tryRC=TRUE, multithread=FALSE)
  taxa = addSpecies(taxa, refFasta = sptrainingset, tryRC=TRUE)
  unname(head(taxa))
  taxa=as.data.frame(taxa)
  taxa=as.matrix(taxa)
  

  if (tree=="TRUE"){
    seqs <- getSequences(seqtab)
    names(seqs) <- seqs # This propagates to the tip labels of the tree
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
  
    phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phangAlign)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit = pml(treeNJ, data=phangAlign)
    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
    
    merged = phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),  tax_table(taxa),phy_tree(fitGTR$tree))
  } else{
    merged = phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),  tax_table(taxa))
}


  
  merged@sam_data=sample_data(as.data.frame(sample_names(merged)))
  sample_names(sample_data(merged))=sample_names(merged)
  sample_data(merged)[ , 2] <- sample_data(merged)[ ,1]
  colnames(merged@sam_data)=c("names", "dummy")
  merged = subset_taxa(merged, Kingdom== "Bacteria")
  saveRDS(merged, paste0(loopnr, "merged.rds"))
 
  getN = function(x) sum(getUniques(x))
  track = cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab), sample_sums(merged))
  colnames(track) = c("Original Reads", "Filtered Reads", "Denoised Reads",  "tabled", "Non-chimeric Reads", "BacterialReads")
  rownames(track) = sample.names
  track1= track[, -7]
  write.table(track1, paste0(loopnr,"tracking_table.txt"), sep="\t")
 
  #write.table(sort(sample_sums(merged)), "reads.per.sample.txt", sep="\t")
  print (paste(c("forward_read_names=", forward_read_names,
                 "seed=",seed,
                 "maxeef=", maxeef,
                 "trainingset=", trainingset, 
                 "truncLenF=", truncLenF, 
                 "trainingset=", trainingset,
                 "sptrainingset=", sptrainingset, 
                 "trimleft=", trimleft, 
                 "tree=", tree, 
                 "multithread=", multithread)))
}


tests=seq(from = 60, to = 200, by = 10)

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

######Part 1 is done! ######

sink()

