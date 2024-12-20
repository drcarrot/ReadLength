
lapply(c('ggfortify', 'psych','msa', 'ggplot2','picante', 'RColorBrewer','phangorn', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn', 'DECIPHER', 'hillR'), require,character.only=TRUE) #add as necessary

print("A log of processing is saved as log.txt")
sink(file="log.analyses.txt")

files <- list.files(pattern="tracking_table.txt", recursive=TRUE)

tracking_all= NULL

#merge all tracking tables
for (f in 1:length(files)) {
   dat <- read.table(files[f], header=T, sep="\t", na.strings="", row.names=NULL)
   dat$forward.trim=paste0(gsub(".tracking_table.txt", "", files[f]),"0")
   tracking_all <- rbind(tracking_all, dat)
}

write.table(tracking_all, "tracking_all.txt")

# Get a list of all files in the directory that end in "merged.RDS"
file_list <- list.files(pattern = "merged.rds")

####rarefaction_depth list#######
rare_depth <- matrix(ncol=3)
colnames(rare_depth)<-c("samplename", "depth", "ntaxa_OG")
# Loop through each file in the list
for (file in file_list) {
  # Load the phyloseq object
  physeq <- readRDS(file)
  # Extract the sample name
  sample_name <- paste0(gsub(".merged.rds", "", file),"0")
  ####RAREFACTION####
  rare_depth <- rbind(rare_depth, c(sample_name, min(sample_sums(physeq)), ntaxa(physeq)))
}

rare_depth <- rare_depth[-1,]
write.table(rare_depth, "rare.depth.txt")

rdepth=read.table("rare.depth.txt")

# Create an empty matrix to store the results
alpha.results_matrix <- matrix(ncol = 4)
colnames(alpha.results_matrix)<-c("samplename", "richness", "simpson", "FaithsPD")
rarefaction_depth <- matrix(ncol=3)
colnames(rarefaction_depth)<-c("samplename", "depth", "ntaxa")
unclassifieds <- matrix(ncol = 7)
colnames(unclassifieds)<-c("samplename", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
classifieds <- matrix(ncol = 7)
colnames(classifieds)<-c("samplename", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
percent_matrix <- matrix(ncol = 7)
colnames(percent_matrix)<-c("samplename", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
horn_list <- list()
sorensen_list <- list()
wunifrac_list <- list()
unifrac_list <- list()
alpha.results_matrix2 <- matrix(ncol = 4)
colnames(alpha.results_matrix2)<-c("samplename", "richness",  "simpsons", "FaithsPD")
beta.results_matrix <- matrix(ncol = 9)
colnames(beta.results_matrix)<-c("samplename", "hornR",  "hornP", "SorensenR", "SorensenP", "wUnifracR", "wUnifracP", "UnifracR", "UnifracP")
variance.beta.b <- matrix(ncol = 3)
colnames(variance.beta.b)<-c("samplename", "pre",  "post")
variance.beta.s <- matrix(ncol = 3)
colnames(variance.beta.s)<-c("samplename", "pre",  "post")
variance.beta.w <- matrix(ncol = 3)
colnames(variance.beta.w)<-c("samplename", "pre",  "post")
variance.beta.u <- matrix(ncol = 3)
colnames(variance.beta.u)<-c("samplename", "pre",  "post")

#Add metadata
meta=read.table("Meta.txt", header=TRUE, row.names=1)
row.names(meta)=gsub(".fastq.gz", "", row.names(meta))
meta$Time_since_dist[meta$Time_since_dist < 0] <- -1


# Loop through each file in the list
for (file in file_list) {
  # Load the phyloseq object
  physeq <- readRDS(file)
  sample_names(physeq)=gsub(".fastq.gz", "", sample_names(physeq))
  meta=meta[match(sample_names(physeq), row.names(meta)),]
  physeq@sam_data=sample_data(meta)
  physeq=prune_samples(physeq@sam_data$Time_since_dist<3, physeq)
  #physeq@sam_data$Time_since_dist=as.factor(physeq@sam_data$Time_since_dist)

  
  # Extract the sample name
  sample_name <- paste0(gsub(".merged.rds", "", file),"0")

  ####RAREFACTION####
  physeq=rarefy_even_depth(physeq, rngseed = 1, sample.size = (min(rdepth[,2])))
  rarefaction_depth <- rbind(rarefaction_depth, c(sample_name, min(rdepth[,2]), ntaxa(physeq)))
  
  ###DISTANCE MATRICES####
  # Compute horn distances
  distances.b <- vegdist(physeq@otu_table@.Data, method="horn")
  Rb=(adonis2(distances.b~physeq@sam_data$Time_since_dist))$R2[1]
  Pb=(adonis2(distances.b~physeq@sam_data$Time_since_dist))$'Pr(>F)'[1]
  disp.b=betadisper(distances.b, physeq@sam_data$Time_since_dist)$group.distances
  
  distances.s <- vegdist(physeq@otu_table@.Data, binary = TRUE)
  Rs=(adonis2(distances.s~physeq@sam_data$Time_since_dist))$R2[1]
  Ps=(adonis2(distances.s~physeq@sam_data$Time_since_dist))$'Pr(>F)'[1]
  disp.s=betadisper(distances.s, physeq@sam_data$Time_since_dist)$group.distances

if (is.null(physeq@phy_tree)) {
  Rw=NA
  Pw=NA
  disp.w=NA
  } else {
  distances.w <- UniFrac(physeq, weighted=TRUE)
  Rw=(adonis2(distances.w~physeq@sam_data$Time_since_dist))$R2[1]
  Pw=(adonis2(distances.w~physeq@sam_data$Time_since_dist))$'Pr(>F)'[1]
  disp.w=betadisper(distances.w, physeq@sam_data$Time_since_dist)$group.distances
}

if (is.null(physeq@phy_tree)) {
  Ru=NA
  Pu=NA
  disp.w=NA
  } else {
  distances.u <- UniFrac(physeq, weighted=FALSE)
  Ru=(adonis2(distances.w~physeq@sam_data$Time_since_dist))$R2[1]
  Pu=(adonis2(distances.w~physeq@sam_data$Time_since_dist))$'Pr(>F)'[1]
  disp.u=betadisper(distances.u, physeq@sam_data$Time_since_dist)$group.distances
}
  
  beta.results_matrix  <- rbind(beta.results_matrix , c(sample_name, Rb, Pb, Rs, Ps, Rw, Pw, Ru, Pu))
  variance.beta.b<- rbind(variance.beta.b ,c(sample_name, disp.b))
  variance.beta.s<- rbind(variance.beta.s ,c(sample_name, disp.s))
  variance.beta.w<- rbind(variance.beta.w ,c(sample_name, disp.w))
  variance.beta.u<- rbind(variance.beta.u ,c(sample_name, disp.u))

  # Add the distance matrix to the list
  horn_list[[file]] <- distances.b
  sorensen_list[[file]]<-distances.s
  wunifrac_list[[file]]<-distances.w
  unifrac_list[[file]]<-distances.u


  
  ####ALPHA#######
  # Calculate mean richness and inverse Simpson diversity in all samples
  richness <-(hill_taxa(physeq@otu_table, q = 0))
  inv_simpson <-(hill_taxa(physeq@otu_table, q = 2))
  
if (is.null(physeq@phy_tree)) {
  Faith=NA
  physeq@sam_data$Faith=NA
  } else {
   Faith<- pd(physeq@otu_table,physeq@phy_tree, include.root=FALSE)$PD
   physeq@sam_data$Faith <- pd(physeq@otu_table,physeq@phy_tree, include.root=FALSE)$PD
}
  physeq@sam_data$richness <- (hill_taxa(physeq@otu_table, q = 0))
  physeq@sam_data$inv_simpson <- (hill_taxa(physeq@otu_table, q = 2))
  
  if (is.null(physeq@phy_tree)) {
  faithw=NA
 } else {
    faithw=wilcox.test(physeq@sam_data$Faith ~ physeq@sam_data$Time_since_dist)$statistic
}

  richw=wilcox.test(physeq@sam_data$richness ~ physeq@sam_data$Time_since_dist)$statistic
  invw=wilcox.test(physeq@sam_data$inv_simpson ~ physeq@sam_data$Time_since_dist)$statistic

  # Add the results for undisturbed samples to the matrix
  physeq=prune_samples(physeq@sam_data$Time_since_dist<1, physeq)
  alpha.results_matrix <- rbind(alpha.results_matrix, as.matrix(data.frame(sample_name, physeq@sam_data$richness ,  physeq@sam_data$inv_simpson,  physeq@sam_data$Faith)))
  
  alpha.results_matrix2  <- rbind(alpha.results_matrix2, c(sample_name, richw, invw, faithw))
  
  ####CLASSIFIEDS####
  # Count the number of non-NA's at each taxonomic level
  counts <- c(sample_name, sum(!is.na(physeq@tax_table[,1])),sum(!is.na(physeq@tax_table[,2])), sum(!is.na(physeq@tax_table[,3])), sum(!is.na(physeq@tax_table[,4])), sum(!is.na(physeq@tax_table[,5])), sum(!is.na(physeq@tax_table[,6])))
  
  # Add the results to the matrix
  classifieds <- rbind(classifieds, counts)
  
  ####UNCLASSIFIEDS####
  # Count the number of NA's at each taxonomic level
  na_counts <- c(sample_name, sum(is.na(physeq@tax_table[,1])),sum(is.na(physeq@tax_table[,2])), sum(is.na(physeq@tax_table[,3])), sum(is.na(physeq@tax_table[,4])), sum(is.na(physeq@tax_table[,5])), sum(is.na(physeq@tax_table[,6])))
  
  # Add the results to the matrix
  unclassifieds <- rbind(unclassifieds, na_counts)
  
  # Count the number of unclassified taxa at each taxonomic level
  unclassified_counts <- c(sample_name, 
			if (all(is.na(physeq@tax_table[,1]))) {result <- 1
			} else {1-mean(sample_sums(subset_taxa(physeq, Kingdom!="NA"))/min(sample_sums(physeq)))},
                        if (all(is.na(physeq@tax_table[,2]))) {result <- 1
			} else {1-mean(sample_sums(subset_taxa(physeq, Phylum!="NA"))/min(sample_sums(physeq)))},
                        if (all(is.na(physeq@tax_table[,3]))) {result <- 1
			} else {1-mean(sample_sums(subset_taxa(physeq, Class!="NA"))/min(sample_sums(physeq)))},
                        if (all(is.na(physeq@tax_table[,4]))) {result <- 1
			} else {1-mean(sample_sums(subset_taxa(physeq, Order!="NA"))/min(sample_sums(physeq)))},
                        if (all(is.na(physeq@tax_table[,5]))) {result <- 1
			} else {1-mean(sample_sums(subset_taxa(physeq, Family!="NA"))/min(sample_sums(physeq)))},
                        if (all(is.na(physeq@tax_table[,6]))) {result <- 1
			} else {1-mean(sample_sums(subset_taxa(physeq, Genus!="NA"))/min(sample_sums(physeq)))})
         
  # Add the results to the matrix
  percent_matrix <- rbind(percent_matrix, unclassified_counts)
  }

# Remove the first row of the matrix (which is empty)
rarefaction_depth <- rarefaction_depth[-1,]

# Remove the first row of the matrix (which is empty)
alpha.results_matrix <- alpha.results_matrix[-1,]
alpha.results_matrix2 <- alpha.results_matrix2[-1,]
unclassifieds <- unclassifieds[-1,]
percent_matrix <- percent_matrix[-1,]
classifieds <- classifieds[-1,]
unclassifieds <-unclassifieds[-1,]
percent_matrix <- percent_matrix[-1,]
beta.results_matrix <- beta.results_matrix[-1,]
variance.beta.b <- variance.beta.b[-1,]
variance.beta.s <- variance.beta.s[-1,]
variance.beta.w <- variance.beta.w[-1,]
variance.beta.u <- variance.beta.u[-1,]


write.table(rarefaction_depth, "rarefaction.depth2.TXT")
write.table(alpha.results_matrix, "alpha.results.TXT")
write.table(alpha.results_matrix2, "alpha.results2.TXT")
write.table(unclassifieds, "unclassifieds2.TXT")
write.table(classifieds, "classifieds2.TXT")
write.table(percent_matrix, "unclassified.percent2.TXT")
saveRDS(horn_list, "horn.list.RDS")
saveRDS(sorensen_list, "sorensen.list.RDS")
saveRDS(wunifrac_list, "wunifrac.list.RDS")
saveRDS(unifrac_list, "unifrac.list.RDS")
write.table(beta.results_matrix, "beta.results.diff.TXT")
write.table(variance.beta.b, "variance.beta.b.TXT")
write.table(variance.beta.s, "variance.beta.s.TXT")
write.table(variance.beta.w, "variance.beta.w.TXT")
write.table(variance.beta.u, "variance.beta.u.TXT")

sink()

