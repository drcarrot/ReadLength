library(phyloseq)
library(hillR)
library(vegan)

####rarefaction_depths#######

rare_depth <- matrix(ncol=2)
colnames(rare_depth)<-c("samplename", "depth")
# Loop through each file in the list
for (file in file_list) {
  # Load the phyloseq object
  physeq <- readRDS(file)
  # Extract the sample name
  sample_name <- paste0(gsub(".merged.rds", "", file),"0")
  ####RAREFACTION####
  physeq=rarefy_even_depth(physeq, rngseed = 1, sample.size = min(sample_sums(physeq)))
  rare_depth <- rbind(rare_depth, c(sample_name, min(sample_sums(physeq)))) }

rare_depth <- rare_depth[-1,]
write.table(rare_depth, "rare.depth.txt")

rdepth=read.table("rare.depth.txt")

# Create an empty matrix to store the results
alpha.results_matrix <- matrix(ncol = 3)
colnames(alpha.results_matrix)<-c("samplename", "richness", "simpson")
rarefaction_depth <- matrix(ncol=3)
colnames(rarefaction_depth)<-c("samplename", "depth", "ntaxa")
unclassifieds <- matrix(ncol = 7)
colnames(unclassifieds)<-c("samplename", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
classifieds <- matrix(ncol = 7)
colnames(classifieds)<-c("samplename", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
percent_matrix <- matrix(ncol = 7)
colnames(percent_matrix)<-c("samplename", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
bray_list <- list()
sorensen_list <- list()

#Add metadata
meta=read.table("Meta.txt", header=TRUE)
meta$Time_since_dist[meta$Time_since_dist < 0] <- -1

# Get a list of all files in the directory that end in "merged.RDS"
file_list <- list.files(pattern = "merged.rds")

# Loop through each file in the list
for (file in file_list) {
  # Load the phyloseq object
  physeq <- readRDS(file)
  
  # Extract the sample name
  sample_name <- paste0(gsub(".merged.rds", "", file),"0")

  ####RAREFACTION####
  physeq=rarefy_even_depth(physeq, rngseed = 1, sample.size = min(min(rdepth[,2])))
  rarefaction_depth <- rbind(rarefaction_depth, c(sample_name, min(rdepth[,2]), ntaxa(physeq)))
  
 ####DISTANCE MATRICES####
  # Compute Bray-Curtis distances
  distances.b <- vegdist(physeq@otu_table@.Data)
  distances.s <- vegdist(physeq@otu_table@.Data, binary = TRUE)
  
  # Add the distance matrix to the list
  bray_list[[file]] <- distances.b
  sorensen_list[[file]]<-distances.s

 meta=meta[match(physeq@sam_data$names, row.names(meta)),]
  physeq@sam_data=sample_data(meta)
  physeq=prune_samples(physeq@sam_data$Time_since_dist<1, physeq)

  ####ALPHA#######
  # Calculate mean richness and inverse Simpson diversity
  richness <-(hill_taxa(physeq@otu_table, q = 0))
  inv_simpson <-(hill_taxa(physeq@otu_table, q = 2))
  
  # Add the results to the matrix
  alpha.results_matrix <- rbind(alpha.results_matrix, as.matrix(data.frame(sample_name, richness, inv_simpson)))
  
  ####CLASSIFIEDS####
  # Count the number of non-NA's at each taxonomic level
  counts <- c(sample_name, sum(!is.na(physeq@tax_table[,1])),sum(!is.na(physeq@tax_table[,2])), sum(!is.na(physeq@tax_table[3])), sum(!is.na(physeq@tax_table[,4])), sum(!is.na(physeq@tax_table[,5])), sum(!is.na(physeq@tax_table[,6])))
  
  # Add the results to the matrix
  classifieds <- rbind(classifieds, counts)
  

  ####UNCLASSIFIEDS####
  # Count the number of NA's at each taxonomic level
  na_counts <- c(sample_name, sum(is.na(physeq@tax_table[,1])),sum(is.na(physeq@tax_table[,2])), sum(is.na(physeq@tax_table[,3])), sum(is.na(physeq@tax_table[,4])), sum(is.na(physeq@tax_table[,5])), sum(is.na(physeq@tax_table[,6])))
  
  # Add the results to the matrix
  unclassifieds <- rbind(unclassifieds, na_counts)
  

  # Count the number of unclassified taxa at each taxonomic level
  unclassified_counts <- c(sample_name, 
                           1-mean(sample_sums(subset_taxa(physeq, Kingdom!="NA"))/min(sample_sums(physeq))),
                           1-mean(sample_sums(subset_taxa(physeq, Phylum!="NA"))/min(sample_sums(physeq))),
                           1-mean(sample_sums(subset_taxa(physeq, Class!="NA"))/min(sample_sums(physeq))),
                           1-mean(sample_sums(subset_taxa(physeq, Order!="NA"))/min(sample_sums(physeq))),
                           1-mean(sample_sums(subset_taxa(physeq, Family!="NA"))/min(sample_sums(physeq))),
                           1-mean(sample_sums(subset_taxa(physeq, Genus!="NA"))/min(sample_sums(physeq))))
         
  # Add the results to the matrix
  percent_matrix <- rbind(percent_matrix, unclassified_counts)
  
 
}

# Remove the first row of the matrix (which is empty)
rarefaction_depth <- rarefaction_depth[-1,]

# Remove the first row of the matrix (which is empty)
alpha.results_matrix <- alpha.results_matrix[-1,]

# Remove the first row of the matrix (which is empty)
unclassifieds <- unclassifieds[-1,]
percent_matrix <- percent_matrix[-1,]
classifieds <- classifieds[-1,]
row.names(classifieds) <-classifieds[,1]
row.names(unclassifieds) <-unclassifieds[,1]
row.names(percent_matrix) <- percent_matrix[,1]

write.table(rarefaction_depth, "rarefaction.depth2.txt")
write.table(alpha.results_matrix, "alpha.results2.txt")
write.table(unclassifieds, "unclassifieds2.txt")
write.table(classifieds, "classifieds2.txt")
write.table(percent_matrix, "unclassified.percent2.txt")
saveRDS(bray_list, "bray.list2.rds")
saveRDS(sorensen_list, "sorensen.list2.rds")
