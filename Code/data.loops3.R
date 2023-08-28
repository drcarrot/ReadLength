library(phyloseq)
library(hillR)
library(vegan)

rdepth=read.table("rarefaction.depth.txt")

# Create an empty matrix to store the results
alpha.results_matrix <- matrix(ncol = 3)
colnames(alpha.results_matrix)<-c("samplename", "richness",  "simpsons")

beta.results_matrix <- matrix(ncol = 5)
colnames(beta.results_matrix)<-c("samplename", "BrayR",  "BrayP", "SorensenR", "SorensenP")

# Get a list of all files in the directory that end in "merged.RDS"
file_list <- list.files(pattern = "merged.rds")

#Add metadata
meta=read.table("Meta.txt", header=TRUE)
meta$Time_since_dist[meta$Time_since_dist < 0] <- -1

# Loop through each file in the list
for (file in file_list) {
  # Load the phyloseq object
  physeq <- readRDS(file)
  
  # Extract the sample name
  sample_name <- paste0(gsub(".merged.rds", "", file),"0")
  
  ####RAREFACTION####
  physeq=rarefy_even_depth(physeq, rngseed = 1, sample.size = min(min(rdepth[,2])))
  meta=meta[match(physeq@sam_data$names, row.names(meta)),]
  physeq@sam_data=sample_data(meta)
  physeq=prune_samples(physeq@sam_data$Time_since_dist<2, physeq)
  
  ####ALPHA#######
  # Calculate mean richness and inverse Simpson diversity
  physeq@sam_data$richness <- (hill_taxa(physeq@otu_table, q = 0))
  physeq@sam_data$inv_simpson <- (hill_taxa(physeq@otu_table, q = 2))

  richp=wilcox.test(physeq@sam_data$richness ~ physeq@sam_data$Time_since_dist)$p.value
  invp=wilcox.test(physeq@sam_data$inv_simpson ~ physeq@sam_data$Time_since_dist)$p.value
 
alpha.results_matrix  <- rbind(alpha.results_matrix , c(sample_name, richp, invp))
 
  
####DISTANCE MATRICES####
  # Compute Bray-Curtis distances
  distances.b <- vegdist(physeq@otu_table@.Data)
 Rb=(adonis2(distances.b~physeq@sam_data$Time_since_dist))$R2[1]
 Pb=(adonis2(distances.b~physeq@sam_data$Time_since_dist))$'Pr(>F)'[1]

  distances.s <- vegdist(physeq@otu_table@.Data, binary = TRUE)
 Rs=(adonis2(distances.s~physeq@sam_data$Time_since_dist))$R2[1]
 Ps=(adonis2(distances.s~physeq@sam_data$Time_since_dist))$'Pr(>F)'[1]

beta.results_matrix  <- rbind(beta.results_matrix , c(sample_name, Rb, Pb, Rs, Ps))

}

# Remove the first row of the matrix (which is empty)
alpha.results_matrix <- alpha.results_matrix[-1,]
beta.results_matrix <- beta.results_matrix[-1,]

write.table(alpha.results_matrix, "alpha.results.diff.txt")
write.table(beta.results_matrix, "beta.results.diff.txt")
