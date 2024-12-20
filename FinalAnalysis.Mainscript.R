setwd(dir = "../Revision/")

lapply(c('ggfortify','Kendall', 'colorspace', 'scales','msa','dada2', 'tidyr','ggplot2','picante', 'RColorBrewer','phangorn', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn', 'DECIPHER', 'hillR'), require,character.only=TRUE) #add as necessary



key=data.frame(Study=c("Jurburg_Pig" , "Kennedy" , "Venkataraman", "Dong", "Overholt", "Qian", "Fuentes", "Jurburg_microcosm",  "Kruistuim" ), ID= c(1,2,3,4,5,6,7,8,9), Environment=c("Animal", "Animal", "Animal", "Water", "Water","Water","Soil", "Soil", "Soil"))
key$Environment=as.factor(key$Environment)
key$Color=c("#FF9966","#FF6633", "#CC3300","#99CCFF", "#6699CC", "#003366","#CCFF99","#99CC66","#336600")

key$Environment=factor(key$Environment, levels = c("Animal", "Water", "Soil"))
key$ID=factor(key$ID, levels = c("1", "2", "3","4","5","6","7","8","9"))
#animal_colors <- c("#FF9966","#FF6633", "#CC3300")  # 3 shades of peach
#water_colors <- c("#99CCFF", "#6699CC", "#003366")# 3 shades of blue
#soil_colors <- c ("#CCFF99","#99CC66","#336600")# 3 shades of green

# Combine colors based on the Environment
#key$Color=0
#key[which(key$Environment=="Water"),]$Color=water_colors
#key[which(key$Environment=="Animal"),]$Color=animal_colors
#key[which(key$Environment=="Soil"),]$Color=soil_colors

#Add rarefaction depth to each study
key$rarefaction <- NA
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name
  subdir <- key$Study[i]
  
  # Construct the file path for rare.depth.txt
  file_path <- file.path(subdir, "rare.depth.txt")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the rare.depth.txt file as a data frame
    raredepth <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Calculate the minimum depth
    min_depth <- min(raredepth$depth, na.rm = TRUE)
    
    # Add the minimum depth to the rarefaction column in the key data frame
    key$rarefaction[i] <- min_depth
  } else {
    warning(paste("File not found:", file_path))
  }
}



write.table(key, "key.txt")
####################################

# Initialize an empty data frame for the master table
tracking.master <- data.frame()

# Loop through each subdirectory in key$Study
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and corresponding metadata from key
  subdir <- key$Study[i]
  study_id <- key$ID[i]
  study_color <- key$Color[i]
  study_environment <- key$Environment[i]
  
  # Construct the file path for tracking_all.txt
  file_path <- file.path(subdir, "tracking_all.txt")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the tracking_table.txt file as a data frame
    tracking_table <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Add metadata columns to the tracking_table
    tracking_table$Study <- subdir
    tracking_table$ID <- study_id
    tracking_table$Color <- study_color
    tracking_table$Environment <- study_environment
    
    # Append to the master table
    tracking.master <- rbind(tracking.master, tracking_table)
  } else {
    warning(paste("File not found:", file_path))
  }
}

write.table(tracking.master, "tracking.master.txt")

tracking.master.=tracking.master%>% group_by(forward.trim, Study,  Environment, ID) %>%
  summarize(mean=mean(tabled/Original.Reads*100),
            sd=sd(100*tabled/Original.Reads))

results.processing <- tracking.master%>%
  group_by(ID) %>%
  summarise(
    Kendall_tau = Kendall(tabled/Original.Reads*100, forward.trim)$tau,
    p_value = Kendall(tabled/Original.Reads*100, forward.trim)$sl
  )


F1.1=ggplot(tracking.master., aes(x=forward.trim, y=mean, color=ID))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color=ID),width = 0.2, alpha=0.5)+
  xlab("Read Length")+
  scale_color_manual(values = key$Color)+
  ylab("% Preserved Reads")+
  facet_wrap(~Environment, nrow = 1)+
  theme_minimal()

ggsave(
  filename = "F1.1.pdf", 
  plot = F1.1,          
  width = 8.27,               # A4 width in inches (210mm)
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
)

####################Number of Bacterial reads
ggplot(tracking.master, aes(x=forward.trim, y=100*BacterialReads/Original.Reads, color=ID))+
  geom_point()+
  xlab("Read Length")+
  ylab("% Preserved Reads")+
  facet_wrap(~Environment, nrow = 1)+
  scale_color_manual(values = key$Color)+
  theme_minimal()

ggsave(
  filename = "S3.pdf",  #
  plot = S3,           
  width = 8.27,               # A4 width in inches (210mm)
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
)

##############rarefaction depth##########
# Initialize an empty data frame for the master table
rare.master <- data.frame()

# Loop through each subdirectory in key$Study
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and corresponding metadata from key
  subdir <- key$Study[i]
  study_id <- key$ID[i]
  study_color <- key$Color[i]
  study_environment <- key$Environment[i]
  
  # Construct the file path for rare.depth.TXT
  file_path <- file.path(subdir, "rare.depth.TXT")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the rare.depth.TXT file as a data frame
    rare.depth <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Add metadata columns to the tracking_table
    rare.depth$Study <- subdir
    rare.depth$ID <- study_id
    rare.depth$Color <- study_color
    rare.depth$Environment <- study_environment
    
    # Append to the master table
    rare.master <- rbind(rare.master, rare.depth)
  } else {
    warning(paste("File not found:", file_path))
  }
}

write.table(rare.master, "rare.master.txt")

results.processing2 <- rare.master%>%
  group_by(ID) %>%
  summarise(
    Kendall_tau = Kendall(depth, samplename)$tau,
    p_value = Kendall(depth, samplename)$sl
  )

F1.2=ggplot(rare.master, aes(x=samplename, y=depth, color=ID))+
  geom_point()+
  geom_line()+
  xlab("Read Length")+
  scale_color_manual(values = key$Color)+
  ylab("Max. rarefaction depth")+
  facet_wrap(~Environment, nrow = 1)+
  theme_minimal()

ggsave(
  filename = "F1.2.pdf", 
  plot = F1.2,          
  width = 8.27,               # A4 width in inches (210mm)
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
)


################# percent unclassified reads
# Initialize an empty data frame for the master table
unclassified.percent.master <- data.frame()

# Loop through each subdirectory in key$Study
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and file path
  subdir <- key$Study[i]
  file_path <- file.path(subdir, "unclassified.percent2.TXT")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the unclassified.percent2.TXT file as a data frame
    unclassified_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, row.names=NULL)
    
    # Add columns for Study, ID, Color, and Environment
    unclassified_data <- unclassified_data %>%
      mutate(
        Study = key$Study[i],
        ID = key$ID[i],
        Color = key$Color[i],
        Environment = key$Environment[i]
      )
    
    # Combine with the master table
    unclassified.percent.master <- bind_rows(unclassified.percent.master, unclassified_data)
  } else {
    warning(paste("File not found:", file_path))
  }
}


write.table(unclassified.percent.master, "unclassified.percent.master.txt")

unclassified.percent.master. =unclassified.percent.master %>%  pivot_longer(cols = c(Kingdom, Phylum, Class,Order, Family,Genus), 
                                                                            names_to = "Taxonomy", 
                                                                            values_to = "Proportion.unclassified")


unclassified.percent.master.$Taxonomy=as.factor(unclassified.percent.master.$Taxonomy)

unclassified.percent.master.$Taxonomy =  factor(unclassified.percent.master.$Taxonomy, levels=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'))

results.taxonomy <- unclassified.percent.master.%>%
  group_by(ID, Taxonomy) %>%
  summarise(
    Kendall_tau = Kendall(Proportion.unclassified, samplename)$tau,
    p_value = Kendall(Proportion.unclassified, samplename)$sl
  )

results.taxonomy=results.taxonomy[which(results.taxonomy$Taxonomy!="Kingdom"),]

results.taxonomy2 <- unclassified.percent.master.%>%
  group_by(Environment, Taxonomy, samplename) %>%
  summarise(mean(Proportion.unclassified)
  )

F2=ggplot(unclassified.percent.master., aes(x=as.numeric(samplename), y=Proportion.unclassified*100, color=Taxonomy))+
  geom_line()+  
  xlab("Read Length")+
  ylab("Unclassified taxa (%)")+
  viridis::scale_color_viridis(discrete=TRUE)+
  theme_minimal()+
  facet_wrap(~factor(unclassified.percent.master.$ID, levels = c("1", "4", "7","2","5","8","3","6","9")), ncol=3)+
  theme(legend.position = "none")
ggsave(
      filename = "F2.pdf",  
      plot = F2,          
      width = 8.27,               # A4 width in inches (210mm)
      height = 8.27,
      dpi = 300,                  # High resolution (dots per inch)
      units = "in"                # Specify dimensions in inches
    )
    

##################Figure 2.1 Alpha Diversity
# Initialize an empty data frame for the master table
alpha <- data.frame()

# Loop through each subdirectory in key$Study
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and file path
  subdir <- key$Study[i]
  file_path <- file.path(subdir, "alpha.results.TXT")
  
  # Check if the file exists
  if (file.exists(file_path)) {
       alpha_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
    
    # Add columns for Study, ID, Color, and Environment
    alpha_data <- alpha_data %>%
      mutate(
        Study = key$Study[i],
        ID = key$ID[i],
        Color = key$Color[i],
        Environment = key$Environment[i]
      )
    
    # Combine with the master table
    alpha <- bind_rows(alpha, alpha_data)
  } else {
    warning(paste("File not found:", file_path))
  }
}


write.table(alpha, "alpha.master.txt")

alpha. =alpha %>% pivot_longer(cols = c(richness, simpson, FaithsPD), 
                               names_to = "index", 
                               values_to = "alpha")
alpha.$index=as.factor(alpha.$index)

alpha.$index=factor(alpha.$index, levels = c("richness", "simpson", "FaithsPD"))

alpha..=alpha.%>% group_by(samplename, Study, index, Environment, ID) %>%
  summarize(mean=mean(alpha),
            sd=sd(alpha)  )

results.alpha <- alpha..%>%
  group_by(Environment, ID,index) %>%
  summarise( Kendall_tau = Kendall(mean, samplename)$tau,
             p_value = Kendall(mean, samplename)$sl
  )

F3=ggplot(alpha.., aes(x=as.numeric(samplename), y=mean, color=ID))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color=ID),width = 0.2, alpha=0.5)+
  xlab("Read Length")+
  ylab("alpha diversity")+
  scale_color_manual(values = key$Color)+
  theme_minimal()+
  facet_wrap(index~Environment, ncol=3, scales = "free_y")

ggsave(
  filename = "F3.pdf",  
  plot = F3,          
  width = 8.27,               # A4 width in inches (210mm)
  height = 8.27,
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
)


##################Figure 2.2 Alpha Diversity
# Initialize an empty data frame for the master table
alpha2 <- data.frame()

# Loop through each subdirectory in key$Study
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and file path
  subdir <- key$Study[i]
  file_path <- file.path(subdir, "alpha.results2.TXT")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the unclassified.percent2.TXT file as a data frame
    alpha_data2 <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Add columns for Study, ID, Color, and Environment
    alpha_data2 <- alpha_data2 %>%
      mutate(
        Study = key$Study[i],
        ID = key$ID[i],
        Color = key$Color[i],
        Environment = key$Environment[i]
      )
    
    # Combine with the master table
    alpha2 <- bind_rows(alpha2, alpha_data2)
  } else {
    warning(paste("File not found:", file_path))
  }
}


write.table(alpha2, "alpha.master2.txt")

alpha2. =alpha2 %>%   pivot_longer(cols = c(richness, simpsons, FaithsPD), 
                                   names_to = "index", 
                                   values_to = "wstatistic")
alpha2.$alpha=as.factor(alpha2.$index)

results.alpha2 <- alpha2.%>%
  group_by(Environment, ID,index) %>%
  summarise( Kendall_tau = Kendall(wstatistic, samplename)$tau,
             p_value = Kendall(wstatistic, samplename)$sl
  )

S2=ggplot(alpha2., aes(x=as.numeric(samplename), y=wstatistic, color=ID))+
  geom_line()+  
  xlab("Read Length")+
  ylab("w-statistic")+
  scale_color_manual(values = key$Color)+
  theme_minimal()+
  facet_wrap(index~Environment, ncol=3, scales = "free_y")
ggsave(
  filename = "S2.pdf",  
  plot = S2,          
  width = 8.27,               # A4 width in inches (210mm)
  height = 8.27,
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
)


##################permanovas in Beta Diversity
# Initialize an empty data frame for the master table
betadiff.master <- data.frame()

# Loop through each subdirectory in key$Study
for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and file path
  subdir <- key$Study[i]
  file_path <- file.path(subdir, "beta.results.diff.TXT")
  
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the unclassified.percent2.TXT file as a data frame
    betadiff <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Add columns for Study, ID, Color, and Environment
    betadiff <- betadiff %>%
      mutate(
        Study = key$Study[i],
        ID = key$ID[i],
        Color = key$Color[i],
        Environment = key$Environment[i]
      )
    
    # Combine with the master table
    betadiff.master <- bind_rows(betadiff.master, betadiff)
  } else {
    warning(paste("File not found:", file_path))
  }
}

write.table(betadiff.master, "betadiff.master.txt")

 betadiff.master. =betadiff.master %>%   pivot_longer(cols = c(hornR, SorensenR,  wUnifracR, UnifracR), 
                                                      names_to = "index", 
                                                      values_to = "R2")
 betadiff.master.$index=as.factor(betadiff.master.$index)
 
 betadiff.master.$index=factor(betadiff.master.$index, levels = c("hornR", "SorensenR", "wUnifracR", "UnifracR"))
 
 F5=ggplot(betadiff.master., aes(x=as.numeric(samplename), y=R2, color=ID))+
   geom_point()+
   geom_line()+
   xlab("Read Length")+
   ylab("R2")+
   scale_color_manual(values = key$Color)+
   theme_minimal()+
   facet_wrap(index~Environment, ncol=3)
 ggsave(
   filename = "F5.pdf",  
   plot = F5,          
   width = 8.27,               # A4 width in inches (210mm)
   height = 8.27,
   dpi = 300,                  # High resolution (dots per inch)
   units = "in"                # Specify dimensions in inches
 )

 #Now, plot the p values
 
 betadiff.master.. =betadiff.master %>%   pivot_longer(cols = c(hornP, SorensenP,  wUnifracP, UnifracP), 
                                                       names_to = "index", 
                                                       values_to = "P")
 betadiff.master..$index=as.factor(betadiff.master..$index)
 
 betadiff.master..$index=factor(betadiff.master..$index, levels = c("hornP", "SorensenP", "wUnifracP", "UnifracP"))
 
 results.betadiff <-betadiff.master.%>%
   group_by(Environment, ID,index) %>%
   summarise( Kendall_tau = Kendall(R2, samplename)$tau,
              p_value = Kendall(R2, samplename)$sl
   )
 
 results.betadiff2 <-betadiff.master..%>%
   group_by(Environment, ID,index) %>%
   summarise( Kendall_tau = Kendall(P, samplename)$tau,
              p_value = Kendall(P, samplename)$sl
   )
 
 S4=ggplot(betadiff.master.., aes(x=as.numeric(samplename), y=P, color=ID))+
   geom_point()+
   geom_line()+
   xlab("Read Length")+
   ylab("P")+
   scale_color_manual(values = key$Color)+
   theme_minimal()+
   facet_wrap(index~Environment, ncol=3)+
   scale_y_log10() 
 ggsave(
   filename = "S4.pdf",  
   plot = S4,          
   width = 8.27,               # A4 width in inches (210mm)
   height = 8.27,
   dpi = 300,                  # High resolution (dots per inch)
   units = "in"                # Specify dimensions in inches
 ) 

 
 ##############-Variance in beta diversity###########3
 ##################Figure 2.1 Alpha Diversity
 # Initialize an empty data frame for the master table
 betavar.b <- data.frame()
 betavar.s <- data.frame()
 betavar.u <- data.frame()
 betavar.w <- data.frame()
 
 #MORISITA
 # Loop through each subdirectory in key$Study
 for (i in seq_len(nrow(key))) {
   # Get the subdirectory name and file path
   subdir <- key$Study[i]
   file_path <- file.path(subdir, "variance.beta.b.TXT")
   
   # Check if the file exists
   if (file.exists(file_path)) {
     betavar <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
     
     # Add columns for Study, ID, Color, and Environment
     betavar <- betavar %>%
       mutate(
         Study = key$Study[i],
         ID = key$ID[i],
         Color = key$Color[i],
         Environment = key$Environment[i],
         index="Morisita"
       )
     
     # Combine with the master table
     betavar.b <- bind_rows(betavar.b, betavar)
   } else {
     warning(paste("File not found:", file_path))
   }
 }
 
 #SORENSEN
 # Loop through each subdirectory in key$Study
 for (i in seq_len(nrow(key))) {
   # Get the subdirectory name and file path
   subdir <- key$Study[i]
   file_path <- file.path(subdir, "variance.beta.s.TXT")
   
   # Check if the file exists
   if (file.exists(file_path)) {
     betavar <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
     
     # Add columns for Study, ID, Color, and Environment
     betavar <- betavar %>%
       mutate(
         Study = key$Study[i],
         ID = key$ID[i],
         Color = key$Color[i],
         Environment = key$Environment[i],
         index="Sorensen"
       )
     
     # Combine with the master table
     betavar.s <- bind_rows(betavar.s, betavar)
   } else {
     warning(paste("File not found:", file_path))
   }
 }
 
 #wUnifrac
 # Loop through each subdirectory in key$Study
 for (i in seq_len(nrow(key))) {
   # Get the subdirectory name and file path
   subdir <- key$Study[i]
   file_path <- file.path(subdir, "variance.beta.w.TXT")
   
   # Check if the file exists
   if (file.exists(file_path)) {
     betavar <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
     
     # Add columns for Study, ID, Color, and Environment
     betavar <- betavar %>%
       mutate(
         Study = key$Study[i],
         ID = key$ID[i],
         Color = key$Color[i],
         Environment = key$Environment[i],
         index="wUnifrac"
       )
     
     # Combine with the master table
     betavar.w <- bind_rows(betavar.w, betavar)
   } else {
     warning(paste("File not found:", file_path))
   }
 }
 
 #Unifrac
 # Loop through each subdirectory in key$Study
 for (i in seq_len(nrow(key))) {
   # Get the subdirectory name and file path
   subdir <- key$Study[i]
   file_path <- file.path(subdir, "variance.beta.u.TXT")
   
   # Check if the file exists
   if (file.exists(file_path)) {
     betavar <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
     
     # Add columns for Study, ID, Color, and Environment
     betavar <- betavar %>%
       mutate(
         Study = key$Study[i],
         ID = key$ID[i],
         Color = key$Color[i],
         Environment = key$Environment[i],
         index="Unifrac"
       )
     
     # Combine with the master table
     betavar.u <- bind_rows(betavar.u, betavar)
   } else {
     warning(paste("File not found:", file_path))
   }
 }
 
 betavar.master=rbind(betavar.b,betavar.s,betavar.u,betavar.w)
 
 write.table(betavar.master, "betavar.master.txt")
 
betavar.master$index=as.factor(betavar.master$index)
 
betavar.master$index=factor(betavar.master$index, levels = c("Morisita", "Sorensen", "wUnifrac", "Unifrac"))
 

F4=ggplot(betavar.master, aes(x=as.numeric(samplename), y=pre, color=ID))+
  geom_point()+
  geom_line()+
  xlab("Read Length")+
  ylab("variance")+
  scale_color_manual(values = key$Color)+
  theme_minimal()+
  facet_wrap(index~Environment, ncol=3)+
  scale_y_log10() 
ggsave(
  filename = "F4.pdf",  
  plot = F4,          
  width = 8.27, 
  height = 8.27,  # A4 width in inches (210mm)
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
) 

#Here we look at variance post disturbance/variance pre disturbance
S3=ggplot(betavar.master, aes(x=as.numeric(samplename), y=post/pre, color=ID))+
  geom_point()+
  geom_line()+
  xlab("Read Length")+
  ylab("variance")+
  scale_color_manual(values = key$Color)+
  theme_minimal()+
  facet_wrap(index~Environment, ncol=3)+
  scale_y_log10() 
ggsave(
  filename = "S3.pdf",  
  plot = S3,          
  width = 8.27,               # A4 width in inches (210mm)
  height = 8.27,
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
) 


results.betavar <-betavar.master%>%
  group_by(Environment, ID,index) %>%
  summarise( Kendall_tau = Kendall(pre, samplename)$tau,
             p_value = Kendall(pre, samplename)$sl
  )

#####Figure 6 MANTEL TESTS

mantels=data.frame()

for (i in seq_len(nrow(key))) {
  # Get the subdirectory name and file path
  subdir <- key$Study[i]
  file_path.s <- file.path(subdir, "sorensen.list.RDS")
  file_path.h <- file.path(subdir, "horn.list.RDS")
  file_path.w <- file.path(subdir, "wunifrac.list.RDS")
  file_path.u <- file.path(subdir, "unifrac.list.RDS")

  horns <- readRDS(file_path.h)
  sorensens <- readRDS(file_path.s)
  unifracs <- readRDS(file_path.u)
  wunifracs <- readRDS(file_path.w)

  names(horns)=paste0(gsub(".merged.rds", "", names(horns)),"0")
  horn.mantel= matrix(ncol = 7, nrow = length(horns))
  for (j in 1:length(horns)) {
    colnames(horn.mantel)<-c("samplename", "r", "p", "index", "ID", "Color", "Environment")
    horn.mantel[j,1]=names(horns[j])
    horn.mantel[j,2]=mantel(horns[[j]], horns[[as.character(max(as.numeric(names(horns))))]], method="pearson")$statistic
    horn.mantel[j,3]=mantel(horns[[j]], horns[[as.character(max(as.numeric(names(horns))))]], method="pearson")$signif
    horn.mantel[j,4]="Morisita-Horn"
    horn.mantel[j,5]= key$ID[i]
    horn.mantel[j,6]= key$Color[i]
    horn.mantel[j,7]=key$Environment[i]
  }
  names(sorensens)=paste0(gsub(".merged.rds", "", names(sorensens)),"0")
  sorensen.mantel= matrix(ncol = 7, nrow = length(sorensens))
  for (j in 1:length(sorensens)) {
    colnames(sorensen.mantel)<-c("samplename", "r", "p", "index", "ID", "Color", "Environment")
    sorensen.mantel[j,1]=names(sorensens[j])
    sorensen.mantel[j,2]=mantel(sorensens[[j]], sorensens[[as.character(max(as.numeric(names(sorensens))))]], method="pearson")$statistic
    sorensen.mantel[j,3]=mantel(sorensens[[j]], sorensens[[as.character(max(as.numeric(names(sorensens))))]], method="pearson")$signif
    sorensen.mantel[j,4]="sorensen"
    sorensen.mantel[j,5]= key$ID[i]
    sorensen.mantel[j,6]= key$Color[i]
    sorensen.mantel[j,7]=key$Environment[i]
  }
  names(unifracs)=paste0(gsub(".merged.rds", "", names(unifracs)),"0")
  unifrac.mantel= matrix(ncol = 7, nrow = length(unifracs))
  for (j in 1:length(unifracs)) {
    colnames(unifrac.mantel)<-c("samplename", "r", "p", "index", "ID", "Color", "Environment")
    unifrac.mantel[j,1]=names(unifracs[j])
    unifrac.mantel[j,2]=mantel(unifracs[[j]], unifracs[[as.character(max(as.numeric(names(unifracs))))]], method="pearson")$statistic
    unifrac.mantel[j,3]=mantel(unifracs[[j]], unifracs[[as.character(max(as.numeric(names(unifracs))))]], method="pearson")$signif
    unifrac.mantel[j,4]="Unifrac"
    unifrac.mantel[j,5]= key$ID[i]
    unifrac.mantel[j,6]= key$Color[i]
    unifrac.mantel[j,7]=key$Environment[i]
  }
  names(wunifracs)=paste0(gsub(".merged.rds", "", names(wunifracs)),"0")
  wunifrac.mantel= matrix(ncol = 7, nrow = length(wunifracs))
  for (j in 1:length(wunifracs)) {
    colnames(wunifrac.mantel)<-c("samplename", "r", "p", "index", "ID", "Color", "Environment")
    wunifrac.mantel[j,1]=names(wunifracs[j])
    wunifrac.mantel[j,2]=mantel(wunifracs[[j]], wunifracs[[as.character(max(as.numeric(names(wunifracs))))]], method="pearson")$statistic
    wunifrac.mantel[j,3]=mantel(wunifracs[[j]], wunifracs[[as.character(max(as.numeric(names(wunifracs))))]], method="pearson")$signif
    wunifrac.mantel[j,4]="Weighted unifrac"
    wunifrac.mantel[j,5]= key$ID[i]
    wunifrac.mantel[j,6]= key$Color[i]
    wunifrac.mantel[j,7]=key$Environment[i]
  }
  # Combine with the master table
 mantels<- bind_rows(mantels, as.data.frame(horn.mantel),
                     as.data.frame(sorensen.mantel), 
                     as.data.frame(unifrac.mantel), 
                     as.data.frame(wunifrac.mantel))
}

mantels$r=as.numeric(mantels$r)
mantels$p=as.numeric(mantels$p)
write.table(mantels, "mantels.txt")

F6=ggplot(mantels, aes(x=as.numeric(samplename), y=r, color=ID))+
  geom_point()+
  geom_line()+
  xlab("Read Length")+
  ylab("Mantel's R")+
  scale_color_manual(values = key$Color)+
  theme_minimal()+
  facet_wrap(index~Environment, ncol=3)

ggsave(
  filename = "F6.pdf",  
  plot = F6,          
  width = 8.27,               # A4 width in inches (210mm)
  height=8.27,
  dpi = 300,                  # High resolution (dots per inch)
  units = "in"                # Specify dimensions in inches
) 

results.mantel <-mantels%>%
  group_by(Environment, ID,index) %>%
  summarise( Kendall_tau = Kendall(r, samplename)$tau,
             p_value = Kendall(r, samplename)$sl
  )

results.mantel2 <-mantels%>%
  group_by(Environment, ID,index) %>%
  summarise( Kendall_tau = Kendall(p, samplename)$tau,
             p_value = Kendall(p, samplename)$sl
  )
############FIGURE S1####################
##Note: this only works if there is one fastq.gz file per directory!
# subdirs <-paste0("./",key$Study)
# 
# # Initialize a list to store plots
# all_plots <- vector("list", length = nrow(key))
# 
# # Loop through each subdirectory
# for (i in seq_len(nrow(key))) {
#   subdir <- subdirs[i]
#   id <- key$ID[i]
#   fastq_files <- list.files(path = subdir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
#   all_plots[[id]] <- plotQualityProfile(fastq_files) + 
#       ggtitle(paste0("ID: ", id)) +
#       theme(plot.title = element_text(size = 10))
#   }
# 
# # Combine the plots into a grid with the desired layout
# grid_plot <- plot_grid(plotlist = all_plots[order(key$Environment)], ncol = 3, align = "v")


