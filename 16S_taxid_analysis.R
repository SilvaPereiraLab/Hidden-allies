#DATA ANALYSIS OF AVS TABLE OBTAINED FROM DADA2 PIPELINE

###
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(DECIPHER); packageVersion("DECIPHER")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(msa); packageVersion("msa")
library(phangorn); packageVersion("phangorn")
library(ggtree); packageVersion("ggtree")
library(scales); packageVersion("scales")
library(ape)
library(tidytree)
library(aplot)
library(tidyverse)
library(Hmisc)
library(tidyr)
library(tibble)
library(ggcorrplot)
library(mia)
library(decontam); packageVersion("decontam")
library(scater)
library(ggsignif) #plot significance level to plot
library(gridExtra)
library(ggtreeExtra)
library(microbiome)
library(dendsort)
library(vegan)
library(pheatmap)
library(ggpubr)
###

set.seed(10)
#Importing dataframes
tax_table <- read.csv("~\\SequencingV3V4\\R_analysis\\all_raw\\tax_table.csv", header=TRUE, row.names="ASV", sep =";")
refseq <- read.csv("~\\SequencingV3V4\\R_analysis\\all_raw\\refseq.csv", header=TRUE, row.names="ASV", sep =";")
abund_table <- read.csv("~\\SequencingV3V4\\R_analysis\\all_raw\\abund_table_all.csv", header=TRUE, row.names="ID", sep =";")
sam_data <- read.csv("~\\SequencingV3V4\\R_analysis\\all_raw\\sam_data.csv", header=TRUE, row.names="Code", sep =";")


#Sequences information 
dna <- Biostrings::DNAStringSet(refseq$Sequence)
names(dna) <- row.names(refseq)

#Construction of phylogenetic tree using phangron package
mult <- msa(dna, method="ClustalW", type="dna", order="input") #alignment of sequences
phang.align <- as.phyDat(mult, type="DNA")

#Model test
#mt <- modelTest(phang.align, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
                #control = pml.control(trace = 0))
#infer a phylogenetic tree with best found model
fitGTRGI <- pml_bb(phang.align, model="GTR+G(4)+I")
plot(fitGTRGI)
#mid-point-rooted maximum likelihood phylogenetic tree
midpoint_fitGTRGI = midpoint(fitGTRGI$tree)
plot(midpoint_fitGTRGI)
add.scale.bar()



#Creating matrix with all information to use in phyloseq 
asv_info <- phyloseq(otu_table(as.matrix(abund_table), taxa_are_rows = FALSE), 
                     sample_data(sam_data), 
                     tax_table(as.matrix(tax_table)),
                     phy_tree(midpoint_fitGTRGI))
asv_info <- merge_phyloseq(asv_info, dna)
asv_info

#Save dataframe as RData
save(asv_info, file="16S_antibio_taxid_table_dada2_tree_final_all_abs_sequences.RData")
load("~\\SequencingV3V4\\R_analysis\\all_raw\\16S_antibio_taxid_table_dada2_tree_final_all.RData")


#Change sample data 
sample_data(asv_info) <- sam_data
#Clean ASV table names
tax_table(asv_info)[, colnames(tax_table(asv_info))] <- gsub(tax_table(asv_info)[, colnames(tax_table(asv_info))], 
                                                 pattern = "Methylobacterium-Methylorubrum", replacement = "MM")
tax_table(asv_info)[, colnames(tax_table(asv_info))] <- gsub(tax_table(asv_info)[, colnames(tax_table(asv_info))], 
                                                             pattern = "Burkholderia-Caballeronia-Paraburkholderia", replacement = "BCP")




#####################################################Clean samples#####################
#Delete the ASVs present in NC from samples
asv_NC = subset_samples(asv_info, Type == 'NC')
asv_NC = filter_taxa(asv_NC, function(x) sum(x) > 0, TRUE)

#Remove Taxa of NC in isolates samples 
prune_negatives = function(physeq, negs, samps) {
  negs.n1 = prune_taxa(taxa_sums(negs)>=1, negs) 
  samps.n1 = prune_taxa(taxa_sums(samps)>=1, samps) 
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(negs.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samps.n1))
}

asv_clean = prune_negatives(asv_info, asv_NC, asv_info)


#Subset conditions:
#Positive cotrol
asv_PC = subset_samples(asv_clean, Type == 'PC')
asv_PC = filter_taxa(asv_PC, function(x) sum(x) > 0, TRUE)
#Aspergillus fimugatus isolates
asv_Afu = subset_samples(asv_clean, Type == 'Afu Isolate')
asv_Afu = filter_taxa(asv_Afu, function(x) sum(x) > 0, TRUE)
#Other Aspergillus isolates
asv_no_Afu = subset_samples(asv_clean, Type == 'Other Isolate')
asv_no_Afu = filter_taxa(asv_no_Afu, function(x) sum(x) > 0, TRUE)
#all Aspergillus and PC
asv_asp = subset_samples(asv_clean, Type %in% c('Afu Isolate', 'Other Isolate', 'PC'))
asv_asp = filter_taxa(asv_asp, function(x) sum(x) > 0, TRUE)

############################PLOTS##############################

############################MAIN FIGURES#######################

#Fumigatus Samples 
#Group OTUs by specie classification
afu_genus <- tax_glom(asv_Afu, taxrank="genus", NArm = TRUE)
afu_genus = filter_taxa(afu_genus, function(x) sum(x) > 0, TRUE)
###Look to 30 most abundant ASVs
top40 <- names(sort(taxa_sums(afu_genus), decreasing=TRUE))[1:40]
asv.top40 <- prune_taxa(top40, afu_genus)

# Generate a custom color palette with 40 diverse colors
library(RColorBrewer)
n <- 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_select = sample(col_vector, n) #select randomly the colors
#pie(rep(1,n), col=sample(col_vector, n))
#Different colors (divergent)
div_col_pals <- brewer.pal.info[brewer.pal.info$category == 'div', ]
col_vector <- unlist(mapply(brewer.pal, div_col_pals$maxcolors, rownames(div_col_pals)))
col_vector_select <- sample(col_vector, n)
#pie(rep(1,n), col = col_vector_select) 

colors = c("#B2ABD2","#9970AB","#35978F","#BABABA","#BF812D","#F7F7F7","#053061",
                    "#1B7837","#D53E4F","#74ADD1","#878787","#A50026","#FDDBC7",
                    "#E6F5D0","#E08214","#7F3B08","#2D004B","#01665E","#F1B6DA",
                    "#FEE08B","#00441B","#F46D43","#4393C3","#FEE0B6","#DE77AE",
                    "#D8DAEB","#313695","#C2A5CF","#FEE08B","#F7F7a4","#A6DBA0",
                    "#D73027","#67001F","#543039","#80CDC1","#543008","#C7EAE5",
                    "#FDAE61","#1A9850","#762A83", "#D9F0D3")
                    

# Create the bar plot using the custom color palette
afu_genus_bar <- plot_bar(afu_genus, fill = "genus", x = 'Name') + 
  geom_bar(aes(fill = genus), stat = "identity", 
           position = "stack", width=1, colour=NA,size=0) +
  # Use the custom color palette
  scale_color_manual(values = colors, na.value = "#CACFD2") + 
  scale_fill_manual(values = colors, na.value = "#CACFD2") +
  theme(legend.key.size = unit(0.5, "cm"),
        #axis.text.x=element_blank(),
        #strip.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size =10, color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 11, face = 'bold', vjust = 2),
        axis.title.x = element_text(size = 11, face = 'bold', vjust = -3),
        legend.title = element_text(size = 11, face = "bold", vjust = 1.5),
        legend.text = element_text(size = 12)) + 
  ylab('Relative Abundance (%)') +
  guides(fill = guide_legend(ncol = 2))

afu_genus_bar


#Core heatmap
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/core-microbiota.html

#Use the format_to_besthit function from microbiomeutilities to get the best classification of the ASVs
afu_genus.f <- aggregate_taxa(afu_genus, "genus")
afu_genus.f <- subset_taxa(afu_genus.f, genus!="Unknown")
# Core with compositionals:
prevalences <- c(0, 1)
detections = c(0.1, 0.5, 1, 5, 15, 30)


p.core <- plot_core(afu_genus.f, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(8, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = 0.75) + 
  xlab("Detection Threshold 
       (Relative Abundance (%))")
#asthetics 
p.core + theme(axis.text.y = element_text(size=10, color = 'black'),
               axis.text.x = element_text(size = 12, color = 'black', angle = 90, hjust=1,vjust=0.3),
               axis.title.y = element_text(size=11, face = 'bold'),
               legend.text = element_text(size =10, colour = 'black'), 
               legend.title = element_text(size = 11, colour = 'black', face = 'bold')) +
  scale_fill_gradientn(name = "Prevalence", colors = rev(brewer.pal(8, "Spectral")),
                       breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  coord_flip() +
  scale_y_discrete(position = 'left')


###############################SUPP FIGURES#############################
#Negative control
#Group OTUs by specie classification
nc_genus <- tax_glom(asv_NC, taxrank="genus", NArm = FALSE)
sample_data(nc_genus) <- data.frame(NC_type = c('Amplification', 'Amplification', 'Extraction', 'Extraction'),
                                    Sample = c('16SCN1', '16SCN2','16SCN3','16SCN4'),
                                    row.names = 'Sample')

tax_table(nc_genus)[, colnames(tax_table(nc_genus))] <- gsub(tax_table(nc_genus)[, colnames(tax_table(nc_genus))], 
                                                             pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", replacement = "ANPR")

##Heatmap
nc_genus_df <- psmelt(nc_genus)

nc_genus_heatmap <- ggplot(nc_genus, aes(x = Sample, y=genus, fill = Abundance)) +
  geom_tile(color = "#283747",lwd = 0.1,linetype = 1)  + 
  #facet_wrap(~ NC_type, scales = "free") + 
  scale_fill_gradientn(name = "Relative Abundance (%)", 
                       colours= c("#FDFEFE", "#7FB3D5", "#A9DFBF", "#F9E79F"), 
                       na.value = "transparent", breaks=c(0,5,15,25,35,45),
                       values = c(0, rescale(5, from = range(nc_genus_df$Abundance)), 1))+
  theme(strip.background=element_rect(fill="#ebebeb"), 
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold", vjust = 1.8),
        legend.key.size = unit(0.8, "cm"),
        axis.text.x = element_text(size =11, color = "black", angle = 90),
        axis.text.y = element_text(size= 11, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=11))

nc_genus_heatmap

##Bar plot
n <- 20
#Different colors (divergent)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_select = sample(col_vector, n, replace = FALSE)

colors = c("#1A9850","#FDAE61","#92C5DE","#878787","#A6DBA0","#D6604D","#C7EAE5",
                    "#DFC27D","#FDDBC7","#ABD9E9","#5AAE61","#7FBC41","#D9F0D3",
                    "#FFFFBF","#DE77AE","#9970AB","#003C30","#F4A582","#74ADD1","#D1E5F0")

nc_genus_bar = plot_bar(nc_genus, fill="genus") + 
  geom_bar(stat = "identity", position = "stack", 
           width=0.9, colour=NA,size=0) + facet_wrap(~ NC_type, scales = "free") + 
  scale_color_manual(values = col_vector_select, na.value = "#CACFD2") + 
  scale_fill_manual(values = col_vector_select, na.value = "#CACFD2") +
  theme(legend.key.size = unit(0.5, "cm"),
        strip.text = element_text(size = 11, face = 'bold'),
        axis.text.x = element_text(size =10, color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 11, face = 'bold', vjust = 2),
        axis.title.x = element_text(size = 11, face = 'bold', vjust = -1),
        legend.title = element_text(size = 11, face = "bold", vjust = 1.5),
        legend.text = element_text(size = 12)) + 
  ylab('Relative Abundance (%)') +
  guides(fill=guide_legend(ncol=1))

nc_genus_bar



####Other Aspergillus
#select 3 fumigatus isolates to compare

asperg <- prune_samples(sample_names(asv_asp) %in% c("16SAf01", "16SAf02", "16SAf03", "16SAfla04", "16SAn05", "16SAt25"), asv_asp)

#Group OTUs by genus classification (deleting NA)
asp_genus <- tax_glom(asperg, taxrank="genus", NArm = TRUE)


# Generate a custom color palette with 50 diverse colors

colors = c("#B2ABD2","#9970AB","#35978F","#BABABA","#BF812D","#F7F7F7","#053061",
                    "#1B7837","#D53E4F","#74ADD1","#878787","#A50026","#FDDBC7",
                    "#E6F5D0","#E08214","#7F3B08","#2D004B","#01665E","#F1B6DA",
                    "#FEE08B","#00441B","#F46D43","#4393C3","#FEE0B6","#DE77AE",
                    "#D8DAEB","#313695","#C2A5CF","#FEE08B","#F7F7a4","#A6DBA0",
                    "#D73027","#67001F","#543039","#80CDC1","#543008","#C7EAE5",
                    "#FDAE61","#1A9850","#762A83", "#D9F0D3")
                    

# Create the bar plot using the custom color palette
asp_genus_bar <- plot_bar(asp_genus, fill = "genus", x = 'Name') + 
  geom_bar(aes(fill = genus), stat = "identity", 
           position = "stack", width=1, colour=NA,size=0) + facet_wrap(~Type, scales = "free") +
  ylim(0, 75) +
  # Use the custom color palette
  scale_color_manual(values = colors, na.value = "#CACFD2") + 
  scale_fill_manual(values = colors, na.value = "#CACFD2") +
  theme(legend.key.size = unit(0.5, "cm"),
        #axis.text.x=element_blank(),
        #strip.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size =10, color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 11, face = 'bold', vjust = 2),
        axis.title.x = element_text(size = 11, face = 'bold', vjust = -3),
        legend.title = element_text(size = 11, face = "bold", vjust = 1.5),
        legend.text = element_text(size = 12)) + 
  ylab('Relative Abundance (%)') +
  guides(fill = guide_legend(ncol = 2))

asp_genus_bar

#####Core bacteriome##############
###Tree of the core bacteriome 


asv_core = subset_taxa(asv_Afu, genus %in% c("Gemmata", "BCP", "Bradyrhizobium", "Bryobacter", "Puia", "Ralstonia", "Edaphobacter",
                                               "Caulobacter", "Sediminibacterium", "MM", "Nevskia", "Acidibacter", "Reyranella", "Acidipila", 
                                               "Afipia", "Asinibacterium", "Granulicella", "Fimbriimonas", "Labrys"))

##Getting info to add
asv_core_df <- psmelt(asv_core)
#write.csv(asv_core_df, "~\\SequencingV3V4\\R_analysis\\all_raw\\asv_core_df.csv")

# Filter out rows where Abundance = 0
asv_core_df_mean <- asv_core_df %>%
  filter(Abundance != 0)

# Remove the columns Sample, Name, and Type
asv_core_df_mean <- select(asv_core_df_mean, -Sample, -Name, -Type)

asv_core_df_mean <- asv_core_df_mean %>%
  group_by(OTU) %>%
  mutate(Mean_Abundance = mean(Abundance), Sample_Count = n())

# Remove the columns Sample, Name, and Type
asv_core_df_mean <- select(asv_core_df_mean, -Abundance)
# Remove duplicated rows while keeping the first occurrence
asv_core_df_mean <- distinct(asv_core_df_mean, .keep_all = TRUE)

asv_core_df_mean <- data.frame(asv_core_df_mean, row.names = "OTU")

asv_core_genus <- asv_core_df_mean[,"genus", drop = FALSE] 

##Creating tree plot
#create the tree
core_log10 <- transform(asv_core, 'log10')


core_tree = ggtree(asv_core, ladderize = FALSE, size = 1, layout = "circular") +
  geom_tiplab(align =TRUE, size = 0) +
  geom_treescale(fontsize=3, linesize=1, offset=0.4)
  

core_tree <- rotate_tree(core_tree, -90)
print(core_tree)

colors = c("#B2ABD2","#9970AB","#35978F","#67001F","#BF812D","#D9F0D3","#053061",
                    "#1B7837","#D53E4F","#74ADD1","#878787","#A50026","#FDDBC7",
                    "#F1B6DA","#FEE08B","#FDAE61","#2D004B","#01665E","#D8DAEB")

genus_core = gheatmap(core_tree, asv_core_genus, offset=0, width=.08,
                 colnames_angle=0, colnames_offset_y = 0) +
  theme(legend.text = element_text(size = 11, color = 'black'),
        legend.key.size = unit(0.5, "cm"))+
  guides(fill = guide_legend(ncol = 1)) + 
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(ncol = 1))
  
genus_core


melt_simple <- psmelt(core_z) %>%
  select(OTU, val=Abundance)

tree_log10 = genus_core + geom_fruit(data=melt_simple, geom=geom_boxplot,
                       mapping = aes(y=OTU, x=val, group=label,fill=genus), 
                       pwidth = 0.8, size=.2, outlier.size=0.8, outlier.stroke=0.1, offset = 0.4,
                       outlier.shape=21, axis.params=list(axis = "log10", text.size  = 2,
                                                          hjust = 1,vjust  = 0.5, nbreak = 2,),
                       grid.params=list()) 


    
core_tree_log10 <- tree_log10 + 
  scale_fill_manual(values = colors) +
  theme(legend.position = "none")

core_tree_log10



########################################CORRELATIONAL ANALYSIS##########################################
######Heatmap

###Uisng all bacteriome
#Updating the sample data with virulance and MICs
sam_data_afu <- read.csv("~\\SequencingV3V4\\R_analysis\\all_raw\\CSV original\\sam_data_afu.csv", header=TRUE, row.names="Code", sep =";")
sample_data(asv_Afu) <- sam_data_afu

#group to genus level
asv_Afu_genus <- tax_glom(asv_Afu, taxrank="genus", NArm = TRUE)
#Transform data to log10 score
asv_Afu_log10  <- transform(asv_Afu_genus, 'log10')

#Transform data to a dataframe 
asv_Afu_log10_df <- psmelt(asv_Afu_log10) #transform to dataframe
asv_Afu_log10_df_na = replace(asv_Afu_log10_df, asv_Afu_log10_df == 0, NA)


#dataframe for distance
asv_Afu_log10_df_adundance_na = asv_Afu_log10_df_na[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)

asv_Afu_log10_df_adundance_na = data.frame(asv_Afu_log10_df_adundance_na , row.names = "Name")
#dataframe for plot
asv_Afu_log10_df_adundance = asv_Afu_log10_df[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)

asv_Afu_log10_df_adundance = data.frame(asv_Afu_log10_df_adundance , row.names = "Name")

#Hierarchical clustering using Bray-Curtis distances

cluster_sapmles <- hclust(vegdist(asv_Afu_log10_df_adundance, method="bray", na.rm = TRUE))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
cluster_sapmles <- sort_hclust(cluster_sapmles)
plot(cluster_sapmles, main = "Sorted Dendrogram", xlab = "", sub = "")
#check number of possible clusters by heaigth 
clusters = cutree(cluster_sapmles, k = NULL, h = 0.2)
max(clusters)
#Get the isolates virulance and resistance data
isolates_data <- sam_data_afu[, (ncol(sam_data_afu) - 4):ncol(sam_data_afu)]
isolates_data <- data.frame(isolates_data, row.names = "Name")

#Define the colors
my_colour <- list(Infection.capacity = c("High" = "#4A235A", "Medium-high" = "#7D3C98", 
                                         "Medium-low" = "#BB8FCE", "Low" = "#EBDEF0"),
  Amphotericin.B = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"),
  Voriconazole = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"),
  Posaconazole = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"))



#The plot
pheatmap(asv_Afu_log10_df_adundance_na, na_col = "transparent", 
         cluster_cols=FALSE, cluster_rows = cluster_sapmles,
         annotation_row = isolates_data, 
         color = colorRampPalette(c("#EBF5FB", "#7FB3D5", "#A9DFBF", "#F9E79F", "#F5CBA7", "#EB984E"))(100),
         annotation_colors = my_colour,
         fontsize=10.5, border_color= "white",
         cutree_rows = 11)



###Excluding the core bacteriome
# Create a vector of genera to exclude
core_bacteriome <- c("Gemmata", "BCP", "Bradyrhizobium", "Bryobacter", "Puia", "Ralstonia", "Edaphobacter",
                       "Caulobacter", "Sediminibacterium", "MM", "Nevskia", "Acidibacter", "Reyranella", "Acidipila", 
                       "Afipia", "Asinibacterium", "Granulicella", "Fimbriimonas", "Labrys")
#Create vector with all genera
tax_table <- data.frame(asv_Afu@tax_table)
genus_names <- unique(tax_table$genus)
#Exclude core from all genera
# Exclude values from list1 that are in list2
genera_exc_core <- genus_names[!genus_names %in% core_bacteriome]

#phyloseq without core
asv_exc_core <-  subset_taxa(asv_Afu, genus %in% genera_exc_core)

####Correlation
#group to genus level
asv_Afu_genus <- tax_glom(asv_exc_core, taxrank="genus", NArm = TRUE)
#Transform data to log10
asv_Afu_log10  <- transform(asv_Afu_genus, 'log10')

#Transform data to a dataframe 
asv_Afu_log10_df <- psmelt(asv_Afu_log10) #transform to dataframe
asv_Afu_log10_df_na = replace(asv_Afu_log10_df, asv_Afu_log10_df == 0, NA)


#dataframe for distance
asv_Afu_log10_df_adundance_na = asv_Afu_log10_df_na[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)

asv_Afu_log10_df_adundance_na = data.frame(asv_Afu_log10_df_adundance_na , row.names = "Name")
#dataframe for plot
asv_Afu_log10_df_adundance = asv_Afu_log10_df[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)

asv_Afu_log10_df_adundance = data.frame(asv_Afu_log10_df_adundance , row.names = "Name")

#Hierarchical clustering using Bray-Curtis distances

cluster_sapmles <- hclust(vegdist(asv_Afu_log10_df_adundance, method="bray", na.rm = TRUE))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
cluster_sapmles <- sort_hclust(cluster_sapmles)
plot(cluster_sapmles, main = "Sorted Dendrogram", xlab = "", sub = "")

clusters = cutree(cluster_sapmles, k = NULL, h = 0.4)
max(clusters)
#Get the isolates virulance and resistance data
isolates_data <- sam_data_afu[, (ncol(sam_data_afu) - 4):ncol(sam_data_afu)]
isolates_data <- data.frame(isolates_data, row.names = "Name")

#Define the colors
my_colour <- list(Infection.capacity = c("High" = "#4A235A", "Medium-high" = "#7D3C98", 
                                         "Medium-low" = "#BB8FCE", "Low" = "#EBDEF0"),
                  Amphotericin.B = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"),
                  Voriconazole = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"),
                  Posaconazole = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"))



#The plot
pheatmap(asv_Afu_log10_df_adundance_na, na_col = "transparent", 
         cluster_cols=FALSE, cluster_rows = cluster_sapmles,
         annotation_row = isolates_data, 
         color = colorRampPalette(c("#EBF5FB", "#7FB3D5", "#A9DFBF", "#F9E79F", "#F5CBA7", "#EB984E"))(100),
         annotation_colors = my_colour,
         fontsize=10.5, border_color= "white")


#####Spearman’s correlation analysis 
##Creating the two dataframes
#group to genus level
asv_Afu_genus <- tax_glom(asv_Afu, taxrank="genus", NArm = TRUE)

#Transform data to a dataframe 
asv_Afu_df <- psmelt(asv_Afu_genus) #transform to dataframe

#abundance dataframe
asv_Afu_df_adundance = asv_Afu_df[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)


#traits dataframe
traits_df <- read.csv("~\\SequencingV3V4\\R_analysis\\all_raw\\CSV original\\mic_vir.csv", header=TRUE, sep =";")

#Join the two dataframes
corr_df  <- merge(asv_Afu_df_adundance, traits_df, by = "Name")
corr_df <- data.frame(corr_df, row.names = "Name")

corr_df <- corr_df %>%
  mutate_if(is.character, as.numeric)

#Calcul of the correlation matrix and p-value
corr_spearman <- rcorr(as.matrix(corr_df), type = "spearman") #dataframe needs to be transformed to matrix
#clean correlation matrix
corr_matrix <- corr_spearman$r #correlation matrix
corr_matrix <- corr_matrix[, !(names(corr_df) %in% c("Infection.capacity","Amphotericin.B","Voriconazole","Posaconazole"))]
corr_matrix <- corr_matrix[c("Infection.capacity","Amphotericin.B","Voriconazole","Posaconazole"), ]
#clean p-value matrix
p_value <- corr_spearman$P #p-value matrix
p_value <- p_value[, !(names(corr_df) %in% c("Infection.capacity","Amphotericin.B","Voriconazole","Posaconazole"))]
p_value <- p_value[c("Infection.capacity","Amphotericin.B","Voriconazole","Posaconazole"), ]
#clean p-value matrix


corrplot(corr_matrix, method = "color", p.mat = p_value, 
         sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", 
         tl.col = 'black', pch.cex = 0.9, 
         col = colorRampPalette(c("#1F618D","#6D9EC1", "white", "#D98880", "#922B21"))(100),
         tl.cex = 0.95, cl.pos = "b", cl.ratio = 0.6) 


#####Alpha diversity vs infection capacity

##Alpha and Beta diversities 
#To use mia package the phyloseq object need to me converted 
asv_tse <- makeTreeSummarizedExperimentFromPhyloseq(asv_Afu) #convert phyloseq to TSE
asv_tse

#Calculate alpha diversity
asv_tse <- mia::estimateRichness(asv_tse, abund_values = "counts", index = "observed", name="richness")
asv_tse  <- mia::estimateDiversity(asv_tse,index = "shannon", name = "shannon")
asv_tse  <- mia::estimateEvenness(asv_tse , abund_values = "counts", index="pielou")
#head(colData(asv_tse)) #check dataframe
#save as data frame
alpha_df <- data.frame(colData(asv_tse))
alpha_df$Infection.capacity <- factor(alpha_df$Infection.capacity  , levels=c("High", "Medium-high","Medium-low", "Low"))

# Define the variables
variables <- c("High", "Medium-high","Medium-low", "Low")
# Generate all possible pairwise combinations
comparisons <- combn(variables, 2, simplify = FALSE)

#plot boxplots
richness_plt <- ggplot(alpha_df, aes(x=Infection.capacity, y=richness, fill=Infection.capacity)) + 
  geom_boxplot(add = "jitter") +
  scale_fill_manual(values=c("High" = "#512E5F", "Medium-high" = "#884EA0", 
                               "Medium-low" = "#C39BD3", "Low" = "#EBDEF0")) +
  #stat_compare_means(comparisons = comparisons, label = "p.signif", hide.ns = T) +
  stat_compare_means(label.y = 75) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 11, face = 'bold', vjust = 2),
        axis.title.x = element_text(size = 11, face = 'bold', vjust = -3),
        #legend.title = element_text(size = 11, face = "bold", vjust = 1.5),
        #legend.text = element_text(size = 12),
        legend.position = "none") + 
  ylab("Richness") + xlab("")
  
richness_plt

#plot boxplots
shannon_plt <- ggplot(alpha_df, aes(x=Infection.capacity, y=shannon, fill=Infection.capacity)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("High" = "#512E5F", "Medium-high" = "#884EA0", 
                             "Medium-low" = "#C39BD3", "Low" = "#EBDEF0")) +
  #stat_compare_means(comparisons = comparisons, label = "p.signif", hide.ns = F) +
  stat_compare_means(label.y = 3.4) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 11, face = 'bold', vjust = 2),
        axis.title.x = element_text(size = 11, face = 'bold', vjust = -3),
        #legend.title = element_text(size = 11, face = "bold", vjust = 1.5),
        #legend.text = element_text(size = 12),
        legend.position = "none") + 
  ylab("Shannon Index") + xlab("")

shannon_plt

#plot boxplots
pielou_plt <- ggplot(alpha_df, aes(x=Infection.capacity, y=pielou, fill=Infection.capacity)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("High" = "#512E5F", "Medium-high" = "#884EA0", 
                             "Medium-low" = "#C39BD3", "Low" = "#EBDEF0")) +
  #stat_compare_means(comparisons = comparisons, label = "p.signif", hide.ns = F) +
  stat_compare_means(label.y = 0.84) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 11, face = 'bold', vjust = 2),
        axis.title.x = element_text(size = 11, face = 'bold', vjust = -3)) + 
  ylab("Evenness") + xlab("")

pielou_plt


aplha_div_plot <- grid.arrange(shannon_plt, richness_plt, pielou_plt, nrow=3)
aplha_div_plot

