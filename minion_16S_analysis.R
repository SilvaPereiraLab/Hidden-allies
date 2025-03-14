####16S-full length sequencing data analysis 
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
library(scater)
library(ggsignif) #plot significance level to plot
library(gridExtra)
library(RColorBrewer)
library(microbiome)
library(vegan)
library(dendsort)
library(ggtreeExtra)
library(ggnewscale)
###


#Load phyloseq object
load("D:\\MinION\\New_Analysis\\all_phyloseq.RData")

genus <- unique(data.frame(all_info@tax_table)$genus)
genus_bact <- genus[!genus %in% c("Ralstonia", "Methylobacterium-Methylorubrum","Burkholderia-Caballeronia-Paraburkholderia")]
#Separete into main groups 
#Ralstonia, Methylobcter, Brukholderia
all_ralstonia <-  subset_taxa(all_info, genus == "Ralstonia")
#phylogenetic agglomeration
h<- 0.5  #define the desired 
agg_ralstonia <- tip_glom(all_ralstonia, h = h)

all_mm <-  subset_taxa(all_info, genus == "Methylobacterium-Methylorubrum")
#phylogenetic agglomeration
h<- 0.3  #define the desired 
agg_mm <- tip_glom(all_mm, h = h)

all_bcp <-  subset_taxa(all_info, genus == "Burkholderia-Caballeronia-Paraburkholderia")
#phylogenetic agglomeration
h<- 0.25  #define the desired 
agg_bcp <- tip_glom(all_bcp, h = h)

#mixed bacteria genus
other_mix <- c("Novosphingobium", "Shewanella","Serratia","Pantoea",
               "Bacillus", "Facklamia", "Planomicrobium", "Aquabacterium",
               "Pelomonas","Delftia","Comamonas")
all_other_mix <- subset_taxa(all_info, genus %in% other_mix)
h <- 0.13
agg_other_mix <- tip_glom(all_other_mix, h = h)

#other not mixed
other <- genus_bact[!genus_bact %in% other_mix]

all_other <- subset_taxa(all_info, genus %in% other)
h <- 0.2
agg_other <- tip_glom(all_other, h = h)

########Join all agg#######
#create datafarmes 
tax_table_01 <- data.frame(agg_ralstonia@tax_table)
abund_table_01 <- data.frame(agg_ralstonia@otu_table)
seq_agg_01 <- data.frame(agg_ralstonia@refseq)
colnames(seq_agg_01)[colnames(seq_agg_01) == "agg_ralstonia.refseq"] <- "Seq"

tax_table_02 <- data.frame(agg_bcp@tax_table)
abund_table_02 <- data.frame(agg_bcp@otu_table)
seq_agg_02 <- data.frame(agg_bcp@refseq)
colnames(seq_agg_02)[colnames(seq_agg_02) == "agg_bcp.refseq"] <- "Seq"

tax_table_03 <- data.frame(agg_mm@tax_table)
abund_table_03 <- data.frame(agg_mm@otu_table)
seq_agg_03 <- data.frame(agg_mm@refseq)
colnames(seq_agg_03)[colnames(seq_agg_03) == "agg_mm.refseq"] <- "Seq"

tax_table_04 <- data.frame(agg_other_mix@tax_table)
abund_table_04 <- data.frame(agg_other_mix@otu_table)
seq_agg_04 <- data.frame(agg_other_mix@refseq)
colnames(seq_agg_04)[colnames(seq_agg_04) == "agg_other_mix.refseq"] <- "Seq"

tax_table_05 <- data.frame(agg_other@tax_table)
abund_table_05 <- data.frame(agg_other@otu_table)
seq_agg_05 <- data.frame(agg_other@refseq)
colnames(seq_agg_05)[colnames(seq_agg_05) == "agg_other.refseq"] <- "Seq"

#join tables
tax_table_all <- rbind(tax_table_03, tax_table_04, tax_table_05, tax_table_01,tax_table_02)

seq_all <- rbind(seq_agg_03, seq_agg_04, seq_agg_05, seq_agg_01, seq_agg_02)


abund_table_all <- bind_rows(abund_table_03, abund_table_04, abund_table_05, abund_table_01,
                             abund_table_02)


abund_table_all[is.na(abund_table_all)] <- 0

#Creating phyloseq object for all

#Sequences information 
dna <- Biostrings::DNAStringSet(seq_all$Seq)
names(dna) <- row.names(seq_all)

all_info_agg <- phyloseq(otu_table(as.matrix(abund_table_all), taxa_are_rows = TRUE),
                     tax_table(as.matrix(tax_table_all)),
                     refseq(dna))

all_info_agg <- merge_phyloseq(all_info_agg, dna)
all_info_agg

#Construction of phylogenetic tree using phangron package
mults_all <- msa(all_info_agg@refseq, method="ClustalW", type="dna", order="input") #alignment of sequences
phang.align <- as.phyDat(mults_all, type="DNA")

#Model test
#mt <- modelTest(phang.align, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
#control = pml.control(trace = 0))

#infer a phylogenetic tree with best found model
fitGTRGI_all <- pml_bb(phang.align, model="GTR+G(4)+I")

plot(fitGTRGI_all)
#mid-point-rooted maximum likelihood phylogenetic tree
midpoint_fitGTRGI_all = midpoint(fitGTRGI_all$tree)
plot(midpoint_fitGTRGI_all)
add.scale.bar()

#Adding tree info to phyloseq object
all_info_agg <- merge_phyloseq(all_info_agg,  phy_tree(midpoint_fitGTRGI_all))
#add sample data 
sam_data <- read.csv("D:\\MinION\\New_Analysis\\sam_data.csv", header=TRUE, row.names="Code", sep =";")
sample_data(all_info_agg) <- sam_data

#Save 
save(all_info_agg, file="D:\\MinION\\New_Analysis\\all_agg_phyloseq.RData")
load("D:\\MinION\\New_Analysis\\all_agg_phyloseq.RData")


#Clean ASV table names
tax_table(all_info_agg)[, colnames(tax_table(all_info_agg))] <- gsub(tax_table(all_info_agg)[, colnames(tax_table(all_info_agg))], 
                                                             pattern = "Methylobacterium-Methylorubrum", replacement = "MM")
tax_table(all_info_agg)[, colnames(tax_table(all_info_agg))] <- gsub(tax_table(all_info_agg)[, colnames(tax_table(all_info_agg))], 
                                                             pattern = "Burkholderia-Caballeronia-Paraburkholderia", replacement = "BCP")




#########Data analysis#######
#Core heatmap
#Use the format_to_besthit function from microbiomeutilities to get the best classification of the ASVs
all_genus <- aggregate_taxa(all_info_agg, "genus")
all_genus <- subset_taxa(all_genus, genus!="Unknown")
# Core with compositionals:
prevalences <- c(0, 1)
detections = c(0.1, 0.5, 1, 5)


p.core <- plot_core(all_genus, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(8, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = 0.75) + 
  xlab("Detection Threshold 
       (Relative Abundance (%))")
#asthetics 
p.core + theme(axis.text.y = element_text(size=10, color = 'black'),
               axis.text.x = element_text(size = 11, color = 'black', angle = 90, hjust=1,vjust=0.3),
               axis.title.y = element_text(size=10, face = 'bold'),
               legend.text = element_text(size =10, colour = 'black'), 
               legend.title = element_text(size = 11, colour = 'black', face = 'bold')) +
  scale_fill_gradientn(name = "Prevalence", colors = rev(brewer.pal(8, "Spectral")),
                       breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  coord_flip() +
  scale_y_discrete(position = 'left')

#Tree with abundances 

##Getting info to add
all_agg_df <- psmelt(all_info_agg)
# Filter out rows where Abundance = 0
all_agg_df <- all_agg_df %>% filter(Abundance != 0)

##Dataframe for log10 abundance
all_abund_genus_df <- select(all_agg_df, OTU, Abundance, Name)
all_abund_genus_df <- all_abund_genus_df %>% mutate(Abundance = log10(Abundance)) #log10 abundance transformation
all_abund_genus_df  <- distinct(all_abund_genus_df , .keep_all = TRUE) #remove duplicates
all_abund_genus_df <- all_abund_genus_df %>% pivot_wider(names_from = Name, values_from = Abundance)
all_abund_genus_df <- data.frame(all_abund_genus_df, row.names = "OTU")

##Dataframe for genus
all_agg_genus_df <- select(all_agg_df, OTU, genus)
all_agg_genus_df  <- distinct(all_agg_genus_df , .keep_all = TRUE) #remove duplicates
all_agg_genus_df <- data.frame(all_agg_genus_df, row.names = "OTU")



all_tree = ggtree(all_info_agg, ladderize = FALSE, size = 0.8, layout = "circular") +
  geom_tiplab(aes(label=genus), align =TRUE, size =3) +
  geom_treescale(fontsize=3, linesize=1, offset=0.4)




#all_tree <- rotate_tree(all_tree, -90)
print(all_tree)



abund_tree = gheatmap(all_tree, all_abund_genus_df, offset=0, width=.4,
         colnames_angle=0, colnames_offset_y = 0) +
  theme(legend.text = element_text(size = 10, color = 'black'),
        legend.key.size = unit(0.6, "cm"))+
  scale_fill_gradientn(colours = c("#EBF5FB", "#7FB3D5","#A9DFBF", "#F9E79F", "#F5CBA7", "#EB984E"),
                       na.value = "white", name = "Relative Abundance (log10)")
abund_tree

all_genus_abund <- abund_tree + new_scale_fill() 



colors = c("#B2ABD2","#35978F","#9970AB","#BABABA","#BF812D","#F7F7F7","#053061",
                    "#1B7837","#D53E4F","#74ADD1","#878787","#A50026","#FDDBC7",
                    "#E6F5D0","#E08214","#7F3B08","#2D004B","#01665E","#F1B6DA",
                    "#FEE08B","#00441B","#F46D43","#4393C3","#FEE0B6","#DE77AE",
                    "#C2A5CF","#313695","#D8DAEB","#FEE08B","#A6DBA0", "#F7F7a4",
                    "#D73027","#67001F","#543039","#80CDC1","#543008","#C7EAE5",
                    "#FDAE61","#1A9850","#762A83", "#D9F0D3", "#87CEEB", "#FF6347", 
                    "#B0E0E6", "#FFA07A", "#C1FFC1", "#D8BFD8", "#FFEC8B", "#98FB98","#ADD8E6")
                
                    

all_genus = gheatmap(all_genus_abund, all_agg_genus_df, offset=-1, width=0.95,
                     colnames_angle=0, colnames_offset_y = 0) +
  theme(legend.text = element_text(size = 9, color = 'black'),
        legend.key.size = unit(0.4, "cm")) + 
  scale_fill_manual(values = colors) + 
  guides(fill = guide_legend(ncol = 2))

all_genus



all_genus = gheatmap(all_tree, all_agg_genus_df, offset=-0.05, width=.05,
                     colnames_angle=0, colnames_offset_y = 0) +
  theme(legend.text = element_text(size = 9, color = 'black'),
        legend.key.size = unit(0.4, "cm")) + 
  scale_fill_manual(values = colors) + 
  guides(fill = guide_legend(ncol = 2))

all_genus

all_genus_abund <- all_genus + new_scale_fill() 


gheatmap(all_genus_abund, all_abund_genus_df, offset=0.01, width=.4,
          colnames_angle=0, colnames_offset_y = 0) +
theme(legend.text = element_text(size = 10, color = 'black'),
      legend.key.size = unit(0.6, "cm"))+
scale_fill_gradientn(colours = c("#EBF5FB", "#7FB3D5","#A9DFBF", "#F9E79F", "#F5CBA7", "#EB984E"),
                      na.value = "white", name = "Relative Abundance (log10)")


#Presence/absence heatmap
pa_genus_df <- select(all_agg_df, Name, genus)
pa_genus_df  <- distinct(pa_genus_df , .keep_all = TRUE) #remove duplicates

# Use dplyr to create a new data frame with counts
counts <- pa_genus_df %>%
  group_by(genus) %>%
  summarise(count = n_distinct(Name))


counts_pa_genus_df <-merge(pa_genus_df, counts, by = "genus", all= T)

#dataframe for distance
counts_pa_genus_df = counts_pa_genus_df[, c("Name", "count", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = count)

counts_pa_genus_df[is.na(counts_pa_genus_df)] <- 0

counts_pa_genus_df = data.frame(counts_pa_genus_df , row.names = "Name")


#The plot
pheatmap(counts_pa_genus_df, na_col = "transparent",
         color = c("transparent", "#EBF5FB", "#7FB3D5","#7FB3D5","#7FB3D5","#7FB3D5","#7FB3D5","#7FB3D5","#7FB3D5", "#154360"),
         fontsize=10.5, border_color= "white", cluster_rows=FALSE)


colors = c("#B2ABD2","#9970AB","#35978F","#67001F","#BF812D","#D9F0D3","#053061",
                    "#1B7837","#D53E4F","#74ADD1","#878787","#A50026","#FDDBC7",
                    "#F1B6DA","#FEE08B","#FDAE61","#2D004B","#01665E","#D8DAEB")
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




#Heat map with correlation
#group to genus level
all_genus <- tax_glom(all_info_agg, taxrank="genus", NArm = TRUE)
#Transform data to log10 score
all_log10  <- transform(all_genus, 'log10')

#Transform data to a dataframe 
all_log10_df <- psmelt(all_log10) #transform to dataframe
all_log10_df_na = replace(all_log10_df, all_log10_df == 0, NA)

#dataframe for distance
all_df_adundance_na = all_log10_df_na[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)

all_df_adundance_na = data.frame(all_df_adundance_na , row.names = "Name")
#dataframe for plot
all_log10_df_adundance = all_log10_df[, c("Name", "Abundance", "genus")] %>% 
  pivot_wider(names_from = genus, values_from = Abundance)

all_log10_df_adundance = data.frame(all_log10_df_adundance , row.names = "Name")

#Hierarchical clustering using Bray-Curtis distances

cluster_samples <- hclust(vegdist(all_log10_df_adundance, method="bray", na.rm = TRUE))
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

cluster_sapmles <- sort_hclust(cluster_samples)
plot(cluster_sapmles, main = "Sorted Dendrogram", xlab = "", sub = "")
#check number of possible clusters by heaigth 
clusters = cutree(cluster_sapmles, k = NULL, h = 0.3)
max(clusters)
#Get the isolates virulance and resistance data
isolates_data <- sam_data[, (ncol(sam_data) - 4):ncol(sam_data)]
isolates_data <- data.frame(isolates_data, row.names = "Name")

#Define the colors
my_colour <- list(Infection.capacity = c("High" = "#4A235A", "Medium-high" = "#7D3C98", 
                                         "Medium-low" = "#BB8FCE", "Low" = "#EBDEF0"),
                  Amphotericin.B = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"),
                  Voriconazole = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"),
                  Posaconazole = c("Resistant" = "#C0392B", "Susceptible" = "#16A085"))



#The plot
pheatmap(all_df_adundance_na, na_col = "transparent", 
         cluster_cols=FALSE, cluster_rows = cluster_sapmles,
         annotation_row = isolates_data, 
         color = colorRampPalette(c("#EBF5FB", "#7FB3D5", "#A9DFBF", "#F9E79F", "#F5CBA7", "#EB984E"))(100),
         annotation_colors = my_colour,
         fontsize=10.5, border_color= "white",
         cutree_rows = 3)



ggtree(all_info_agg, ladderize = FALSE, size = 0.8, open.angle=10) + 
  geom_text2(aes(subset=!isTip, label=label), hjust=1.5, size=2, vjust = 1) +
  #geom_text(aes(label=node), hjust=-.3, size =2) +
  geom_treescale(fontsize=3, linesize=1, offset=0.2) +
  geom_tiplab(aes(label=genus), align =TRUE, size =3.5)





plot_bar(all_info_agg, fill="genus") + 
  geom_bar(aes(color=genus), stat="identity", position="stack") 

bar_absolut


any(is.na(data.frame(agg_all@tax_table)$genus))