#DATA ANALYSIS OF AVS TABLE OBTAINED FROM DADA2 PIPELINE

#Mycelium ciprofloxacin treatment

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
library(ggtreeExtra)
library(microbiome)
###
set.seed(10)
#Importing dataframes
tax_table <- read.csv("~/Mycelium_ciprofloxacin_treatment/Data/CSV_final_files/tax_table.csv", header=TRUE, row.names="ASV", sep =";")
refseq <- read.csv("~/Mycelium_ciprofloxacin_treatment/Data/CSV_final_files/refseq.csv", header=TRUE, row.names="ASV", sep =";")
abund_table <- read.csv("~/Mycelium_ciprofloxacin_treatment/Data/CSV_final_files/abund_table.csv", header=TRUE, row.names="Isolate", sep =";")
sam_data <- read.csv("~/Mycelium_ciprofloxacin_treatment/Data/CSV_final_files/sam_data.csv", header=TRUE, row.names="ID", sep =";")
abund_table_rev <- read.csv("~/Mycelium_ciprofloxacin_treatment/Data/CSV_final_files/abund_table_rev.csv", header=TRUE, row.names="ASV", sep =";")
#Sequences information 
dna <- Biostrings::DNAStringSet(refseq$sequence)
names(dna) <- row.names(refseq)

#Construction of phylogenetic tree using phangron package
mult <- msa(dna, method="ClustalW", type="dna", order="input") #alignment of sequences
phang.align <- as.phyDat(mult, type="DNA")

#Model test
mt <- modelTest(phang.align, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
                control = pml.control(trace = 0))
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
save(asv_info, file="16S_antibio_taxid_table_dada2_tree_final.RData")
load("~/Mycelium_ciprofloxacin_treatment/Data/16S_antibio_taxid_table_dada2_tree_final.RData")


#Convert abundances to relative abundances 
asv_info.relativ <- transform_sample_counts(asv_info, function(OTU) OTU/sum(OTU)*100)

###################################################
###Look to all ASVs
##Heatmap at class level total ASVs
#phyloseq object need to be summed at class level
class_pobject <- tax_glom(asv_info.relativ, taxrank="class")
class_df <- psmelt(class_pobject) #transform to dataframe

#plot heatmap
heatmap_class <- ggplot(class_df, aes(x=Antibiotic, y=class, fill = Abundance)) +
  geom_tile(color = "#283747",lwd = 0.1,linetype = 1) + facet_wrap(~Isolate) +
  scale_fill_gradientn(name = "Relative Abundance (%)", 
                       colours= c("#FDFEFE", "#7FB3D5", "#A9DFBF", "#F9E79F"), 
                       na.value = "transparent",
                       values = c(0, rescale(5, from = range(class_df$Abundance)), 1))+
  theme(strip.background=element_rect(fill="#ebebeb"), 
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold", vjust = 1.8),
        legend.key.size = unit(0.8, "cm"),
        axis.text.x = element_text(size =11, color = "black"),
        axis.text.y = element_text(size= 11, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=13))+
  coord_fixed() #makes tiles in square shape 
  
heatmap_class
#save with height:650 


##Correlation coefficients between the samples (using dataframe: abund_table_rev)  
corr_pearson <- rcorr(as.matrix(abund_table_rev)) #dataframe needs to be transformed to matrix
corr_matrix <- corr_pearson$r #correlation matrix
p_value <- corr_pearson$P #p-value matrix 

corrplot(corr_matrix, type = "lower", order = "hclust",diag=FALSE, method= "color",
         p.mat = p_value, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig")


corr_pearson_plot <- ggcorrplot(corr_matrix, hc.order = TRUE, type = "lower", 
           outline.col = "white", legend.title = "Pearson Correlation",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#D98880"),
           p.mat = p_value) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold", vjust = 2))
corr_pearson_plot
#save with height:472

##Alpha and Beta diversities 
#To use mia package the phyloseq object need to me converted 
asv_tse <- makeTreeSummarizedExperimentFromPhyloseq(asv_info.relativ) #convert phyloseq to TSE
asv_tse

#Calculate observed richness ()
asv_tse <- mia::estimateRichness(asv_tse, abund_values = "counts", index = "observed", name="observed")
#calculate diversity: 
asv_tse  <- mia::estimateDiversity(asv_tse , abund_values = "counts",index = "shannon", name = "shannon")
asv_tse  <- mia::estimateDominance(asv_tse , abund_values = "counts", index="relative")
#Calculate eveness
asv_tse  <- mia::estimateEvenness(asv_tse , abund_values = "counts", index="pielou")
head(colData(asv_tse)) #check dataframe

plot_shannon <- plotColData(asv_tse, "shannon", "Antibiotic", colour_by = "Isolate") +
  ylab(expression("Shannon Index")) + 
  theme_gray()

plot_evenness <- plotColData(asv_tse, "pielou", "Antibiotic", colour_by = "Isolate") +
  theme_gray() + 
  ylab(expression("Evenness Index"))

plot_relative <- plotColData(asv_tse, "relative", "Antibiotic", colour_by = "Isolate") +
  theme_gray() + 
  ylab(expression("Dominance Index"))


aplha_div_plot <- grid.arrange(plot_shannon, plot_evenness, plot_relative, ncol=3)
aplha_div_plot

#Checking p-values using non-parametric Kolmogorov-Smirnov test for two-group comparisons
#Create data frame from the collected data
colData_df <- as.data.frame(colData(asv_tse))

# Split the values by group for each calculated index
spl_pielou <- split(df$pielou, df$Antibiotic) #evenness index
spl_shannon <- split(df$shannon, df$Antibiotic) #diversity index
spl_relative <- split(df$relative, df$Antibiotic) #dominance index


# Kolmogorov-Smironv test comparing the two groups: Antibiotic - yes and now
pvalue_pielou <- ks.test(spl_pielou$no, spl_pielou$yes)$p.value
pvalue_shannon <- ks.test(spl_shannon$no, spl_shannon$yes)$p.value
pvalue_relative <- ks.test(spl_relative$no, spl_relative$yes)$p.value
pvalue_list <- c(pvalue_pielou,pvalue_shannon,pvalue_relative)
# Adjust the p-value using Benjamini & Hochberg method
padj <- p.adjust(pvalue_list, method ="BH")
padj

Index <- c("evenness", "diversity", "domimance")
Signif <- c("NS", "NS", "NS")

pvalue_df <- data.frame(Index, pvalue_list, Signif)


###################################################


###Look to 100 most abundant ASVs
top100 <- names(sort(taxa_sums(asv_info.relativ), decreasing=TRUE))[1:100]
asv.top100 <- prune_taxa(top100, asv_info.relativ)



##Creating tree plot
#First step is to know witch node correspond to witch class, to define later in dataframe
tree_nodes = ggtree(asv.top100, ladderize = FALSE, size = 1, aes(color=class)) + 
  #geom_point(aes(color = class),na.rm=TRUE, size = 2) +
  geom_text2(aes(subset=!isTip, label=label), hjust=-2, size=0.8) +
  geom_treescale(fontsize=3, linesize=0.8, offset=-2) +
  geom_text(aes(label=node), hjust=-.3, size =2) +
  geom_tiplab(align =TRUE, size = 2) +
  geom_nodelab(aes(label=label), hjust=-.05, size=0.1) 
print(tree_nodes)

#Define the dataframe with the class nodes to hilight 
nodes_top100 <- c(172,174,134,185,102,121,184,149,90,74,181,34)
class_node_top100 <- c("Acidobacteriae", "Actinobacteria", "Alphaproteobacteria","Bacilli", 
                       "Bacteroidia", "Clostridia", "Cyanobacteriia", "Gammaproteobacteria", 
                       "Oligoflexia", "Polyangia", "Vampirivibrionia", "Verrucomicrobiae")
nodes_data_top100 <- data.frame(node=nodes_top100, type=class_node_top100)

#Retrieve just tree info from phyloseq object
tree_info_top100 = asv.top100@phy_tree
#To plot the tree with the the heatmap, we need to transform phyloseq object to dataframe (selection columns of interest)
top100_df <- psmelt(asv.top100) %>% select(label = OTU, Abundance, Colony, Isolate)

#plot the tree with highlighted class
tree_top100_class <- ggtree(tree_info_top100, ladderize = FALSE, size = 0.8) + 
  geom_text2(aes(subset=!isTip, label=label), hjust=1.5, size=2, vjust = 1) +
  #geom_text(aes(label=node), hjust=-.3, size =2) +
  #geom_treescale(fontsize=3, linesize=1, offset=-5) +
  geom_tiplab(align =TRUE, size =2) +
  geom_hilight(data=nodes_data_top100, aes(node=node, fill=type)) +
  hexpand(.01)+
  theme(legend.text = element_text(size =13), 
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.7, "cm")) +
  guides(fill=guide_legend(title="Class")) 
tree_top100_class


##Creating heatmap of relative abundance
heatmap_top100 <- ggplot(top100_df, aes(x=Antibiotic, y=label, fill = Abundance)) +
  geom_tile(color = "#283747",lwd = 0.1,linetype = 1) + facet_wrap(~Isolate) +
  theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + 
  scale_fill_gradientn(name = "Relative Abundance (%)", 
                       colours= c("#FDFEFE","#7FB3D5", "#2471A3", "#A9DFBF", "#F9E79F", "#FAD7A0"), 
                       na.value = "transparent",breaks=c(0,5,10,15,20,25), 
                       values = c(0, rescale(6, from = range(top100_df$Abundance)), 1)) + #makes color scale discontinuous
  theme(strip.background=element_rect(fill="#ebebeb"), 
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold", vjust = 1.5),
        legend.key.size = unit(0.8, "cm"),
        legend.direction="horizontal",
        axis.text.x = element_text(size =11, color = "black"),
        axis.title.x = element_text(size=13)) +
  coord_fixed(ratio = 0.25)
heatmap_top100

#Combine the two plots
heatmap_top100 %>% insert_left(tree_top100_class)
#100 top bacteria correspond to >90% of relative abundance 
top100_df %>% group_by(Isolate, Colony) %>% summarise(abundance = sum(Abundance))


#################
##Beta-diversity
#UniFrac metric incorporates phylogenic information into beta-diversity analysis


# PCoA plot using the unweighted UniFrac as distance
uwunifrac_dist = phyloseq::distance(asv_info.relativ, method="unifrac", weighted=F)
un.ordination = ordinate(asv_info.relativ, method="PCoA", distance=uwunifrac_dist)

#plot
pcoa_uw <- plot_ordination(asv_info.relativ, un.ordination, color = "Antibiotic") + geom_point(size = 2.8) + 
  xlim(-0.5,0.5) + ylim(-0.5, 0.5) + coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#7FB3D5", "#F1948A")) +
  theme(axis.text.x = element_text(size =10, color = "black"),
        axis.text.y = element_text(size =10, color = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold", vjust = 1.5))
#PERMANOVA
permanova_uw = adonis(uwunifrac_dist ~ Antibiotic, data = metadf)
permanova_uw$aov.tab

# PCoA plot using the weighted UniFrac as distance
wunifrac_dist = phyloseq::distance(asv_info.relativ, method="wunifrac")
w.ordination = ordinate(asv_info.relativ, method="PCoA", distance=wunifrac_dist)

#plot
pcoa_w <- plot_ordination(asv_info.relativ, w.ordination, color = "Antibiotic") + geom_point(size = 2.8) + 
  xlim(-0.5,0.5) + ylim(-0.5, 0.5) + coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#7FB3D5", "#F1948A")) +
  theme(axis.text.x = element_text(size =10, color = "black"),
        axis.text.y = element_text(size =10, color = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold", vjust = 1.5))
pcoa_w
#PERMANOVA
permanova_w = adonis(wunifrac_dist ~ Antibiotic, data = metadf)
permanova_w$aov.tab


cowplot::plot_grid(pcoa_w, pcoa_uw, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))


#Generate distances
ord_unifrac <- ordinate(asv_info.relativ, method = "PCoA", 
                        distance = "wunifrac", grouping_column = 'Antibiotic') 
ord_unifrac_un <- ordinate(asv_info.relativ, method = "PCoA", distance = "unifrac")  




pcoa_uw <- plot_ordination(asv_info.relativ, ord_unifrac_un, color = "Antibiotic") + geom_point(size = 2.8) +
  xlim(-0.5,0.5) + ylim(-0.5, 0.5) + coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#7FB3D5", "#F1948A")) +
  theme(axis.text.x = element_text(size =10, color = "black"),
        axis.text.y = element_text(size =10, color = "black"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold", vjust = 1.5))


cowplot::plot_grid(pcoa_w, pcoa_uw, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))


metadf <- data.frame(sample_data(asv_info.relativ))
unifrac.dist <- UniFrac(asv_info.relativ, 
                        weighted = TRUE, 
                        normalized = FALSE,  
                        parallel = FALSE, 
                        fast = TRUE)
unifrac.dist
permanova <- adonis(unifrac.dist ~ Antibiotic, data = metadf)
permanova

anosim <- anosim(unifrac.dist, grouping =)


###############################
####Analysis of ASVs for paper

otu_table <- as.data.frame(asv_info.relativ@otu_table) %>%
  mutate_all(~ if_else(. < 0.01, 0, 1))
# Extract the names of columns for Af1, Af2, and Af3
af_names <- c("Af1", "Af2", "Af3")

# Extract the names of columns for Af1_cip, Af2_cip, and Af3_cip
af_cip_names <- paste0(af_names, "_cip")

# Count the number of columns where the value is 1 in either Af1_cip, Af2_cip, or Af3_cip
# and 0 in either Af1, Af2, or Af3
count <- sum(colSums(otu_table[af_cip_names, ]) > 0 & colSums(otu_table[af_names,]) == 0)

sum(colSums(otu_table[af_cip_names, ]) == 0 & colSums(otu_table[af_names,]) > 0)



library(ggVennDiagram)

column_sums_cip <- colSums(otu_table[af_cip_names, ])
asv_cip <- names(column_sums_cip[column_sums_cip > 0])


column_sums <- colSums(otu_table[af_names, ])
asv_nocip <- names(column_sums[column_sums > 0])

ggVennDiagram(list(asv_cip, asv_nocip), 
              category.names = c("yes", "no")) +
              ggplot2::scale_fill_gradient(low='#CD6155', high = '#5499C7', mid = '#EBDEF0') + coord_flip()



common_asv <- intersect(asv_cip, asv_nocip)

result <- setdiff(asv_cip, common_asv)

df = as.data.frame(asv_info.relativ@otu_table)[, result]

comun_cip = df[af_cip_names, ]
comun_nocip = df[af_names, ]

rowSums(comun_cip)
rowSums(comun_nocip)

colSums(comun_cip) /3
colSums(comun_nocip) /3

# Function to count unique values
count_unique <- function(column) {
  length(unique(column))
}

# Apply the function to each column
unique_counts <- sapply(as.data.frame(asv_info.relativ@tax_table), count_unique)


rowSums(as.data.frame(asv.top100@otu_table))
rowSums(as.data.frame(asv_info.relativ@otu_table))








