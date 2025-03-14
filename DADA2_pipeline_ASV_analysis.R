#DADA2 PIPELINE FOR ILLUMINA 16S SEQUENCING ANALYSIS


#call installed packages and checks the version
library(devtools)
library(dada2); packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(DECIPHER); packageVersion("DECIPHER")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(MicEco)

(require(stats))

#IMPORTANT NOTES TO USE DADA2 
#1. sample are splited into individual per-sample fastq files
#2. primers, adapters, linkers have been removed (done in the script below) 
#3. forward and reverse fastq files contain reads in matched order (done in the script bellow)


#Set working directory (is better to create a folder just with seq files)
path = "C:\\SequencingV3V4\\R_analysis\\all_raw"
list.files(path) #important to check if all files are present

#Generate matched lists of the forward and reverse read files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

#Extract samples names
sample.names <- sapply(strsplit(basename(fnFs), "-N"), `[`, 1)


#Inspect the quality profiles of the reads
plotQualityProfile(fnFs[1:6]) #forward
plotQualityProfile(fnRs[1:6]) #reverse


##"pre-filter" the sequences 
#remove Ns (ambiguous bases), remove primers 
fnFs.pre_filt <- file.path(path, "pre_filt", paste0(sample.names, "_F_prefilt.fastq")) # Put filtered files in pre_filtered subdirectory
fnRs.pre_filt <- file.path(path, "pre_filt", paste0(sample.names, "_R_prefilt.fastq"))
filterAndTrim(fnFs, fnFs.pre_filt, maxN = 0, multithread = FALSE, trimLeft = 19, trimRight=20)
filterAndTrim(fnRs, fnRs.pre_filt, maxN = 0, multithread = FALSE, trimLeft = 19, trimRight=20)

##Filter and Trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#maxEE parameter sets the maximum number of "expected errors" allowed in a read (may need to be changed depending on the samples)
#truncLen parameter trims end of the red to n bp
out <- filterAndTrim(fnFs.pre_filt, filtFs, fnRs.pre_filt, filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs=TRUE)
head(out)

#Inspect the quality profiles of the reads
plotQualityProfile(filtFs[1:6]) #forward
plotQualityProfile(filtRs[1:6])

##Learn the Error Rates
#Estimation of the error rates and inference of sample composition using parametric error model
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
#Visualization of the the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Sample Inference
#Applying the core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaFs[1]
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaRs[1]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #(n samples, n ASVs)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

##Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##Assign taxonomy using IdTaxa machine learning
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/PhD/Sequencing_analysis_22/Mycelium_ciprofloxacin_treatment/Data/SILVA_SSU_r138_2019.RData") #CHECK IF FILE PATH IS CORRECT
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#creating phyloseq object 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxid))

#add DNA info
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
#add samples metadata
sam_data <- read.csv("C:\\SequencingV3V4\\R_analysis\\all_raw\\sam_data.csv", header=T, row.names="Code",sep =";")

ps@sam_data <- sample_data(sam_data)
#Removing singletons
abundance_matrix <- as.matrix(otu_table(ps))
# Replace values lower than 0.005 with 0
abundance_matrix[abundance_matrix < 1] <- 0
# Update the abundance matrix in the phyloseq object
otu_table(ps) <- as.matrix(abundance_matrix)

#Removing ASVs representing less than 0.005% of the total high-quality sequences

#asv_info.relativ <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU)*100)
#asv_info.relativ_0.005 <- ps_prune(asv_info.relativ, min.abundance = 0.005)
asv_info_filt = filter_taxa(ps, function(x) sum(x) > 1, TRUE)


# Now, 'asv_info_subset' contains the ASVs whose relative abundance is higher than 0.005%.

#Save dataframe as RData 
save(asv_info_filt, file="16S_taxid_table_output_dada2_v3v4_all.RData")

#Save each separate ps dataframe as csv file - tax_table, refseq, abund_table, sam_table
write.csv(asv_info_filt@tax_table,"C:\\SequencingV3V4\\R_analysis\\all_raw\\tax_table.csv", row.names = TRUE)
write.csv(asv_info_filt@refseq,"C:\\SequencingV3V4\\R_analysis\\all_raw\\refseq.csv", row.names = TRUE)
write.csv(asv_info_filt@otu_table,"C:\\SequencingV3V4\\R_analysis\\all_raw\\abund_table_all.csv", row.names = TRUE)