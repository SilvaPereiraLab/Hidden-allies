#DADA2 PIPELINE FOR ILLUMINA 16S SEQUENCING ANALYSIS

#

#call installed packages and checks the version
library(dada2); packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(DECIPHER); packageVersion("DECIPHER")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

(require(stats))

#IMPORTANT NOTES TO USE DADA2 
#1. sample are splited into individual per-sample fastq files
#2. primers, adapters, linkers have been removed (done in the script below) 
#3. forward and reverse fastq files contain reads in matched order (done in the script bellow)


#Set working directory, is the folder were all files for analysis are stored (is better to create a folder just with seq files)
path = "C:/Mycelium_ciprofloxacin_treatment/Data/Sequences"
list.files(path) #important to check if all files are present

#Generate matched lists of the forward and reverse read files
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))
#Extract samples names
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)


##"pre-filter" the sequences 
#remove Ns (ambiguous bases), remove primers 
fnFs.pre_filt <- file.path(path, "pre_filt", basename(fnFs)) # Put filtered files in pre_filtered subdirectory
fnRs.pre_filt <- file.path(path, "pre_filt", basename(fnRs))
filterAndTrim(fnFs, fnFs.pre_filt, maxN = 0, multithread = FALSE, trimLeft = 19, trimRight=20)
filterAndTrim(fnRs, fnRs.pre_filt, maxN = 0, multithread = FALSE, trimLeft = 19, trimRight=20)

##Check if primers can still be found in the sequences (sanity check)
#primers sequences
FWD <- "GTGYCAGCMGCCGCGGTAA"  #515F
REV <- "GGACTACNVGGGTWTCTAAT" #806R
#function to create all possible orientation of the primers 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
#applying the function to the primers to generate primers orients
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
#function to 
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
#checking the primer presence in the sequences using created functions
#change fnFs.pre_filt[[n]] to select sequence (n)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.pre_filt[2]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.pre_filt[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.pre_filt[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.pre_filt[[2]]))
#Inspect the quality profiles of the reads
plotQualityProfile(fnFs.pre_filt[1:6]) #forward
plotQualityProfile(fnRs.pre_filt[1:6]) #reverse

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
dadaFs[1:6]
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[1:6]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 30)
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
load("~/Mycelium_ciprofloxacin_treatment/Data/SILVA_SSU_r138_2019.RData") #CHECK IF FILE PATH IS CORRECT
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


##DATA TREATMNET
#Creating a dataframe with information of samples
samples.out <- rownames(seqtab.nochim)
isolate <- c("Af_CI.01", "Af_CI.01","Af_CI.02","Af_CI.02","Af_CI.03","Af_CI.03")
treatment <- c("yes", "no", "yes", "no", "yes", "no")
samdf <- data.frame(Isolate=isolate, Treatment=treatment)
rownames(samdf) <- samples.out
#
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxid))
#
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- names(refseq)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Save dataframe as RData 
save(ps, file="16S_taxid_table_output_dada2.RData")
load("~/Mycelium_ciprofloxacin_treatment/Data/16S_taxid_table_output_dada2.RData")
#to load this dataframe 

#Save each separate ps dataframe as csv file - tax_table, refseq, abund_table, sam_table
write.csv(ps@tax_table,"~/Mycelium_ciprofloxacin_treatment/Data/\\tax_table.csv", row.names = TRUE)
write.csv(ps@refseq,"~/Mycelium_ciprofloxacin_treatment/Data/\\refseq.csv", row.names = TRUE)
write.csv(ps@otu_table,"~/Mycelium_ciprofloxacin_treatment/Data/\\abund_table.csv", row.names = TRUE)
write.csv(ps@sam_data,"~/Mycelium_ciprofloxacin_treatment/Data/\\sam_data.csv", row.names = TRUE)







