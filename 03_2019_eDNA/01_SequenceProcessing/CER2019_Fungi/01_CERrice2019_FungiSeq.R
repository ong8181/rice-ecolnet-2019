####
#### Rice field samples (from 2019.5.20-9.16, 9 plots)
#### No.1 CMR-011 Fungi sequence analysis by DADA2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("00_SessionInfo")
dir.create("01_CERrice2019_FungiSeqOut")

# Load library and functions
library(dada2); packageVersion("dada2") # 1.16.0, 2021.1.13
library(ShortRead); packageVersion("ShortRead") # 1.46.0, 2021.1.13
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.13
source("../functions/F01_HelperFunctions.R")

# Load sequence reads
# CMR-011-Fungi
path <- "/01_SequenceProcessing/seqdata_all/CMR-011-FUN"
fnFs <- sort(list.files(path, pattern=".forward.fastq", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern=".reverse.fastq", full.names = T)) # Reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
#sample_names <- sapply(strsplit(fnFs, "_"), `[`, 5)

# Identify primers
FWD <- "CTHGGTCATTTAGAGGAASTAA" # ITS1FKYO1
REV <- "TTYRCTRCGTTCTTCATC" # ITS2KYO2
FWD_orients <- AllOrients(FWD)#; AllOrients(FWD_rc)
REV_orients <- AllOrients(REV)#; AllOrients(REV_rc)

# Pre-filtering to remove Ns
fnFs_filtN <- file.path(path, "01_filtN", str_sub(basename(fnFs), start = 1, end = -4)) # Put N-filterd files in filtN/ subdirectory
fnRs_filtN <- file.path(path, "01_filtN", str_sub(basename(fnRs), start = 1, end = -4))
filterAndTrim(fnFs, fnFs_filtN, fnRs, fnRs_filtN, maxN = 0, multithread = TRUE, compress = FALSE)
# Compress fastq
pigz <- "~/miniconda3/bin/pigz"
#system2(pigz, args = "--version")
system2(pigz, args = c("-p 72", paste0(file.path(path, "01_filtN"), "/*.fastq")))
fnFs_filtN <- file.path(path, "01_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs_filtN <- file.path(path, "01_filtN", basename(fnRs))

# Identify primers
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs_filtN[[4]]),
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs_filtN[[4]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs_filtN[[4]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs_filtN[[4]]))

# Tell the path to cutadapt command
cutadapt <- "~/miniconda3/bin/cutadapt"
system2(cutadapt, args = "--version")

# Remove primers
path_cut <- file.path(path, "02_cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)
fnFs_cut <- file.path(path_cut, basename(fnFs))
fnRs_cut <- file.path(path_cut, basename(fnRs))

FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1_flags <- paste("-g", FWD, "-a", REV_RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2_flags <- paste("-G", REV, "-A", FWD_RC)

# Run Cutadapt
# Cutadapt commands may crash from R console. In that case, excute it from terminal.
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-j 72", # Multithred option
                             R1_flags, R2_flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs_cut[i], "-p", fnRs_cut[i], # output files
                             fnFs_filtN[i], fnRs_filtN[i])) # input files
}

# Sanity check
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs_cut[[4]]), 
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs_cut[[4]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs_cut[[4]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs_cut[[4]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path_cut, pattern = ".forward.fastq", full.names = TRUE))
cutRs <- sort(list.files(path_cut, pattern = ".reverse.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
sample_names <- unname(sapply(cutFs, function(x) strsplit(basename(x), "_")[[1]][2]))
head(sample_names)
# Visualize quality# Visualize quality
#plotQualityProfile(fnFs[101:103])
#plotQualityProfile(fnRs[101:103])

# Perform quality filtering
filtFs <- file.path(path, "03_filtered", paste0(sample_names, "_F_filt.fastq"))
filtRs <- file.path(path, "03_filtered", paste0(sample_names, "_R_filt.fastq"))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     # Output sequences of Claident already trimmed Ns and primers
                     #trimLeft = c(0,0),
                     trimRight = c(0, 50), # Visual inspection suggests the end of reverse reads is low-quality
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=F, minLen = 20,
                     compress=FALSE, multithread=TRUE) # On Windows set multithread=FALSE
# Compress files
system2(pigz, args = c("-p 72", paste0(file.path(path, "03_filtered"), "/*.fastq")))
filtFs <- file.path(path, "03_filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "03_filtered", paste0(sample_names, "_R_filt.fastq.gz"))
head(out)

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample_names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(path, "03_filtered", paste0(sample_names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "03_filtered", paste0(sample_names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
#min_nbases <- 5e+09
min_nbases <- 300 * sum(out[,2])
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
#plotErrors(errF, nominalQ = T)
#plotErrors(errR, nominalQ = T)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[out[,2]>0 & out[,1]>0]
names(derepRs) <- sample_names[out[,2]>0 & out[,1]>0]

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dadaFs[[1]]

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang = TRUE, minOverlap = 20, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab)

# Track reads thourhg the pipeline
out2 <- out[out[,2]>0 & out[,1]>0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2]>0 & out[,1]>0]
head(track)

# Taxa output for claident tax assginment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">Taxa%05d", 1:length(seqs)), seqs)), ncol=1)
write.table(seqs_out, "01_CERrice2019_FungiSeqOut/FungiASV_seqs.fa", col.names = FALSE, row.names = FALSE, quote = FALSE)

save.image("01_CERrice2019_FungiSeqOut/01_CERrice2019_FungiSeq_AllOut.RData")

# Reduce size
rm(dadaFs)
rm(dadaRs)
rm(derepFs)
rm(derepRs)

save(list = ls(all.names = TRUE),
     file = "01_CERrice2019_FungiSeqOut/01_CERrice2019_FungiSeqOut.RData")

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_FungiSeq_%s.txt", substr(Sys.time(), 1, 10)))

