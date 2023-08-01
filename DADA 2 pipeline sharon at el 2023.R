library(dada2)
library(ShortRead)
library(Biostrings)


#### RAW SEQUENCE PROCESSING TO GENERATE ASVs (using dada2), done according to: : https://benjjneb.github.io/dada2/ITS_workflow.html####

path <- name the folder with the sequences 

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- add the path to cutadapt
system2(cutadapt, args = "--version") 

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[c(1:25)])

pdf("plotQualityProfile_F.pdf", height = 60, width = 60)
plotQualityProfile(cutFs[c(1:25)])
dev.off()

plotQualityProfile(cutRs[c(1:25)])

pdf("plotQualityProfile_R.pdf", height = 60, width = 60)
plotQualityProfile(cutRs[c(1:25)])
dev.off()

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out)

filtFs <- file.path(list.files(file.path(path.cut, "filtered"),pattern = "_R1_001.fastq.gz", full.names = TRUE))
filtRs <- file.path(list.files(file.path(path.cut, "filtered"),pattern = "_R2_001.fastq.gz", full.names = TRUE))

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

pdf("plotErrorsF.pdf", height = 20, width = 20)
plotErrors(errF, nominalQ = TRUE)
dev.off()
pdf("plotErrorsR.pdf", height = 20, width = 20)
plotErrors(errR, nominalQ = TRUE)
dev.off()

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

## using onlt forward sequences## 

seqtabAll.F <- makeSequenceTable(dadaFs)
table(nchar(getSequences(seqtabAll.F)))

seqtabNoC.F <- removeBimeraDenovo(seqtabAll.F)
write.csv(seqtabNoC.F,"seqtabNoC.F.csv")

briefToSeq.F <- colnames(seqtabNoC.F)
names(briefToSeq.F) <- paste0("ASV", seq(ncol(seqtabNoC.F))) # ASV1, ASV2, ...
seqtab.brief.F <- seqtabNoC.F
colnames(seqtab.brief.F) <- names(briefToSeq.F)
briefToSeq.F["ASV1"]
write.csv(seqtab.brief.F,"seqtabNoC2.F.csv")

fastaRef.UNITEF <- "C:/Users/user/Desktop/Sequences2021/sh_general_release_dynamic_s_02.02.2019_dev.fasta"
taxidbp50.F.UNITEF <- assignTaxonomy(seqtabNoC.F, refFasta = fastaRef.UNITEF, minBoot = 50, outputBootstraps = TRUE,  tryRC = TRUE, multithread = TRUE)
fastaRef.UNITEA <- "C:/Users/user/Desktop/Sequences2021/sh_general_release_dynamic_s_all_02.02.2019_dev.fasta"
taxidbp50.F.UNITEA <- assignTaxonomy(seqtabNoC.F, refFasta = fastaRef.UNITEA, minBoot = 50, outputBootstraps = TRUE,  tryRC = TRUE, multithread = TRUE)

## export taxonomy tables## 

rownames(taxidbp50.F.UNITEF[[1]]) <- names(briefToSeq.F)
rownames(taxidbp50.F.UNITEF[[2]]) <- names(briefToSeq.F)
rownames(taxidbp50.F.UNITEA[[1]]) <- names(briefToSeq.F)
rownames(taxidbp50.F.UNITEA[[2]]) <- names(briefToSeq.F)
write.csv(taxidbp50.F.UNITEF[[1]],"taxidbp50.F.UNITEFa.csv")
write.csv(taxidbp50.F.UNITEF[[2]],"taxidbp50.F.UNITEFb.csv")
write.csv(taxidbp50.F.UNITEA[[1]],"taxidbp50.F.UNITEAa.csv")
write.csv(taxidbp50.F.UNITEA[[2]],"taxidbp50.F.UNITEAb.csv")

taxidbp50.dfm1 <- as.data.frame(taxidbp50.F.UNITEA[[1]]) 
taxidbp50.dfm2 <- as.data.frame(taxidbp50.F.UNITEA[[2]])
dim(taxidbp50.dfm1);taxidbp50.dfm1[1,];rownames(taxidbp50.dfm1)[1]
dim(taxidbp50.dfm2);taxidbp50.dfm2[1,];rownames(taxidbp50.dfm2)[1]

# remove non fungal ASVs## 

nonFungalASVs <- c("ASV163", "ASV369", ....)
allASVs <- rownames(taxidbp50.dfm1)
goodAsvs <- allASVs[!(allASVs %in% nonFungalASVs)]
taxidbp50.flt <- taxidbp50.dfm1[goodAsvs, ]

seqtab.brief.flt <- seqtab.brief.F[,goodAsvs]
dim(seqtab.brief.F)
dim(seqtab.brief.flt)
write.csv(seqtab.brief.flt,"seqtab.brief.flt.csv")

## arrange tax table ## 

tax <- taxidbp50.flt

tax[] <- lapply(tax, gsub, pattern='s__', replacement='')
tax[] <- lapply(tax, gsub, pattern='g__', replacement='')
tax[] <- lapply(tax, gsub, pattern='f__', replacement='')
tax[] <- lapply(tax, gsub, pattern='o__', replacement='')
tax[] <- lapply(tax, gsub, pattern='c__', replacement='')
tax[] <- lapply(tax, gsub, pattern='p__', replacement='')
tax[] <- lapply(tax, gsub, pattern='k__', replacement='')

tax.clean <- tax
tax.clean[is.na(tax.clean)] <- ""
tax.clean[c(1:10),]

tax$asvs <- rownames(tax)

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,1] == ""){
    tax.clean$Kingdom[i] <- paste("Fungi", tax$asvs[i], sep = "_")
    tax.clean[i, 1:7] <- tax.clean$Kingdom[i]
  } else if (tax.clean[i,2] == ""){
    tax.clean$Phylum[i] <- paste(tax.clean$Kingdom[i], tax$asvs[i], sep = "_")
    tax.clean[i, 2:7] <- tax.clean$Phylum[i]
  } else if (tax.clean[i,3] == ""){
    tax.clean$Class[i] <- paste(tax.clean$Phylum[i], tax$asvs[i], sep = "_")
    tax.clean[i, 3:7] <- tax.clean$Class[i]
  } else if (tax.clean[i,4] == ""){
    tax.clean$Order[i] <- paste(tax.clean$Class[i], tax$asvs[i], sep = "_")
    tax.clean[i, 4:7] <- tax.clean$Order[i]
  } else if (tax.clean[i,5] == ""){
    tax.clean$Family[i] <- paste(tax.clean$Order[i], tax$asvs[i], sep = "_")
    tax.clean[i, 5:7] <- tax.clean$Family[i]
  } else if (tax.clean[i,6] == ""){
    tax.clean$Genus[i] <- paste(tax.clean$Family[i], tax$asvs[i], sep = "_")
    tax.clean[i, 6:7] <- tax.clean$Genus[i]
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax$asvs[i], sep = "_")
  } else {
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = "_")
  }
}
tax.clean[c(1:10),]

write.csv(tax.clean,"taxidbp50clean.F.UNITEAa.csv")
