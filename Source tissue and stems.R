library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(tidyr)
library(tidyverse)   
library(DESeq2) 
library(picante)
library(reshape2)
library(rstatix)
library(ggthemes)
library(ComplexHeatmap)
library(ggvenn)

## uploading taxonomy files from DADA2 ##

samdf <- read.csv("metadata.csv",header = TRUE)
rownames(samdf) <- samdf$SampleID
seqtab.brief.flt <- read.csv("seqtab.brief.flt.csv", header = TRUE, row.names = 1)
rownames(seqtab.brief.flt) <- samdf$SampleID
tax.final <- read.csv("taxidbp50clean.F.UNITEAa.csv", row.names = 1, header = TRUE)

## creating phyloseq object ##

ps <- phyloseq(otu_table(seqtab.brief.flt, taxa_are_rows=FALSE), sample_data(samdf), tax_table(as.matrix(tax.final)))
ps2 <- tax_glom(ps, "Species", NArm = TRUE)
ps3 <- prune_taxa(taxa_sums(ps2) > 2, ps2) #remove singletons and doubletons
ps3 <- prune_samples(rowSums(otu_table(ps)) > 100, ps3) #remove samples with less than 100 seqs

## abundace filtering according to tissues ## 

## Seed ##

ps.seed <- subset_samples(ps3, Sample_type == "Seed", drop_zeroes = TRUE)
ps.seed <- prune_taxa(taxa_sums(ps.seed) > 0, ps.seed)

psra.seed <- transform_sample_counts(ps.seed, function(x){x / sum(x)}) # ra transforms, but just for filtering, not for subsequent analysis
ps5.seed <- filter_taxa(psra.seed, function(x) sum(x) > .005, TRUE) # filtering

keepseed<- taxa_names(ps5.seed)
ps.seed.counts <- prune_taxa(keepseed, ps.seed)

## Stem ##

ps.stem <- subset_samples(ps3, Sample_type == "Stem", drop_zeroes = TRUE)
ps.stem <- prune_taxa(taxa_sums(ps.stem) > 0, ps.stem)

psra.stem <- transform_sample_counts(ps.stem, function(x){x / sum(x)}) # ra transforms, but just for filtering, not for subsequent analysis
ps5.stem <- filter_taxa(psra.stem, function(x) sum(x) > .005, TRUE) # filtering

keepstem<- taxa_names(ps5.stem)
ps.stem.counts <- prune_taxa(keepstem, ps.stem)

## Embryo ##

ps.embryo <- subset_samples(ps3, Sample_type == "Embryo", drop_zeroes = TRUE)
ps.embryo <- prune_taxa(taxa_sums(ps.embryo) > 0, ps.embryo)

psra.embryo <- transform_sample_counts(ps.embryo, function(x){x / sum(x)}) # ra transforms, but just for filtering, not for subsequent analysis
ps5.embryo <- filter_taxa(psra.embryo, function(x) sum(x) > .005, TRUE) # filtering

keepembryo<- taxa_names(ps5.embryo)
ps.embryo.counts <- prune_taxa(keepembryo, ps.embryo)

## Callus ##

ps.callus <- subset_samples(ps3, Sample_type == "Callus", drop_zeroes = TRUE)
ps.callus <- prune_taxa(taxa_sums(ps.callus) > 0, ps.callus)

psra.callus <- transform_sample_counts(ps.callus, function(x){x / sum(x)}) # ra transforms, but just for filtering, not for subsequent analysis
ps5.callus <- filter_taxa(psra.callus, function(x) sum(x) > .005, TRUE) # filtering

keepcallus<- taxa_names(ps5.callus)
ps.callus.counts <- prune_taxa(keepcallus, ps.callus)

### merging the data ##

ps.merge.clean <- merge_phyloseq(ps5.stem, ps5.seed, ps5.embryo, ps5.callus)

ps.merge.counts <- merge_phyloseq(ps.seed.counts, ps.stem.counts, ps.embryo.counts, ps.callus.counts)

## transforming to Hillenger data ##

ps.merge.helt <- ps.merge.clean
otu_table(ps.merge.helt) <-otu_table(decostand(otu_table(ps.merge.clean), method = "hellinger"), taxa_are_rows=FALSE)

## alpha diversity analysis ##

alpha <- data.frame(
  "Observed" = estimate_richness(ps.merge.counts, measures = "Observed"),
  "Shannon" = estimate_richness(ps.merge.helt, measures = "Shannon"),  
  "Sample_type" = sample_data(ps.merge.counts)$Sample_type, 
  "Origin" = sample_data(ps.merge.counts)$Origin, 
  "Treatment" = sample_data(ps.merge.counts)$Treatment,
  "Time_point" = sample_data(ps.merge.counts)$Time_point)
head(alpha)

mean_summery_alpha <- alpha %>%
  group_by(Origin, Treatment) %>%
  dplyr::summarise(mean_observed = mean(Observed),
                   mean_shannon = mean(Shannon), .groups = 'keep', 
                   sd_observed = sd(Observed), sd_shannon = sd(Shannon))

T0_Alpha <- read.csv("Alpha_T0.csv", header = TRUE)
T1_Alpha <- read.csv("Alpha_T1.csv", header = TRUE)
seed_Alpha <- read.csv("Alpha_seed.csv", header = TRUE)
seed_fungicide_Alpha <- read.csv("Alpha_seed_fungicide.csv", header = TRUE)
embryo_Alpha <- read.csv("Alpha_embryo.csv", header = TRUE)
callus_Alpha <- read.csv("Alpha_callus.csv", header = TRUE)

kruskal.test(Observed ~ Treatment, T0_Alpha)
pairwise.wilcox.test(T0_Alpha$Observed, T0_Alpha$Treatment, p.adjust.method="BH")
kruskal.test(Shannon ~ Treatment, T0_Alpha)
pairwise.wilcox.test(T0_Alpha$Shannon, T0_Alpha$Treatment, p.adjust.method="BH")

kruskal.test(Observed ~ Treatment, T1_Alpha)
pairwise.wilcox.test(T1_Alpha$Observed, T1_Alpha$Treatment, p.adjust.method="BH")
kruskal.test(Shannon ~ Treatment, T1_Alpha)
pairwise.wilcox.test(T1_Alpha$Shannon, T1_Alpha$Treatment, p.adjust.method="BH")

wilcox.test(Observed ~ Treatment, data = seed_Alpha, p.adjust.method="BH")
wilcox.test(Shannon ~ Treatment, data = seed_Alpha, p.adjust.method="BH") 

kruskal.test(Observed ~ Treatment, seed_fungicide_Alpha)
pairwise.wilcox.test(seed_fungicide_Alpha$Observed, seed_fungicide_Alpha$Treatment, p.adjust.method="BH")
kruskal.test(Shannon ~ Treatment, seed_fungicide_Alpha)
pairwise.wilcox.test(seed_fungicide_Alpha$Shannon, seed_fungicide_Alpha$Treatment, p.adjust.method="BH")

kruskal.test(Observed ~ Treatment, embryo_Alpha)
pairwise.wilcox.test(embryo_Alpha$Observed, embryo_Alpha$Treatment, p.adjust.method="BH")
kruskal.test(Shannon ~ Treatment, embryo_Alpha)
pairwise.wilcox.test(embryo_Alpha$Shannon, embryo_Alpha$Treatment, p.adjust.method="BH")

kruskal.test(Observed ~ Treatment, callus_Alpha)
pairwise.wilcox.test(callus_Alpha$Observed, callus_Alpha$Treatment, p.adjust.method="BH")
kruskal.test(Shannon ~ Treatment, callus_Alpha)
pairwise.wilcox.test(callus_Alpha$Shannon, callus_Alpha$Treatment, p.adjust.method="BH")

## beta diversity ##

ordps.bray <- ordinate(ps.merge.helt, method = "PCoA", distance = "bray")

## PCOA by tissue ##

ps.origin <- subset_samples(ps.merge.helt, Tissue == "Origin", drop_zeroes = TRUE)
ordps.ps.origin <- ordinate(ps.origin , method = "PCoA", distance = "bray")
plot_ordination(ps.origin, ordps.ps.origin, color = "Treatment") +  stat_ellipse(aes(group = Treatment)

ps.stem2 <- subset_samples(ps.merge.helt, Tissue == "Stem", drop_zeroes = TRUE)
ordps.ps.stem2 <- ordinate(ps.stem2 , method = "PCoA", distance = "bray")
plot_ordination(ps.stem2, ordps.ps.stem2, color = "Treatment") +  coord_equal() +  stat_ellipse(aes(group = Treatment))

## pairwise adonis according to: https://github.com/pmartinezarbizu/pairwiseAdonis Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4 ##

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='BH',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      } 
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
    F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
} 

O.T0 <- subset_samples(ps.merge.helt, Time_point == "T0", drop_zeroes = TRUE) ## T0- source tissues ##
O.T0.bray <- phyloseq::distance(O.T0, method = "bray")
sampledf.bray.O.T0 <- data.frame(sample_data(O.T0))
adonis(O.T0.bray ~ Treatment, data = sampledf.bray.O.T0)
pairwise.adonis(O.T0.bray, sample_data(O.T0)$Treatment)

O.T1 <- subset_samples(ps.merge.helt, Time_point == "T1", drop_zeroes = TRUE) ## T1- stems ##
O.T1.bray <- phyloseq::distance(O.T1, method = "bray")
sampledf.bray.O.T1 <- data.frame(sample_data(O.T1))
adonis(O.T1.bray ~ Treatment, data = sampledf.bray.O.T1)
pairwise.adonis(O.T1.bray, sample_data(O.T1)$Treatment)

## diferential abundance origin tissue ##

# Loop t-test for all taxa
results <- data.frame()
for (i in 1:length(taxa_names(ps.origin))){
  
  asv <- taxa_names(ps.origin)[i]
  
  Control_seed_helt <- mean(otu_table(subset_samples(ps.origin, Treatment == "Control_seed"))[,i])
  Seed_fungicides_helt <- mean(otu_table(subset_samples(ps.origin, Treatment == "Seed_fungicides"))[,i])
  Embryo_helt <- mean(otu_table(subset_samples(ps.origin, Treatment == "Embryo"))[,i])
  Callus_helt <- mean(otu_table(subset_samples(ps.origin, Treatment == "Callus"))[,i])
  
  pval <- kruskal.test(as.vector(otu_table(ps.origin)[,i]) ~ sample_data(ps.origin)$Treatment)$p.value
  results <- rbind(results, data.frame(OTU = asv, 
                                       Control_seed_Mean_HELT = Control_seed_helt, 
                                       Seed_fungicides_Mean_HELT = Seed_fungicides_helt,
                                       Embryo_Mean_HELT = Embryo_helt,
                                       Callus_Mean_HELT = Callus_helt,
                                       P = pval))
}


# Add taxa designations to table
resultsorigin <- inner_join(results, rownames_to_column(data.frame(tax_table(ps.origin)), "OTU"), by = "OTU")
resultsorigin

# Computing FDR corrected p-values and filtering ## 
kruskal_resultsorigin <- resultsorigin %>%
  arrange(P) %>%
  mutate(BH_FDR = p.adjust(P, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, P, BH_FDR, everything())
write.csv(kruskal_resultsorigin, "kruskalorigin.csv")

kruskal_resultsorigin_noBH <- resultsorigin %>%
  arrange(P) %>%
  mutate(BH_FDR = p.adjust(P, "BH")) %>%
  dplyr::select(OTU, P, BH_FDR, everything())
write.csv(kruskal_resultsorigin_noBH, "kruskaloriginnoBHfilter.csv")

## pairwise wilcox ##
ps.Cladosporium_sphaerospermum <- subset_taxa(ps.origin, Species == "Cladosporium_sphaerospermum")
melt.ps.Cladosporium_sphaerospermum <- psmelt(ps.Cladosporium_sphaerospermum)
pairwise.wilcox.test(melt.ps.Cladosporium_sphaerospermum$Abundance,melt.ps.Cladosporium_sphaerospermum$Treatment, p.adjust.method="BH")

ps.Mycosphaerella_tassiana <- subset_taxa(ps.origin, Species == "Mycosphaerella_tassiana")
melt.ps.Mycosphaerella_tassiana <- psmelt(ps.Mycosphaerella_tassiana)
pairwise.wilcox.test(melt.ps.Mycosphaerella_tassiana$Abundance,melt.ps.Mycosphaerella_tassiana$Treatment, p.adjust.method="BH")

ps.Alternaria_infectoria <- subset_taxa(ps.origin, Species == "Alternaria_infectoria")
melt.ps.Alternaria_infectoria <- psmelt(ps.Alternaria_infectoria)
pairwise.wilcox.test(melt.ps.Alternaria_infectoria$Abundance,melt.ps.Alternaria_infectoria$Treatment, p.adjust.method="BH")

ps.Cladosporium_halotolerans <- subset_taxa(ps.origin, Species == "Cladosporium_halotolerans")
melt.ps.Cladosporium_halotolerans <- psmelt(ps.Cladosporium_halotolerans)
pairwise.wilcox.test(melt.ps.Cladosporium_halotolerans$Abundance,melt.ps.Cladosporium_halotolerans$Treatment, p.adjust.method="BH")

## diferential abundance stems ##

# Loop t-test for all taxa
resultsstems <- data.frame()
for (i in 1:length(taxa_names(ps.stem2))){
  
  asv <- taxa_names(ps.stem2)[i]
  
  Callus_stem_MS_fungicides_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Callus_stem_MS_fungicides"))[,i])
  Control_stem_MS_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Control_stem_MS"))[,i])
  Callus_stem_MS_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Callus_stem_MS"))[,i])
  Embryo_stem_MS_fungicides_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Embryo_stem_MS_fungicides"))[,i])
  Embryo_stem_MS_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Embryo_stem_MS"))[,i])
  Seed_fungicides_stem_MS_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Seed_fungicides_stem_MS"))[,i])
  Seed_fungicides_stem_MS_fungicides_helt <- mean(otu_table(subset_samples(ps.stem2, Treatment == "Seed_fungicides_stem_MS_fungicides"))[,i])
  
  pval <- kruskal.test(as.vector(otu_table(ps.stem2)[,i]) ~ sample_data(ps.stem2)$Treatment)$p.value
  resultsstems <- rbind(resultsstems, data.frame(OTU = asv, 
                                                 Callus_stem_MS_fungicides_Mean = Callus_stem_MS_fungicides_helt, 
                                                 Control_stem_MS_Mean = Control_stem_MS_helt,
                                                 Callus_stem_MS_Mean = Callus_stem_MS_helt,
                                                 Embryo_stem_MS_fungicides_Mean = Embryo_stem_MS_fungicides_helt,
                                                 Embryo_stem_MS_Mean = Embryo_stem_MS_helt,
                                                 Seed_fungicides_stem_MS_Mean = Seed_fungicides_stem_MS_helt,
                                                 Seed_fungicides_stem_MS_fungicides_Mean = Seed_fungicides_stem_MS_fungicides_helt,
                                                 P = pval))
}


# Add taxa designations to table
resultsstems <- inner_join(resultsstems, rownames_to_column(data.frame(tax_table(ps.stems2)), "OTU"), by = "OTU")
resultsstems

# Computing FDR corrected p-values and filtering ## 
kruskal_resultsstems <- resultsstems %>%
  arrange(P) %>%
  mutate(BH_FDR = p.adjust(P, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, P, BH_FDR, everything())
write.csv(kruskal_resultsstems, "kruskalstems.csv")

kruskal_resultsstems_noBH <- resultsstems %>%
  arrange(P) %>%
  mutate(BH_FDR = p.adjust(P, "BH")) %>%
  dplyr::select(OTU, P, BH_FDR, everything())
write.csv(kruskal_resultsstems_noBH, "kruskalstemoBHfilter.csv")


## evenness analysis ##

data_evenness <- data.frame(
  evenness(ps.merge.helt, index = "pielou", zeroes = TRUE, detection = 0), "Treatment" = sample_data(ps.merge.helt)$Treatment, "Origin" = sample_data(ps.merge.helt)$Origin)
head(data_evenness)
write.csv(data_evenness, "Evenness_first_run.csv")

data_evenness %>%
  gather(key = metric, value = value, c("pielou")) %>%
  mutate(metric = factor(metric, levels = c("pielou"))) %>%
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "", y = "")

data_evenness_stems %>%
  group_by(Treatment) %>%
  dplyr::summarise(mean_pielou = mean(pielou), 
                   sd_pielou = sd(pielou))

evenness_origin <- read.csv("evenness_origin.csv", header = TRUE)
evenness_stems <- read.csv("evenness_stems.csv", header = TRUE)

kruskal.test(pielou ~ Treatment, evenness_origin)
pairwise.wilcox.test(evenness_origin$pielou, evenness_origin$Treatment, p.adjust.method="BH")

kruskal.test(pielou ~ Treatment, evenness_stems)
pairwise.wilcox.test(evenness_stems$pielou, evenness_stems$Treatment, p.adjust.method="BH")


data_evenness_stems <- data.frame(
  evenness(ps.stem2, index = "pielou", zeroes = TRUE, detection = 0), "Treatment" = sample_data(ps.stem2)$Treatment)
head(data_evenness)
kruskal.test(pielou ~ Treatment, data_evenness_stems)
pairwise.wilcox.test(data_evenness_stems$pielou, data_evenness_stems$Treatment, p.adjust.method="BH")

## Venn diagram and upsetplot##

Controlseed <- list("Alternaria_infectoria",
                    "Parastagonospora_ASV4",
                    "Alternaria_alternata",
                    "Stemphylium_ASV6",
                    "Tremellales_ASV267",
                    "Auricularia_auricula-judae",
                    "Coriolopsis_gallica",
                    "Peniophora_versiformis")
ControlstemMS <- list("Alternaria_infectoria",
                      "Alternaria_alternata",
                      "Stemphylium_ASV6",
                      "Cladosporium_sphaerospermum",
                      "Mycosphaerella_tassiana",
                      "Sporobolomyces_roseus",
                      "Cladosporium_cladosporioides",
                      "Phaeosphaeria_juncicola",
                      "Cladosporium_halotolerans",
                      "Coriolopsis_gallica",
                      "Blumeria_graminis",
                      "Nigrospora_oryzae",
                      "Symmetrospora_coprosmae",
                      "Cladosporium_ramotenellum",
                      "Wallemia_ASV73",
                      "Debaryomyces_hansenii",
                      "Cladosporium_ASV8",
                      "Dioszegia_hungarica",
                      "Erythrobasidium_hasegawianum",
                      "Preussia_flanaganii",
                      "Aureobasidium_pullulans",
                      "Cladosporium_exasperatum",
                      "Saxophila_tyrrhenica",
                      "Coprinellus_micaceus",
                      "Peniophora_ASV55",
                      "Cryptomarasmius_corbariensis",
                      "Acremonium_sclerotigenum",
                      "Lepista_sordida",
                      "Penicillium_astrolabium",
                      "Botrytis_cinerea",
                      "Rhodotorula_graminis",
                      "Phaeosphaeria_ASV107",
                      "Coniosporium_apollinis",
                      "Knufia_perforans",
                      "Crocicreas_ASV92",
                      "Phialocephala_fluminis",
                      "Cladosporium_dominicanum",
                      "Cystobasidium_lysinophilum",
                      "Arthrinium_marii",
                      "Schizophyllum_commune",
                      "Wallemia_hederae",
                      "Irpex_lacteus",
                      "Vishniacozyma_carnescens",
                      "Neodevriesiaceae_ASV93",
                      "Septoriella_phragmitis",
                      "Coniothyrium_sidae")
Seedfungicides <- list("Alternaria_infectoria",
                       "Alternaria_alternata",
                       "Cladosporium_cladosporioides",
                       "Cladosporium_ASV8",
                       "Schizophyllum_commune",
                       "Sporobolomyces_roseus",
                       "Stemphylium_ASV6",
                       "Coriolopsis_gallica")
SeedfungicidesstemMS <- list("Alternaria_infectoria",
                             "Cladosporium_sphaerospermum",
                             "Mycosphaerella_tassiana",
                             "Sporobolomyces_roseus",
                             "Cladosporium_cladosporioides",
                             "Alternaria_alternata",
                             "Preussia_flanaganii",
                             "Cladosporium_ASV8",
                             "Cladosporium_ramotenellum",
                             "Spiroplana_centripeta",
                             "Coprinellus_micaceus",
                             "Filobasidium_ASV78",
                             "Symmetrospora_coprosmae",
                             "Candida_heveicola",
                             "Cutaneotrichosporon_curvatus",
                             "Peniophora_ASV164",
                             "Cystobasidium_lysinophilum",
                             "Botrytis_cinerea",
                             "Naganishia_uzbekistanensis",
                             "Acremonium_sclerotigenum",
                             "Blumeria_graminis",
                             "Nectria_pseudotrichia",
                             "Peniophora_ASV195",
                             "Stemphylium_ASV6",
                             "Candida_sake",
                             "Arthrinium_marii",
                             "Alternaria_terricola",
                             "Hypoxylon_begae",
                             "Wallemia_ASV73",
                             "Cryptomarasmius_corbariensis",
                             "Rhodotorula_mucilaginosa",
                             "Devriesia_sardiniae",
                             "Peniophora_ASV34",
                             "Ramularia_eucalypti",
                             "Cladosporium_ASV650",
                             "Debaryomyces_hansenii",
                             "Beauveria_bassiana",
                             "Coriolopsis_gallica",
                             "Torula_hollandica",
                             "Neodevriesiaceae_ASV93",
                             "Wallemia_sebi",
                             "Septoriella_phragmitis",
                             "Capnodiales_ASV486",
                             "Pestalotiopsis_microspora",
                             "Penicillium_astrolabium",
                             "Candida_orthopsilosis",
                             "Vishniacozyma_taibaiensis",
                             "Erysiphe_pisi",
                             "Sarocladium_bactrocephalum",
                             "Amylosporus_campbellii",
                             "Cladosporium_halotolerans",
                             "Idriella_rara",
                             "Aspergillus_intermedius",
                             "Gibberella_intricans",
                             "Trichothecium_roseum",
                             "Wallemia_hederae",
                             "Vishniacozyma_carnescens",
                             "Penicillium_kongii")
SeedfungicidesstemMSfungicides <- list("Alternaria_infectoria",
                                       "Alternaria_alternata",
                                       "Cladosporium_sphaerospermum",
                                       "Botrytis_ASV45",
                                       "Mycosphaerella_tassiana",
                                       "Neodevriesia_ASV321",
                                       "Sporobolomyces_roseus",
                                       "Cladosporium_ramotenellum",
                                       "Stemphylium_ASV127",
                                       "Cladosporium_cladosporioides",
                                       "Alternaria_chlamydosporigena",
                                       "Parengyodontium_album",
                                       "Amyloporia_sinuosa",
                                       "Phaeosphaeria_juncicola",
                                       "Rigidoporus_vinctus",
                                       "Toxicocladosporium_strelitziae",
                                       "Sarocladium_kiliense",
                                       "Peniophora_albobadia",
                                       "Phaeosphaeria_ASV107",
                                       "Cladosporium_exasperatum",
                                       "Venturia_atriseda",
                                       "Stemphylium_ASV6",
                                       "Blumeria_graminis",
                                       "Cladosporium_ASV8",
                                       "Inocutis_tamaricis",
                                       "Coniothyrium_ASV405",
                                       "Peniophora_ASV55",
                                       "Hortaea_werneckii",
                                       "Alternaria_oregonensis",
                                       "Acremonium_polychromum",
                                       "Arthrinium_marii",
                                       "Nectria_pseudotrichia",
                                       "Aspergillus_intermedius",
                                       "Penicillium_astrolabium")
Embryo <- list("Alternaria_alternata",
               "Sporobolomyces_roseus",
               "Cladosporium_sphaerospermum",
               "Mycosphaerella_tassiana",
               "Aureobasidium_pullulans",
               "Candida_tropicalis",
               "Cladosporium_halotolerans",
               "Debaryomyces_hansenii",
               "Cladosporium_ASV8",
               "Alternaria_infectoria",
               "Cladosporium_cladosporioides",
               "Stemphylium_ASV6",
               "Ramularia_eucalypti",
               "Coriolopsis_gallica",
               "Cladosporium_ramotenellum",
               "Symmetrospora_coprosmae",
               "Antrodia_serialiformis",
               "Cladosporium_dominicanum",
               "Cryptococcus_laurentii",
               "Botrytis_cinerea",
               "Peniophora_ASV189",
               "Filobasidium_chernovii",
               "Naganishia_diffluens",
               "Toxicocladosporium_irritans",
               "Sarocladium_implicatum",
               "Ganoderma_adspersum",
               "Golovinomyces_ASV178",
               "Chaetomium_angustispirale",
               "Neosetophoma_ASV118",
               "Preussia_flanaganii",
               "Peniophora_ASV86",
               "Sporobolomyces_salicinus",
               "Cladosporium_exasperatum",
               "Penicillium_astrolabium",
               "Uncobasidium_ASV354",
               "Alternaria_lolii",
               "Filobasidium_oeirense",
               "Parastagonospora_ASV4",
               "Wallemia_ASV73",
               "Alternaria_terricola",
               "Nectria_pseudotrichia",
               "Acremonium_sclerotigenum",
               "Bipolaris_axonopicola",
               "Coprinellus_ASV363",
               "Pseudoidium_neolycopersici",
               "Verticillium_leptobactrum",
               "Hyaloseta_nolinae",
               "Botrytis_caroliniana",
               "Pyrenophora_japonica",
               "Fungi_ASV441",
               "Septoriella_phragmitis",
               "Fungi_ASV521",
               "Stemphylium_vesicarium",
               "Filobasidium_wieringae",
               "Aspergillus_ruber",
               "Verticillium_dahliae",
               "Hypomyces_ASV522")
EmbryostemMS <- list("Alternaria_infectoria",
                     "Cladosporium_cladosporioides",
                     "Mycosphaerella_tassiana",
                     "Cladosporium_sphaerospermum",
                     "Cladosporium_ramotenellum",
                     "Gibberella_intricans",
                     "Cladosporium_ASV8",
                     "Sporobolomyces_roseus",
                     "Debaryomyces_hansenii",
                     "Alternaria_alternata",
                     "Botrytis_cinerea",
                     "Nectria_pseudotrichia",
                     "Filobasidium_chernovii",
                     "Aureobasidium_pullulans",
                     "Epicoccum_nigrum",
                     "Blumeria_graminis",
                     "Filobasidium_ASV78",
                     "Peniophora_ASV55",
                     "Cladosporium_exasperatum",
                     "Schizophyllum_commune",
                     "Wallemia_hederae",
                     "Parengyodontium_album",
                     "Botrytis_caroliniana",
                     "Trichaptum_biforme",
                     "Zymoseptoria_brevis",
                     "Candida_heveicola",
                     "Dioszegia_ASV346",
                     "Coniosporium_apollinis",
                     "Coniothyrium_sidae",
                     "Phaeosphaeria_ASV107",
                     "Stemphylium_ASV6",
                     "Filobasidium_ASV25",
                     "Wallemia_ASV73",
                     "Pterulaceae_ASV229",
                     "Cystofilobasidium_capitatum",
                     "Neoascochyta_europaea",
                     "Cladosporium_halotolerans",
                     "Simocybe_haustellaris",
                     "Cyberlindnera_jadinii",
                     "Cryptomarasmius_corbariensis",
                     "Xylariales_ASV314",
                     "Lecanicillium_ASV261",
                     "Neodevriesiaceae_ASV93",
                     "Golovinomyces_ASV178",
                     "Coprinellus_micaceus",
                     "Erysiphe_polygoni",
                     "Peniophora_nuda",
                     "Monodictys_castaneae",
                     "Cystobasidium_slooffiae",
                     "Meruliaceae_ASV245",
                     "Alternaria_oregonensis",
                     "Alternaria_chlamydosporigena",
                     "Aspergillus_intermedius",
                     "Stemphylium_ASV127",
                     "Tritirachium_oryzae",
                     "Stemphylium_vesicarium",
                     "Trichothecium_roseum",
                     "Sporobolomyces_ASV450",
                     "Scytalidium_circinatum",
                     "Erythrobasidium_hasegawianum",
                     "Penicillium_astrolabium",
                     "Leptospora_ASV165",
                     "Peniophora_ASV546",
                     "Arthrinium_marii",
                     "Nigrospora_oryzae",
                     "Hyphodermella_rosae",
                     "Zopfiella_lundqvistii",
                     "Torula_ASV213",
                     "Knufia_ASV505",
                     "Chaetomium_madrasense",
                     "Phaeosphaeria_juncicola",
                     "Issatchenkia_orientalis",
                     "Pseudallescheria_boydii",
                     "Vishniacozyma_heimaeyensis",
                     "Stagonospora_ASV589",
                     "Acremonium_sclerotigenum",
                     "Cladosporium_dominicanum",
                     "Buckleyzyma_aurantiaca",
                     "Vishniacozyma_carnescens",
                     "Aspergillus_versicolor",
                     "Chaetomium_angustispirale",
                     "Buckleyzyma_phyllomatis")
EmbryostemMSfungicides <- list("Alternaria_infectoria",
                               "Cladosporium_sphaerospermum",
                               "Alternaria_alternata",
                               "Cladosporium_cladosporioides",
                               "Mycosphaerella_tassiana",
                               "Sporobolomyces_roseus",
                               "Botrytis_ASV45",
                               "Aureobasidium_pullulans",
                               "Cladosporium_halotolerans",
                               "Schizophyllum_commune",
                               "Botrytis_cinerea",
                               "Chaetomium_angustispirale",
                               "Didymella_urticicola",
                               "Coprinopsis_foetidellus",
                               "Cladosporium_exasperatum",
                               "Candida_sake",
                               "Blumeria_graminis",
                               "Monodictys_castaneae",
                               "Phaeosphaeria_juncicola",
                               "Irpex_lacteus",
                               "Phaeosphaeria_ASV107",
                               "Wallemia_hederae",
                               "Sarocladium_kiliense",
                               "Coriolopsis_gallica",
                               "Golovinomyces_ASV178",
                               "Naganishia_diffluens",
                               "Neodevriesia_ASV366",
                               "Cryptococcus_saitoi",
                               "Ramularia_eucalypti",
                               "Filobasidium_chernovii",
                               "Beauveria_bassiana",
                               "Coprinopsis_ASV367",
                               "Meruliaceae_ASV454",
                               "Cladosporium_dominicanum",
                               "Symmetrospora_coprosmae",
                               "Acremonium_charticola",
                               "Cladosporium_ramotenellum",
                               "Cutaneotrichosporon_curvatus",
                               "Wallemia_ASV88",
                               "Parengyodontium_album",
                               "Lophiostoma_ASV565",
                               "Devriesia_ASV497",
                               "Wallemia_sebi",
                               "Saxophila_tyrrhenica",
                               "Discosia_pseudoartocreas",
                               "Alternaria_lolii",
                               "Mycosphaerellaceae_ASV424",
                               "Zopfiella_marina",
                               "Zopfiella_lundqvistii",
                               "Cladosporium_ASV8",
                               "Antrodia_xantha",
                               "Septoriella_phragmitis",
                               "Myriostoma_coliforme",
                               "Dioszegia_hungarica",
                               "Wallemia_ASV276",
                               "Stemphylium_vesicarium",
                               "Verticillium_leptobactrum",
                               "Penicillium_kongii")
Callus <- list("Cladosporium_sphaerospermum",
               "Mycosphaerella_tassiana",
               "Cladosporium_halotolerans",
               "Cladosporium_cladosporioides",
               "Alternaria_infectoria",
               "Alternaria_alternata",
               "Stemphylium_ASV6",
               "Aureobasidium_pullulans",
               "Sporobolomyces_roseus",
               "Nectria_pseudotrichia",
               "Sporormiaceae_ASV68",
               "Botrytis_cinerea",
               "Debaryomyces_hansenii",
               "Sporobolomyces_symmetricus",
               "Schizophyllum_commune",
               "Buckleyzyma_phyllomatis",
               "Penicillium_astrolabium",
               "Coriolopsis_gallica",
               "Filobasidium_oeirense",
               "Cladosporium_ramotenellum",
               "Cladosporium_exasperatum",
               "Cladosporium_ASV8",
               "Hortaea_werneckii",
               "Noosia_banksiae",
               "Cryptomarasmius_corbariensis",
               "Psathyrella_ASV137",
               "Arthrinium_marii",
               "Fusarium_proliferatum",
               "Peniophora_ASV34",
               "Cladosporium_dominicanum",
               "Paradendryphiella_arenariae",
               "Pilatoporus_bondartsevae",
               "Pyrenochaeta_ASV201",
               "Rhodotorula_diobovata",
               "Sporormiaceae_ASV286",
               "Cistella_ASV269",
               "Leucosporidium_yakuticum",
               "Protomyces_inouyei",
               "Pseudogymnoascus_appendiculatus",
               "Alutaceodontia_alutacea",
               "Parastagonospora_ASV4",
               "Wallemia_ASV88",
               "Coprinellus_micaceus",
               "Candida_orthopsilosis",
               "Vishniacozyma_victoriae",
               "Kodamaea_ohmeri",
               "Volvopluteus_gloiocephalus",
               "Lentinus_arcularius",
               "Xylodon_sambuci",
               "Coprinopsis_ASV344",
               "Tritirachium_oryzae",
               "Trimmatostroma_salinum",
               "Buckleyzyma_aurantiaca",
               "Pichia_membranifaciens",
               "Periglandula_ASV240",
               "Acremonium_sclerotigenum",
               "Bartalinia_robillardoides",
               "Stemphylium_ASV320",
               "Kwoniella_ASV323",
               "Leptospora_ASV165",
               "Wallemia_ASV276",
               "Ramularia_eucalypti",
               "Aspergillus_versicolor",
               "Alternaria_terricola",
               "Libertasomyces_myopori",
               "Sarocladium_implicatum",
               "Alternaria_chlamydosporigena",
               "Dipodascaceae_ASV540")
CallusstemMS <- list("Mycosphaerella_tassiana",
                     "Cladosporium_cladosporioides",
                     "Filobasidium_ASV25",
                     "Cladosporium_sphaerospermum",
                     "Sporobolomyces_roseus",
                     "Coriolopsis_gallica",
                     "Cladosporium_ASV8",
                     "Cladosporium_halotolerans",
                     "Alternaria_alternata",
                     "Filobasidium_oeirense",
                     "Cladosporium_ramotenellum",
                     "Blumeria_graminis",
                     "Peniophora_albobadia",
                     "Alternaria_infectoria",
                     "Symmetrospora_foliicola",
                     "Naganishia_diffluens",
                     "Aureobasidium_pullulans",
                     "Naganishia_albida",
                     "Acremonium_charticola",
                     "Auriculariales_ASV53",
                     "Candida_tropicalis",
                     "Lepista_sordida",
                     "Cryptomarasmius_corbariensis",
                     "Coprinopsis_ASV122",
                     "Vishniacozyma_heimaeyensis",
                     "Coprinellus_ASV123",
                     "Steccherinaceae_ASV158",
                     "Botrytis_ASV45",
                     "Stemphylium_ASV6",
                     "Schizophyllum_commune",
                     "Cystobasidium_lysinophilum",
                     "Torula_ASV213",
                     "Omphalotus_olearius",
                     "Phaeosphaeria_ASV107",
                     "Irpex_lacteus",
                     "Tubaria_conspersa",
                     "Sistotrema_brinkmannii",
                     "Peniophora_nuda",
                     "Torula_ASV101",
                     "Verticillium_leptobactrum",
                     "Ceratobasidium_ASV202",
                     "Neosetophoma_ASV118",
                     "Inocutis_tamaricis",
                     "Toxicocladosporium_strelitziae",
                     "Peniophora_ASV34",
                     "Psathyrella_trinitatensis",
                     "Vishniacozyma_victoriae",
                     "Phaeosphaeria_juncicola",
                     "Coprinellus_micaceus",
                     "Lentinus_arcularius",
                     "Wallemia_ASV73",
                     "Filobasidium_ASV138",
                     "Alternaria_chlamydosporigena",
                     "Neodevriesiaceae_ASV93",
                     "Dioszegia_hungarica",
                     "Cystidiodontia_laminifera",
                     "Wallemia_hederae",
                     "Phanerochaete_ASV235",
                     "Fusarium_proliferatum",
                     "Acremonium_sclerotigenum",
                     "Botrytis_cinerea",
                     "Issatchenkia_orientalis",
                     "Leptospora_ASV165",
                     "Ramularia_eucalypti",
                     "Rhodotorula_graminis",
                     "Cladosporium_exasperatum",
                     "Toxicocladosporium_irritans",
                     "Stachybotrys_chartarum",
                     "Dekkera_anomala",
                     "Kodamaea_ohmeri",
                     "Aspergillus_intermedius",
                     "Sistotrema_oblongisporum",
                     "Cystobasidium_pinicola",
                     "Xylariales_ASV314",
                     "Leucosporidium_yakuticum",
                     "Arthrinium_marii",
                     "Cladosporium_dominicanum",
                     "Nectria_pseudotrichia",
                     "Candida_sake",
                     "Debaryomyces_hansenii",
                     "Hortaea_werneckii",
                     "Wallemia_ASV88",
                     "Symmetrospora_coprosmae",
                     "Alternaria_terricola",
                     "Vishniacozyma_carnescens",
                     "Penicillium_kongii",
                     "Hyphodermella_rosae",
                     "Penicillium_astrolabium")
CallusstemMSfungicides <- list("Alternaria_alternata",
                               "Mycosphaerella_tassiana",
                               "Sporobolomyces_roseus",
                               "Cladosporium_ramotenellum",
                               "Cladosporium_sphaerospermum",
                               "Cladosporium_cladosporioides",
                               "Wallemia_hederae",
                               "Alternaria_infectoria",
                               "Coriolopsis_gallica",
                               "Cladosporium_ASV8",
                               "Gibberella_intricans",
                               "Stemphylium_ASV6",
                               "Phaeosphaeria_juncicola",
                               "Naganishia_diffluens",
                               "Preussia_flanaganii",
                               "Clonostachys_miodochialis",
                               "Knufia_perforans",
                               "Sistotremastrum_guttuliferum",
                               "Coprinellus_micaceus",
                               "Botrytis_cinerea",
                               "Blumeria_graminis",
                               "Preussia_persica",
                               "Cladosporium_halotolerans",
                               "Vishniacozyma_carnescens",
                               "Parastagonospora_ASV4",
                               "Buckleyzyma_phyllomatis",
                               "Phialocephala_fluminis",
                               "Aureobasidium_pullulans",
                               "Acrodontium_crateriforme",
                               "Epicoccum_nigrum",
                               "Symmetrospora_coprosmae",
                               "Crocicreas_ASV92",
                               "Trichosporonaceae_ASV550",
                               "Plectosphaerellaceae_ASV125",
                               "Cladosporium_dominicanum",
                               "Coprinellus_radians",
                               "Wallemia_ASV73",
                               "Auriculariales_ASV53",
                               "Holtermanniella_wattica",
                               "Schizophyllum_commune",
                               "Buckleyzyma_aurantiaca",
                               "Amyloporia_sinuosa",
                               "Rhodotorula_mucilaginosa",
                               "Candida_heveicola",
                               "Hansfordia_ASV607",
                               "Arthrinium_marii",
                               "Fusariella_ASV200",
                               "Cladosporium_exasperatum",
                               "Zygosporium_ASV128",
                               "Neoerysiphe_nevoi",
                               "Candida_railenensis",
                               "Capnobotryella_ASV368",
                               "Sarocladium_implicatum",
                               "Ramularia_eucalypti",
                               "Aspergillus_intermedius",
                               "Stemphylium_vesicarium",
                               "Auricularia_auricula-judae",
                               "Porodisculus_pendulus",
                               "Candida_sake",
                               "Saxophila_tyrrhenica",
                               "Fonsecaea_ASV389",
                               "Peniophora_ASV407",
                               "Neodevriesia_ASV321",
                               "Alternaria_terricola",
                               "Golovinomyces_montagnei",
                               "Aspergillus_versicolor",
                               "Trichaptum_biforme",
                               "Toxicocladosporium_irritans",
                               "Phaeosphaeria_ASV107",
                               "Penicillium_kongii",
                               "Nectria_pseudotrichia",
                               "Irpex_lacteus")

stems = list(ControlstemMS = c("Alternaria_infectoria",
                               "Alternaria_alternata",
                               "Stemphylium_ASV6",
                               "Cladosporium_sphaerospermum",
                               "Mycosphaerella_tassiana",
                               "Sporobolomyces_roseus",
                               "Cladosporium_cladosporioides",
                               "Phaeosphaeria_juncicola",
                               "Cladosporium_halotolerans",
                               "Coriolopsis_gallica",
                               "Blumeria_graminis",
                               "Nigrospora_oryzae",
                               "Symmetrospora_coprosmae",
                               "Cladosporium_ramotenellum",
                               "Wallemia_ASV73",
                               "Debaryomyces_hansenii",
                               "Cladosporium_ASV8",
                               "Dioszegia_hungarica",
                               "Erythrobasidium_hasegawianum",
                               "Preussia_flanaganii",
                               "Aureobasidium_pullulans",
                               "Cladosporium_exasperatum",
                               "Saxophila_tyrrhenica",
                               "Coprinellus_micaceus",
                               "Peniophora_ASV55",
                               "Cryptomarasmius_corbariensis",
                               "Acremonium_sclerotigenum",
                               "Lepista_sordida",
                               "Penicillium_astrolabium",
                               "Botrytis_cinerea",
                               "Rhodotorula_graminis",
                               "Phaeosphaeria_ASV107",
                               "Coniosporium_apollinis",
                               "Knufia_perforans",
                               "Crocicreas_ASV92",
                               "Phialocephala_fluminis",
                               "Cladosporium_dominicanum",
                               "Cystobasidium_lysinophilum",
                               "Arthrinium_marii",
                               "Schizophyllum_commune",
                               "Wallemia_hederae",
                               "Irpex_lacteus",
                               "Vishniacozyma_carnescens",
                               "Neodevriesiaceae_ASV93",
                               "Septoriella_phragmitis",
                               "Coniothyrium_sidae"), 
             SeedfungicidesstemMS = c("Alternaria_infectoria",
                                      "Cladosporium_sphaerospermum",
                                      "Mycosphaerella_tassiana",
                                      "Sporobolomyces_roseus",
                                      "Cladosporium_cladosporioides",
                                      "Alternaria_alternata",
                                      "Preussia_flanaganii",
                                      "Cladosporium_ASV8",
                                      "Cladosporium_ramotenellum",
                                      "Spiroplana_centripeta",
                                      "Coprinellus_micaceus",
                                      "Filobasidium_ASV78",
                                      "Symmetrospora_coprosmae",
                                      "Candida_heveicola",
                                      "Cutaneotrichosporon_curvatus",
                                      "Peniophora_ASV164",
                                      "Cystobasidium_lysinophilum",
                                      "Botrytis_cinerea",
                                      "Naganishia_uzbekistanensis",
                                      "Acremonium_sclerotigenum",
                                      "Blumeria_graminis",
                                      "Nectria_pseudotrichia",
                                      "Peniophora_ASV195",
                                      "Stemphylium_ASV6",
                                      "Candida_sake",
                                      "Arthrinium_marii",
                                      "Alternaria_terricola",
                                      "Hypoxylon_begae",
                                      "Wallemia_ASV73",
                                      "Cryptomarasmius_corbariensis",
                                      "Rhodotorula_mucilaginosa",
                                      "Devriesia_sardiniae",
                                      "Peniophora_ASV34",
                                      "Ramularia_eucalypti",
                                      "Cladosporium_ASV650",
                                      "Debaryomyces_hansenii",
                                      "Beauveria_bassiana",
                                      "Coriolopsis_gallica",
                                      "Torula_hollandica",
                                      "Neodevriesiaceae_ASV93",
                                      "Wallemia_sebi",
                                      "Septoriella_phragmitis",
                                      "Capnodiales_ASV486",
                                      "Pestalotiopsis_microspora",
                                      "Penicillium_astrolabium",
                                      "Candida_orthopsilosis",
                                      "Vishniacozyma_taibaiensis",
                                      "Erysiphe_pisi",
                                      "Sarocladium_bactrocephalum",
                                      "Amylosporus_campbellii",
                                      "Cladosporium_halotolerans",
                                      "Idriella_rara",
                                      "Aspergillus_intermedius",
                                      "Gibberella_intricans",
                                      "Trichothecium_roseum",
                                      "Wallemia_hederae",
                                      "Vishniacozyma_carnescens",
                                      "Penicillium_kongii"), 
             SeedfungicidesstemMSfungicides = c("Alternaria_infectoria",
                                                "Alternaria_alternata",
                                                "Cladosporium_sphaerospermum",
                                                "Botrytis_ASV45",
                                                "Mycosphaerella_tassiana",
                                                "Neodevriesia_ASV321",
                                                "Sporobolomyces_roseus",
                                                "Cladosporium_ramotenellum",
                                                "Stemphylium_ASV127",
                                                "Cladosporium_cladosporioides",
                                                "Alternaria_chlamydosporigena",
                                                "Parengyodontium_album",
                                                "Amyloporia_sinuosa",
                                                "Phaeosphaeria_juncicola",
                                                "Rigidoporus_vinctus",
                                                "Toxicocladosporium_strelitziae",
                                                "Sarocladium_kiliense",
                                                "Peniophora_albobadia",
                                                "Phaeosphaeria_ASV107",
                                                "Cladosporium_exasperatum",
                                                "Venturia_atriseda",
                                                "Stemphylium_ASV6",
                                                "Blumeria_graminis",
                                                "Cladosporium_ASV8",
                                                "Inocutis_tamaricis",
                                                "Coniothyrium_ASV405",
                                                "Peniophora_ASV55",
                                                "Hortaea_werneckii",
                                                "Alternaria_oregonensis",
                                                "Acremonium_polychromum",
                                                "Arthrinium_marii",
                                                "Nectria_pseudotrichia",
                                                "Aspergillus_intermedius",
                                                "Penicillium_astrolabium"), 
             EmbryostemMS = c("Alternaria_infectoria",
                              "Cladosporium_cladosporioides",
                              "Mycosphaerella_tassiana",
                              "Cladosporium_sphaerospermum",
                              "Cladosporium_ramotenellum",
                              "Gibberella_intricans",
                              "Cladosporium_ASV8",
                              "Sporobolomyces_roseus",
                              "Debaryomyces_hansenii",
                              "Alternaria_alternata",
                              "Botrytis_cinerea",
                              "Nectria_pseudotrichia",
                              "Filobasidium_chernovii",
                              "Aureobasidium_pullulans",
                              "Epicoccum_nigrum",
                              "Blumeria_graminis",
                              "Filobasidium_ASV78",
                              "Peniophora_ASV55",
                              "Cladosporium_exasperatum",
                              "Schizophyllum_commune",
                              "Wallemia_hederae",
                              "Parengyodontium_album",
                              "Botrytis_caroliniana",
                              "Trichaptum_biforme",
                              "Zymoseptoria_brevis",
                              "Candida_heveicola",
                              "Dioszegia_ASV346",
                              "Coniosporium_apollinis",
                              "Coniothyrium_sidae",
                              "Phaeosphaeria_ASV107",
                              "Stemphylium_ASV6",
                              "Filobasidium_ASV25",
                              "Wallemia_ASV73",
                              "Pterulaceae_ASV229",
                              "Cystofilobasidium_capitatum",
                              "Neoascochyta_europaea",
                              "Cladosporium_halotolerans",
                              "Simocybe_haustellaris",
                              "Cyberlindnera_jadinii",
                              "Cryptomarasmius_corbariensis",
                              "Xylariales_ASV314",
                              "Lecanicillium_ASV261",
                              "Neodevriesiaceae_ASV93",
                              "Golovinomyces_ASV178",
                              "Coprinellus_micaceus",
                              "Erysiphe_polygoni",
                              "Peniophora_nuda",
                              "Monodictys_castaneae",
                              "Cystobasidium_slooffiae",
                              "Meruliaceae_ASV245",
                              "Alternaria_oregonensis",
                              "Alternaria_chlamydosporigena",
                              "Aspergillus_intermedius",
                              "Stemphylium_ASV127",
                              "Tritirachium_oryzae",
                              "Stemphylium_vesicarium",
                              "Trichothecium_roseum",
                              "Sporobolomyces_ASV450",
                              "Scytalidium_circinatum",
                              "Erythrobasidium_hasegawianum",
                              "Penicillium_astrolabium",
                              "Leptospora_ASV165",
                              "Peniophora_ASV546",
                              "Arthrinium_marii",
                              "Nigrospora_oryzae",
                              "Hyphodermella_rosae",
                              "Zopfiella_lundqvistii",
                              "Torula_ASV213",
                              "Knufia_ASV505",
                              "Chaetomium_madrasense",
                              "Phaeosphaeria_juncicola",
                              "Issatchenkia_orientalis",
                              "Pseudallescheria_boydii",
                              "Vishniacozyma_heimaeyensis",
                              "Stagonospora_ASV589",
                              "Acremonium_sclerotigenum",
                              "Cladosporium_dominicanum",
                              "Buckleyzyma_aurantiaca",
                              "Vishniacozyma_carnescens",
                              "Aspergillus_versicolor",
                              "Chaetomium_angustispirale",
                              "Buckleyzyma_phyllomatis"), 
             EmbryostemMSfungicides = c("Alternaria_infectoria",
                                        "Cladosporium_sphaerospermum",
                                        "Alternaria_alternata",
                                        "Cladosporium_cladosporioides",
                                        "Mycosphaerella_tassiana",
                                        "Sporobolomyces_roseus",
                                        "Botrytis_ASV45",
                                        "Aureobasidium_pullulans",
                                        "Cladosporium_halotolerans",
                                        "Schizophyllum_commune",
                                        "Botrytis_cinerea",
                                        "Chaetomium_angustispirale",
                                        "Didymella_urticicola",
                                        "Coprinopsis_foetidellus",
                                        "Cladosporium_exasperatum",
                                        "Candida_sake",
                                        "Blumeria_graminis",
                                        "Monodictys_castaneae",
                                        "Phaeosphaeria_juncicola",
                                        "Irpex_lacteus",
                                        "Phaeosphaeria_ASV107",
                                        "Wallemia_hederae",
                                        "Sarocladium_kiliense",
                                        "Coriolopsis_gallica",
                                        "Golovinomyces_ASV178",
                                        "Naganishia_diffluens",
                                        "Neodevriesia_ASV366",
                                        "Cryptococcus_saitoi",
                                        "Ramularia_eucalypti",
                                        "Filobasidium_chernovii",
                                        "Beauveria_bassiana",
                                        "Coprinopsis_ASV367",
                                        "Meruliaceae_ASV454",
                                        "Cladosporium_dominicanum",
                                        "Symmetrospora_coprosmae",
                                        "Acremonium_charticola",
                                        "Cladosporium_ramotenellum",
                                        "Cutaneotrichosporon_curvatus",
                                        "Wallemia_ASV88",
                                        "Parengyodontium_album",
                                        "Lophiostoma_ASV565",
                                        "Devriesia_ASV497",
                                        "Wallemia_sebi",
                                        "Saxophila_tyrrhenica",
                                        "Discosia_pseudoartocreas",
                                        "Alternaria_lolii",
                                        "Mycosphaerellaceae_ASV424",
                                        "Zopfiella_marina",
                                        "Zopfiella_lundqvistii",
                                        "Cladosporium_ASV8",
                                        "Antrodia_xantha",
                                        "Septoriella_phragmitis",
                                        "Myriostoma_coliforme",
                                        "Dioszegia_hungarica",
                                        "Wallemia_ASV276",
                                        "Stemphylium_vesicarium",
                                        "Verticillium_leptobactrum",
                                        "Penicillium_kongii"), 
             CallusstemMS = c("Mycosphaerella_tassiana",
                              "Cladosporium_cladosporioides",
                              "Filobasidium_ASV25",
                              "Cladosporium_sphaerospermum",
                              "Sporobolomyces_roseus",
                              "Coriolopsis_gallica",
                              "Cladosporium_ASV8",
                              "Cladosporium_halotolerans",
                              "Alternaria_alternata",
                              "Filobasidium_oeirense",
                              "Cladosporium_ramotenellum",
                              "Blumeria_graminis",
                              "Peniophora_albobadia",
                              "Alternaria_infectoria",
                              "Symmetrospora_foliicola",
                              "Naganishia_diffluens",
                              "Aureobasidium_pullulans",
                              "Naganishia_albida",
                              "Acremonium_charticola",
                              "Auriculariales_ASV53",
                              "Candida_tropicalis",
                              "Lepista_sordida",
                              "Cryptomarasmius_corbariensis",
                              "Coprinopsis_ASV122",
                              "Vishniacozyma_heimaeyensis",
                              "Coprinellus_ASV123",
                              "Steccherinaceae_ASV158",
                              "Botrytis_ASV45",
                              "Stemphylium_ASV6",
                              "Schizophyllum_commune",
                              "Cystobasidium_lysinophilum",
                              "Torula_ASV213",
                              "Omphalotus_olearius",
                              "Phaeosphaeria_ASV107",
                              "Irpex_lacteus",
                              "Tubaria_conspersa",
                              "Sistotrema_brinkmannii",
                              "Peniophora_nuda",
                              "Torula_ASV101",
                              "Verticillium_leptobactrum",
                              "Ceratobasidium_ASV202",
                              "Neosetophoma_ASV118",
                              "Inocutis_tamaricis",
                              "Toxicocladosporium_strelitziae",
                              "Peniophora_ASV34",
                              "Psathyrella_trinitatensis",
                              "Vishniacozyma_victoriae",
                              "Phaeosphaeria_juncicola",
                              "Coprinellus_micaceus",
                              "Lentinus_arcularius",
                              "Wallemia_ASV73",
                              "Filobasidium_ASV138",
                              "Alternaria_chlamydosporigena",
                              "Neodevriesiaceae_ASV93",
                              "Dioszegia_hungarica",
                              "Cystidiodontia_laminifera",
                              "Wallemia_hederae",
                              "Phanerochaete_ASV235",
                              "Fusarium_proliferatum",
                              "Acremonium_sclerotigenum",
                              "Botrytis_cinerea",
                              "Issatchenkia_orientalis",
                              "Leptospora_ASV165",
                              "Ramularia_eucalypti",
                              "Rhodotorula_graminis",
                              "Cladosporium_exasperatum",
                              "Toxicocladosporium_irritans",
                              "Stachybotrys_chartarum",
                              "Dekkera_anomala",
                              "Kodamaea_ohmeri",
                              "Aspergillus_intermedius",
                              "Sistotrema_oblongisporum",
                              "Cystobasidium_pinicola",
                              "Xylariales_ASV314",
                              "Leucosporidium_yakuticum",
                              "Arthrinium_marii",
                              "Cladosporium_dominicanum",
                              "Nectria_pseudotrichia",
                              "Candida_sake",
                              "Debaryomyces_hansenii",
                              "Hortaea_werneckii",
                              "Wallemia_ASV88",
                              "Symmetrospora_coprosmae",
                              "Alternaria_terricola",
                              "Vishniacozyma_carnescens",
                              "Penicillium_kongii",
                              "Hyphodermella_rosae",
                              "Penicillium_astrolabium"), CallusstemMSfungicides = c("Alternaria_alternata",
                                                                                     "Mycosphaerella_tassiana",
                                                                                     "Sporobolomyces_roseus",
                                                                                     "Cladosporium_ramotenellum",
                                                                                     "Cladosporium_sphaerospermum",
                                                                                     "Cladosporium_cladosporioides",
                                                                                     "Wallemia_hederae",
                                                                                     "Alternaria_infectoria",
                                                                                     "Coriolopsis_gallica",
                                                                                     "Cladosporium_ASV8",
                                                                                     "Gibberella_intricans",
                                                                                     "Stemphylium_ASV6",
                                                                                     "Phaeosphaeria_juncicola",
                                                                                     "Naganishia_diffluens",
                                                                                     "Preussia_flanaganii",
                                                                                     "Clonostachys_miodochialis",
                                                                                     "Knufia_perforans",
                                                                                     "Sistotremastrum_guttuliferum",
                                                                                     "Coprinellus_micaceus",
                                                                                     "Botrytis_cinerea",
                                                                                     "Blumeria_graminis",
                                                                                     "Preussia_persica",
                                                                                     "Cladosporium_halotolerans",
                                                                                     "Vishniacozyma_carnescens",
                                                                                     "Parastagonospora_ASV4",
                                                                                     "Buckleyzyma_phyllomatis",
                                                                                     "Phialocephala_fluminis",
                                                                                     "Aureobasidium_pullulans",
                                                                                     "Acrodontium_crateriforme",
                                                                                     "Epicoccum_nigrum",
                                                                                     "Symmetrospora_coprosmae",
                                                                                     "Crocicreas_ASV92",
                                                                                     "Trichosporonaceae_ASV550",
                                                                                     "Plectosphaerellaceae_ASV125",
                                                                                     "Cladosporium_dominicanum",
                                                                                     "Coprinellus_radians",
                                                                                     "Wallemia_ASV73",
                                                                                     "Auriculariales_ASV53",
                                                                                     "Holtermanniella_wattica",
                                                                                     "Schizophyllum_commune",
                                                                                     "Buckleyzyma_aurantiaca",
                                                                                     "Amyloporia_sinuosa",
                                                                                     "Rhodotorula_mucilaginosa",
                                                                                     "Candida_heveicola",
                                                                                     "Hansfordia_ASV607",
                                                                                     "Arthrinium_marii",
                                                                                     "Fusariella_ASV200",
                                                                                     "Cladosporium_exasperatum",
                                                                                     "Zygosporium_ASV128",
                                                                                     "Neoerysiphe_nevoi",
                                                                                     "Candida_railenensis",
                                                                                     "Capnobotryella_ASV368",
                                                                                     "Sarocladium_implicatum",
                                                                                     "Ramularia_eucalypti",
                                                                                     "Aspergillus_intermedius",
                                                                                     "Stemphylium_vesicarium",
                                                                                     "Auricularia_auricula-judae",
                                                                                     "Porodisculus_pendulus",
                                                                                     "Candida_sake",
                                                                                     "Saxophila_tyrrhenica",
                                                                                     "Fonsecaea_ASV389",
                                                                                     "Peniophora_ASV407",
                                                                                     "Neodevriesia_ASV321",
                                                                                     "Alternaria_terricola",
                                                                                     "Golovinomyces_montagnei",
                                                                                     "Aspergillus_versicolor",
                                                                                     "Trichaptum_biforme",
                                                                                     "Toxicocladosporium_irritans",
                                                                                     "Phaeosphaeria_ASV107",
                                                                                     "Penicillium_kongii",
                                                                                     "Nectria_pseudotrichia",
                                                                                     "Irpex_lacteus"))

names(stems) <- c("Control MS", "Seed fungicides MS","Seed fungicides MS fungicides", "Embryo MS","Embryo MS fungicides", "Callus MS", "Callus MS fungicides")
list_to_matrix(stems)
All_matrix = make_comb_mat(stems)
All_matrix     
UpSet(All_matrix, pt_size = unit(5, "mm"), lwd = 3, 
      top_annotation = upset_top_annotation(All_matrix, add_numbers = TRUE),
      right_annotation = upset_right_annotation(All_matrix, add_numbers = TRUE), comb_order = order(comb_size(All_matrix)))

originVENN = list(Embryo, Controlseed, Seedfungicides, Callus)
names(originVENN) <- c( "Embryo", "Control seed","Seed fungicides", "Callus")
ggvenn(
  originVENN, 
  fill_color = c("#6B990F", "#666666", "#990F0F", "#0F6B99"),
  stroke_size = 1, set_name_size = 6.5, text_size = 10, show_percentage = FALSE
)