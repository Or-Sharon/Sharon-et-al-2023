library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(tidyr)
library(rstatix)
library(ampvis2)

## uploading taxonomy files from DADA2 ##

samdf <- read.csv("metadataCLEANSEEDS.csv",header = TRUE)
rownames(samdf) <- samdf$SampleID
seqtab.brief.flt <- read.csv("seqtab.brief.flt.csv", header = TRUE, row.names = 1)
rownames(seqtab.brief.flt) <- samdf$SampleID
tax.final <- read.csv("taxidbp50clean.F.UNITEAa.csv", row.names = 1, header = TRUE)

## creating phyloseq object ## 

ps <- phyloseq(otu_table(seqtab.brief.flt, taxa_are_rows=FALSE), sample_data(samdf), tax_table(as.matrix(tax.final)))

ps2 <- tax_glom(ps, "Species", NArm = TRUE)

ps3 <- prune_taxa(taxa_sums(ps2) > 2, ps2) #remove singletons and doubletons
ps3 <- prune_samples(rowSums(otu_table(ps)) > 100, ps3) #remove samples with less than 100 seqs

psra <- transform_sample_counts(ps3, function(x){x / sum(x)}) # ra transforms, but just for filtering, not for subsequent analysis
ps5 <- filter_taxa(psra, function(x) sum(x) > .005, TRUE) # filtering

keep<- taxa_names(ps5)
ps.counts <- prune_taxa(keep, ps3)

## transforming to Hillenger data ##

ps.helt <- ps5
otu_table(ps.helt) <-otu_table(decostand(otu_table(ps5), method = "hellinger"), taxa_are_rows=FALSE)

## alpha diversity ##

alpha <- data.frame(
  "Observed" = estimate_richness(ps.counts, measures = "Observed"),
  "Shannon" = estimate_richness(ps.helt, measures = "Shannon"),  
  "Origin" = sample_data(ps.counts)$Origin, 
  "Treatment" = sample_data(ps.counts)$Treatment)

head(alpha)

mean_summery_alpha <- alpha %>%
  group_by(Origin, Treatment) %>%
  dplyr::summarise(mean_observed = mean(Observed),
                   mean_shannon = mean(Shannon), .groups = 'keep', 
                   sd_observed = sd(Observed), sd_shannon = sd(Shannon))

kruskal.test(Observed ~ Treatment, alpha)
pairwise.wilcox.test(alpha$Observed, alpha$Treatment, p.adjust.method="BH")
kruskal.test(Shannon ~ Treatment, alpha)
pairwise.wilcox.test(alpha$Shannon, alpha$Treatment, p.adjust.method="BH")

## beta diversity ##

ordps.bray <- ordinate(ps.helt, method = "PCoA", distance = "bray")

plot_ordination(ps.helt, ordps.bray, color = "Treatment") + coord_equal()  

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

seed.bray <- phyloseq::distance(ps.helt, method = "bray")
sampledf.bray <- data.frame(sample_data(ps.helt))
adonis(seed.bray ~ Treatment, data = sampledf.bray)
pairwise.adonis(seed.bray, sample_data(ps.helt)$Treatment)

## upset plot ##

All <- list(control = c("Alternaria_infectoria",
                        "Alternaria_metachromatica",
                        "Parastagonospora_ASV8",
                        "Alternaria_alternata",
                        "Cladosporium_sphaerospermum",
                        "Mycosphaerella_tassiana",
                        "Cladosporium_dominicanum",
                        "Coriolopsis_gallica",
                        "Phaeosphaeria_juncicola",
                        "Cladosporium_grevilleae",
                        "Wallemiales_ASV52",
                        "Cladosporium_cladosporioides",
                        "Cystidiodontia_laminifera",
                        "Aspergillus_versicolor",
                        "Sporobolomyces_roseus",
                        "Pichia_ASV79",
                        "Dioszegia_hungarica",
                        "Cladosporium_halotolerans",
                        "Agaricus_subrufescens",
                        "Cladosporium_ramotenellum",
                        "Schizophyllum_commune",
                        "Neocatenulostroma_ASV157",
                        "Amyloporia_sinuosa",
                        "Zopfiella_ASV77",
                        "Cystofilobasidium_macerans",
                        "Bjerkandera_adusta",
                        "Acremonium_sclerotigenum",
                        "Sarocladium_implicatum",
                        "Penicillium_astrolabium",
                        "Debaryomyces_hansenii",
                        "Rhodosporidiobolus_fluvialis",
                        "Irpex_lacteus",
                        "Naganishia_diffluens",
                        "Stereum_complicatum",
                        "Symmetrospora_coprosmae",
                        "Aureobasidium_pullulans",
                        "Toxicocladosporium_irritans",
                        "Lepista_sordida",
                        "Fusarium_proliferatum",
                        "Trichoderma_longibrachiatum"), seedfungicideMS = c("Acremonium_sclerotigenum",
                                                                            "Penicillium_citrinum",
                                                                            "Cladosporium_grevilleae",
                                                                            "Cystobasidium_lysinophilum",
                                                                            "Cladosporium_halotolerans",
                                                                            "Bjerkandera_adusta",
                                                                            "Alternaria_alternata",
                                                                            "Acremonium_persicinum",
                                                                            "Filobasidium_ASV131",
                                                                            "Oliveonia_pauxilla",
                                                                            "Peniophora_albobadia",
                                                                            "Magnaporthe_oryzae",
                                                                            "Lecanicillium_ASV232",
                                                                            "Moesziomyces_aphidis",
                                                                            "Geosmithia_ASV292",
                                                                            "Nalanthamala_vermoesenii",
                                                                            "Rhodotorula_diobovata",
                                                                            "Cladosporium_dominicanum",
                                                                            "Alternaria_infectoria",
                                                                            "Schizophyllum_commune",
                                                                            "Cystidiodontia_laminifera",
                                                                            "Debaryomyces_hansenii",
                                                                            "Fuscoporia_contigua",
                                                                            "Coriolopsis_gallica",
                                                                            "Aspergillus_versicolor",
                                                                            "Candida_parapsilosis",
                                                                            "Auriculariales_ASV62",
                                                                            "Meruliaceae_ASV49",
                                                                            "Cladosporium_sphaerospermum",
                                                                            "Penicillium_astrolabium",
                                                                            "Cladosporium_cladosporioides",
                                                                            "Gymnopus_dryophilus"), seedfungicideMSfungicide = c("Acremonium_sclerotigenum",
                                                                                                                                 "Penicillium_citrinum",
                                                                                                                                 "Cladosporium_cladosporioides",
                                                                                                                                 "Alternaria_alternata",
                                                                                                                                 "Cladosporium_sphaerospermum",
                                                                                                                                 "Cladosporium_halotolerans",
                                                                                                                                 "Cladosporium_grevilleae",
                                                                                                                                 "Peniophora_albobadia",
                                                                                                                                 "Coriolopsis_gallica",
                                                                                                                                 "Sporobolomyces_salicinus",
                                                                                                                                 "Orbilia_serpentina",
                                                                                                                                 "Cryptomarasmius_corbariensis",
                                                                                                                                 "Trichaptum_biforme",
                                                                                                                                 "Hortaea_werneckii",
                                                                                                                                 "Cladosporium_dominicanum",
                                                                                                                                 "Verticillium_leptobactrum",
                                                                                                                                 "Peniophora_nuda",
                                                                                                                                 "Sarocladium_implicatum",
                                                                                                                                 "Henningsomyces_puber",
                                                                                                                                 "Filobasidium_ASV60",
                                                                                                                                 "Amylosporus_ASV164",
                                                                                                                                 "Rhodosporidiobolus_fluvialis",
                                                                                                                                 "Ascomycota_ASV78",
                                                                                                                                 "Ceriporia_viridans",
                                                                                                                                 "Aureobasidium_pullulans",
                                                                                                                                 "Alternaria_infectoria",
                                                                                                                                 "Toxicocladosporium_irritans",
                                                                                                                                 "Stereum_ostrea",
                                                                                                                                 "Clitocybe_ASV182",
                                                                                                                                 "Mycoacia_uda",
                                                                                                                                 "Meruliaceae_ASV49",
                                                                                                                                 "Amyloporia_sinuosa",
                                                                                                                                 "Beauveria_bassiana",
                                                                                                                                 "Vishniacozyma_heimaeyensis",
                                                                                                                                 "Coprinopsis_ASV212",
                                                                                                                                 "Eupenidiella_venezuelensis",
                                                                                                                                 "Tubaria_conspersa",
                                                                                                                                 "Cladosporium_ASV268",
                                                                                                                                 "Ceratobasidiaceae_ASV281",
                                                                                                                                 "Gymnopus_dryophilus",
                                                                                                                                 "Gibberella_intricans",
                                                                                                                                 "Pleurotus_pulmonarius",
                                                                                                                                 "Mycosphaerella_tassiana",
                                                                                                                                 "Aspergillus_olivicola",
                                                                                                                                 "Trichaptum_abietinum",
                                                                                                                                 "Stemphylium_ASV95",
                                                                                                                                 "Fusarium_penzigii",
                                                                                                                                 "Stereum_complicatum",
                                                                                                                                 "Cladosporium_ramotenellum",
                                                                                                                                 "Stemphylium_ASV401",
                                                                                                                                 "Cyberlindnera_jadinii",
                                                                                                                                 "Trichoderma_longibrachiatum",
                                                                                                                                 "Moesziomyces_aphidis",
                                                                                                                                 "Auriculariales_ASV62",
                                                                                                                                 "Cladosporium_exasperatu"), embryoms = c("Acremonium_sclerotigenum",
                                                                                                                                                                          "Cladosporium_sphaerospermum",
                                                                                                                                                                          "Cladosporium_cladosporioides",
                                                                                                                                                                          "Schizophyllum_commune",
                                                                                                                                                                          "Mycosphaerella_tassiana",
                                                                                                                                                                          "Coriolopsis_gallica",
                                                                                                                                                                          "Alternaria_alternata",
                                                                                                                                                                          "Debaryomyces_hansenii",
                                                                                                                                                                          "Cladosporium_halotolerans",
                                                                                                                                                                          "Phaeosphaeria_juncicola",
                                                                                                                                                                          "Wallemia_ASV29",
                                                                                                                                                                          "Cladosporium_dominicanum",
                                                                                                                                                                          "Cryptomarasmius_corbariensis",
                                                                                                                                                                          "Cladosporium_exasperatum",
                                                                                                                                                                          "Cladosporium_grevilleae",
                                                                                                                                                                          "Verticillium_leptobactrum",
                                                                                                                                                                          "Vishniacozyma_heimaeyensis",
                                                                                                                                                                          "Aureobasidium_pullulans",
                                                                                                                                                                          "Dioszegia_hungarica",
                                                                                                                                                                          "Burgoa_ASV89",
                                                                                                                                                                          "Tubaria_conspersa",
                                                                                                                                                                          "Fusarium_proliferatum",
                                                                                                                                                                          "Cladosporium_ramotenellum",
                                                                                                                                                                          "Toxicocladosporium_irritans",
                                                                                                                                                                          "Penicillium_astrolabium",
                                                                                                                                                                          "Rhodotorula_mucilaginosa",
                                                                                                                                                                          "Alternaria_mimicula",
                                                                                                                                                                          "Alternaria_terricola",
                                                                                                                                                                          "Lepista_sordida",
                                                                                                                                                                          "Alternaria_infectoria",
                                                                                                                                                                          "Sarocladium_implicatum",
                                                                                                                                                                          "Rhodotorula_diobovata",
                                                                                                                                                                          "Peniophora_albobadia",
                                                                                                                                                                          "Peniophora_ASV142",
                                                                                                                                                                          "Botrytis_cinerea",
                                                                                                                                                                          "Ganoderma_adspersum",
                                                                                                                                                                          "Phanerochaete_tuberculata",
                                                                                                                                                                          "Stemphylium_ASV95",
                                                                                                                                                                          "Debaryomycetaceae_ASV242",
                                                                                                                                                                          "Hyphopichia_burtonii",
                                                                                                                                                                          "Basidiomycota_ASV192",
                                                                                                                                                                          "Filobasidium_ASV60",
                                                                                                                                                                          "Pseudopithomyces_chartarum",
                                                                                                                                                                          "Peniophora_ASV16",
                                                                                                                                                                          "Hortaea_werneckii",
                                                                                                                                                                          "Aspergillus_ruber",
                                                                                                                                                                          "Alternaria_obclavata",
                                                                                                                                                                          "Cordycipitaceae_ASV442",
                                                                                                                                                                          "Irpex_lacteus",
                                                                                                                                                                          "Bulleribasidium_ASV373",
                                                                                                                                                                          "Trimmatostroma_salinum",
                                                                                                                                                                          "Trichoderma_longibrachiatum"),
            embryomsfungicide = c("Cladosporium_sphaerospermum",
                                  "Coriolopsis_gallica",
                                  "Sporobolomyces_roseus",
                                  "Schizophyllum_commune",
                                  "Acremonium_sclerotigenum",
                                  "Cladosporium_halotolerans",
                                  "Cyberlindnera_jadinii",
                                  "Alternaria_infectoria",
                                  "Cladosporium_cladosporioides",
                                  "Naganishia_diffluens",
                                  "Candida_tropicalis",
                                  "Gibberella_intricans",
                                  "Phanerochaete_ASV84",
                                  "Gloeophyllum_abietinum",
                                  "Cladosporium_grevilleae",
                                  "Simocybe_haustellaris",
                                  "Auriculariales_ASV62",
                                  "Alternaria_alternata",
                                  "Stemphylium_ASV95",
                                  "Neodevriesiaceae_ASV219",
                                  "Sistotrema_brinkmannii",
                                  "Peniophora_ASV16",
                                  "Clitocybe_dealbata",
                                  "Byssomerulius_corium",
                                  "Peniophora_albobadia",
                                  "Ceriporia_viridans",
                                  "Cladosporium_ramotenellum",
                                  "Antrodia_xantha",
                                  "Phaeosphaeria_juncicola",
                                  "Capnodiales_ASV313",
                                  "Candida_parapsilosis",
                                  "Cladosporium_exasperatum",
                                  "Xylariaceae_ASV251",
                                  "Zygosporium_ASV253",
                                  "Dioszegia_ASV46",
                                  "Nigrospora_oryzae",
                                  "Tubaria_conspersa",
                                  "Symmetrospora_coprosmae",
                                  "Xepicula_ASV309",
                                  "Filobasidium_ASV60",
                                  "Nectriaceae_ASV371",
                                  "Naganishia_ASV193",
                                  "Agaricus_punjabensis",
                                  "Cladosporium_dominicanum",
                                  "Mycosphaerella_tassiana",
                                  "Fuscoporia_contigua",
                                  "Crustoderma_marianum",
                                  "Peniophora_ASV368",
                                  "Parengyodontium_album",
                                  "Peniophora_ASV406",
                                  "Meruliaceae_ASV49",
                                  "Thanatephorus_cucumeris",
                                  "Coprinellus_radians",
                                  "Rhizoctonia_fusispora",
                                  "Wallemiales_ASV356",
                                  "Rhodotorula_mucilaginosa",
                                  "Alternaria_metachromatica",
                                  "Penicillium_astrolabium",
                                  "Pleurotaceae_ASV358",
                                  "Basidiomycota_ASV386",
                                  "Knufia_epidermidis",
                                  "Cantharellales_ASV372",
                                  "Peniophora_ASV408",
                                  "Fomitopsis_ASV357",
                                  "Fusarium_proliferatum",
                                  "Acremonium_hyalinulum",
                                  "Ganoderma_resinaceum",
                                  "Henningsomyces_puber",
                                  "Trichoderma_longibrachiatum",
                                  "Sarocladium_implicatum"), 
            callusms = c("Cladosporium_halotolerans",
                         "Cladosporium_cladosporioides",
                         "Botrytis_cinerea",
                         "Stypella_variabilis",
                         "Acremonium_sclerotigenum",
                         "Amyloporia_sinuosa",
                         "Radulomyces_molaris",
                         "Pyrenochaetopsis_leptospora",
                         "Ceriporia_viridans",
                         "Cladosporium_sphaerospermum",
                         "Sistotrema_brinkmannii",
                         "Peniophora_ASV16",
                         "Rhodosporidiobolus_odoratus",
                         "Auricularia_auricula-judae",
                         "Sporobolomyces_roseus",
                         "Dioszegia_ASV46",
                         "Geosmithia_ASV257",
                         "Coriolopsis_gallica",
                         "Zopfiella_ASV34",
                         "Schizophyllum_commune",
                         "Aspergillus_westerdijkiae",
                         "Cladosporium_grevilleae",
                         "Typhula_uncialis",
                         "Trametes_hirsuta",
                         "Alternaria_infectoria",
                         "Peniophora_ASV100",
                         "Steccherinaceae_ASV447",
                         "Erysiphe_cruciferarum",
                         "Byssomerulius_corium",
                         "Wallemia_ASV475",
                         "Alternaria_alternata",
                         "Trichoderma_longibrachiatum",
                         "Cladosporium_exasperatum"), callusmsfungicide = c("Cladosporium_cladosporioides",
                                                                            "Cladosporium_grevilleae",
                                                                            "Alternaria_alternata",
                                                                            "Aureobasidium_pullulans",
                                                                            "Cladosporium_halotolerans",
                                                                            "Cladosporium_sphaerospermum",
                                                                            "Coprinopsis_ASV222",
                                                                            "Sporobolomyces_roseus",
                                                                            "Nigrospora_oryzae",
                                                                            "Henningsomyces_puber",
                                                                            "Wallemia_ASV29",
                                                                            "Wallemia_ASV35",
                                                                            "Peniophora_ASV16",
                                                                            "Auriculariales_ASV62",
                                                                            "Zygosaccharomyces_rouxii",
                                                                            "Coriolopsis_gallica",
                                                                            "Zopfiella_ASV34",
                                                                            "Sarocladium_implicatum",
                                                                            "Dictyosporiaceae_ASV113",
                                                                            "Peniophora_laxitexta",
                                                                            "Clitocybe_dealbata",
                                                                            "Meruliaceae_ASV49",
                                                                            "Leptospora_ASV68",
                                                                            "Schizophyllum_commune",
                                                                            "Neocamarosporium_ASV179",
                                                                            "Acremonium_sclerotigenum",
                                                                            "Cladosporium_dominicanum",
                                                                            "Gymnopus_dryophilus",
                                                                            "Mycosphaerella_tassiana",
                                                                            "Phanerochaete_ASV84",
                                                                            "Verticillium_leptobactrum",
                                                                            "Pluteus_cinereofuscus",
                                                                            "Dothideaceae_ASV146",
                                                                            "Auriculariales_ASV274",
                                                                            "Candida_parapsilosis",
                                                                            "Trichoderma_lixii",
                                                                            "Basidiobolomycota_ASV133",
                                                                            "Uncobasidium_ASV203",
                                                                            "Wallemia_ASV204",
                                                                            "Peniophora_albobadia",
                                                                            "Phanerochaete_ASV240",
                                                                            "Hemimycena_ochrogaleata",
                                                                            "Cyberlindnera_jadinii",
                                                                            "Wallemia_ASV374",
                                                                            "Penicillium_citrinum",
                                                                            "Aspergillus_versicolor",
                                                                            "Coprinopsis_episcopalis",
                                                                            "Trichothecium_roseum",
                                                                            "Trichoderma_longibrachiatum",
                                                                            "Moesziomyces_aphidis",
                                                                            "Alternaria_infectoria"))
names(All) <- c("Control seed"," F1 seed fungicide MS","F1 seed fungicide MS fungicides","F1 embryo MS", "F1 embryo MS fungicides", "F1 callus MS", "F1 callus MS fungicides")
list_to_matrix(All)
All_matrix = make_comb_mat(All)
All_matrix     
UpSet(All_matrix, pt_size = unit(5, "mm"),lwd = 3,
      top_annotation = upset_top_annotation(All_matrix, add_numbers = TRUE),
      right_annotation = upset_right_annotation(All_matrix, add_numbers = TRUE), comb_order = order(comb_size(All_matrix)))


## heatmap realtive abundace ##

#require the devtools package to source gist
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

#go converty

ampvis2_obj <- phyloseq_to_ampvis2(ps.helt)

amp_heatmap(
  data = ampvis2_obj,
  tax_aggregate = "Species", measure = "mean", 
  group_by = "Treatment",  tax_show = 10, round = 3, min_abundance = 0.001, 
  plot_values = TRUE, normalise = FALSE, plot_colorscale = "sqrt") 


## differential abundace ##

# Loop t-test for all taxa
results <- data.frame()
for (i in 1:length(taxa_names(ps.helt))){
  
  asv <- taxa_names(ps.helt)[i]
  
  Control_seed_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "Control_seed"))[,i])
  F1_embryo_MS_fungicides_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "F1_embryo_MS_fungicides"))[,i])
  F1_embryo_MS_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "F1_embryo_MS"))[,i])
  F1_seed_fungicide_MS_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "F1_seed_fungicide_MS"))[,i])
  F1_seed_fungicide_MS_fungicides_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "F1_seed_fungicide_MS_fungicides"))[,i])
  F1_callus_MS_fungicides_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "F1_callus_MS_fungicides"))[,i])
  F1_callus_MS_helt <- mean(otu_table(subset_samples(ps.helt, Treatment == "F1_callus_MS"))[,i])
  
  pval <- kruskal.test(as.vector(otu_table(ps.helt)[,i]) ~ sample_data(ps.helt)$Treatment)$p.value
  results <- rbind(results, data.frame(OTU = asv, 
                                       Control_seed_Mean = Control_seed_helt, 
                                       F1_embryo_MS_fungicides_Mean = F1_embryo_MS_fungicides_helt,
                                       F1_embryo_MS_Mean = F1_embryo_MS_helt,
                                       F1_seed_fungicide_MS_Mean = F1_seed_fungicide_MS_helt,
                                       F1_seed_fungicide_MS_fungicides_Mean = F1_seed_fungicide_MS_fungicides_helt,
                                       F1_callus_MS_fungicides_Mean = F1_callus_MS_fungicides_helt,
                                       F1_callus_MS_Mean = F1_callus_MS_helt,
                                       P = pval))
}


# Add taxa designations to table
results <- inner_join(results, rownames_to_column(data.frame(tax_table(ps.helt)), "OTU"), by = "OTU")
results

# Computing FDR corrected p-values and filtering ## 
kruskal_results <- results %>%
  arrange(P) %>%
  mutate(BH_FDR = p.adjust(P, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, P, BH_FDR, everything())
write.csv(kruskal_results, "kruskalF1seeds.csv")

kruskal_results_noBH <- results %>%
  arrange(P) %>%
  mutate(BH_FDR = p.adjust(P, "BH")) %>%
  dplyr::select(OTU, P, BH_FDR, everything())
write.csv(kruskal_results_noBH, "kruskalF1seedsoBHfilter.csv")

## evenness analysis ## 

data_evenness <- data.frame(
  evenness(ps.helt, index = "pielou", zeroes = TRUE, detection = 0), "Treatment" = sample_data(ps.helt)$Treatment)
head(data_evenness)
write.csv(data_evenness, "Evenness_second_run.csv")

data_evenness %>%
  group_by(Treatment) %>%
  dplyr::summarise(mean_pielou = mean(pielou), 
                   sd_pielou = sd(pielou))
kruskal.test(pielou ~ Treatment, data_evenness)
pairwise.wilcox.test(data_evenness$pielou, data_evenness$Treatment, p.adjust.method="BH")
