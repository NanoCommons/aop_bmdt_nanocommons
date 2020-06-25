###########################################################################
## AOP-dependent and AOP-independent NM grouping
## Based on the workflows described in
## 1. Labib, S. et al. Nano-risk Science: application of toxicogenomics in 
## an adverse outcome pathway framework for risk assessment of multi-walled
## carbon nanotubes. Particle and Fibre Toxicology 13, 15 (2016).
## 2. Halappanavar, S. et al. Ranking of nanomaterial potency to induce 
## pathway perturbations associated with lung responses. NanoImpact 14, 
## 100158 (2019).
## 3. Rahman, L., Wu, D., Johnston, M., Williams, A. & Halappanavar, S. 
## Toxicogenomics analysis of mouse lung responses following exposure to 
## titanium dioxide nanomaterials reveal their disease potential at high 
## doses. Mutagenesis 32, 59-76 (2017).
## Author: Irene Liampa (@irliampa, irini.liampa@gmail.com)
## Code created : September 2019
## Last Revision : June 2020
##
###########################################################################

# Import libraries --------------------------------------------------------

library(GEOquery)
library(readxl)

library(maanova)
library(AnnotationDbi)
library(mgug4122a.db)

library(enrichR)
library(CTDquerier)

library(tidyverse)
library(stringr)
library(reshape2)

# Download data -----------------------------------------------------------

gseRahman16 <- getGEO("GSE81570")
#save(gseRahman16, file = "gseRahman16.RData") 

exprRahman16 <- exprs(gseRahman16[[1]])
phenoRahman16 <- pData(gseRahman16[[1]])
sampleNames(gseRahman16[[1]])
#experimentData(gseRahman16[[1]])

phenoRahman16$Array <- sapply(phenoRahman16$title, str_match, "\\[[^\\]\\[]*]")
phenoRahman16$Array <- str_replace_all(phenoRahman16$Array, "\\[|\\]", "")
phenoRahman16$Array <- as.factor(phenoRahman16$Array)

# to download separate datasets from the SuperSeries
# run separateSetsDownload.R

# Basic QC ----------------------------------------------------------------

# take 40 random samples to check dataset status
pdf("random_boxplots.pdf")
boxplot(exprRahman16[, sample(1:ncol(exprRahman16), size = 40)], pch = 16, cex = .9)
title(main = "Rahman et al, 2016 GSE81570 Dataset")
dev.off()

# Data are normalized!!!

# Number of samples in the complete datasets
nRahman16 <- ncol(exprRahman16)
nRahman16 #233

# Clean phenotype data
phenoRahman16 <- phenoRahman16[, c(1, 2, 8, 13:15, 17, 24, 27, 36, 48:56)]
write_csv(phenoRahman16, path = "phenoRahman16.csv")

# Annotation object
gpl <- featureData(gseRahman16[[1]])

# Reform dataset 
phenoRahman16$Dose <- phenoRahman16$`dosage:ch1`
phenoRahman16$Time <- phenoRahman16$`time point:ch1`

# remove the "control" sample, since missing info
phenoRahman16 <- phenoRahman16[!(phenoRahman16$`dosage:ch1` == "control"), ]
exprRahman16 <- exprRahman16[, rownames(phenoRahman16)]

phenoRahman16$Dose <- ifelse(phenoRahman16$Dose == "0ug", "negative", 
                             (ifelse(phenoRahman16$Dose == "54ug", "middle", 
                              ifelse(phenoRahman16$Dose == "162ug" , "high", "veryHigh"))))
phenoRahman16$Dose <- factor(phenoRahman16$Dose, levels = c("negative", "middle", "high", 
                                                            "veryHigh"))
phenoRahman16$nDose <- ifelse(phenoRahman16$Dose == "negative", 0, 
                              (ifelse(phenoRahman16$Dose == "middle", 54, 
                               ifelse(phenoRahman16$Dose == "high" , 162, 486))))
phenoRahman16$fDose <- factor(phenoRahman16$nDose)

phenoRahman16$t <- ifelse(phenoRahman16$Time == "1d post-exposure", "1", "28")
phenoRahman16$t <- as.numeric(as.character(phenoRahman16$t))
phenoRahman16$Time <- factor(phenoRahman16$Time)

phenoRahman16$NM <- phenoRahman16$`exposed to:ch1`
phenoRahman16$NM <- ifelse(phenoRahman16$NM == "anatase TiO2NPs of 8 nm diameter", "anatase8", 
                   ifelse(phenoRahman16$NM == "anatase TiO2NPs of 20nm diameter", "anatase20", 
                          ifelse(phenoRahman16$NM == "anatase TiO2NPs of 300 nm diameter", 
                                 "anatase300", ifelse(phenoRahman16$NM == "mix rutile/anatase 
                                        TiO2NPs of 20 nm diameter" , "mix20", ifelse(phenoRahman16$NM
                                                == "rutile hydrophilic TiO2NPs of 20 nm diameter",
                                               "rutileHydrophilic20", "rutileHydrophobic20")))))
phenoRahman16$NM <- factor(phenoRahman16$NM, levels = c( "anatase8", "anatase20", "anatase300", 
                               "mix20", "rutileHydrophilic20", "rutileHydrophobic20"))
phenoRahman16$NMc <- factor(phenoRahman16$NM, levels = c("control", levels(phenoRahman16$NM)))
phenoRahman16[phenoRahman16$Dose == "negative", ncol(phenoRahman16)] <- "control"

phenoRahman16 <- phenoRahman16[, (ncol(phenoRahman16)-7):ncol(phenoRahman16)]

rm(gseRahman16, nRahman16)

save.image("rah16Reformed.RData")


# maanova -----------------------------------------------------------------

rm(gseRahman16)

pheno <- phenoRahman16
pheno$Slide <- c(1:nrow(pheno))
pheno$nArray <- ifelse(pheno$Array == "set 1", "1", ifelse(pheno$Array == "set 2", "2", 
                                                           ifelse(pheno$Array == "set 3", "3", ifelse(pheno$Array == "set 4", "4", 
                                                           ifelse(pheno$Array == "set 5", "5", "6")))))
pheno$nArray <- as.numeric(pheno$nArray)
pheno$Sample <- ifelse(pheno$NMc == "control", 0, ifelse(pheno$NMc == "anatase8", 1, 
                                                         ifelse(pheno$NMc == "anatase20", 2, ifelse(pheno$NMc == "anatase300", 3,
                                                         ifelse(pheno$NMc == "mix20", 4, ifelse(pheno$NMc == "rutileHydrophilic20", 5, 6))))))

dataRahman <- read.madata(datafile = exprRahman16, designfile = pheno, header = T)

# Differential Expression Analysis

# Model descriptions/rationale
# fixed effects 
modnFull <- fitmaanova(dataRahman, ~ nArray + nDose + t + NMc, random = ~ nArray, 
                       covariate = ~ nDose + t, verbose = TRUE)

modnFull2 <- fitmaanova(dataRahman, ~ nArray + NMc, random = ~ 1, 
                       covariate = ~  nArray, verbose = TRUE)  # this is what Halappanavar et al., 2019 report

# Results
testnFull <- matest(dataRahman, modnFull, term = "NMc",
                    shuffle.method = "resid", n.perm = 100, verbose = TRUE)
testnFullAdj <- adjPval(testnFull , method = "jsFDR")
resultnFullAdj <- summarytable(testnFullAdj)
resultnFullAdj
summarytable(testnFullAdj, outfile='maanovaRahman.csv') 
significantSet = summarytable(testnFullAdj, method ='Pvalperm', 
                              test=c('F1','Fs'), whichTest='Fs.Pvalperm', 
                              threshold = 0.05, outfile='maanovaRahmanSignificant.csv')

testnFull2 <- matest(dataRahman, modnFull2, term = "NMc",
                    shuffle.method = "resid", n.perm = 100, verbose = TRUE)
testnFullAdj2 <- adjPval(testnFull2, method = "jsFDR")
resultnFullAdj2 <- summarytable(testnFullAdj2)
resultnFullAdj2
summarytable(testnFullAdj2, outfile='maanovaRahman.csv') 
significantSet2 = summarytable(testnFullAdj2, method ='Pvalperm', 
                              test=c('F1','Fs'), whichTest='Fs.Pvalperm', 
                              threshold = 0.05, outfile='maanovaRahmanSignificant.csv') 
save.image('maanovaRahman.RData')

# GSEA --------------------------------------------------------------------

significantResults <- resultnFullAdj[rownames(significantSet), ]

# Filter by foldchange
significantResults <- significantResults[which(significantResults[, 1] >= 0.5), ]
significantResults <- cbind(significantResults, rownames(significantResults))
colnames(significantResults)[6] <- 'PROBEID'

# Annotate
keytypes(mgug4122a.db)
keys <- keytypes(mgug4122a.db)[c(3, 6, 10, 20, 23)]

keysAnova <- significantResults[, 6]
annotated <- AnnotationDbi::select(mgug4122a.db, keys = keysAnova,
                                               columns = keys, keytype = 'PROBEID')
symbolsAnova <- unique(annotated$SYMBOL)
symbolsAnova <- toupper(na.omit(symbolsAnova))
symbols <- as.vector(symbolsAnova)

# Merge information
significantResults <- merge(significantResults, annotated,
                                  by.x=0, by.y="PROBEID")

# Write results
write.csv(significantResults, file = 'significantResultsAnova.csv')
write.csv(symbolsAnova, file = 'symbolsAnova.csv')


# GSEA --------------------------------------------------------------------

dbs <- listEnrichrDbs()
dbs.go <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
dbs.pathways <- c("KEGG_2019_Mouse", "KEGG_2019_Human", "WikiPathways_2019_Mouse", "WikiPathways_2019_Human")
dbs.diseases <- c("OMIM_Disease", "DisGeNET")
dbs.phenotypes <- "MGI_Mammalian_Phenotype_Level_4_2019"
dbs.tissues <- c("Tissue_Protein_Expression_from_ProteomicsDB", "Tissue_Protein_Expression_from_Human_Proteome_Map", 
                 "Jensen_TISSUES", "ARCHS4_Tissues")

go <- enrichr(symbols, dbs.go)
pathways <- enrichr(symbols, dbs.pathways)
diseases <- enrichr(symbols, dbs.diseases)
phenotypes <- enrichr(symbols, dbs.phenotypes)
tissues <- enrichr(symbols, dbs.tissues)


# Filter by gene count ----------------------------------------------------

go$GO_Biological_Process_2018$GeneCount <- sapply(strsplit(as.character(go$GO_Biological_Process_2018$Overlap),'/'), "[", 1)
bp5 <- go$GO_Biological_Process_2018[go$GO_Biological_Process_2018$GeneCount > 4, ]

pathways$KEGG_2019_Mouse$GeneCount <- sapply(strsplit(as.character(pathways$KEGG_2019_Mouse$Overlap),'/'), "[", 1)
keggMm5 <- pathways$KEGG_2019_Mouse[pathways$KEGG_2019_Mouse$GeneCount > 4, ]

pathways$WikiPathways_2019_Mouse$GeneCount <- sapply(strsplit(as.character(pathways$WikiPathways_2019_Mouse$Overlap),'/'), "[", 1)
wikiMm5 <- pathways$WikiPathways_2019_Mouse[pathways$WikiPathways_2019_Mouse$GeneCount > 4, ]

save.image('RahmanAnovaGSEA.RData')

# Center Intensities ------------------------------------------------------

# get log2 intensities (ratios) for the DEGs
signIntensities <- exprRahman16[significantResults$Row.names, ]

# group average using MEAN INTENSITY

# 1d
controlSamples1d <- rownames(pheno[which(pheno$NMc == "control" & 
                                                   pheno$Time == "1d post-exposure"), ])

caseSamples1d <- rownames(pheno[which(pheno$NMc != "control" & 
                                        pheno$Time == "1d post-exposure"), ])
  
# 28d
controlSamples28d <- rownames(pheno[which(pheno$NMc == "control" & 
                                           pheno$Time == "28d post-exposure"), ])

caseSamples28d <- rownames(pheno[which(pheno$NMc != "control" & 
                                        pheno$Time == "28d post-exposure"), ])  

controlIntensities1d <- signIntensities[, controlSamples1d]
caseIntensities1d <- signIntensities[, caseSamples1d]

controlIntensities28d <- signIntensities[, controlSamples28d]
caseIntensities28d <- signIntensities[, caseSamples28d]

controlMeans1d <- apply(controlIntensities1d, MARGIN = 1, mean)
controlMeans28d <- apply(controlIntensities28d, MARGIN = 1, mean)

# Center log2 transformed intensities using median of controls (mean subtraction)
controlCent1d <- (controlIntensities1d - controlMeans1d)
caseCent1d <- (caseIntensities1d - controlMeans1d)

controlCent28d <- (controlIntensities28d - controlMeans28d)
caseCent28d <- (caseIntensities28d - controlMeans28d)

# Samples
Samples1d <- rownames(phenoRahman16[which(phenoRahman16$Time == "1d post-exposure"), ])
controlSamples1d <- rownames(phenoRahman16[which(phenoRahman16$NMc == "control" & 
                                                    phenoRahman16$Time == "1d post-exposure"), ])
anatase8Samples1d <- rownames(phenoRahman16[which(phenoRahman16$NM =="anatase8" & 
                                                     phenoRahman16$Time == "1d post-exposure"), ])
anatase20Samples1d <- rownames(phenoRahman16[which(phenoRahman16$NM =="anatase20" & 
                                                      phenoRahman16$Time == "1d post-exposure"), ])
anatase300Samples1d <- rownames(phenoRahman16[which(phenoRahman16$NM =="anatase300" & 
                                                       phenoRahman16$Time == "1d post-exposure"), ])
mix20Samples1d <- rownames(phenoRahman16[which(phenoRahman16$NM =="mix20" & 
                                                  phenoRahman16$Time == "1d post-exposure"), ])
rutileHydrophilic20Samples1d <- rownames(phenoRahman16[which(phenoRahman16$NM 
                                                              =="rutileHydrophilic20" & 
                                                                phenoRahman16$Time == "1d post-exposure"), ])
rutileHydrophobic20Samples1d <- rownames(phenoRahman16[which(phenoRahman16$NM 
                                                              =="rutileHydrophobic20" & 
                                                                phenoRahman16$Time == "1d post-exposure"), ])
Samples28d <- rownames(phenoRahman16[which(phenoRahman16$Time == "28d post-exposure"), ])
controlSamples28d <- rownames(phenoRahman16[which(phenoRahman16$NMc == "control" & 
                                                    phenoRahman16$Time == "28d post-exposure"), ])
anatase8Samples28d <- rownames(phenoRahman16[which(phenoRahman16$NM =="anatase8" & 
                                                     phenoRahman16$Time == "28d post-exposure"), ])
anatase20Samples28d <- rownames(phenoRahman16[which(phenoRahman16$NM =="anatase20" & 
                                                      phenoRahman16$Time == "28d post-exposure"), ])
anatase300Samples28d <- rownames(phenoRahman16[which(phenoRahman16$NM =="anatase300" & 
                                                       phenoRahman16$Time == "28d post-exposure"), ])
mix20Samples28d <- rownames(phenoRahman16[which(phenoRahman16$NM =="mix20" & 
                                                  phenoRahman16$Time == "28d post-exposure"), ])
rutileHydrophilic20Samples28d <- rownames(phenoRahman16[which(phenoRahman16$NM 
                                                              =="rutileHydrophilic20" & 
                                                                phenoRahman16$Time == "28d post-exposure"), ])
rutileHydrophobic20Samples28d <- rownames(phenoRahman16[which(phenoRahman16$NM 
                                                              =="rutileHydrophobic20" & 
                                                                phenoRahman16$Time == "28d post-exposure"), ])
rm(annotated, dbs, dbs.diseases, dbs.go, dbs.pathways, dbs.phenotypes, dbs.tissues)

# add dose info
intensities1d <- cbind(controlCent1d, caseCent1d)
intensities28d <- cbind(controlCent28d, caseCent28d)

doseNM1d <- phenoRahman16[Samples1d, c("nDose", "NM")]
intensities1d <- intensities1d[, rownames(doseNM1d)]
DRdata1d <- cbind(doseNM1d, t(intensities1d))

doseNM28d <- phenoRahman16[Samples28d, c("nDose", "NM")]
intensities28d <- intensities28d[, rownames(doseNM28d)]
DRdata28d <- cbind(doseNM28d, t(intensities28d))

# stratify data 
DRanatase81d <- DRdata1d[which(DRdata1d$NM == "anatase8"), -2]
DRanatase201d <- DRdata1d[which(DRdata1d$NM == "anatase20"), -2]
DRanatase3001d <- DRdata1d[which(DRdata1d$NM == "anatase300"), -2]
DRmix201d <- DRdata1d[which(DRdata1d$NM == "mix20"), -2]
DRrutileHydrophilic201d <- DRdata1d[which(DRdata1d$NM == "rutileHydrophilic20"), -2]
DRrutileHydrophobic201d <- DRdata1d[which(DRdata1d$NM == "rutileHydrophobic20"), -2]

DRanatase828d <- DRdata28d[which(DRdata28d$NM == "anatase8"), -2]
DRanatase2028d <- DRdata28d[which(DRdata28d$NM == "anatase20"), -2]
DRanatase30028d <- DRdata28d[which(DRdata28d$NM == "anatase300"), -2]
DRmix2028d <- DRdata28d[which(DRdata28d$NM == "mix20"), -2]
DRrutileHydrophilic2028d <- DRdata28d[which(DRdata28d$NM == "rutileHydrophilic20"), -2]
DRrutileHydrophobic2028d <- DRdata28d[which(DRdata28d$NM == "rutileHydrophobic20"), -2]

write.table(t(DRanatase201d), file = "DRanatase201dT.txt", sep = '\t')
write.table(t(DRanatase2028d), file = "DRanatase2028dT.txt", sep = '\t')

write.table(t(DRanatase81d), file = "DRanatase81dT.txt", sep = '\t')
write.table(t(DRanatase828d), file = "DRanatase828dT.txt", sep = '\t')

write.table(t(DRanatase3001d), file = "DRanatase3001dT.txt", sep = '\t')
write.table(t(DRanatase30028d), file = "DRanatase30028dT.txt", sep = '\t')

write.table(t(DRmix201d), file = "DRmix201dT.txt", sep = '\t')
write.table(t(DRmix2028d), file = "DRmix2028dT.txt", sep = '\t')

write.table(t(DRrutileHydrophilic201d), file = "DRrutileHydrophilic20201dT.txt", sep = '\t')
write.table(t(DRrutileHydrophilic2028d), file = "DRrutileHydrophilic2028dT.txt", sep = '\t')

write.table(t(DRrutileHydrophobic201d), file = "DRrutileHydrophobic201dT.txt", sep = '\t')
write.table(t(DRrutileHydrophobic2028d), file = "DRrutileHydrophobic2028dT.txt", sep = '\t')

save.image("rah16Centered.RData")

###########################################################################
## This is a part to handle the output files for the best models from the 
## BMD analysis of the stratified centered intensities using the 
## BMDExpress2 program for the exposure to 6 TiO2 nanoparticles at day 1 
## p.e. (post exposure) and day 28.
## It implements the Part III of the AOP/BMDt workflow, as it is described 
## in Deliverable Report 6.3 
###########################################################################


# Input -------------------------------------------------------------------

# 1d

# anatase 20nm

# best models only
bestAnatase201d <- read.delim("BestModels_DRanatase201dT_BMD.txt", header = FALSE, skip = 40,
                              row.names = 1, fill = TRUE)
bestAnatase201d <- bestAnatase201d[, -c(2:4)]
colnames(bestAnatase201d) <- bestAnatase201d[1, ]
bestAnatase201d <- bestAnatase201d[-1, ]

# extract BMDL
bestBMDAnatase201d <- bestAnatase201d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDAnatase201d)
bestBMDAnatase201d <- apply(bestBMDAnatase201d, 2, as.numeric)
rownames(bestBMDAnatase201d) <- elements

anatase201dBMDlower <- apply(bestBMDAnatase201d, 1, min)
anatase201dBMDlower <- as.matrix(anatase201dBMDlower)
colnames(anatase201dBMDlower) <- "BMDlower"

# anatase 8nm
# best models only
bestAnatase81d <- read.delim("BestModels_DRanatase81dT_BMD.txt", header = FALSE, skip = 36,
                              row.names = 1, fill = TRUE)
bestAnatase81d <- bestAnatase81d[, -c(2:4)]
colnames(bestAnatase81d) <- bestAnatase81d[1, ]
bestAnatase81d <- bestAnatase81d[-1, ]

# extract BMDL
bestBMDAnatase81d <- bestAnatase81d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDAnatase81d)
bestBMDAnatase81d <- apply(bestBMDAnatase81d, 2, as.numeric)
rownames(bestBMDAnatase81d) <- elements

anatase81dBMDlower <- apply(bestBMDAnatase81d, 1, min)
anatase81dBMDlower <- as.matrix(anatase81dBMDlower)
colnames(anatase81dBMDlower) <- "BMDlower"


# anatase 300nm
# best models only
bestAnatase3001d <- read.delim("BestModels_DRanatase3001dT_BMD.txt", header = FALSE, skip = 36,
                             row.names = 1, fill = TRUE)
bestAnatase3001d <- bestAnatase3001d[, -c(2:4)]
colnames(bestAnatase3001d) <- bestAnatase3001d[1, ]
bestAnatase3001d <- bestAnatase3001d[-1, ]

# extract BMDL
bestBMDAnatase3001d <- bestAnatase3001d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDAnatase3001d)
bestBMDAnatase3001d <- apply(bestBMDAnatase3001d, 2, as.numeric)
rownames(bestBMDAnatase3001d) <- elements

anatase3001dBMDlower <- apply(bestBMDAnatase3001d, 1, min)
anatase3001dBMDlower <- as.matrix(anatase3001dBMDlower)
colnames(anatase3001dBMDlower) <- "BMDlower"


# mix 20nm
# best models only
bestMix201d <- read.delim("BestModels_DRmix201dT_BMD.txt", header = FALSE, skip = 40,
                              row.names = 1, fill = TRUE)
bestMix201d <- bestMix201d[, -c(2:4)]
colnames(bestMix201d) <- bestMix201d[1, ]
bestMix201d <- bestMix201d[-1, ]

# extract BMDL
bestBMDMix201d <- bestMix201d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDMix201d)
bestBMDMix201d <- apply(bestBMDMix201d, 2, as.numeric)
rownames(bestBMDMix201d) <- elements

mix201dBMDlower <- apply(bestBMDMix201d, 1, min)
mix201dBMDlower <- as.matrix(mix201dBMDlower)
colnames(mix201dBMDlower) <- "BMDlower"


# rutile hydrophilic 20nm
# best models only
bestRutileHydrophilic201d <- read.delim("BestModels_DRrutileHydrophilic201dT_BMD.txt", 
                                        header = FALSE, skip = 40,
                                        row.names = 1, fill = TRUE)
colnames(bestRutileHydrophilic201d) <- bestRutileHydrophilic201d[1, ]
bestRutileHydrophilic201d <- bestRutileHydrophilic201d[-1, -c(2:4)]

# extract BMDL
bestBMDrutileHydrophilic201d <- bestRutileHydrophilic201d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDrutileHydrophilic201d)
bestBMDrutileHydrophilic201d <- apply(bestBMDrutileHydrophilic201d, 2, as.numeric)
rownames(bestBMDrutileHydrophilic201d) <- elements

rutileHydrophilic201dBMDlower <- apply(bestBMDrutileHydrophilic201d, 1, min)
rutileHydrophilic201dBMDlower <- as.matrix(rutileHydrophilic201dBMDlower)
colnames(rutileHydrophilic201dBMDlower) <- "BMDlower"


# rutile hydrophobic 20nm
# best models only
bestRutileHydrophobic201d <- read.delim("BestModels_DRrutileHydrophobic201dT_BMD.txt", 
                                        header = FALSE, skip = 36,
                                        fill = TRUE, blank.lines.skip = TRUE)
colnames(bestRutileHydrophobic201d) <- bestRutileHydrophobic201d[1, ]
bestRutileHydrophobic201d <- bestRutileHydrophobic201d[-1, -c(3:5)]
rownames(bestRutileHydrophobic201d) <- bestRutileHydrophobic201d[, 1]

# extract BMDL
bestBMDrutileHydrophobic201d <- bestRutileHydrophobic201d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDrutileHydrophobic201d)
bestBMDrutileHydrophobic201d <- apply(bestBMDrutileHydrophobic201d, 2, as.numeric)
rownames(bestBMDrutileHydrophobic201d) <- elements

rutileHydrophobic201dBMDlower <- apply(bestBMDrutileHydrophobic201d, 1, min)
rutileHydrophobic201dBMDlower <- as.matrix(rutileHydrophobic201dBMDlower)
colnames(rutileHydrophobic201dBMDlower) <- "BMDlower"


# 28d


# anatase 20nm
# best models only
bestAnatase2028d <- read.delim("BestModels_DRanatase2028dT_BMD.txt", header = FALSE, 
                              skip = 34, row.names = 1, fill = TRUE)
colnames(bestAnatase2028d) <- bestAnatase2028d[1, ]
bestAnatase2028d <- bestAnatase2028d[-1, -c(2:4)]

# extract BMDL
bestBMDAnatase2028d <- bestAnatase2028d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDAnatase2028d)
bestBMDAnatase2028d <- apply(bestBMDAnatase2028d, 2, as.numeric)
rownames(bestBMDAnatase2028d) <- elements

anatase2028dBMDlower <- apply(bestBMDAnatase2028d, 1, min)
anatase2028dBMDlower <- as.matrix(anatase2028dBMDlower)
colnames(anatase2028dBMDlower) <- "BMDlower"


# anatase 8nm
# best models only
bestAnatase828d <- read.delim("BestModels_DRanatase828dT_BMD.txt", 
                             header = FALSE, skip = 40,
                             row.names = 1, fill = TRUE)
colnames(bestAnatase828d) <- bestAnatase828d[1, ]
bestAnatase828d <- bestAnatase828d[-1, -c(2:4)]

# extract BMDL
bestBMDAnatase828d <- bestAnatase828d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDAnatase828d)
bestBMDAnatase828d <- apply(bestBMDAnatase828d, 2, as.numeric)
rownames(bestBMDAnatase828d) <- elements

anatase828dBMDlower <- apply(bestBMDAnatase828d, 1, min)
anatase828dBMDlower <- as.matrix(anatase828dBMDlower)
colnames(anatase828dBMDlower) <- "BMDlower"


# anatase 300nm
# best models only
bestAnatase30028d <- read.delim("BestModels_DRanatase30028dT_BMD.txt", 
                               header = FALSE, skip = 36,
                               row.names = 1, fill = TRUE)
colnames(bestAnatase30028d) <- bestAnatase30028d[1, ]
bestAnatase30028d <- bestAnatase30028d[-1, -c(2:4)]


# extract BMDL
bestBMDAnatase30028d <- bestAnatase30028d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDAnatase30028d)
bestBMDAnatase30028d <- apply(bestBMDAnatase30028d, 2, as.numeric)
rownames(bestBMDAnatase30028d) <- elements

anatase30028dBMDlower <- apply(bestBMDAnatase30028d, 1, min)
anatase30028dBMDlower <- as.matrix(anatase30028dBMDlower)
colnames(anatase30028dBMDlower) <- "BMDlower"


# mix 20nm
# best models only
bestMix2028d <- read.delim("BestModels_DRmix2028dT_BMD.txt", 
                          header = FALSE, skip = 35,
                          fill = TRUE)
colnames(bestMix2028d) <- bestMix2028d[1, ]
bestMix2028d <- bestMix2028d[-1, -c(3:5)]
bestMix2028d <- bestMix2028d[which(bestMix2028d$`Probe Id` != 'true' &  bestMix2028d$`Probe Id` != ""),]
rownames(bestMix2028d) <- bestMix2028d[, 1]

# extract BMDL
bestBMDMix2028d <- bestMix2028d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDMix2028d)
bestBMDMix2028d <- apply(bestBMDMix2028d, 2, as.numeric)
rownames(bestBMDMix2028d) <- elements

mix2028dBMDlower <- apply(bestBMDMix2028d, 1, min)
mix2028dBMDlower <- as.matrix(mix2028dBMDlower)
colnames(mix2028dBMDlower) <- "BMDlower"


# rutile hydrophilic 20nm
# best models only
bestRutileHydrophilic2028d <- read.delim("BestModels_DRrutileHydrophilic2028dT_BMD.txt", 
                                        header = FALSE, skip = 40,
                                        row.names = 1, fill = TRUE)
colnames(bestRutileHydrophilic2028d) <- bestRutileHydrophilic2028d[1, ]
bestRutileHydrophilic2028d <- bestRutileHydrophilic2028d[-1, -c(2:4)]

# extract BMDL
bestBMDrutileHydrophilic2028d <- bestRutileHydrophilic2028d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDrutileHydrophilic2028d)
bestBMDrutileHydrophilic2028d <- apply(bestBMDrutileHydrophilic2028d, 2, as.numeric)
rownames(bestBMDrutileHydrophilic2028d) <- elements

rutileHydrophilic2028dBMDlower <- apply(bestBMDrutileHydrophilic2028d, 1, min)
rutileHydrophilic2028dBMDlower <- as.matrix(rutileHydrophilic2028dBMDlower)
colnames(rutileHydrophilic2028dBMDlower) <- "BMDlower"


# rutile hydrophobic 20nm
# best models only
bestRutileHydrophobic2028d <- read.delim("BestModels_DRrutileHydrophobic2028dT_BMD.txt", 
                                        header = FALSE, skip = 37,
                                        fill = TRUE, blank.lines.skip = TRUE)
colnames(bestRutileHydrophobic2028d) <- bestRutileHydrophobic2028d[1, ]
bestRutileHydrophobic2028d <- bestRutileHydrophobic2028d[-1, -c(3:5)]
rownames(bestRutileHydrophobic2028d) <- bestRutileHydrophobic2028d[, 1]

# extract BMDL
bestBMDrutileHydrophobic2028d <- bestRutileHydrophobic2028d[ , c("BMDL", "BMDU")]
elements <- rownames(bestBMDrutileHydrophobic2028d)
bestBMDrutileHydrophobic2028d <- apply(bestBMDrutileHydrophobic2028d, 2, as.numeric)
rownames(bestBMDrutileHydrophobic2028d) <- elements

rutileHydrophobic2028dBMDlower <- apply(bestBMDrutileHydrophobic2028d, 1, min)
rutileHydrophobic2028dBMDlower <- as.matrix(rutileHydrophobic2028dBMDlower)
colnames(rutileHydrophobic2028dBMDlower) <- "BMDlower"

rm(list = ls(pattern = "best"))
rm(elements)

# Filter out out of range BMDL values -------------------------------------
anatase201dBMDlower <- anatase201dBMDlower[which(anatase201dBMDlower[, 1] >= 0  
                                                 & anatase201dBMDlower[, 1] <= 486), ]
anatase81dBMDlower <- anatase81dBMDlower[which(anatase81dBMDlower[, 1] >= 0  
                                                 & anatase81dBMDlower[, 1] <= 486), ]
anatase3001dBMDlower <- anatase3001dBMDlower[which(anatase3001dBMDlower[, 1] >= 0  
                                                 & anatase3001dBMDlower[, 1] <= 486), ]
mix201dBMDlower <- mix201dBMDlower[which(mix201dBMDlower[, 1] >= 0  
                                                 & mix201dBMDlower[, 1] <= 486), ]
rutileHydrophilic201dBMDlower <- rutileHydrophilic201dBMDlower[which(rutileHydrophilic201dBMDlower[, 1] >= 0  
                                                 & rutileHydrophilic201dBMDlower[, 1] <= 486), ]
rutileHydrophobic201dBMDlower <- rutileHydrophobic201dBMDlower[which(rutileHydrophobic201dBMDlower[, 1] >= 0  
                                                 & rutileHydrophobic201dBMDlower[, 1] <= 486), ]

anatase2028dBMDlower <- anatase2028dBMDlower[which(anatase2028dBMDlower[, 1] >= 0  
                                                 & anatase2028dBMDlower[, 1] <= 486), ]
anatase828dBMDlower <- anatase828dBMDlower[which(anatase828dBMDlower[, 1] >= 0  
                                               & anatase828dBMDlower[, 1] <= 486), ]
anatase30028dBMDlower <- anatase30028dBMDlower[which(anatase30028dBMDlower[, 1] >= 0  
                                                   & anatase30028dBMDlower[, 1] <= 486), ]
mix2028dBMDlower <- mix2028dBMDlower[which(mix2028dBMDlower[, 1] >= 0  
                                         & mix2028dBMDlower[, 1] <= 486), ]
rutileHydrophilic2028dBMDlower <- rutileHydrophilic2028dBMDlower[which(rutileHydrophilic2028dBMDlower[, 1] >= 0  
                                                                     & rutileHydrophilic2028dBMDlower[, 1] <= 486), ]
rutileHydrophobic2028dBMDlower <- rutileHydrophobic2028dBMDlower[which(rutileHydrophobic2028dBMDlower[, 1] >= 0  
                                                                     & rutileHydrophobic2028dBMDlower[, 1] <= 486), ]

save.image("NMlowerBMD.RData")


# Annotate BMDts ----------------------------------------------------------

# Annotate BMDL for common probes between NMs
allProbes <- unique(c(names(anatase201dBMDlower), names(anatase81dBMDlower),
                      names(anatase3001dBMDlower), names(mix201dBMDlower),
                      names(rutileHydrophilic201dBMDlower), names(rutileHydrophobic201dBMDlower),
                      names(anatase2028dBMDlower), names(anatase828dBMDlower),
                      names(anatase30028dBMDlower), names(mix2028dBMDlower),
                      names(rutileHydrophilic2028dBMDlower), names(rutileHydrophobic2028dBMDlower)))

keytypes(mgug4122a.db)
keys <- keytypes(mgug4122a.db)[c(3, 6, 10, 20, 23)]

annotate <- AnnotationDbi::select(mgug4122a.db, keys=allProbes, columns=keys, keytype="PROBEID")
symbols <- unique(annotate$SYMBOL)
write.csv(symbols , file = "symbols.csv")

# BMD GSEA ----------------------------------------------------------------

bmdPathways<- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse")

enrichedBMDpathways <-  enrichr(symbols[!(symbols == "NA")], bmdPathways)

enrichedBMDpathways$WikiPathways_2019_Mouse
fibrosisWiki <- enrichedBMDpathways$WikiPathways_2019_Mouse[5, ]

enrichedBMDkegg <- enrichedBMDpathways$KEGG_2019_Mouse
enrichedBMDkegg$Genes_count <- sapply(strsplit(as.character(enrichedBMDkegg$Overlap),'/'), 
                                      "[", 1)
enrichedBMDkegg <- enrichedBMDkegg[enrichedBMDkegg$Genes_count >4, ]

enrichedBMDkegg <- enrichedBMDkegg[order(enrichedBMDkegg$Genes_count, decreasing = T), ]
enrichedBMDkeggCounts10 <- enrichedBMDkegg[1:10, ]
enrichedBMDkegg <- enrichedBMDkegg[order(enrichedBMDkegg$Old.Adjusted.P.value, 
                                         decreasing = F), ]
enrichedBMDkeggPvalue10 <- enrichedBMDkegg[1:10, ]
enrichedBMDkegg <- enrichedBMDkegg[order(enrichedBMDkegg$Combined.Score, decreasing = F), ]
enrichedBMDkeggScore10 <- enrichedBMDkegg[1:10, ]

save.image("NMlowerBMDGsea.RData")
save.image()

# CTD ---------------------------------------------------------------------

# search for LUNG FIBROSIS term to get disease associated gene list 
ctdFibrosis <- CTDquerier::query_ctd_dise(terms = "Pulmonary Fibrosis", verbose = TRUE )
fibrosisGenes <- ctdFibrosis@gene_interactions
fibrosisChemicals <- ctdFibrosis@chemicals_interactions
fibrosisPathways <- ctdFibrosis@kegg

fibrosisKegg <- fibrosisPathways[str_detect(fibrosisPathways[ ,"Pathway.ID"], 
                                            pattern = "KEGG"), ]
fibrosisReactome <- fibrosisPathways[str_detect(fibrosisPathways[ ,"Pathway.ID"], 
                                                pattern = "REACT"), ]
save(fibrosisGenes, fibrosisChemicals, fibrosisPathways, fibrosisKegg, fibrosisReactome, 
     ctdFibrosis, file = "CTDresults.RData")

# BMD for pathways --------------------------------------------------------

keggGenes <- strsplit(enrichedBMDkeggPvalue10$Genes, ';')
names(keggGenes) <- enrichedBMDkeggPvalue10$Term
keggGenesAll <- strsplit(enrichedBMDkegg$Genes, ';')
names(keggGenesAll) <- enrichedBMDkegg$Term

# annotate BMD/BMDL values

# 1d
anatase81d <- cbind(anatase81dBMDlower)
anatase81d <- cbind(rownames(anatase81d), anatase81d)
colnames(anatase81d)[1] <- "PROBEID"
anatase81d <- merge(anatase81d, annotate, by.x=0, by.y="PROBEID")

anatase201d <- cbind(anatase201dBMDlower)
anatase201d <- cbind(rownames(anatase201d), anatase201d)
colnames(anatase201d)[1] <- "PROBEID"
anatase201d <- merge(anatase201d, annotate, by.x=0, by.y="PROBEID")

anatase3001d <- cbind(anatase3001dBMDlower)
anatase3001d <- cbind(rownames(anatase3001d), anatase3001d)
colnames(anatase3001d)[1] <- "PROBEID"
anatase3001d <- merge(anatase3001d, annotate, by.x=0, by.y="PROBEID")

mix201d <- cbind(mix201dBMDlower)
mix201d <- cbind(rownames(mix201d), mix201d)
colnames(mix201d)[1] <- "PROBEID"
mix201d <- merge(mix201d, annotate, by.x=0, by.y="PROBEID")

rutileHydrophilic201d <- cbind(rutileHydrophilic201dBMDlower)
rutileHydrophilic201d <- cbind(rownames(rutileHydrophilic201d), rutileHydrophilic201d)
colnames(rutileHydrophilic201d)[1] <- "PROBEID"
rutileHydrophilic201d <- merge(rutileHydrophilic201d, annotate, by.x=0, by.y="PROBEID")

rutileHydrophobic201d <- cbind(rutileHydrophobic201dBMDlower)
rutileHydrophobic201d <- cbind(rownames(rutileHydrophobic201d), rutileHydrophobic201d)
colnames(rutileHydrophobic201d)[1] <- "PROBEID"
rutileHydrophobic201d <- merge(rutileHydrophobic201d, annotate, by.x=0, by.y="PROBEID")

# 28d
anatase828d <- cbind(anatase828dBMDlower)
anatase828d <- cbind(rownames(anatase828d), anatase828d)
colnames(anatase828d)[1] <- "PROBEID"
anatase828d <- merge(anatase828d, annotate, by.x=0, by.y="PROBEID")

anatase2028d <- cbind(anatase2028dBMDlower)
anatase2028d <- cbind(rownames(anatase2028d), anatase2028d)
colnames(anatase2028d)[1] <- "PROBEID"
anatase2028d <- merge(anatase2028d, annotate, by.x=0, by.y="PROBEID")

anatase30028d <- cbind(anatase30028dBMDlower)
anatase30028d <- cbind(rownames(anatase30028d), anatase30028d)
colnames(anatase30028d)[1] <- "PROBEID"
anatase30028d <- merge(anatase30028d, annotate, by.x=0, by.y="PROBEID")

mix2028d <- cbind(mix2028dBMDlower)
mix2028d <- cbind(rownames(mix2028d), mix2028d)
colnames(mix2028d)[1] <- "PROBEID"
mix2028d <- merge(mix2028d, annotate, by.x=0, by.y="PROBEID")

rutileHydrophilic2028d <- cbind(rutileHydrophilic2028dBMDlower)
rutileHydrophilic2028d <- cbind(rownames(rutileHydrophilic2028d), rutileHydrophilic2028d)
colnames(rutileHydrophilic2028d)[1] <- "PROBEID"
rutileHydrophilic2028d <- merge(rutileHydrophilic2028d, annotate, by.x=0, by.y="PROBEID")

rutileHydrophobic2028d <- cbind(rutileHydrophobic2028dBMDlower)
rutileHydrophobic2028d <- cbind(rownames(rutileHydrophobic2028d), rutileHydrophobic2028d)
colnames(rutileHydrophobic2028d)[1] <- "PROBEID"
rutileHydrophobic2028d <- merge(rutileHydrophobic2028d, annotate, by.x=0, by.y="PROBEID")
save.image()

# Merge with Pathways -----------------------------------------------------

# anatase 8nm 1d ----------------------------------------------------------

keggBMDanatase81d <- list()
keggBMDanatase81dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDanatase81d[[i]] <- anatase81d[which(toupper(anatase81d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDanatase81d[[i]] <- as.numeric(keggBMDanatase81d[[i]])
  keggBMDanatase81dmedian[i] <- median(keggBMDanatase81d[[i]])
}
names(keggBMDanatase81d) <- names(keggGenes)
names(keggBMDanatase81dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDanatase81d <- list()
keggAllBMDanatase81dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDanatase81d[[i]] <- anatase81d[which(toupper(anatase81d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDanatase81d[[i]] <- as.numeric(keggAllBMDanatase81d[[i]])
  keggAllBMDanatase81dmedian[i] <- median(keggAllBMDanatase81d[[i]])
}
names(keggAllBMDanatase81d) <- names(keggGenesAll)
names(keggAllBMDanatase81dmedian) <- names(keggGenesAll)

fibrosisBMDanatase81d <- keggBMDanatase81d[names(keggBMDanatase81d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDanatase81dMedian <- keggBMDanatase81dmedian[names(keggBMDanatase81dmedian)
                                                       %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase81d <- keggAllBMDanatase81d[names(keggAllBMDanatase81d)
                                                 %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase81dMedian <- keggAllBMDanatase81dmedian[names
                                                (fibrosisAllBMDanatase81d)]

fibrosisAllBMDanatase81dMedian <- t(rbind(fibrosisAllBMDanatase81dMedian))
fibrosisAllBMDanatase81dMedian <- rbind(fibrosisAllBMDanatase81dMedian)

names(fibrosisAllBMDanatase81dMedian) <- names(fibrosisAllBMDanatase81d)
fibrosisAllBMDanatase81dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase81dMedian)
colnames(fibrosisAllBMDanatase81dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDanatase81dMedianTop10  <- fibrosisAllBMDanatase81dMedian[order(
                          fibrosisAllBMDanatase81dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDanatase81dMedianTop10 <- fibrosisAllBMDanatase81dMedianTop10[1:10]
fibrosisAllBMDanatase81dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase81dMedianTop10)
colnames(fibrosisAllBMDanatase81dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDanatase81dMedianTop10[, 1])
fibrosisAllBMDanatase81dMedianTop10 <- rbind(fibrosisAllBMDanatase81dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDanatase81dMedianTop10, file = "fibrosisAllBMDanatase81dMedianTop10.txt")


# anatase 20nm 1d ---------------------------------------------------------

keggBMDanatase201d <- list()
keggBMDanatase201dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDanatase201d[[i]] <- anatase201d[which(toupper(anatase201d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDanatase201d[[i]] <- as.numeric(keggBMDanatase201d[[i]])
  keggBMDanatase201dmedian[i] <- median(keggBMDanatase201d[[i]])
}
names(keggBMDanatase201d) <- names(keggGenes)
names(keggBMDanatase201dmedian) <- names(keggGenes)

keggAllBMDanatase201d <- list()
keggAllBMDanatase201dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDanatase201d[[i]] <- anatase201d[which(toupper(anatase201d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDanatase201d[[i]] <- as.numeric(keggAllBMDanatase201d[[i]])
  keggAllBMDanatase201dmedian[i] <- median(keggAllBMDanatase201d[[i]])
}
names(keggAllBMDanatase201d) <- names(keggGenesAll)
names(keggAllBMDanatase201dmedian) <- names(keggGenesAll)

fibrosisBMDanatase201d <- keggBMDanatase201d[names(keggBMDanatase201d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDanatase201dMedian <- keggBMDanatase201dmedian[names(keggBMDanatase201dmedian)
                                                       %in% fibrosisKegg[, 3]]
fibrosisAllBMDanatase201d <- keggAllBMDanatase201d[names(keggAllBMDanatase201d)
                                                 %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase201dMedian <- keggAllBMDanatase201dmedian[names(fibrosisAllBMDanatase201d)]

fibrosisAllBMDanatase201dMedian <- t(rbind(fibrosisAllBMDanatase201dMedian))
fibrosisAllBMDanatase201dMedian <- rbind(fibrosisAllBMDanatase201dMedian)
names(fibrosisAllBMDanatase201dMedian) <- names(fibrosisAllBMDanatase201d)
fibrosisAllBMDanatase201dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase201dMedian)
colnames(fibrosisAllBMDanatase201dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDanatase201dMedianTop10  <- fibrosisAllBMDanatase201dMedian[order(
  fibrosisAllBMDanatase201dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDanatase201dMedianTop10 <- fibrosisAllBMDanatase201dMedianTop10[1:10]
fibrosisAllBMDanatase201dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase201dMedianTop10)
colnames(fibrosisAllBMDanatase201dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDanatase201dMedianTop10[, 1])
fibrosisAllBMDanatase201dMedianTop10 <- rbind(fibrosisAllBMDanatase201dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDanatase201dMedianTop10, file = "fibrosisAllBMDanatase201dMedianTop10.txt")



# anatase 300nm 1d --------------------------------------------------------

keggBMDanatase3001d <- list()
keggBMDanatase3001dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDanatase3001d[[i]] <- anatase3001d[which(toupper(anatase3001d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDanatase3001d[[i]] <- as.numeric(keggBMDanatase3001d[[i]])
  keggBMDanatase3001dmedian[i] <- median(keggBMDanatase3001d[[i]])
}
names(keggBMDanatase3001d) <- names(keggGenes)
names(keggBMDanatase3001dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDanatase3001d <- list()
keggAllBMDanatase3001dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDanatase3001d[[i]] <- anatase3001d[which(toupper(anatase3001d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDanatase3001d[[i]] <- as.numeric(keggAllBMDanatase3001d[[i]])
  keggAllBMDanatase3001dmedian[i] <- median(keggAllBMDanatase3001d[[i]])
}
names(keggAllBMDanatase3001d) <- names(keggGenesAll)
names(keggAllBMDanatase3001dmedian) <- names(keggGenesAll)

fibrosisBMDanatase3001d <- keggBMDanatase3001d[names(keggBMDanatase3001d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDanatase3001dMedian <- keggBMDanatase3001dmedian[names(keggBMDanatase3001dmedian)
                                                       %in% fibrosisKegg[, 3]]
fibrosisAllBMDanatase3001d <- keggAllBMDanatase3001d[names(keggAllBMDanatase3001d)
                                                 %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase3001dMedian <- keggAllBMDanatase3001dmedian[names(fibrosisAllBMDanatase3001d)]

fibrosisAllBMDanatase3001dMedian <- t(rbind(fibrosisAllBMDanatase3001dMedian))
fibrosisAllBMDanatase3001dMedian <- rbind(fibrosisAllBMDanatase3001dMedian)
names(fibrosisAllBMDanatase3001dMedian) <- names(fibrosisAllBMDanatase3001d)
fibrosisAllBMDanatase3001dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase3001dMedian)
colnames(fibrosisAllBMDanatase3001dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDanatase3001dMedianTop10  <- fibrosisAllBMDanatase3001dMedian[order(
  fibrosisAllBMDanatase3001dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDanatase3001dMedianTop10 <- fibrosisAllBMDanatase3001dMedianTop10[1:10]
fibrosisAllBMDanatase3001dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase3001dMedianTop10)
colnames(fibrosisAllBMDanatase3001dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDanatase3001dMedianTop10[, 1])
fibrosisAllBMDanatase3001dMedianTop10 <- rbind(fibrosisAllBMDanatase3001dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDanatase3001dMedianTop10, file = "fibrosisAllBMDanatase3001dMedianTop10.txt")


# mix 20nm 1d -------------------------------------------------------------

keggBMDmix201d <- list()
keggBMDmix201dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDmix201d[[i]] <- mix201d[which(toupper(mix201d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDmix201d[[i]] <- as.numeric(keggBMDmix201d[[i]])
  keggBMDmix201dmedian[i] <- median(keggBMDmix201d[[i]])
}
names(keggBMDmix201d) <- names(keggGenes)
names(keggBMDmix201dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDmix201d <- list()
keggAllBMDmix201dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDmix201d[[i]] <- mix201d[which(toupper(mix201d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDmix201d[[i]] <- as.numeric(keggAllBMDmix201d[[i]])
  keggAllBMDmix201dmedian[i] <- median(keggAllBMDmix201d[[i]])
}
names(keggAllBMDmix201d) <- names(keggGenesAll)
names(keggAllBMDmix201dmedian) <- names(keggGenesAll)

fibrosisBMDmix201d <- keggBMDmix201d[names(keggBMDmix201d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDmix201dMedian <- keggBMDmix201dmedian[names(keggBMDmix201dmedian)
                                                       %in% fibrosisKegg[, 3]]
fibrosisAllBMDmix201d <- keggAllBMDmix201d[names(fibrosisAllBMDanatase3001d)]

fibrosisAllBMDmix201dMedian <- keggAllBMDmix201dmedian[names(keggBMDmix201dmedian)
                                                             %in% fibrosisKegg[, 3]]
fibrosisAllBMDmix201dMedian <- t(rbind(fibrosisAllBMDmix201dMedian))
fibrosisAllBMDmix201dMedian <- rbind(fibrosisAllBMDmix201dMedian)
names(fibrosisAllBMDmix201dMedian) <- names(fibrosisAllBMDmix201d)
fibrosisAllBMDmix201dMedianTop10 <- do.call(rbind, fibrosisAllBMDmix201dMedian)
colnames(fibrosisAllBMDmix201dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDmix201dMedianTop10  <- fibrosisAllBMDmix201dMedian[order(
  fibrosisAllBMDmix201dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDmix201dMedianTop10 <- fibrosisAllBMDmix201dMedianTop10[1:10]
fibrosisAllBMDmix201dMedianTop10 <- do.call(rbind, fibrosisAllBMDmix201dMedianTop10)
colnames(fibrosisAllBMDmix201dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDmix201dMedianTop10[, 1])
fibrosisAllBMDmix201dMedianTop10 <- rbind(fibrosisAllBMDmix201dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDmix201dMedianTop10, file = "fibrosisAllBMDmix201dMedianTop10.txt")


# rutileHydrophilic 20nm 1d -----------------------------------------------


keggBMDrutileHydrophilic201d <- list()
keggBMDrutileHydrophilic201dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDrutileHydrophilic201d[[i]] <- rutileHydrophilic201d[which(toupper(rutileHydrophilic201d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDrutileHydrophilic201d[[i]] <- as.numeric(keggBMDrutileHydrophilic201d[[i]])
  keggBMDrutileHydrophilic201dmedian[i] <- median(keggBMDrutileHydrophilic201d[[i]])
}
names(keggBMDrutileHydrophilic201d) <- names(keggGenes)
names(keggBMDrutileHydrophilic201dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDrutileHydrophilic201d <- list()
keggAllBMDrutileHydrophilic201dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDrutileHydrophilic201d[[i]] <- rutileHydrophilic201d[which(toupper(rutileHydrophilic201d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDrutileHydrophilic201d[[i]] <- as.numeric(keggAllBMDrutileHydrophilic201d[[i]])
  keggAllBMDrutileHydrophilic201dmedian[i] <- median(keggAllBMDrutileHydrophilic201d[[i]])
}
names(keggAllBMDrutileHydrophilic201d) <- names(keggGenesAll)
names(keggAllBMDrutileHydrophilic201dmedian) <- names(keggGenesAll)

fibrosisBMDrutileHydrophilic201d <- keggBMDrutileHydrophilic201d[names(keggBMDrutileHydrophilic201d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDrutileHydrophilic201dMedian <- keggBMDrutileHydrophilic201dmedian[names(keggBMDrutileHydrophilic201dmedian)
                                                       %in% fibrosisKegg[, 3]]
fibrosisAllBMDrutileHydrophilic201d <- keggAllBMDrutileHydrophilic201d[names(fibrosisAllBMDrutileHydrophilic201d)]

fibrosisAllBMDrutileHydrophilic201dMedian <- keggAllBMDrutileHydrophilic201dmedian[names(keggBMDrutileHydrophilic201dmedian)
                                                             %in% fibrosisKegg[, 3]]
fibrosisAllBMDrutileHydrophilic201dMedian <- t(rbind(fibrosisAllBMDrutileHydrophilic201dMedian))
fibrosisAllBMDrutileHydrophilic201dMedian <- rbind(fibrosisAllBMDrutileHydrophilic201dMedian)
names(fibrosisAllBMDrutileHydrophilic201dMedian) <- names(fibrosisAllBMDrutileHydrophilic201d)
fibrosisAllBMDrutileHydrophilic201dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophilic201dMedian)
colnames(fibrosisAllBMDrutileHydrophilic201dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDrutileHydrophilic201dMedianTop10  <- fibrosisAllBMDrutileHydrophilic201dMedian[order(
  fibrosisAllBMDrutileHydrophilic201dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDrutileHydrophilic201dMedianTop10 <- fibrosisAllBMDrutileHydrophilic201dMedianTop10[1:10]
fibrosisAllBMDrutileHydrophilic201dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophilic201dMedianTop10)
colnames(fibrosisAllBMDrutileHydrophilic201dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDrutileHydrophilic201dMedianTop10[, 1])
fibrosisAllBMDrutileHydrophilic201dMedianTop10 <- rbind(fibrosisAllBMDrutileHydrophilic201dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDrutileHydrophilic201dMedianTop10, file = "fibrosisAllBMDrutileHydrophilic201dMedianTop10.txt")


# rutileHydrophobic 20nm 1d -----------------------------------------------

keggBMDrutileHydrophobic201d <- list()
keggBMDrutileHydrophobic201dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDrutileHydrophobic201d[[i]] <- rutileHydrophobic201d[which(toupper(rutileHydrophobic201d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDrutileHydrophobic201d[[i]] <- as.numeric(keggBMDrutileHydrophobic201d[[i]])
  keggBMDrutileHydrophobic201dmedian[i] <- median(keggBMDrutileHydrophobic201d[[i]])
}
names(keggBMDrutileHydrophobic201d) <- names(keggGenes)
names(keggBMDrutileHydrophobic201dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDrutileHydrophobic201d <- list()
keggAllBMDrutileHydrophobic201dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDrutileHydrophobic201d[[i]] <- rutileHydrophobic201d[which(toupper(rutileHydrophobic201d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDrutileHydrophobic201d[[i]] <- as.numeric(keggAllBMDrutileHydrophobic201d[[i]])
  keggAllBMDrutileHydrophobic201dmedian[i] <- median(keggAllBMDrutileHydrophobic201d[[i]])
}
names(keggAllBMDrutileHydrophobic201d) <- names(keggGenesAll)
names(keggAllBMDrutileHydrophobic201dmedian) <- names(keggGenesAll)

fibrosisBMDrutileHydrophobic201d <- keggBMDrutileHydrophobic201d[names(keggBMDrutileHydrophobic201d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDrutileHydrophobic201dMedian <- keggBMDrutileHydrophobic201dmedian[names(keggBMDrutileHydrophobic201dmedian)
                                                       %in% fibrosisKegg[, 3]]
fibrosisAllBMDrutileHydrophobic201d <- keggAllBMDrutileHydrophobic201d[names(keggAllBMDrutileHydrophobic201d)
                                                 %in% fibrosisKegg[, 3]]

fibrosisAllBMDrutileHydrophobic201dMedian <- keggAllBMDrutileHydrophobic201dmedian[
                                          names(fibrosisAllBMDrutileHydrophobic201d)]

fibrosisAllBMDrutileHydrophobic201dMedian <- t(rbind(fibrosisAllBMDrutileHydrophobic201dMedian))
fibrosisAllBMDrutileHydrophobic201dMedian <- rbind(fibrosisAllBMDrutileHydrophobic201dMedian)
names(fibrosisAllBMDrutileHydrophobic201dMedian) <- names(fibrosisAllBMDrutileHydrophobic201d)
fibrosisAllBMDrutileHydrophobic201dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophobic201dMedian)
colnames(fibrosisAllBMDrutileHydrophobic201dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDrutileHydrophobic201dMedianTop10  <- fibrosisAllBMDrutileHydrophobic201dMedian[order(
  fibrosisAllBMDrutileHydrophobic201dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDrutileHydrophobic201dMedianTop10 <- fibrosisAllBMDrutileHydrophobic201dMedianTop10[1:10]
fibrosisAllBMDrutileHydrophobic201dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophobic201dMedianTop10)
colnames(fibrosisAllBMDrutileHydrophobic201dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDrutileHydrophobic201dMedianTop10[, 1])
fibrosisAllBMDrutileHydrophobic201dMedianTop10 <- rbind(fibrosisAllBMDrutileHydrophobic201dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDrutileHydrophobic201dMedianTop10, file = "fibrosisAllBMDrutileHydrophobic201dMedianTop10.txt")


# anatase 8nm 28d ----------------------------------------------------------

keggBMDanatase828d <- list()
keggBMDanatase828dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDanatase828d[[i]] <- anatase828d[which(toupper(anatase828d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDanatase828d[[i]] <- as.numeric(keggBMDanatase828d[[i]])
  keggBMDanatase828dmedian[i] <- median(keggBMDanatase828d[[i]])
}
names(keggBMDanatase828d) <- names(keggGenes)
names(keggBMDanatase828dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDanatase828d <- list()
keggAllBMDanatase828dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDanatase828d[[i]] <- anatase828d[which(toupper(anatase828d$SYMBOL) %in% 
                                                  keggGenesAll[[i]]), 3]
  keggAllBMDanatase828d[[i]] <- as.numeric(keggAllBMDanatase828d[[i]])
  keggAllBMDanatase828dmedian[i] <- median(keggAllBMDanatase828d[[i]])
}
names(keggAllBMDanatase828d) <- names(keggGenesAll)
names(keggAllBMDanatase828dmedian) <- names(keggGenesAll)

fibrosisBMDanatase828d <- keggBMDanatase828d[names(keggBMDanatase828d)
                                           %in% fibrosisKegg[, 3]]

fibrosisBMDanatase828dMedian <- keggBMDanatase828dmedian[names(keggBMDanatase828dmedian)
                                                       %in% fibrosisKegg[, 3]]
fibrosisAllBMDanatase828d <- keggAllBMDanatase828d[names(keggAllBMDanatase828d)
                                                 %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase828dMedian <- keggAllBMDanatase828dmedian[names(fibrosisAllBMDanatase828d)]
fibrosisAllBMDanatase828dMedian <- t(rbind(fibrosisAllBMDanatase828dMedian))
fibrosisAllBMDanatase828dMedian <- rbind(fibrosisAllBMDanatase828dMedian)
names(fibrosisAllBMDanatase828dMedian) <- names(fibrosisAllBMDanatase828d)
fibrosisAllBMDanatase828dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase828dMedian)
colnames(fibrosisAllBMDanatase828dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDanatase828dMedianTop10  <- fibrosisAllBMDanatase828dMedian[order(
  fibrosisAllBMDanatase828dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDanatase828dMedianTop10 <- fibrosisAllBMDanatase828dMedianTop10[1:10]
fibrosisAllBMDanatase828dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase828dMedianTop10)
colnames(fibrosisAllBMDanatase828dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDanatase828dMedianTop10[, 1])
fibrosisAllBMDanatase828dMedianTop10 <- rbind(fibrosisAllBMDanatase828dMedianTop10, 
                                             AOPmedianBMDt)
write.table(fibrosisAllBMDanatase828dMedianTop10, file = "fibrosisAllBMDanatase828dMedianTop10.txt")


# anatase 20nm 28d ---------------------------------------------------------


keggBMDanatase2028d <- list()
keggBMDanatase2028dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDanatase2028d[[i]] <- anatase2028d[which(toupper(anatase2028d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDanatase2028d[[i]] <- as.numeric(keggBMDanatase2028d[[i]])
  keggBMDanatase2028dmedian[i] <- median(keggBMDanatase2028d[[i]])
}
names(keggBMDanatase2028d) <- names(keggGenes)
names(keggBMDanatase2028dmedian) <- names(keggGenes)

keggAllBMDanatase2028d <- list()
keggAllBMDanatase2028dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDanatase2028d[[i]] <- anatase2028d[which(toupper(anatase2028d$SYMBOL) %in% 
                                                    keggGenesAll[[i]]), 3]
  keggAllBMDanatase2028d[[i]] <- as.numeric(keggAllBMDanatase2028d[[i]])
  keggAllBMDanatase2028dmedian[i] <- median(keggAllBMDanatase2028d[[i]])
}
names(keggAllBMDanatase2028d) <- names(keggGenesAll)
names(keggAllBMDanatase2028dmedian) <- names(keggGenesAll)

fibrosisBMDanatase2028d <- keggBMDanatase2028d[names(keggBMDanatase2028d)
                                             %in% fibrosisKegg[, 3]]

fibrosisBMDanatase2028dMedian <- keggBMDanatase2028dmedian[names(keggBMDanatase2028dmedian)
                                                         %in% fibrosisKegg[, 3]]
fibrosisAllBMDanatase2028d <- keggAllBMDanatase2028d[names(keggAllBMDanatase2028d)
                                                   %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase2028dMedian <- keggAllBMDanatase2028dmedian[names(fibrosisAllBMDanatase2028d)]
fibrosisAllBMDanatase2028dMedian <- t(rbind(fibrosisAllBMDanatase2028dMedian))
fibrosisAllBMDanatase2028dMedian <- rbind(fibrosisAllBMDanatase2028dMedian)
names(fibrosisAllBMDanatase2028dMedian) <- names(fibrosisAllBMDanatase2028d)
fibrosisAllBMDanatase2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase2028dMedian)
colnames(fibrosisAllBMDanatase2028dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDanatase2028dMedianTop10  <- fibrosisAllBMDanatase2028dMedian[order(
  fibrosisAllBMDanatase2028dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDanatase2028dMedianTop10 <- fibrosisAllBMDanatase2028dMedianTop10[1:10]
fibrosisAllBMDanatase2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase2028dMedianTop10)
colnames(fibrosisAllBMDanatase2028dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDanatase2028dMedianTop10[, 1])
fibrosisAllBMDanatase2028dMedianTop10 <- rbind(fibrosisAllBMDanatase2028dMedianTop10, 
                                              AOPmedianBMDt)
write.table(fibrosisAllBMDanatase2028dMedianTop10, file = "fibrosisAllBMDanatase2028dMedianTop10.txt")



# anatase 300nm 28d --------------------------------------------------------

keggBMDanatase30028d <- list()
keggBMDanatase30028dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDanatase30028d[[i]] <- anatase30028d[which(toupper(anatase30028d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDanatase30028d[[i]] <- as.numeric(keggBMDanatase30028d[[i]])
  keggBMDanatase30028dmedian[i] <- median(keggBMDanatase30028d[[i]])
}
names(keggBMDanatase30028d) <- names(keggGenes)
names(keggBMDanatase30028dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDanatase30028d <- list()
keggAllBMDanatase30028dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDanatase30028d[[i]] <- anatase30028d[which(toupper(anatase30028d$SYMBOL) %in% 
                                                      keggGenesAll[[i]]), 3]
  keggAllBMDanatase30028d[[i]] <- as.numeric(keggAllBMDanatase30028d[[i]])
  keggAllBMDanatase30028dmedian[i] <- median(keggAllBMDanatase30028d[[i]])
}
names(keggAllBMDanatase30028d) <- names(keggGenesAll)
names(keggAllBMDanatase30028dmedian) <- names(keggGenesAll)

fibrosisBMDanatase30028d <- keggBMDanatase30028d[names(keggBMDanatase30028d)
                                               %in% fibrosisKegg[, 3]]

fibrosisBMDanatase30028dMedian <- keggBMDanatase30028dmedian[names(keggBMDanatase30028dmedian)
                                                           %in% fibrosisKegg[, 3]]
fibrosisAllBMDanatase30028d <- keggAllBMDanatase30028d[names(keggAllBMDanatase30028d)
                                                     %in% fibrosisKegg[, 3]]

fibrosisAllBMDanatase30028dMedian <- keggAllBMDanatase30028dmedian[names(fibrosisAllBMDanatase30028d)]
fibrosisAllBMDanatase30028dMedian <- t(rbind(fibrosisAllBMDanatase30028dMedian))
fibrosisAllBMDanatase30028dMedian <- rbind(fibrosisAllBMDanatase30028dMedian)
names(fibrosisAllBMDanatase30028dMedian) <- names(fibrosisAllBMDanatase30028d)
fibrosisAllBMDanatase30028dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase30028dMedian)
colnames(fibrosisAllBMDanatase30028dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDanatase30028dMedianTop10  <- fibrosisAllBMDanatase30028dMedian[order(
  fibrosisAllBMDanatase30028dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDanatase30028dMedianTop10 <- fibrosisAllBMDanatase30028dMedianTop10[1:10]
fibrosisAllBMDanatase30028dMedianTop10 <- do.call(rbind, fibrosisAllBMDanatase30028dMedianTop10)
colnames(fibrosisAllBMDanatase30028dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDanatase30028dMedianTop10[, 1])
fibrosisAllBMDanatase30028dMedianTop10 <- rbind(fibrosisAllBMDanatase30028dMedianTop10, 
                                               AOPmedianBMDt)
write.table(fibrosisAllBMDanatase30028dMedianTop10, file = "fibrosisAllBMDanatase30028dMedianTop10.txt")


# mix 20nm 28d -------------------------------------------------------------

keggBMDmix2028d <- list()
keggBMDmix2028dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDmix2028d[[i]] <- mix2028d[which(toupper(mix2028d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDmix2028d[[i]] <- as.numeric(keggBMDmix2028d[[i]])
  keggBMDmix2028dmedian[i] <- median(keggBMDmix2028d[[i]])
}
names(keggBMDmix2028d) <- names(keggGenes)
names(keggBMDmix2028dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDmix2028d <- list()
keggAllBMDmix2028dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDmix2028d[[i]] <- mix2028d[which(toupper(mix2028d$SYMBOL) %in% 
                                            keggGenesAll[[i]]), 3]
  keggAllBMDmix2028d[[i]] <- as.numeric(keggAllBMDmix2028d[[i]])
  keggAllBMDmix2028dmedian[i] <- median(keggAllBMDmix2028d[[i]])
}
names(keggAllBMDmix2028d) <- names(keggGenesAll)
names(keggAllBMDmix2028dmedian) <- names(keggGenesAll)

fibrosisBMDmix2028d <- keggBMDmix2028d[names(keggBMDmix2028d)
                                     %in% fibrosisKegg[, 3]]

fibrosisBMDmix2028dMedian <- keggBMDmix2028dmedian[names(keggBMDmix2028dmedian)
                                                 %in% fibrosisKegg[, 3]]
fibrosisAllBMDmix2028d <- keggAllBMDmix2028d[names(keggAllBMDmix2028d)
                                           %in% fibrosisKegg[, 3]]

fibrosisAllBMDmix2028dMedian <- keggAllBMDmix2028dmedian[names(fibrosisAllBMDmix2028d)]
fibrosisAllBMDmix2028dMedian <- t(rbind(fibrosisAllBMDmix2028dMedian))
fibrosisAllBMDmix2028dMedian <- rbind(fibrosisAllBMDmix2028dMedian)
names(fibrosisAllBMDmix2028dMedian) <- names(fibrosisAllBMDmix2028d)
fibrosisAllBMDmix2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDmix2028dMedian)
colnames(fibrosisAllBMDmix2028dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDmix2028dMedianTop10  <- fibrosisAllBMDmix2028dMedian[order(
  fibrosisAllBMDmix2028dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDmix2028dMedianTop10 <- fibrosisAllBMDmix2028dMedianTop10[1:10]
fibrosisAllBMDmix2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDmix2028dMedianTop10)
colnames(fibrosisAllBMDmix2028dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDmix2028dMedianTop10[, 1])
fibrosisAllBMDmix2028dMedianTop10 <- rbind(fibrosisAllBMDmix2028dMedianTop10, 
                                          AOPmedianBMDt)
write.table(fibrosisAllBMDmix2028dMedianTop10, file = "fibrosisAllBMDmix2028dMedianTop10.txt")


# rutileHydrophilic 20nm 28d -----------------------------------------------


keggBMDrutileHydrophilic2028d <- list()
keggBMDrutileHydrophilic2028dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDrutileHydrophilic2028d[[i]] <- rutileHydrophilic2028d[which(toupper(rutileHydrophilic2028d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDrutileHydrophilic2028d[[i]] <- as.numeric(keggBMDrutileHydrophilic2028d[[i]])
  keggBMDrutileHydrophilic2028dmedian[i] <- median(keggBMDrutileHydrophilic2028d[[i]])
}
names(keggBMDrutileHydrophilic2028d) <- names(keggGenes)
names(keggBMDrutileHydrophilic2028dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDrutileHydrophilic2028d <- list()
keggAllBMDrutileHydrophilic2028dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDrutileHydrophilic2028d[[i]] <- rutileHydrophilic2028d[which(toupper(rutileHydrophilic2028d$SYMBOL) %in% 
                                                                        keggGenesAll[[i]]), 3]
  keggAllBMDrutileHydrophilic2028d[[i]] <- as.numeric(keggAllBMDrutileHydrophilic2028d[[i]])
  keggAllBMDrutileHydrophilic2028dmedian[i] <- median(keggAllBMDrutileHydrophilic2028d[[i]])
}
names(keggAllBMDrutileHydrophilic2028d) <- names(keggGenesAll)
names(keggAllBMDrutileHydrophilic2028dmedian) <- names(keggGenesAll)

fibrosisBMDrutileHydrophilic2028d <- keggBMDrutileHydrophilic2028d[names(keggBMDrutileHydrophilic2028d)
                                                                 %in% fibrosisKegg[, 3]]

fibrosisBMDrutileHydrophilic2028dMedian <- keggBMDrutileHydrophilic2028dmedian[names(keggBMDrutileHydrophilic2028dmedian)
                                                                             %in% fibrosisKegg[, 3]]
fibrosisAllBMDrutileHydrophilic2028d <- keggAllBMDrutileHydrophilic2028d[names(keggAllBMDrutileHydrophilic2028d)
                                                                       %in% fibrosisKegg[, 3]]

fibrosisAllBMDrutileHydrophilic2028dMedian <- keggAllBMDrutileHydrophilic2028dmedian[
                                          names(fibrosisAllBMDrutileHydrophilic2028d)]
fibrosisAllBMDrutileHydrophilic2028dMedian <- t(rbind(fibrosisAllBMDrutileHydrophilic2028dMedian))
fibrosisAllBMDrutileHydrophilic2028dMedian <- rbind(fibrosisAllBMDrutileHydrophilic2028dMedian)
names(fibrosisAllBMDrutileHydrophilic2028dMedian) <- names(fibrosisAllBMDrutileHydrophilic2028d)
fibrosisAllBMDrutileHydrophilic2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophilic2028dMedian)
colnames(fibrosisAllBMDrutileHydrophilic2028dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDrutileHydrophilic2028dMedianTop10  <- fibrosisAllBMDrutileHydrophilic2028dMedian[order(
  fibrosisAllBMDrutileHydrophilic2028dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDrutileHydrophilic2028dMedianTop10 <- fibrosisAllBMDrutileHydrophilic2028dMedianTop10[1:10]
fibrosisAllBMDrutileHydrophilic2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophilic2028dMedianTop10)
colnames(fibrosisAllBMDrutileHydrophilic2028dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDrutileHydrophilic2028dMedianTop10[, 1])
fibrosisAllBMDrutileHydrophilic2028dMedianTop10 <- rbind(fibrosisAllBMDrutileHydrophilic2028dMedianTop10, 
                                                        AOPmedianBMDt)
write.table(fibrosisAllBMDrutileHydrophilic2028dMedianTop10, file = "fibrosisAllBMDrutileHydrophilic2028dMedianTop10.txt")


# rutileHydrophobic 20nm 28d -----------------------------------------------

keggBMDrutileHydrophobic2028d <- list()
keggBMDrutileHydrophobic2028dmedian <- list()
for (i in 1:length(keggGenes)){
  keggBMDrutileHydrophobic2028d[[i]] <- rutileHydrophobic2028d[which(toupper(rutileHydrophobic2028d$SYMBOL) %in% keggGenes[[i]]), 3]
  keggBMDrutileHydrophobic2028d[[i]] <- as.numeric(keggBMDrutileHydrophobic2028d[[i]])
  keggBMDrutileHydrophobic2028dmedian[i] <- median(keggBMDrutileHydrophobic2028d[[i]])
}
names(keggBMDrutileHydrophobic2028d) <- names(keggGenes)
names(keggBMDrutileHydrophobic2028dmedian) <- names(keggGenes)

# Merge with Pathways
keggAllBMDrutileHydrophobic2028d <- list()
keggAllBMDrutileHydrophobic2028dmedian <- list()
for (i in 1:length(keggGenesAll)){
  keggAllBMDrutileHydrophobic2028d[[i]] <- rutileHydrophobic2028d[which(toupper(rutileHydrophobic2028d$SYMBOL) %in% 
                                                                        keggGenesAll[[i]]), 3]
  keggAllBMDrutileHydrophobic2028d[[i]] <- as.numeric(keggAllBMDrutileHydrophobic2028d[[i]])
  keggAllBMDrutileHydrophobic2028dmedian[i] <- median(keggAllBMDrutileHydrophobic2028d[[i]])
}
names(keggAllBMDrutileHydrophobic2028d) <- names(keggGenesAll)
names(keggAllBMDrutileHydrophobic2028dmedian) <- names(keggGenesAll)

fibrosisBMDrutileHydrophobic2028d <- keggBMDrutileHydrophobic2028d[names(keggBMDrutileHydrophobic2028d)
                                                                 %in% fibrosisKegg[, 3]]

fibrosisBMDrutileHydrophobic2028dMedian <- keggBMDrutileHydrophobic2028dmedian[names(keggBMDrutileHydrophobic2028dmedian)
                                                                             %in% fibrosisKegg[, 3]]
fibrosisAllBMDrutileHydrophobic2028d <- keggAllBMDrutileHydrophobic2028d[names(keggAllBMDrutileHydrophobic2028d)
                                                                       %in% fibrosisKegg[, 3]]

fibrosisAllBMDrutileHydrophobic2028dMedian <- keggAllBMDrutileHydrophobic2028dmedian[
                                              names(fibrosisAllBMDrutileHydrophobic2028d)]
fibrosisAllBMDrutileHydrophobic2028dMedian <- t(rbind(fibrosisAllBMDrutileHydrophobic2028dMedian))
fibrosisAllBMDrutileHydrophobic2028dMedian <- rbind(fibrosisAllBMDrutileHydrophobic2028dMedian)
names(fibrosisAllBMDrutileHydrophobic2028dMedian) <- names(fibrosisAllBMDrutileHydrophobic2028d)
fibrosisAllBMDrutileHydrophobic2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophobic2028dMedian)
colnames(fibrosisAllBMDrutileHydrophobic2028dMedianTop10) <- "medianAOPpathwaysBMDt"  
fibrosisAllBMDrutileHydrophobic2028dMedianTop10  <- fibrosisAllBMDrutileHydrophobic2028dMedian[order(
  fibrosisAllBMDrutileHydrophobic2028dMedianTop10[, 1], decreasing = FALSE), ]

fibrosisAllBMDrutileHydrophobic2028dMedianTop10 <- fibrosisAllBMDrutileHydrophobic2028dMedianTop10[1:10]
fibrosisAllBMDrutileHydrophobic2028dMedianTop10 <- do.call(rbind, fibrosisAllBMDrutileHydrophobic2028dMedianTop10)
colnames(fibrosisAllBMDrutileHydrophobic2028dMedianTop10) <- "medianAOPpathwaysBMDt"
AOPmedianBMDt <- median(fibrosisAllBMDrutileHydrophobic2028dMedianTop10[, 1])
fibrosisAllBMDrutileHydrophobic2028dMedianTop10 <- rbind(fibrosisAllBMDrutileHydrophobic2028dMedianTop10, 
                                                        AOPmedianBMDt)
write.table(fibrosisAllBMDrutileHydrophobic2028dMedianTop10, file = "fibrosisAllBMDrutileHydrophobic2028dMedianTop10.txt")


save.image("rah16Bmdl.RData")



