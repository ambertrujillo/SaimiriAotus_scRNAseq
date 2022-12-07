#########################################################################
############## Zinbwave to make lots of zeros useful ####################
#########################################################################

load("RBC_clustering.Rdata")

# Format metadata
parasite_load <- FetchData(object = RBC.combined, vars = c("parasite_load")) # gets list of all cells and their SRR numbers
species <- FetchData(object = RBC.combined, vars = c("species")) # gets list of all cells and their SRR numbers
individual <- FetchData(object = RBC.combined, vars = c("individual"))

SRRs <- FetchData(object = RBC.combined, vars = c("orig.ident")) # gets list of all cells and their SRR numbers

parasite_info <- read.csv("parasite_info.csv")

chloroquine_resistance = subset(parasite_info, select=c(chloroquine_resistance))
vivax_strain = subset(parasite_info, select=c(vivax_strain))

# Load data for zinbwave
zin_RBC <- SummarizedExperiment(list(counts=as.matrix(GetAssayData(object = RBC.combined, slot = "counts"))), colData=c(parasite_load, species, individual, chloroquine_resistance, vivax_strain))

zin_RBC

#class: SummarizedExperiment 
#dim: 16294 924 
#metadata(0):
#assays(1): counts
#rownames(16294): ENSSBOG00000026468 ENSSBOG00000029669 ...
#  ENSSBOG00000035683 ENSSBOG00000018438
#rowData names(0):
#colnames(924): AAGTCTGGTGTGTGCC.1_1 ACACTGACACCGGAAA.1_1 ...
#  TGGCCAGTCAACCAAC.1_10 TTAGGCAAGGGCATGT.1_10
#colData names(5): parasite_load species individual
#  chloroquine_resistance vivax_strain

# Filter lowloy expressed genes, genes that do not have at least 1 reads in at least 1 samples
filter <- rowSums(assay(zin_RBC)>1)>1
table(filter)
#FALSE  TRUE (We have 527 genes)
#16105   189 

zin_RBC <- zin_RBC[filter,]

# Identify the 100 most variable genes, which will be the input for zinbwave. 
#assay(Cluster9) %>% log1p %>% rowVars -> vars
#names(vars) <- rownames(Cluster9)
#vars <- sort(vars, decreasing = TRUE)
#head(vars)
#             Aotus-BEST1 Aotus-ENSANAG00000026505              Aotus-FBXW8 
#               2.8992824                1.4142759                1.3896312 
#             Aotus-ERCC6 Aotus-ENSANAG00000029701              Aotus-GNPTG 
#               1.3045575                1.1114966                0.9995429 

#Cluster9 <- Cluster9[names(vars)[1:100],]

# Remove cells that have 0 gene counts, otherwise it will throw: "Sample 3 has only 0 counts!"
filterCols <- colSums(assay(zin_RBC)) > 0
zin_RBC <- zin_RBC[,filterCols]

assayNames(zin_RBC)[1] <- "counts" #To avoid needing to specify which assay we should use for the zinbwave workflow

# ZINB-WaVE
#zinb <- zinbFit(Cluster9, K=0)
#table(colData(Cluster9)$orig.ident)
#SRR11008270 SRR11008271 SRR11008274 SRR11008275 SRR11008277 SRR11008278 
#         23           9          14           8          18          10 

#table(rowSums(assay(Cluster9)) == 0)
#which(rowSums(assay(Cluster9)) == 0)


# Using strict filter >5)>5 = there are < 100 genes and filtering based on the 100 most variable genes doesn't work
# Including all genes, not just 100 most variable = Sample 8 has only 0 counts! (same error different sample)
# Use zinbFit first, throws same error on zinbFit step


# getting the data ready for edgeR

zinb <- zinbFit(zin_RBC, K=2, epsilon=1000)
RBC_zinb <- zinbwave(zin_RBC, fitted_model = zinb, K = 6, epsilon=1000, observationalWeights = TRUE)

# View resulting W matrix

W <- reducedDim(RBC_zinb)

pdf(file="Plots/cluster_plot_allgenes.pdf")
data.frame(W, bio=colData(zin_RBC)$species,
           parasite_load=colData(zin_RBC)$parasite_load) %>%
    ggplot(aes(W1, W2, colour=parasite_load, shape=bio)) + geom_point() + 
    scale_color_gradient(low="blue", high="red") + theme_classic()
dev.off()

pdf(file="Plots/cluster_plot_allgenes_individuals.pdf")
data.frame(W, bio=colData(zin_RBC)$species,
           individual=colData(zin_RBC)$individual) %>%
    ggplot(aes(W1, W2, colour=individual, shape=bio)) + geom_point() + 
    scale_color_discrete() + theme_classic()
dev.off()

pdf(file="Plots/cluster_plot_allgenes_parasite_info.pdf")
data.frame(W, bio=colData(zin_RBC)$species,
           strain=colData(zin_RBC)$vivax_strain) %>%
    ggplot(aes(W1, W2, colour=strain, shape=bio)) + geom_point() + 
    scale_color_discrete() + theme_classic()
dev.off()


library(Rtsne)
tsne_data <- Rtsne(W, pca = FALSE, perplexity=10, max_iter=5000, check_duplicates = FALSE)

pdf(file="Plots/tsne_plot_allgenes.pdf")
data.frame(Dim1=tsne_data$Y[,1], Dim2=tsne_data$Y[,2], 
           bio=colData(zin_RBC)$species,
           parasite_load=colData(zin_RBC)$parasite_load) %>%
    ggplot(aes(Dim1, Dim2, colour=parasite_load, shape=bio)) + geom_point() + 
    scale_color_gradient(low="blue", high="red") + theme_classic() +
    ggtitle("t-SNE")
dev.off()   

# Differential Expression
#calculate weights
weights <- assay(RBC_zinb, "weights")

library(edgeR)

dge <- DGEList(assay(RBC_zinb))
dge <- calcNormFactors(dge)

design <- model.matrix(~parasite_load + species + parasite_load:species, data = colData(zin_RBC))
dge$weights <- weights

dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)

#Considering each coefficient
lrt_pl <- glmWeightedF(fit, coef = "parasite_load")
lrt_species <- glmWeightedF(fit, coef = "speciesSaimiri")
lrt_plSpecies <- glmWeightedF(fit, coef = "parasite_load:speciesSaimiri")

tt_pl = topTags(lrt_pl, n=Inf, adjust.method="BH", p.value=1) # infinity tags for functional enrichment
tt_species = topTags(lrt_species, n=Inf, adjust.method="BH", p.value=1) # infinity tags for functional enrichment
tt_plSpecies = topTags(lrt_plSpecies, n=Inf, adjust.method="BH", p.value=1) # infinity tags for functional enrichment

tt_pl = data.frame(tt_pl)
tt_species = data.frame(tt_species)
tt_plSpecies = data.frame(tt_plSpecies)

sig.genes.pl = subset(tt_pl, PValue < 0.05, select=c(logFC, PValue, FDR))
dn.genes.pl = subset(sig.genes.pl, logFC < 0 , select=c(logFC, PValue, FDR))
up.genes.pl = subset(sig.genes.pl, logFC > 0 , select=c(logFC, PValue, FDR))

sig.genes.species = subset(tt_species, PValue < 0.05, select=c(logFC, PValue, FDR))
aotus.genes.species = subset(sig.genes.species, logFC < 0 , select=c(logFC, PValue, FDR))
saimiri.genes.species = subset(sig.genes.species, logFC > 0 , select=c(logFC, PValue, FDR))

sig.genes.plSpecies = subset(tt_plSpecies, PValue < 0.05, select=c(logFC, PValue, FDR))
dn.genes.plSpecies = subset(sig.genes.plSpecies, logFC < 0 , select=c(logFC, PValue, FDR))
up.genes.plSpecies = subset(sig.genes.plSpecies, logFC > 0 , select=c(logFC, PValue, FDR))

#coef = parasite_load
write.table(sig.genes.pl,
            file=paste0("coef_pl/sig_DE_genes_pl", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(dn.genes.pl,
            file=paste0("coef_pl/sig_dn_DE_genes_pl", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(up.genes.pl,
            file=paste0("coef_pl/sig_up_DE_genes_pl", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

#coef = species
write.table(sig.genes.species,
            file=paste0("coef_species/sig_DE_genes_species", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(aotus.genes.species,
            file=paste0("coef_species/sig_aotus_DE_genes_species", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(saimiri.genes.species,
            file=paste0("coef_species/sig_saimiri_DE_genes_species", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

#coef = plSpecies
write.table(sig.genes.plSpecies,
            file=paste0("coef_plSpecies/sig_DE_genes_plSpecies", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(dn.genes.plSpecies,
            file=paste0("coef_plSpecies/sig_dn_DE_genes_plSpecies", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)

write.table(up.genes.plSpecies,
            file=paste0("coef_plSpecies/sig_up_DE_genes_plSpecies", suffix=".txt"),
            row.names=TRUE, col.names=TRUE, quote=FALSE)


save.image("after_zinbwave.Rdata")
