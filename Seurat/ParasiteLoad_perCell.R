load("Orthologous_matrices.Rdata")

# Seperate host and pathogen data
# Aotus 
#SRR11008270
Aotus_subset.SRR11008270 <- rownames(SRR11008270.data)[grep("^Aotus", rownames(SRR11008270.data))]
head(Aotus_subset.SRR11008270, 10)
Aotus_PL.SRR11008270 <- SRR11008270.data[Aotus_subset.SRR11008270, ]
Aotus_PL.SRR11008270 <- data.frame(colSums(Aotus_PL.SRR11008270))
Aotus_PL.SRR11008270$cells <- rownames(Aotus_PL.SRR11008270)

Vivax_subset.SRR11008270 <- rownames(SRR11008270.data)[grep("^Vivax", rownames(SRR11008270.data))]
Vivax_PL.SRR11008270 <- SRR11008270.data[Vivax_subset.SRR11008270, ]
Vivax_PL.SRR11008270 <- data.frame(colSums(Vivax_PL.SRR11008270))
Vivax_PL.SRR11008270$cells <- rownames(Vivax_PL.SRR11008270)

PL.SRR11008270 <- merge(Aotus_PL.SRR11008270, Vivax_PL.SRR11008270, by="cells")

PL.SRR11008270$parasite_load_SRR11008270 = PL.SRR11008270$colSums.Vivax_PL.SRR11008270. / PL.SRR11008270$colSums.Aotus_PL.SRR11008270.

write.csv(PL.SRR11008270, file="../Zinbwave/PL.SRR11008270", quote=FALSE)

#SRR11008271
Aotus_subset.SRR11008271 <- rownames(SRR11008271.data)[grep("^Aotus", rownames(SRR11008271.data))]
Aotus_PL.SRR11008271 <- SRR11008271.data[Aotus_subset.SRR11008271, ]
Aotus_PL.SRR11008271 <- data.frame(colSums(Aotus_PL.SRR11008271))
Aotus_PL.SRR11008271$cells <- rownames(Aotus_PL.SRR11008271)

Vivax_subset.SRR11008271 <- rownames(SRR11008271.data)[grep("^Vivax", rownames(SRR11008271.data))]
Vivax_PL.SRR11008271 <- SRR11008271.data[Vivax_subset.SRR11008271, ]
Vivax_PL.SRR11008271 <- data.frame(colSums(Vivax_PL.SRR11008271))
Vivax_PL.SRR11008271$cells <- rownames(Vivax_PL.SRR11008271)

PL.SRR11008271 <- merge(Aotus_PL.SRR11008271, Vivax_PL.SRR11008271, by="cells")
PL.SRR11008271$parasite_load_SRR11008271 = PL.SRR11008271$colSums.Vivax_PL.SRR11008271. / PL.SRR11008271$colSums.Aotus_PL.SRR11008271.

write.csv(PL.SRR11008271, file="../Zinbwave/PL.SRR11008271", quote=FALSE)

#SRR11008274
Aotus_subset.SRR11008274 <- rownames(SRR11008274.data)[grep("^Aotus", rownames(SRR11008274.data))]
Aotus_PL.SRR11008274 <- SRR11008274.data[Aotus_subset.SRR11008274, ]
Aotus_PL.SRR11008274 <- data.frame(colSums(Aotus_PL.SRR11008274))
Aotus_PL.SRR11008274$cells <- rownames(Aotus_PL.SRR11008274)

Vivax_subset.SRR11008274 <- rownames(SRR11008274.data)[grep("^Vivax", rownames(SRR11008274.data))]
Vivax_PL.SRR11008274 <- SRR11008274.data[Vivax_subset.SRR11008274, ]
Vivax_PL.SRR11008274 <- data.frame(colSums(Vivax_PL.SRR11008274))

Vivax_PL.SRR11008274$cells <- rownames(Vivax_PL.SRR11008274)

PL.SRR11008274 <- merge(Aotus_PL.SRR11008274, Vivax_PL.SRR11008274, by="cells")
PL.SRR11008274$parasite_load_SRR11008274 = PL.SRR11008274$colSums.Vivax_PL.SRR11008274. / PL.SRR11008274$colSums.Aotus_PL.SRR11008274.

write.csv(PL.SRR11008274, file="../Zinbwave/PL.SRR11008274", quote=FALSE)

#SRR11008275
Aotus_subset.SRR11008275 <- rownames(SRR11008275.data)[grep("^Aotus", rownames(SRR11008275.data))]
Aotus_PL.SRR11008275 <- SRR11008275.data[Aotus_subset.SRR11008275, ]
Aotus_PL.SRR11008275 <- data.frame(colSums(Aotus_PL.SRR11008275))
Aotus_PL.SRR11008275$cells <- rownames(Aotus_PL.SRR11008275)

Vivax_subset.SRR11008275 <- rownames(SRR11008275.data)[grep("^Vivax", rownames(SRR11008275.data))]
Vivax_PL.SRR11008275 <- SRR11008275.data[Vivax_subset.SRR11008275, ]
Vivax_PL.SRR11008275 <- data.frame(colSums(Vivax_PL.SRR11008275))
Vivax_PL.SRR11008275$cells <- rownames(Vivax_PL.SRR11008275)

PL.SRR11008275 <- merge(Aotus_PL.SRR11008275, Vivax_PL.SRR11008275, by="cells")
PL.SRR11008275$parasite_load_SRR11008275 = PL.SRR11008275$colSums.Vivax_PL.SRR11008275. / PL.SRR11008275$colSums.Aotus_PL.SRR11008275.

write.csv(PL.SRR11008275, file="../Zinbwave/PL.SRR11008275", quote=FALSE)

#SRR11008277
Aotus_subset.SRR11008277 <- rownames(SRR11008277.data)[grep("^Aotus", rownames(SRR11008277.data))]
Aotus_PL.SRR11008277 <- SRR11008277.data[Aotus_subset.SRR11008277, ]
Aotus_PL.SRR11008277 <- data.frame(colSums(Aotus_PL.SRR11008277))
Aotus_PL.SRR11008277$cells <- rownames(Aotus_PL.SRR11008277)

Vivax_subset.SRR11008277 <- rownames(SRR11008277.data)[grep("^Vivax", rownames(SRR11008277.data))]
Vivax_PL.SRR11008277 <- SRR11008277.data[Vivax_subset.SRR11008277, ]
Vivax_PL.SRR11008277 <- data.frame(colSums(Vivax_PL.SRR11008277))
Vivax_PL.SRR11008277$cells <- rownames(Vivax_PL.SRR11008277)

PL.SRR11008277 <- merge(Aotus_PL.SRR11008277, Vivax_PL.SRR11008277, by="cells")
PL.SRR11008277$parasite_load_SRR11008277 = PL.SRR11008277$colSums.Vivax_PL.SRR11008277. / PL.SRR11008277$colSums.Aotus_PL.SRR11008277.

write.csv(PL.SRR11008277, file="../Zinbwave/PL.SRR11008277", quote=FALSE)

#SRR11008278
Aotus_subset.SRR11008278 <- rownames(SRR11008278.data)[grep("^Aotus", rownames(SRR11008278.data))]
Aotus_PL.SRR11008278 <- SRR11008278.data[Aotus_subset.SRR11008278, ]
Aotus_PL.SRR11008278 <- data.frame(colSums(Aotus_PL.SRR11008278))
Aotus_PL.SRR11008278$cells <- rownames(Aotus_PL.SRR11008278)

Vivax_subset.SRR11008278 <- rownames(SRR11008278.data)[grep("^Vivax", rownames(SRR11008278.data))]
Vivax_PL.SRR11008278 <- SRR11008278.data[Vivax_subset.SRR11008278, ]
Vivax_PL.SRR11008278 <- data.frame(colSums(Vivax_PL.SRR11008278))
Vivax_PL.SRR11008278$cells <- rownames(Vivax_PL.SRR11008278)

PL.SRR11008278 <- merge(Aotus_PL.SRR11008278, Vivax_PL.SRR11008278, by="cells")
PL.SRR11008278$parasite_load_SRR11008278 = PL.SRR11008278$colSums.Vivax_PL.SRR11008278. / PL.SRR11008278$colSums.Aotus_PL.SRR11008278.

write.csv(PL.SRR11008278, file="../Zinbwave/PL.SRR11008278", quote=FALSE)

# Saimiri
#SRR11008269
Saimiri_subset.SRR11008269 <- rownames(SRR11008269.data)[grep("^Saimiri", rownames(SRR11008269.data))]
Saimiri_PL.SRR11008269 <- SRR11008269.data[Saimiri_subset.SRR11008269, ]
Saimiri_PL.SRR11008269 <- data.frame(colSums(Saimiri_PL.SRR11008269))
Saimiri_PL.SRR11008269$cells <- rownames(Saimiri_PL.SRR11008269)

Vivax_subset.SRR11008269 <- rownames(SRR11008269.data)[grep("^Vivax", rownames(SRR11008269.data))]
Vivax_PL.SRR11008269 <- SRR11008269.data[Vivax_subset.SRR11008269, ]
Vivax_PL.SRR11008269 <- data.frame(colSums(Vivax_PL.SRR11008269))
Vivax_PL.SRR11008269$cells <- rownames(Vivax_PL.SRR11008269)

PL.SRR11008269 <- merge(Saimiri_PL.SRR11008269, Vivax_PL.SRR11008269, by="cells")
PL.SRR11008269$parasite_load_SRR11008269 = PL.SRR11008269$colSums.Vivax_PL.SRR11008269. / PL.SRR11008269$colSums.Saimiri_PL.SRR11008269.

write.csv(PL.SRR11008269, file="../Zinbwave/PL.SRR11008269", quote=FALSE)

#SRR11008272
Saimiri_subset.SRR11008272 <- rownames(SRR11008272.data)[grep("^Saimiri", rownames(SRR11008272.data))]
Saimiri_PL.SRR11008272 <- SRR11008272.data[Saimiri_subset.SRR11008272, ]
Saimiri_PL.SRR11008272 <- data.frame(colSums(Saimiri_PL.SRR11008272))
Saimiri_PL.SRR11008272$cells <- rownames(Saimiri_PL.SRR11008272)

Vivax_subset.SRR11008272 <- rownames(SRR11008272.data)[grep("^Vivax", rownames(SRR11008272.data))]
Vivax_PL.SRR11008272 <- SRR11008272.data[Vivax_subset.SRR11008272, ]
Vivax_PL.SRR11008272 <- data.frame(colSums(Vivax_PL.SRR11008272))
Vivax_PL.SRR11008272$cells <- rownames(Vivax_PL.SRR11008272)

PL.SRR11008272 <- merge(Saimiri_PL.SRR11008272, Vivax_PL.SRR11008272, by="cells")
PL.SRR11008272$parasite_load_SRR11008272 = PL.SRR11008272$colSums.Vivax_PL.SRR11008272. / PL.SRR11008272$colSums.Saimiri_PL.SRR11008272.

write.csv(PL.SRR11008272, file="../Zinbwave/PL.SRR11008272", quote=FALSE)

#SRR11008273
Saimiri_subset.SRR11008273 <- rownames(SRR11008273.data)[grep("^Saimiri", rownames(SRR11008273.data))]
Saimiri_PL.SRR11008273 <- SRR11008273.data[Saimiri_subset.SRR11008273, ]
Saimiri_PL.SRR11008273 <- data.frame(colSums(Saimiri_PL.SRR11008273))
Saimiri_PL.SRR11008273$cells <- rownames(Saimiri_PL.SRR11008273)

Vivax_subset.SRR11008273 <- rownames(SRR11008273.data)[grep("^Vivax", rownames(SRR11008273.data))]
Vivax_PL.SRR11008273 <- SRR11008273.data[Vivax_subset.SRR11008273, ]
Vivax_PL.SRR11008273 <- data.frame(colSums(Vivax_PL.SRR11008273))
Vivax_PL.SRR11008273$cells <- rownames(Vivax_PL.SRR11008273)

PL.SRR11008273 <- merge(Saimiri_PL.SRR11008273, Vivax_PL.SRR11008273, by="cells")
PL.SRR11008273$parasite_load_SRR11008273 = PL.SRR11008273$colSums.Vivax_PL.SRR11008273. / PL.SRR11008273$colSums.Saimiri_PL.SRR11008273.

write.csv(PL.SRR11008273, file="../Zinbwave/PL.SRR11008273", quote=FALSE)

#SRR11008276
Saimiri_subset.SRR11008276 <- rownames(SRR11008276.data)[grep("^Saimiri", rownames(SRR11008276.data))]
Saimiri_PL.SRR11008276 <- SRR11008276.data[Saimiri_subset.SRR11008276, ]
Saimiri_PL.SRR11008276 <- data.frame(colSums(Saimiri_PL.SRR11008276))
Saimiri_PL.SRR11008276$cells <- rownames(Saimiri_PL.SRR11008276)

Vivax_subset.SRR11008276 <- rownames(SRR11008276.data)[grep("^Vivax", rownames(SRR11008276.data))]
Vivax_PL.SRR11008276 <- SRR11008276.data[Vivax_subset.SRR11008276, ]
Vivax_PL.SRR11008276 <- data.frame(colSums(Vivax_PL.SRR11008276))
Vivax_PL.SRR11008276$cells <- rownames(Vivax_PL.SRR11008276)

PL.SRR11008276 <- merge(Saimiri_PL.SRR11008276, Vivax_PL.SRR11008276, by="cells")
PL.SRR11008276$parasite_load_SRR11008276 = PL.SRR11008276$colSums.Vivax_PL.SRR11008276. / PL.SRR11008276$colSums.Saimiri_PL.SRR11008276.

write.csv(PL.SRR11008276, file="../Zinbwave/PL.SRR11008276", quote=FALSE)

##### make Seurat object for Log normalization
Aotus_data.SRR11008270.log <- CreateSeuratObject(counts = Aotus_data.SRR11008270, min.cells = 0, min.genes = 0, project = "SRR11008270")
Aotus_data.SRR11008270.log <- AddMetaData(object = Aotus_data.SRR11008270.log, metadata = PL.SRR11008270$parasite_load_SRR11008270, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Aotus_data.SRR11008270.log <- AddMetaData(object = Aotus_data.SRR11008270.log, metadata = PL.SRR11008270$colSums.Vivax_PL.SRR11008270., col.name = "parasite_reads")

Aotus_data.SRR11008271.log <- CreateSeuratObject(counts = Aotus_data.SRR11008271, min.cells = 0, min.genes = 0, project = "SRR11008271")
Aotus_data.SRR11008271.log <- AddMetaData(object = Aotus_data.SRR11008271.log, metadata = PL.SRR11008271$parasite_load_SRR11008271, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Aotus_data.SRR11008271.log <- AddMetaData(object = Aotus_data.SRR11008271.log, metadata = PL.SRR11008271$colSums.Vivax_PL.SRR11008271., col.name = "parasite_reads")

Aotus_data.SRR11008274.log <- CreateSeuratObject(counts = Aotus_data.SRR11008274, min.cells = 0, min.genes = 0, project = "SRR11008274")
Aotus_data.SRR11008274.log <- AddMetaData(object = Aotus_data.SRR11008274.log, metadata = PL.SRR11008274$parasite_load_SRR11008274, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Aotus_data.SRR11008274.log <- AddMetaData(object = Aotus_data.SRR11008274.log, metadata = PL.SRR11008274$colSums.Vivax_PL.SRR11008274., col.name = "parasite_reads")

Aotus_data.SRR11008275.log <- CreateSeuratObject(counts = Aotus_data.SRR11008275, min.cells = 0, min.genes = 0, project = "SRR11008275")
Aotus_data.SRR11008275.log <- AddMetaData(object = Aotus_data.SRR11008275.log, metadata = PL.SRR11008275$parasite_load_SRR11008275, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Aotus_data.SRR11008275.log <- AddMetaData(object = Aotus_data.SRR11008275.log, metadata = PL.SRR11008275$colSums.Vivax_PL.SRR11008275., col.name = "parasite_reads")

Aotus_data.SRR11008277.log <- CreateSeuratObject(counts = Aotus_data.SRR11008277, min.cells = 0, min.genes = 0, project = "SRR11008277")
Aotus_data.SRR11008277.log <- AddMetaData(object = Aotus_data.SRR11008277.log, metadata = PL.SRR11008277$parasite_load_SRR11008277, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Aotus_data.SRR11008277.log <- AddMetaData(object = Aotus_data.SRR11008277.log, metadata = PL.SRR11008277$colSums.Vivax_PL.SRR11008277., col.name = "parasite_reads")

Aotus_data.SRR11008278.log <- CreateSeuratObject(counts = Aotus_data.SRR11008278, min.cells = 0, min.genes = 0, project = "SRR11008278")
Aotus_data.SRR11008278.log <- AddMetaData(object = Aotus_data.SRR11008278.log, metadata = PL.SRR11008278$parasite_load_SRR11008278, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Aotus_data.SRR11008278.log <- AddMetaData(object = Aotus_data.SRR11008278.log, metadata = PL.SRR11008278$colSums.Vivax_PL.SRR11008278., col.name = "parasite_reads")

Saimiri_data.SRR11008269.log <- CreateSeuratObject(counts = Saimiri_data.SRR11008269, min.cells = 0, min.genes = 0, project = "SRR11008269")
Saimiri_data.SRR11008269.log <- AddMetaData(object = Saimiri_data.SRR11008269.log, metadata = PL.SRR11008269$parasite_load_SRR11008269, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Saimiri_data.SRR11008269.log <- AddMetaData(object = Saimiri_data.SRR11008269.log, metadata = PL.SRR11008269$colSums.Vivax_PL.SRR11008269., col.name = "parasite_reads")

Saimiri_data.SRR11008272.log <- CreateSeuratObject(counts = Saimiri_data.SRR11008272, min.cells = 0, min.genes = 0, project = "SRR11008272")
Saimiri_data.SRR11008272.log <- AddMetaData(object = Saimiri_data.SRR11008272.log, metadata = PL.SRR11008272$parasite_load_SRR11008272, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Saimiri_data.SRR11008272.log <- AddMetaData(object = Saimiri_data.SRR11008272.log, metadata = PL.SRR11008272$colSums.Vivax_PL.SRR11008272., col.name = "parasite_reads")

Saimiri_data.SRR11008273.log <- CreateSeuratObject(counts = Saimiri_data.SRR11008273, min.cells = 0, min.genes = 0, project = "SRR11008273")
Saimiri_data.SRR11008273.log <- AddMetaData(object = Saimiri_data.SRR11008273.log, metadata = PL.SRR11008273$parasite_load_SRR11008273, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Saimiri_data.SRR11008273.log <- AddMetaData(object = Saimiri_data.SRR11008273.log, metadata = PL.SRR11008273$colSums.Vivax_PL.SRR11008273., col.name = "parasite_reads")

Saimiri_data.SRR11008276.log <- CreateSeuratObject(counts = Saimiri_data.SRR11008276, min.cells = 0, min.genes = 0, project = "SRR11008276")
Saimiri_data.SRR11008276.log <- AddMetaData(object = Saimiri_data.SRR11008276.log, metadata = PL.SRR11008276$parasite_load_SRR11008276, col.name = "parasite_load") ### I have made sure the cell types match up between datasets
Saimiri_data.SRR11008276.log <- AddMetaData(object = Saimiri_data.SRR11008276.log, metadata = PL.SRR11008276$colSums.Vivax_PL.SRR11008276., col.name = "parasite_reads")


#processing Seurat - with lognormalization
Aotus_data.SRR11008270.log <- NormalizeData(Aotus_data.SRR11008270.log, normalization.method = "LogNormalize", scale.factor = 10000)
Aotus_data.SRR11008270.log <- FindVariableFeatures(Aotus_data.SRR11008270.log, selection.method = "vst", nfeatures = 2000)
SRR11008270.Aotus.scaled <- ScaleData(Aotus_data.SRR11008270.log)
SRR11008270.Aotus.scaled[["individual"]] <- "SRR11008270"
SRR11008270.Aotus.scaled[["species"]] <- "Aotus"

Aotus_data.SRR11008271.log <- NormalizeData(Aotus_data.SRR11008271.log, normalization.method = "LogNormalize", scale.factor = 10000)
Aotus_data.SRR11008271.log <- FindVariableFeatures(Aotus_data.SRR11008271.log, selection.method = "vst", nfeatures = 2000)
SRR11008271.Aotus.scaled <- ScaleData(Aotus_data.SRR11008271.log)
SRR11008271.Aotus.scaled[["individual"]] <- "SRR11008271"
SRR11008271.Aotus.scaled[["species"]] <- "Aotus"

Aotus_data.SRR11008274.log  <- NormalizeData(Aotus_data.SRR11008274.log, normalization.method = "LogNormalize", scale.factor = 10000)
Aotus_data.SRR11008274.log <- FindVariableFeatures(Aotus_data.SRR11008274.log, selection.method = "vst", nfeatures = 2000)
SRR11008274.Aotus.scaled <- ScaleData(Aotus_data.SRR11008274.log)
SRR11008274.Aotus.scaled[["individual"]] <- "SRR11008274"
SRR11008274.Aotus.scaled[["species"]] <- "Aotus"

Aotus_data.SRR11008275.log <- NormalizeData(Aotus_data.SRR11008275.log, normalization.method = "LogNormalize", scale.factor = 10000)
Aotus_data.SRR11008275.log <- FindVariableFeatures(Aotus_data.SRR11008275.log, selection.method = "vst", nfeatures = 2000)
SRR11008275.Aotus.scaled <- ScaleData(Aotus_data.SRR11008275.log)
SRR11008275.Aotus.scaled[["individual"]] <- "SRR11008275"
SRR11008275.Aotus.scaled[["species"]] <- "Aotus"

Aotus_data.SRR11008277.log <- NormalizeData(Aotus_data.SRR11008277.log, normalization.method = "LogNormalize", scale.factor = 10000)
Aotus_data.SRR11008277.log <- FindVariableFeatures(Aotus_data.SRR11008277.log, selection.method = "vst", nfeatures = 2000)
SRR11008277.Aotus.scaled <- ScaleData(Aotus_data.SRR11008277.log)
SRR11008277.Aotus.scaled[["individual"]] <- "SRR11008277"
SRR11008277.Aotus.scaled[["species"]] <- "Aotus"

Aotus_data.SRR11008278.log <- NormalizeData(Aotus_data.SRR11008278.log, normalization.method = "LogNormalize", scale.factor = 10000)
Aotus_data.SRR11008278.log <- FindVariableFeatures(Aotus_data.SRR11008278.log, selection.method = "vst", nfeatures = 2000)
SRR11008278.Aotus.scaled <- ScaleData(Aotus_data.SRR11008278.log)
SRR11008278.Aotus.scaled[["individual"]] <- "SRR11008278"
SRR11008278.Aotus.scaled[["species"]] <- "Aotus"

Saimiri_data.SRR11008269.log <- NormalizeData(Saimiri_data.SRR11008269.log, normalization.method = "LogNormalize", scale.factor = 10000)
Saimiri_data.SRR11008269.log <- FindVariableFeatures(Saimiri_data.SRR11008269.log, selection.method = "vst", nfeatures = 2000)
SRR11008269.Saimiri.scaled <- ScaleData(Saimiri_data.SRR11008269.log)
SRR11008269.Saimiri.scaled[["individual"]] <- "SRR11008269"
SRR11008269.Saimiri.scaled[["species"]] <- "Saimiri"

Saimiri_data.SRR11008272.log <- NormalizeData(Saimiri_data.SRR11008272.log, normalization.method = "LogNormalize", scale.factor = 10000)
Saimiri_data.SRR11008272.log <- FindVariableFeatures(Saimiri_data.SRR11008272.log, selection.method = "vst", nfeatures = 2000)
SRR11008272.Saimiri.scaled <- ScaleData(Saimiri_data.SRR11008272.log)
SRR11008272.Saimiri.scaled[["individual"]] <- "SRR11008272"
SRR11008272.Saimiri.scaled[["species"]] <- "Saimiri"

Saimiri_data.SRR11008273.log <- NormalizeData(Saimiri_data.SRR11008273.log, normalization.method = "LogNormalize", scale.factor = 10000)
Saimiri_data.SRR11008273.log <- FindVariableFeatures(Saimiri_data.SRR11008273.log, selection.method = "vst", nfeatures = 2000)
SRR11008273.Saimiri.scaled <- ScaleData(Saimiri_data.SRR11008273.log)
SRR11008273.Saimiri.scaled[["individual"]] <- "SRR11008273"
SRR11008273.Saimiri.scaled[["species"]] <- "Saimiri"

Saimiri_data.SRR11008276.log <- NormalizeData(Saimiri_data.SRR11008276.log, normalization.method = "LogNormalize", scale.factor = 10000)
Saimiri_data.SRR11008276.log <- FindVariableFeatures(Saimiri_data.SRR11008276.log, selection.method = "vst", nfeatures = 2000)
SRR11008276.Saimiri.scaled <- ScaleData(Saimiri_data.SRR11008276.log)
SRR11008276.Saimiri.scaled[["individual"]] <- "SRR11008276"
SRR11008276.Saimiri.scaled[["species"]] <- "Saimiri"

#Batch correction: canonical correlation analysis (CCA) using Seurat

# The first piece of code will identify variable genes that are highly variable in at least 2/4 datasets. We will use these variable genes in our batch correction.
ob.list.log <- list(SRR11008270.Aotus.scaled, SRR11008271.Aotus.scaled, SRR11008274.Aotus.scaled, SRR11008275.Aotus.scaled, SRR11008277.Aotus.scaled, SRR11008278.Aotus.scaled, SRR11008269.Saimiri.scaled, SRR11008272.Saimiri.scaled, SRR11008273.Saimiri.scaled, SRR11008276.Saimiri.scaled)


# Identify anchors on the datasets, commonly shared variable genes across samples, 
# and integrate samples.
anchors.log <- FindIntegrationAnchors(object.list = ob.list.log, anchor.features = 2000, dims = 1:30)

combined.log <- IntegrateData(anchorset = anchors.log, dims = 1:30)

#set default assay to being "correct" values
DefaultAssay(combined.log) <- "integrated"

# Run the standard workflow for visualization and clustering.
# The integrated data object only stores the commonly shared variable genes.
combined.log <- ScaleData(combined.log, do.center = T, do.scale = F)

combined.log <- RunPCA(combined.log, npcs = 40, ndims.print = 1:5, nfeatures.print = 5)

pdf(file = "Plots/PCA_batchCorrected_individual_log.pdf")
DimPlot(combined.log, dims = c(1, 2), reduction = "pca", split.by = "individual") + ggtitle("LogNormalization")
dev.off()

pdf(file = "Plots/PCA_batchCorrected_species_log.pdf")
DimPlot(combined.log, dims = c(1, 2), reduction = "pca", split.by = "species") + ggtitle("LogNormalization")
dev.off()

# Clustering. Choose the dimensional reduction type to use and the number of aligned 
# canonical correlation vectors to use.
combined.log <- FindNeighbors(combined.log, reduction = "pca", dims = 1:20, k.param = 20)

combined.log <- FindClusters(combined.log, resolution = 0.8, algorithm = 1, random.seed = 100)

# UMAP. Choose the dimensional reduction type to use and the number of aligned 
# canonical correlation vectors to use.
combined.log  <- RunUMAP(combined.log , dims = 1:30, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

# After data integration, use the original expression data in all visualization and DE tests.
# The integrated data must not be used in DE tests as it violates assumptions of independence in DE tests!
DefaultAssay(combined.log) <- "RNA"  

# Visualize the Louvain clustering and the batches on the UMAP. 
# Remember, the clustering is stored in @meta.data in column seurat_clusters 
# and the technology is stored in the column tech. Remember you can also use DimPlot.
pdf(file="Plots/UMAP_batchCorrected_individual_species_log.pdf", width=900/72, height=478/72)
p1 <- DimPlot(combined.log, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("LogNormalization")
p2 <- DimPlot(combined.log, reduction = "umap", group.by = "individual", label = FALSE)
p3 <- DimPlot(combined.log, reduction = "umap", group.by = "species", label = TRUE)
p1 + p2 + p3 
dev.off()

save.image(file="batch_corrected_log.Rdata")
