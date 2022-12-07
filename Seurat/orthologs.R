# orthologous reciprocal gene lists for saimiri (N=19757 for all below matrices) and aotus (N=21456 for all below matrices)
#Saimiri orthologs
all_aotus_genes = rownames(Aotus_data.SRR11008270)
all_aotus_genes = str_replace_all(all_aotus_genes,"Aotus_", "")

aotus_to_saimiri = gorth(
  all_aotus_genes,
  source_organism = "anancymaae",
  target_organism = "sbboliviensis",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

aotus_to_saimiri = subset(aotus_to_saimiri, select=c(input, input_ensg, ortholog_ensg))
# 18174     2
unique_saimiri_orthologs = aotus_to_saimiri[!duplicated(aotus_to_saimiri$ortholog_ensg),] # doesnt matter which one you do first (ortholog_ensg or input_ensg) you always end up with the same number of genes at the end
# 17528     2
unique_saimiri_orthologs = unique_saimiri_orthologs[!duplicated(unique_saimiri_orthologs$input_ensg),]
# 16964     2

#Aotus orthologs
all_saimiri_genes = rownames(Saimiri_data.SRR11008269)
all_saimiri_genes = str_replace_all(all_saimiri_genes,"Saimiri_", "")

saimiri_to_aotus = gorth(
  all_saimiri_genes,
  source_organism = "sbboliviensis",
  target_organism = "anancymaae",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

saimiri_to_aotus = subset(saimiri_to_aotus, select=c(input, input_ensg, ortholog_ensg))
# 17926     2
unique_aotus_orthologs = saimiri_to_aotus[!duplicated(saimiri_to_aotus$ortholog_ensg),]
# 17250     2
unique_aotus_orthologs = unique_aotus_orthologs[!duplicated(unique_aotus_orthologs$input_ensg),]
# 16755     2

# combine ortholog lists
unique_orthologs = merge(unique_aotus_orthologs, unique_saimiri_orthologs, 
    by.x="input_ensg", by.y="ortholog_ensg", suffix=c("aotus", "saimiri"))

unique_orthologs = unique_orthologs[unique_orthologs$ortholog_ensg == unique_orthologs$input_ensgsaimiri, ] 

write.table(unique_orthologs,
            file=paste0("/scratch/aet359/scRNAseq/combined_analysis/unique_orthologs", suffix=".txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE)

# Aotus - using "unique_saimiri_orthologs" this list started with Aotus genes as the input

##########SRR11008270
SRR11008270.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/aotus/Aotus_cell_ranger_70/outs/filtered_feature_bc_matrix")

Aotus_subset.SRR11008270 <- rownames(SRR11008270.data)[grep("^Aotus", rownames(SRR11008270.data))]
head(Aotus_subset.SRR11008270, 10)
Aotus_data.SRR11008270 <- SRR11008270.data[Aotus_subset.SRR11008270, ]

# changing rownames to orthologs
Aotus_data.SRR11008270 = data.frame(Aotus_data.SRR11008270)
Aotus_data.SRR11008270$inputsaimiri = rownames(Aotus_data.SRR11008270)
Aotus_data.SRR11008270$inputsaimiri = str_replace_all(Aotus_data.SRR11008270$inputsaimiri, "Aotus_", "")
Aotus_data.SRR11008270 = full_join(Aotus_data.SRR11008270, unique_orthologs, by="inputsaimiri")
Aotus_data.SRR11008270 = drop_na(Aotus_data.SRR11008270)

rownames(Aotus_data.SRR11008270) = Aotus_data.SRR11008270$input_ensg

Aotus_data.SRR11008270 = subset(Aotus_data.SRR11008270, select=-c(inputaotus, input_ensg, ortholog_ensg, inputsaimiri, input_ensgsaimiri))

##########SRR11008271
SRR11008271.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/aotus/Aotus_cell_ranger_71/outs/filtered_feature_bc_matrix")

Aotus_subset.SRR11008271 <- rownames(SRR11008271.data)[grep("^Aotus", rownames(SRR11008271.data))]
head(Aotus_subset.SRR11008271, 10)
Aotus_data.SRR11008271 <- SRR11008271.data[Aotus_subset.SRR11008271, ]

# changing rownames to orthologs
Aotus_data.SRR11008271 = data.frame(Aotus_data.SRR11008271)
Aotus_data.SRR11008271$inputsaimiri = rownames(Aotus_data.SRR11008271)
Aotus_data.SRR11008271$inputsaimiri = str_replace_all(Aotus_data.SRR11008271$inputsaimiri, "Aotus_", "")
Aotus_data.SRR11008271 = full_join(Aotus_data.SRR11008271, unique_orthologs, by="inputsaimiri")
Aotus_data.SRR11008271 = drop_na(Aotus_data.SRR11008271)

rownames(Aotus_data.SRR11008271) = Aotus_data.SRR11008271$input_ensg

Aotus_data.SRR11008271 = subset(Aotus_data.SRR11008271, select=-c(inputaotus, input_ensg, ortholog_ensg, inputsaimiri, input_ensgsaimiri))


##########SRR11008274
SRR11008274.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/aotus/Aotus_cell_ranger_74/outs/filtered_feature_bc_matrix")

Aotus_subset.SRR11008274 <- rownames(SRR11008274.data)[grep("^Aotus", rownames(SRR11008274.data))]
head(Aotus_subset.SRR11008274, 10)
Aotus_data.SRR11008274 <- SRR11008274.data[Aotus_subset.SRR11008274, ]

# changing rownames to orthologs
Aotus_data.SRR11008274 = data.frame(Aotus_data.SRR11008274)
Aotus_data.SRR11008274$inputsaimiri = rownames(Aotus_data.SRR11008274)
Aotus_data.SRR11008274$inputsaimiri = str_replace_all(Aotus_data.SRR11008274$inputsaimiri, "Aotus_", "")
Aotus_data.SRR11008274 = full_join(Aotus_data.SRR11008274, unique_orthologs, by="inputsaimiri")
Aotus_data.SRR11008274 = drop_na(Aotus_data.SRR11008274)

rownames(Aotus_data.SRR11008274) = Aotus_data.SRR11008274$input_ensg

Aotus_data.SRR11008274 = subset(Aotus_data.SRR11008274, select=-c(inputaotus, input_ensg, ortholog_ensg, inputsaimiri, input_ensgsaimiri))

##########SRR11008275
SRR11008275.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/aotus/Aotus_cell_ranger_75/outs/filtered_feature_bc_matrix")

Aotus_subset.SRR11008275 <- rownames(SRR11008275.data)[grep("^Aotus", rownames(SRR11008275.data))]
head(Aotus_subset.SRR11008275, 10)
Aotus_data.SRR11008275 <- SRR11008275.data[Aotus_subset.SRR11008275, ]

# changing rownames to orthologs
Aotus_data.SRR11008275 = data.frame(Aotus_data.SRR11008275)
Aotus_data.SRR11008275$inputsaimiri = rownames(Aotus_data.SRR11008275)
Aotus_data.SRR11008275$inputsaimiri = str_replace_all(Aotus_data.SRR11008275$inputsaimiri, "Aotus_", "")
Aotus_data.SRR11008275 = full_join(Aotus_data.SRR11008275, unique_orthologs, by="inputsaimiri")
Aotus_data.SRR11008275 = drop_na(Aotus_data.SRR11008275)

rownames(Aotus_data.SRR11008275) = Aotus_data.SRR11008275$input_ensg

Aotus_data.SRR11008275 = subset(Aotus_data.SRR11008275, select=-c(inputaotus, input_ensg, ortholog_ensg, inputsaimiri, input_ensgsaimiri))

##########SRR11008277
SRR11008277.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/aotus/Aotus_cell_ranger_77/outs/filtered_feature_bc_matrix")

Aotus_subset.SRR11008277 <- rownames(SRR11008277.data)[grep("^Aotus", rownames(SRR11008277.data))]
head(Aotus_subset.SRR11008277, 10)
Aotus_data.SRR11008277 <- SRR11008277.data[Aotus_subset.SRR11008277, ]

# changing rownames to orthologs
Aotus_data.SRR11008277 = data.frame(Aotus_data.SRR11008277)
Aotus_data.SRR11008277$inputsaimiri = rownames(Aotus_data.SRR11008277)
Aotus_data.SRR11008277$inputsaimiri = str_replace_all(Aotus_data.SRR11008277$inputsaimiri, "Aotus_", "")
Aotus_data.SRR11008277 = full_join(Aotus_data.SRR11008277, unique_orthologs, by="inputsaimiri")
Aotus_data.SRR11008277 = drop_na(Aotus_data.SRR11008277)

rownames(Aotus_data.SRR11008277) = Aotus_data.SRR11008277$input_ensg

Aotus_data.SRR11008277 = subset(Aotus_data.SRR11008277, select=-c(inputaotus, input_ensg, ortholog_ensg, inputsaimiri, input_ensgsaimiri))

##########SRR11008278
SRR11008278.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/aotus/Aotus_cell_ranger_78/outs/filtered_feature_bc_matrix")

Aotus_subset.SRR11008278 <- rownames(SRR11008278.data)[grep("^Aotus", rownames(SRR11008278.data))]
head(Aotus_subset.SRR11008278, 10)
Aotus_data.SRR11008278 <- SRR11008278.data[Aotus_subset.SRR11008278, ]

# changing rownames to orthologs
Aotus_data.SRR11008278 = data.frame(Aotus_data.SRR11008278)
Aotus_data.SRR11008278$inputsaimiri = rownames(Aotus_data.SRR11008278)
Aotus_data.SRR11008278$inputsaimiri = str_replace_all(Aotus_data.SRR11008278$inputsaimiri, "Aotus_", "")
Aotus_data.SRR11008278 = full_join(Aotus_data.SRR11008278, unique_orthologs, by="inputsaimiri")
Aotus_data.SRR11008278 = drop_na(Aotus_data.SRR11008278)

rownames(Aotus_data.SRR11008278) = Aotus_data.SRR11008278$input_ensg

Aotus_data.SRR11008278 = subset(Aotus_data.SRR11008278, select=-c(inputaotus, input_ensg, ortholog_ensg, inputsaimiri, input_ensgsaimiri))

# Saimiri - using "unique_aotus_orthologs" this list started with Saimiri genes as the input

###########SRR11008269
SRR11008269.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/saimiri/Saimiri_cell_ranger_69/outs/filtered_feature_bc_matrix")

Saimiri_subset.SRR11008269 <- rownames(SRR11008269.data)[grep("^Saimiri", rownames(SRR11008269.data))]
head(Saimiri_subset.SRR11008269, 10)
Saimiri_data.SRR11008269 <- SRR11008269.data[Saimiri_subset.SRR11008269, ]


# changing rownames to orthologs
Saimiri_data.SRR11008269 = data.frame(Saimiri_data.SRR11008269)
Saimiri_data.SRR11008269$inputaotus = rownames(Saimiri_data.SRR11008269)
Saimiri_data.SRR11008269$inputaotus = str_replace_all(Saimiri_data.SRR11008269$inputaotus, "Saimiri_", "")
Saimiri_data.SRR11008269 = full_join(Saimiri_data.SRR11008269, unique_orthologs, by="inputaotus")
Saimiri_data.SRR11008269 = drop_na(Saimiri_data.SRR11008269)

rownames(Saimiri_data.SRR11008269) = Saimiri_data.SRR11008269$input_ensg

Saimiri_data.SRR11008269 = subset(Saimiri_data.SRR11008269, select=-c(inputsaimiri, input_ensg, inputaotus, ortholog_ensg, input_ensgsaimiri))

###########SRR11008272
SRR11008272.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/saimiri/Saimiri_cell_ranger_72/outs/filtered_feature_bc_matrix")

Saimiri_subset.SRR11008272 <- rownames(SRR11008272.data)[grep("^Saimiri", rownames(SRR11008272.data))]
head(Saimiri_subset.SRR11008272, 10)
Saimiri_data.SRR11008272 <- SRR11008272.data[Saimiri_subset.SRR11008272, ]

# changing rownames to orthologs
Saimiri_data.SRR11008272 = data.frame(Saimiri_data.SRR11008272)
Saimiri_data.SRR11008272$inputaotus = rownames(Saimiri_data.SRR11008272)
Saimiri_data.SRR11008272$inputaotus = str_replace_all(Saimiri_data.SRR11008272$inputaotus, "Saimiri_", "")
Saimiri_data.SRR11008272 = full_join(Saimiri_data.SRR11008272, unique_orthologs, by="inputaotus")
Saimiri_data.SRR11008272 = drop_na(Saimiri_data.SRR11008272)

rownames(Saimiri_data.SRR11008272) = Saimiri_data.SRR11008272$input_ensg

Saimiri_data.SRR11008272 = subset(Saimiri_data.SRR11008272, select=-c(inputsaimiri, input_ensg, inputaotus, ortholog_ensg, input_ensgsaimiri))

###########SRR11008273
SRR11008273.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/saimiri/Saimiri_cell_ranger_73/outs/filtered_feature_bc_matrix")

Saimiri_subset.SRR11008273 <- rownames(SRR11008273.data)[grep("^Saimiri", rownames(SRR11008273.data))]
head(Saimiri_subset.SRR11008273, 10)
Saimiri_data.SRR11008273 <- SRR11008273.data[Saimiri_subset.SRR11008273, ]

# changing rownames to orthologs
Saimiri_data.SRR11008273 = data.frame(Saimiri_data.SRR11008273)
Saimiri_data.SRR11008273$inputaotus = rownames(Saimiri_data.SRR11008273)
Saimiri_data.SRR11008273$inputaotus = str_replace_all(Saimiri_data.SRR11008273$inputaotus, "Saimiri_", "")
Saimiri_data.SRR11008273 = full_join(Saimiri_data.SRR11008273, unique_orthologs, by="inputaotus")
Saimiri_data.SRR11008273 = drop_na(Saimiri_data.SRR11008273)

rownames(Saimiri_data.SRR11008273) = Saimiri_data.SRR11008273$input_ensg

Saimiri_data.SRR11008273 = subset(Saimiri_data.SRR11008273, select=-c(inputsaimiri, input_ensg, inputaotus, ortholog_ensg, input_ensgsaimiri))

###########SRR11008276
SRR11008276.data <- Read10X(data.dir = "/scratch/aet359/scRNAseq/cellranger/saimiri/Saimiri_cell_ranger_76/outs/filtered_feature_bc_matrix")

Saimiri_subset.SRR11008276 <- rownames(SRR11008276.data)[grep("^Saimiri", rownames(SRR11008276.data))]
head(Saimiri_subset.SRR11008276, 10)
Saimiri_data.SRR11008276 <- SRR11008276.data[Saimiri_subset.SRR11008276, ]

# changing rownames to orthologs
Saimiri_data.SRR11008276 = data.frame(Saimiri_data.SRR11008276)
Saimiri_data.SRR11008276$inputaotus = rownames(Saimiri_data.SRR11008276)
Saimiri_data.SRR11008276$inputaotus = str_replace_all(Saimiri_data.SRR11008276$inputaotus, "Saimiri_", "")
Saimiri_data.SRR11008276 = full_join(Saimiri_data.SRR11008276, unique_orthologs, by="inputaotus")
Saimiri_data.SRR11008276 = drop_na(Saimiri_data.SRR11008276)

rownames(Saimiri_data.SRR11008276) = Saimiri_data.SRR11008276$input_ensg

Saimiri_data.SRR11008276 = subset(Saimiri_data.SRR11008276, select=-c(inputsaimiri, input_ensg, inputaotus, ortholog_ensg, input_ensgsaimiri))

# confirm that rownames of saimiri dataset match aotus dataset
# Aotus matches
table(rownames(Aotus_data.SRR11008270) %in% rownames(Aotus_data.SRR11008271))
table(rownames(Aotus_data.SRR11008270) %in% rownames(Aotus_data.SRR11008274))
table(rownames(Aotus_data.SRR11008270) %in% rownames(Aotus_data.SRR11008275))
table(rownames(Aotus_data.SRR11008270) %in% rownames(Aotus_data.SRR11008277))
table(rownames(Aotus_data.SRR11008270) %in% rownames(Aotus_data.SRR11008278))
# Saimiri matches
table(rownames(Saimiri_data.SRR11008269) %in% rownames(Saimiri_data.SRR11008272))
table(rownames(Saimiri_data.SRR11008269) %in% rownames(Saimiri_data.SRR11008273))
table(rownames(Saimiri_data.SRR11008269) %in% rownames(Saimiri_data.SRR11008276))
# Aotus matches Saimiri
table(rownames(Aotus_data.SRR11008270) %in% rownames(Saimiri_data.SRR11008269))

dim(Aotus_data.SRR11008270)
dim(Aotus_data.SRR11008271)
dim(Aotus_data.SRR11008274)
dim(Aotus_data.SRR11008275)
dim(Aotus_data.SRR11008277)
dim(Aotus_data.SRR11008278)
dim(Saimiri_data.SRR11008269)
dim(Saimiri_data.SRR11008272)
dim(Saimiri_data.SRR11008273)
dim(Saimiri_data.SRR11008276)

save.image(file = "Orthologous_matrices.Rdata")
