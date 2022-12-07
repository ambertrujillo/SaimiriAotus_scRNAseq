load("batch_corrected_log.Rdata")

# Find Hemoglobin Genes of interest (GOI
#look for gene names
unique_orthologs[str_detect(unique_orthologs$inputaotus, "HBB"), ] #find Ensemb name
#ENSSBOG00000029907
rownames(combined.log)[grep(pattern = "*ENSSBOG00000029907", x = rownames(x = combined.log), ignore.case = FALSE)] #look for it in the Seurat object

unique_orthologs[str_detect(unique_orthologs$inputaotus, "HBZ"), ] #find Ensemb name
#ENSSBOG00000022580
rownames(combined.log)[grep(pattern = "*ENSSBOG00000022580", x = rownames(x = combined.log), ignore.case = FALSE)] #look for it in the Seurat object

unique_orthologs[str_detect(unique_orthologs$inputaotus, "HBA"), ] #find Ensemb name
unique_orthologs[str_detect(unique_orthologs$inputaotus, "HBA"), ] #find Ensemb name
#Not present 

unique_orthologs[str_detect(unique_orthologs$inputaotus, "ENSSBOG00000025868"), ] #HBA #find Ensemb name
#ENSSBOG00000025868
rownames(combined.log)[grep(pattern = "*ENSSBOG00000025868", x = rownames(x = combined.log), ignore.case = FALSE)] #look for it in the Seurat object

unique_orthologs[str_detect(unique_orthologs$input_ensg, "ENSSBOG00000029161"), ] #ACVR2B #find Ensemb name
unique_orthologs[str_detect(unique_orthologs$input_ensg, "ENSSBOG00000033030"), ] #HEMGN #find Ensemb name


hemoglobin.beta <- c("ENSSBOG00000029907")
hemoglobin.alpha <- c("ENSSBOG00000025868", "ENSSBOG00000022580")
ACVR2B <- c("ENSSBOG00000029161")
HEMGN <- c("ENSSBOG00000033030")

hemoglobin.beta.df = data.frame(hemoglobin.beta)
colnames(hemoglobin.beta.df)[1] = "genes"
hemoglobin.alpha.df = data.frame(hemoglobin.alpha)
colnames(hemoglobin.alpha.df)[1] = "genes"
ACVR2B.df = data.frame(ACVR2B)
colnames(ACVR2B.df)[1] = "genes"
HEMGN.df = data.frame(HEMGN)
colnames(HEMGN.df)[1] = "genes"



hemoglobin <- full_join(hemoglobin.beta.df, hemoglobin.alpha.df, by="genes")
hemoglobin <- as.vector(hemoglobin)

pdf(file="Plots/hemoglobin_MarkerGenes.pdf", width=811/72, height=478/72)
DotPlot(combined.log, features = hemoglobin, cols = c("Spectral"), scale = TRUE) + RotatedAxis() +theme(axis.text.x=element_text(size=10)) + ggtitle("Hemoglobin genes")
dev.off()

# Feature Plots
hemoglobin.GOI <- c("ENSSBOG00000029907", "ENSSBOG00000025868") # No HBZ

pdf(file="Plots/hemoglobin_GOI.pdf")
FeaturePlot(combined.log, features = hemoglobin.GOI)
dev.off()


##### Split out Aotus and Saimiri Clusters #####

# Aotus
DefaultAssay(combined.log) <- "RNA" # Need to use uncorrected values for DE analysis
Idents(combined.log) <- "species"

markers.Aotus <- FindMarkers(combined.log, ident.1 = "Aotus")

write.table(markers.Aotus, file="Aotus/Aotus_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

Aotus <- subset(combined.log, idents = "Aotus")

# Recluster Aotus
DefaultAssay(Aotus) <- "integrated" # Need to use uncorrected values for DE analysis

Aotus.recluster <- FindNeighbors(Aotus, dims = 1:10) 
Aotus.recluster <- FindClusters(Aotus.recluster, resolution = 0.5) 

pdf(file="Plots/Aotus/UMAP_batchCorrected_Aotus_cluster_individual_log.pdf", width=900/72, height=478/72)
p1 <- DimPlot(Aotus.recluster, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(Aotus.recluster, reduction = "umap", group.by = "individual", label = FALSE)
p1 + p2 
dev.off()

pdf(file = "Plots/Aotus/Violin_plot_nCounts_clusters.pdf", compress=TRUE)
VlnPlot(Aotus.recluster, features = c("nCount_RNA", "nFeature_RNA"), ncol=2)
dev.off()

table(Idents(Aotus.recluster))

Idents(Aotus.recluster) <- "species"

pdf(file = "Plots/Aotus/Violin_plot_nCounts_individual.pdf", compress=TRUE)
VlnPlot(Aotus.recluster, features = c("nCount_RNA", "nFeature_RNA"), ncol=2)
dev.off()

Idents(Aotus.recluster) <- "seurat_clusters"

pdf(file="Plots/Aotus/hemoglobin_MarkerGenes.pdf", width=811/72, height=478/72)
DotPlot(Aotus.recluster, features = hemoglobin.GOI, cols = c("Spectral"), scale = TRUE) + RotatedAxis() +theme(axis.text.x=element_text(size=10)) + ggtitle("Hemoglobin")
dev.off()

pdf(file="Plots/Aotus/hemoglobin_GOI.pdf")
FeaturePlot(Aotus.recluster, features = hemoglobin.GOI)
dev.off()

# Saimiri
DefaultAssay(combined.log) <- "RNA" # Need to use uncorrected values for DE analysis
Idents(combined.log) <- "species"

markers.Saimiri <- FindMarkers(combined.log, ident.1 = "Saimiri")

write.table(markers.Saimiri, file="Saimiri/Saimiri_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

Saimiri <- subset(combined.log, idents = "Saimiri")

# Recluster Aotus
DefaultAssay(Saimiri) <- "integrated" # Need to use uncorrected values for DE analysis

Saimiri.recluster <- FindNeighbors(Saimiri, dims = 1:10) 
Saimiri.recluster <- FindClusters(Saimiri.recluster, resolution = 0.5) 

pdf(file="Plots/Saimiri/UMAP_batchCorrected_Saimiri_cluster_individual_log.pdf", width=900/72, height=478/72)
p1 <- DimPlot(Saimiri.recluster, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(Saimiri.recluster, reduction = "umap", group.by = "individual", label = FALSE)
p1 + p2 
dev.off()

pdf(file = "Plots/Saimiri/Violin_plot_nCounts_clusters.pdf", compress=TRUE)
VlnPlot(Saimiri.recluster, features = c("nCount_RNA", "nFeature_RNA"), ncol=2)
dev.off()

table(Idents(Saimiri.recluster))

pdf(file="Plots/Saimiri/hemoglobin_MarkerGenes.pdf", width=811/72, height=478/72)
DotPlot(Saimiri.recluster, features = hemoglobin.GOI, cols = c("Spectral"), scale = TRUE) + RotatedAxis() +theme(axis.text.x=element_text(size=10)) + ggtitle("Hemoglobin")
dev.off()

pdf(file="Plots/Saimiri/hemoglobin_GOI.pdf")
FeaturePlot(Saimiri.recluster, features = hemoglobin.GOI)
dev.off()

##########################################################################
###################### Red Blood cells #################################
##########################################################################

##### Combine RBC clusters for Aotus (8 and 13) and Saimiri (2 and 3) #####

#Extract Aotus cluster 8 and 13
DefaultAssay(Aotus.recluster) <- "integrated" # Need to use uncorrected values for DE analysis

Idents(Aotus.recluster) <- "seurat_clusters"

Aotus_cluster8_13 <- subset(Aotus.recluster, idents = c(8, 13))

#Extract Siamiri RBC clusters, 2 and 3

DefaultAssay(Saimiri.recluster) <- "integrated" # Need to use uncorrected values for DE analysis

Idents(Saimiri.recluster) <- "seurat_clusters"

Saimiri_cluster2_3 <- subset(Saimiri.recluster, idents = c(2, 3))

#Merge them
RBC.combined <- merge(Aotus_cluster8_13, y = Saimiri_cluster2_3)

#Rescale after merging: "You do not need to renormalize your data because LogNormalize only log transform your data on a sample level, and each sample is treated independently. The values wont change after merging. 
#On the other hand, you need to re-scale your data because the scaling transformation takes into account the content (samples) of your dataset. The values will be different after merging." 

#set default assay to being "correct" values
DefaultAssay(RBC.combined) <- "RNA"

# Run the standard workflow for visualization and clustering.
# The integrated data object only stores the commonly shared variable genes.

RBC.combined <- FindVariableFeatures(RBC.combined, selection.method = "vst", nfeatures = 2000)

RBC.combined <- ScaleData(RBC.combined, do.center = T, do.scale = F)

RBC.combined <- RunPCA(RBC.combined, npcs = 40, ndims.print = 1:5, nfeatures.print = 5)

RBC.combined <- FindNeighbors(RBC.combined, reduction = "pca", dims = 1:20, k.param = 20)

RBC.combined <- FindClusters(RBC.combined, resolution = 0.8, algorithm = 1, random.seed = 100)

# UMAP. Choose the dimensional reduction type to use and the number of aligned 
# canonical correlation vectors to use.
RBC.combined  <- RunUMAP(RBC.combined, dims = 1:30, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

save.image(file="/scratch/aet359/scRNAseq/combined_analysis/Seurat/RBC_clustering.Rdata")

pdf(file="/scratch/aet359/scRNAseq/combined_analysis/Seurat/Plots/RBC/Clusters2_3_8_13/UMAP_batchCorrected_clusters_species_individual.pdf", width=900/72, height=478/72)
p1 <- DimPlot(RBC.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(RBC.combined, reduction = "umap", group.by = "species", label = TRUE) + xlab("UMAP 1") + ylab("UMAP2")
p3 <- DimPlot(RBC.combined, reduction = "umap", group.by = "individual", label = TRUE)
p1 + p2 + p3
dev.off()

pdf(file = "/scratch/aet359/scRNAseq/combined_analysis/Seurat/Plots/RBC/Clusters2_3_8_13/Violin_plot_nCounts_nFeature.pdf", compress=TRUE)
VlnPlot(RBC.combined, features = c("nCount_RNA", "nFeature_RNA"), ncol=2)
dev.off()

pdf(file="/scratch/aet359/scRNAseq/combined_analysis/Seurat/Plots/RBC/Clusters2_3_8_13/hemoglobin_GOI.pdf", width=811/72, height=478/72)
FeaturePlot(RBC.combined, features = hemoglobin.GOI)
dev.off()

pdf(file="/scratch/aet359/scRNAseq/combined_analysis/Seurat/Plots/RBC/Clusters2_3_8_13/parasite_load_heatplot.pdf")
FeaturePlot(RBC.combined, features = "parasite_load")
dev.off()

distrib_parasite_load_RBC <- data.frame(RBC.combined[[]])

pdf(file="/scratch/aet359/scRNAseq/combined_analysis/Seurat/Plots/RBC/Clusters2_3_8_13/histogram_parasite_load_all_cells.pdf")
ggplot(distrib_parasite_load_RBC, aes(x=log(parasite_load))) + geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#B3FF66") + labs(x="Parasite load (by cell)") + facet_wrap(~species) + theme_classic()
dev.off()

Idents(RBC.combined) <- "species"
table(Idents(RBC.combined))
#Aotus Saimiri 
#   468     456

distrib_parasite_load_RBC$cells = rownames(distrib_parasite_load_RBC)

pdf(file="/scratch/aet359/scRNAseq/combined_analysis/Seurat/Plots/RBC/Clusters2_3_8_13/barplot_parasite_load_all_cells.pdf")
ggplot(distrib_parasite_load_RBC, aes(cells, parasite_load, color=species)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_blank())
dev.off()

Idents(RBC.combined) <- "seurat_clusters"

All_markers <- FindAllMarkers(object = RBC.combined)

write.table(All_markers, file="cluster_marker_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
