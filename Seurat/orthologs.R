# Unzip files from outs/raw_feature_bc_matrix/ directory

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
