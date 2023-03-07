# R script to map Ensembl gene IDs to gene names for the output of PrediXcan.
# 6th March 2023

args <- commandArgs(trailingOnly = TRUE)

# library(biomaRt)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

filename = args[1]
gene_file = read.table(filename,header = TRUE)

# First 2 columns are FID and IID, take the rest which are the genes
gene_ids = colnames(gene_file)[3:length(colnames(gene_file))]

# Extract only the gene ID, i.e. remove the version number
gene_no_version = sub("(.*)[\\.].*","\\1",gene_ids)


# musc_gene_names_no_version = getBM(attributes = c("ensembl_gene_id","external_gene_name","hgnc_symbol"),filters = c("ensembl_gene_id"),values = musc_gene_no_version,mart = ensembl)

# Map the gene IDs to their names using the given mapping
gene_map = read.table("ensemble2gene-tab.map",header = FALSE)
gene_map_unique = unique(gene_map[,c(4,5)])
rownames(gene_map_unique) = gene_map_unique[,1]
mapped_genes = gene_map_unique[gene_no_version,2]
colnames(gene_file)[3:length(colnames(gene_file))] = mapped_genes

# Extract the columns which have been matched to a name, drop the rest and write to file
final_df = gene_file[,c(colnames(gene_file)[1:2],na.omit(mapped_genes))]
filename_noext= sub("(.*)[\\.].*","\\1",filename)
write.table(final_df,file = paste(filename_noext,"_gene_mapped.tsv",sep=""),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)


