########################  Analysis of the data ############################



# Read the Adjancy Matrix -------------------------------------------------

Adjmat   <- read.table("Data/significant_modules_ldak.txt",
              header=T)


# Read the mapping file ---------------------------------------------------

map_file <- read.table("Data/snp_2_gene.map")

# Get the membership sets from the adjmat ---------------------------------

membership_Sets <- get_membership(Adjmat)

# Get Snps ----------------------------------------------------------------

snp_Lists <- Get_snps(membership_Sets,gene_snp_map)



# Save the files ----------------------------------------------------------
for (i in 1:length(snp_Lists)){
  write.table(snp_Lists[[i]],file=paste0(path,"/",i,".txt"),
              col.names=F,
              row.names=F,
              quote=F)
}



