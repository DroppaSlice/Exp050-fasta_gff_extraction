library(ape)
library(tidyr)
library(stringr)
library(dplyr)
library(microseq)
library(rentrez)

#-----finding and listing gffs recursively-----#
file_names <- list.files(path = "inputs/ncbi-dataset/data/", pattern = "*.gff", recursive = T)
#reading gffs into an R list object
pathname <- "inputs/ncbi-dataset/data/"
gff_list <- vector(mode = "list", length = length(file_names))
for(i in 1:length(file_names)){
  print(paste0("Reading in file ", file_names[i]))
  gff_list[[i]] <- ape::read.gff(file = paste0(pathname, file_names[i]))
}
names(gff_list) <- str_extract(file_names, "^[^/]+")

#-----sourcing our functions from 04.1-----#
source(file = "src/04.1-defining_functions.R")

#-----iterating gff extraction over all pas_genes-----#
gene_list <- c("pasA1", "pasA2", "pasB", "pasC", "pasD1", "pasD2")
 
#-----for loop to create accession number conversion tables-----#
for(gene in gene_list){
  print(paste0("Processing gffs for ", gene))
  #read input metadata
  input_df <- read.csv(file = paste0("outputs/blast_tables/", gene, "_blast_metadata.csv"))
  
  #create output directories
  if(!dir.exists(paste0("outputs/gffs/", gene))) {dir.create(paste0("outputs/gffs/", gene))}
  
  acc_v <- vector(mode = "character", length = nrow(input_df))
  for( i in 1:nrow(input_df)){
    print(paste0("processing sample ", i, " out of ", nrow(input_df)))
    acc_v[i] <- convert_entrez(input = input_df[i, "sacc"])
  }
  input_df$assembly.accession <- acc_v
  write.csv(input_df, file = paste0("outputs/blast_tables/", gene, "_accession_conversion.csv"))
}  

#-----separate for loop to run the subset_gff() function-----#

for(gene in gene_list){
  print(paste0("Processing ", gene))
  input_df_filtered <- read.csv(paste0("outputs/blast_tables/", gene, "_accession_conversion.csv")) %>%
    filter(assembly.accession %in% names(gff_list)) %>%
    mutate(names = paste0(sacc, "_", qseqid, "_", "sstart", "-up10000-down10000"))
  
  gff_10k <- lapply(X = split(input_df_filtered, f = seq(nrow(input_df_filtered))), FUN = subset_gff, gffs = gff_list, upstream = 10000, downstream = 10000)
  names(gff_10k) <- input_df_filtered$names
  
  write.gff.list(gffs = gff_10k, outfilepath = paste0("outputs/gffs/", gene, "/"))
}