library(Biostrings)
library(dplyr)
library(stringr)

#-----reading in data-----#
#reading genomes fastas as DNAstringset
dna_seqs <- readDNAStringSet(filepath = "outputs/pas_blast_genomes/ncbi_pas_genomes.fasta")
#reading metadata from blast search
blast_metadata <- read.csv(file = "outputs/prok_complete_genomes_blast_homologFiltered.csv") %>%
  select(qseqid, sacc, sstart, send, sstrand, length) %>%
  mutate(qseqid = as.factor(qseqid))

#-----sourcing our functions from 04.1-----#
source(file = "src/04.1-defining_functions.R")

#-----creating input dataframes for each of the pas hits-----#
pasA2_df <- blast_metadata %>%
  filter(qseqid == "pasA2") %>%
  filter(length >= 0.75*max(length))

#-----running scripts for pasA2-----#
pasA2_10k <- seq_extract(input = pasA2_df, seqs = dna_seqs, upstream = 10000, downstream = 10000)
writeXStringSet(pasA2_10k, filepath = "outputs/fastas/pasA2_ncbi_hits_10kbpFlanks.fasta")
write_fastas(stringset = pasA2_10k, outfilepath = "outputs/fastas/pasA2/")

#-----iterating scripts over remaining pas genes-----#
pas_genes <- c("pasA1", "pasA2", "pasB", "pasC", "pasD1", "pasD2")

for(gene in pas_genes){
  print(paste0("Processing ", gene))
  input_df <- blast_metadata %>%
    filter(qseqid == gene) %>%
    filter(length >= 0.75*max(length))
  write.csv(input_df, file = paste0("outputs/blast_tables/", gene, "_blast_metadata.csv"))
  
  df_10k <- seq_extract(input = input_df, seqs = dna_seqs, upstream = 10000, downstream = 10000)
  writeXStringSet(df_10k, filepath = paste0("outputs/fastas/", gene, "_ncbi_hits_10kbpFlanks.fasta"))
  write_fastas(stringset = df_10k, outfilepath = paste0("outputs/fastas/", gene, "/"), verbose = F)
}