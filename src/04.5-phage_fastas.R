library(Biostrings)
library(dplyr)
library(stringr)

#-----sourcing our functions from 04.1-----#
source(file = "src/04.1-defining_functions.R")

#-----read phage blast data-----#
phage_blast <- read.csv(file = "inputs/ncbi_phage_blast.csv") %>%
  mutate(qseqid = as.factor(qseqid)) %>%
  #manually add phage names to match gff list
  mutate(phage.name = c("Enterobacteria_phage_YYZ_2008", 
                        "Stx2a_WGPS6",
                        "Stx2a_WGPS8",
                        "Stx2_1717",
                        "Enterobacteria_phage_YYZ_2008",
                        "Enterobacteria_phage_2851",
                        "Escherichia_virus_2B8",
                        "Escherichia_virus_2B8")) %>%
  mutate(names = paste0(qseqid, "_", phage.name, "_10kFlanks"))

#-----reading in phage dna sequences-----#
phage_seqs <- readDNAStringSet(filepath = "inputs/phage_genomes/ncbi_pas_bacteriophages.fasta")
names(phage_seqs) <- str_extract(names(phage_seqs), pattern = "^[:graph:]+")

#-----running the seq_extract() function-----#
for(pas in levels(phage_blast$qseqid)){
  phage_blast_filtered <- phage_blast %>% 
    filter(qseqid == pas)
  res <- seq_extract(input = phage_blast_filtered, seqs = phage_seqs, upstream = 5000, downstream = 5000)
  file.name <- paste0(pas, "_phage_5kFlanks.fasta")
  writeXStringSet(x = res, filepath = paste0("outputs/fastas/", file.name))
}