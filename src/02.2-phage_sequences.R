library(Biostrings)
library(dplyr)
library(stringr)

#-----read phage blast hits-----#
phage_blast <- read.csv(file = "inputs/ncbi_phage_blast.csv") %>%
  mutate(qseqid = as.factor(qseqid))

#----read Xstringset of phage genomes-----#
phage_genomes <- readDNAStringSet(filepath = "inputs/phage_genomes/ncbi_pas_bacteriophages.fasta")
names(phage_genomes) <- str_extract(names(phage_genomes), pattern = "[:graph:]+")

#-----loading and running seq_extract()-----#
#loading from source file
source(file = "src/seq_extract.R")

#running on phage dataset for each pas gene
for(pas in levels(phage_blast$qseqid)){
  res <- seq_extract(input = phage_blast, seqs = phage_genomes, upstream = 10000, downstream = 10000)
  file.name <- paste0(pas, "_phage_10kFlanks.fasta")
  writeXStringSet(x = res, filepath = paste0("outputs/fastas/", file.name))
}
