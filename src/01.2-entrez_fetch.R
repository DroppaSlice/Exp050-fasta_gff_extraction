library(rentrez)
library(readr)
library(dplyr)
library(seqinr)
library(stringr)

#-----reading in blast results-----#
blast_hits <- read_tsv(file = "outputs/prok_complete_genomes_blast.txt",
                       col_names = c("qseqid", "sseqid", "sacc", "sallacc", "pident", "length", "qcovs", "qstart", "qend", "sstart", "send", "sstrand", "evalue", "bitscore")) %>%
#initial filtering step to remove partial hits with <75% coverage of sRNA 
   filter(qcovs >=75)

#-----pas_filter() function to filter pasA and pasD homologs-----#
pas_filter <- function(blast_results, id_col = "qseqid", a.limit = 95, d.limit = 98){
  list <- split(blast_results, f = blast_results[,id_col])
  homologs <- c("pasA1", "pasA2", "pasD1", "pasD2")
  for(gene in homologs){
    #filter pasA1 and pasA2 dfs for >= 95% sequence identity hits
    if(gene == "pasA1" | gene == "pasA2"){
      list[[gene]] <- dplyr::filter(list[[gene]], pident >= 95)
    }
    #filter pasD1 and pasD2 dfs for >= 98% sequence identity hits
    if(gene == "pasD1" | gene == "pasS2"){
      list[[gene]] <- dplyr::filter(list[[gene]], pident >= 98)
    }
  }
  df <- bind_rows(list)
  return(df)
}

#-----filtering of blast results-----#
blast_filtered <- pas_filter(blast_hits) %>%
  select(sacc) %>%
  unique() 
#creating a vector of unique NCBI accessions
unique_acc <- blast_filtered$sacc
write(unique_acc, file = "outputs/ncbi_unique_accessions.txt")

#-----splitting accessions-----#
#rentrez has a pretty small query limit, so we split the accessions into chunks of 50 for processing
acc_split <- split(unique_acc, ceiling(seq_along(unique_acc) / 50))

#-----retrieving sequences-----#
#initialize a list 'seqs' to hold our sequences
seqs <- vector(mode = "list", length = length(unique_acc))
#create an iterator k to keep track of our place in the dataset
k <- 1

#-----loop through the list of accession numbers-----#
#first loop through elements of the list
for(i in 1:length(acc_split)){
  print(paste0("Processing section ", i, " of ", length(acc_split)))
  #extract list element as a vector
  acc_v <- acc_split[[i]]
  #next loop through elements of the vector
  for(j in 1:length(acc_v)){
    print(paste0("Processing accession ", j, " from section ", i))
    seqs[[k]] <- entrez_fetch(db = "nucleotide",
                            id = acc_v[[j]],
                            rettype = "fasta")
    #advance our iterator variable
    k <- k + 1
  }
}

#-----cleaning up sequences-----#
#grabbing the names
get_names <- function(x){
  name <- str_extract(x, pattern = ">.*\\n")
  name <- str_remove(name, pattern = "\\n")
  name <- str_remove(name, pattern = ">")
}
seq_names <- lapply(seqs, FUN = get_names)
#removing the names from the seq list
remove_fasta_name <- function(x){
  y <- str_remove(x, pattern = ">.*\\n")
  return(y)
}
seqs_cleaned <- lapply(seqs, remove_fasta_name)

#-----writing to fasta-----#
write.fasta(sequences = seqs_cleaned, names = seq_names, file.out = "outputs/pas_blast_genomes/ncbi_pas_genomes.fasta")