library(Biostrings)
library(dplyr)
library(stringr)

#-----reading in data-----#
#reading genomes fastas as DNAstringset
dna_seqs <- readDNAStringSet(filepath = "outputs/pas_blast_genomes/ncbi_pas_genomes.fasta")
#reading metadata from blast search
blast_metadata <- read.csv(file = "outputs/prok_complete_genomes_blast_homologFiltered.csv") %>%
  select(qseqid, sacc, sstart, send, sstrand, length)

#-----creating an input dataframe for just pasA1 hits-----#
pasA1_df <- blast_metadata %>%
  filter(qseqid == "pasA1") %>%
  filter(length >= 0.75*(max(length)))

names(dna_seqs) <- str_extract(names(dna_seqs), pattern = "^[^\\.]+")

#-----defining new function seq_extract()-----#
seq_extract <- function(input, seqs, upstream = 0, downstream = 0){
  require(Biostrings)
  #initialize vectors to hold our results and names
  res_list <- vector(mode = "list", length = nrow(input))
  names_v <- vector(mode = "character", length = nrow(input))
  #iterate through the input data frame using a for loop
  for(i in 1:nrow(input)){
    #define the row that is being operated on in each step of the loop
    row <- input[i,]
    #define a few of the variables stored in that row
    accession <- getElement(row, "sacc")
    qseqid <- getElement(row, "qseqid")
    sstart <- getElement(row, "sstart")
    send <- getElement(row, "send")
    strand <- getElement(row, "sstrand")
    #define the specific DNA string based on the matched accession
    dna <- seqs[[accession]]
    upper <- length(dna)
    #define a unique name identifier for the extracted region
    seq_name <- paste0(accession, "_", qseqid, "_", sstart, "-up", upstream,"-down", downstream)
    #-----error handing-----#
    if(is.null(dna)){
      print("Error: accession not found in DNA list")
      stop()
    }
    #extracting sequences using subseq()
    if(strand == "plus"){
      #check ends within bounds, skip if out of bounds
      if(sstart - upstream <= 0 | send + downstream > upper){
        next()
      }
      res_list[[i]] <- subseq(dna, start = (sstart - upstream), end = (send + downstream))
    }
    else if(strand == "minus"){
      #check ends within bounds, skip if out of bounds
      if(send - upstream <=0 | sstart + downstream > upper){
        next()
      }
      seq <- subseq(dna, start = (send - upstream), end = (sstart + downstream))
      res_list[[i]] <- reverseComplement(seq)
    }
    #append name to the names vector
    names_v[i] <- seq_name
  }
  names(res_list) <- names_v
  res_list <- res_list[!sapply(res_list, is.null)]
  #convert list into a more memory-efficient DNAStringSet object
  res_set <- DNAStringSet(res_list)
  return(res_set)
}

#-----running the seq_extract for pasA1-----#
#just the pasA1 hit
pasA1_set <- seq_extract(input = pasA1_df, seqs = dna_seqs)
#1kbp flanks
pasA1_1k <- seq_extract(input = pasA1_df, seqs = dna_seqs, upstream = 1000, downstream = 1000)
#5kbp flanks
pasA1_5k <- seq_extract(input = pasA1_df, seqs = dna_seqs, upstream = 5000, downstream = 5000)
#10kbp flanks
pasA1_10k <- seq_extract(input = pasA1_df, seqs = dna_seqs, upstream = 10000, downstream = 10000)

#-----writing fastas-----#
writeXStringSet(x = pasA1_set, filepath = "outputs/fastas/pasA1_ncbi_hits.fasta")
writeXStringSet(x = pasA1_1k, filepath = "outputs/fastas/pasA1_ncbi_hits_1kbpFlanks.fasta")
writeXStringSet(x = pasA1_5k, filepath = "outputs/fastas/pasA1_ncbi_hits_5kbpFlanks.fasta")
writeXStringSet(x = pasA1_10k, filepath = "outputs/fastas/pasA1_ncbi_hits_10kbpFlanks.fasta")
#individual flanks for pasA1 10k flanks
for(i in 1:length(pasA1_10k)){
  print(paste0("writing file ", i, " out of ", length(pasA1_10k)))
  writeXStringSet(x = pasA1_10k[i], filepath = paste0("outputs/fastas/pasA1/", names(pasA1_10k[i]), ".fasta"))
}
