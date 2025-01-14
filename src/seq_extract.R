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