#-----seq_extract() for subsetting genome sequence based on regions-----#
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
        print(paste0("Region for sequence ", qseqid, "_", accession, " is out of bounds"))
        next()
      }
      res_list[[i]] <- subseq(dna, start = (sstart - upstream), end = (send + downstream))
    }
    else if(strand == "minus"){
      #check ends within bounds, skip if out of bounds
      if(send - upstream <=0 | sstart + downstream > upper){
        print(paste0("Region for sequence ", qseqid, "_", accession, " is out of bounds"))
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

#-----write_fastas() function to write sequences in list to individual fastas
write_fastas <- function(stringset, outfilepath, verbose = T){
  for(i in 1:length(stringset)){
    if(verbose == T){print(paste0("writing file ", i, " out of ", length(stringset)))}
    writeXStringSet(x = stringset[i], filepath = paste0(outfilepath, names(stringset[i]), ".fasta"))
  }
} 

#-----convert_entrez() converts nucleotide accession numbers to assembly accessions-----#
convert_entrez <- function(input, target_db = "assembly", verbose = T){
  require(rentrez)
  e_search <- entrez_search(db = target_db, term = input)
  e_summ <- entrez_summary(db = target_db, id = e_search$ids)
  if(length(e_summ) == 55){
    res <- e_summ$assemblyaccession
  }
  else if(length(e_summ) != 55){
    v <- vector(mode = "character", length = length(e_summ))
    for(i in 1:length(e_summ)){v[i] <- e_summ[[i]]$assemblyaccession}
    v <- sort(v, decreasing = T)
    res <- v[1]
  }
  #added step to deal with accession numbers not found by entrez_search()
  if (is.null(res)){
    res = "accession not found"
  }
  return(res)
}

#-----subset_gff() extracts regions of a gff files based on inputs-----#
subset_gff <- function(input, gffs, 
                       upstream = 0, downstream = 0, 
                       n.acc.name = "sacc", a.acc.name = "assembly.accession",
                       strand.name = "sstrand",
                       verbose = T){
  missing_count <- 0
  n.acc <- input[,n.acc.name]
  a.acc <- input[,a.acc.name]
  strand <- input[,strand.name]
  start <- input[,"sstart"]
  end <- input[,"send"]
  if(verbose == T) {print(paste0("working on accession ", a.acc, " molecule ", n.acc))}
  #establishing region bounds
  if(strand == "plus"){
    lower = start - upstream
    upper = end + downstream
  }
  else if(strand == "minus"){
    lower = end - downstream
    upper = start + upstream
  }
  #accessing the gff table by assembly accession
  if(str_detect(a.acc, ",")){
    v_acc <- unlist(base::strsplit(a.acc, split = ","))
    print(paste0(length(v_acc), " accessions detected"))
    for(acc in v_acc){
      if(is.null(gffs[[acc]])){
        print(paste0(acc, " is not found in gff list, skipping."))
        next
      }
      else if(!is.null(gffs[[acc]])){
        print(paste0(acc, " found in list. Processing..."))
        df <- gffs[[acc]]
      }
    }
  }
  else if(str_detect(a.acc, ",") == F){
    df <- gffs[[a.acc]]
  }
  df <- df %>%
    #filtering based on nucleotide accession
    filter(str_detect(seqid, n.acc)) %>%
    #extracting features within the region
    filter(!(end < lower | start > upper))
  return(df)
}

#-----write.gff.list() function to write individual gffs from a list of named gff tables-----#
write.gff.list <- function(gffs, outfilepath){
  require(microseq)
  for(i in 1:length(gffs)){
    writeGFF(gffs[[i]], paste0(outfilepath, names(gffs[i]), ".gff3"))
  }
}