library(ape)
library(tidyr)
library(stringr)
library(dplyr)
library(microseq)

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

#-----reading in blast metadata table-----#
blast_metadata <- read.csv(file = "outputs/prok_complete_genomes_blast_homologFiltered.csv") %>%
  select(qseqid, sacc, sstart, send, sstrand, length)
#making another dataframe for just pasA1 hits
pasA1_df <- blast_metadata %>%
  filter(qseqid == "pasA1") %>%
  filter(length >= 0.75*(max(length)))
#write the pasA1 table
write.csv(pasA1_df, file = "outputs/blast_tables/pasA1_blast_metadata.csv")

#-----constructing ncbi accession lookup table-----#
convert_entrez <- function(input, target_db = "assembly"){
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
#initialize a storage vector
acc_v <- vector(mode = "character", length = nrow(pasA1_df))

#loop to grab all assembly accessions for pasA1 blast hits
for(i in 642:nrow(pasA1_df)){
  acc_v[i] <- convert_entrez(input = pasA1_df[i,"sacc"])
  print(paste0("Processing accession ", i, " out of ", nrow(pasA1_df)))
}
#append as new variable to pasA1 data frame
pasA1_df$assembly.accession <- acc_v

#-----function to access and extract from gff list using accession no.-----#
subset_gff <- function(input, gffs, 
                       upstream = 0, downstream = 0, 
                       n.acc.name = "sacc", a.acc.name = "assembly.accession",
                       strand.name = "sstrand"){
  missing_count <- 0
  n.acc <- input[,n.acc.name]
  a.acc <- input[,a.acc.name]
  strand <- input[,strand.name]
  start <- input[,"sstart"]
  end <- input[,"send"]
  print(paste0("working on accession ", a.acc, " molecule ", n.acc))
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

#-----running gff subset function for pasA1 hits-----#
#filter out entries that do not have a matching gff
pasA1_df_filtered <- pasA1_df %>%
  filter(assembly.accession %in% names(gff_list)) %>%
  mutate(names = paste0(sacc, "_", qseqid, "_", "sstart", "-up10000-down10000"))

#extracting gff list for 10k flanks
pasA1_gffs_10k <- lapply(X = split(pasA1_df_filtered, f = seq(nrow(pasA1_df_filtered))), FUN = subset_gff, gffs = gff_list, upstream = 10000, downstream = 10000)
names(pasA1_gffs_10k) <- pasA1_df_filtered$names

#-----writing out gff files-----#
for(i in 1:length(pasA1_gffs_10k)){
  print(paste0("writing file ", i, " out of ", length(pasA1_gffs_10k)))
  writeGFF(pasA1_gffs_10k[[i]], paste0("outputs/gffs/pasA1/", names(pasA1_gffs_10k[i]), ".gff3"))
}
