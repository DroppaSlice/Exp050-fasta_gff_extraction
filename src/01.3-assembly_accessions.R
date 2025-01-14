library(readr)
library(rentrez)

#-----read and parse accessions-----#
nuc_accessions <- read.table(file = "outputs/unique_hits.txt")
nuc_accessions <- nuc_accessions$V1

#-----script to retrieve assembly using entrez-----#
get_acc <- function(x){
  require(rentrez)
  search_res <- entrez_search(db = "assembly", term = x)
  summ <- entrez_summary(db = "assembly", id = search_res$ids)
  if(length(summ) == 55){res <- summ$assemblyaccession}
  else if(length(summ) != 55){
    v <- vector(mode = "character", length = length(summ))
    for(i in 1:length(summ)){v[i] <- summ[[i]]$assemblyaccession}
    v <- sort(v, decreasing = T)
    res <- v[1]
  }
  return(res)
}
#testing
test_acc <- get_acc(nuc_accessions[1])

#-----running on full list-----#
#initialize an empty list
ass_accessions <- vector(mode = "list", length = length(nuc_accessions))
#loop through vector of nucleotide accessions
for(i in 1:length(nuc_accessions)){
  print(paste0("Processing entry ", i, " of ", length(nuc_accessions)))
  ass_accessions[i] <- get_acc(nuc_accessions[i])
}
#convert results fo a vector
accessions_v <- unlist(ass_accessions) %>%
  unique()

#-----writing accession numbers to file-----#
write(x = accessions_v, file = "outputs/assembly_accessions.txt")