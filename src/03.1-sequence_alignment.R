library(DECIPHER)
library(Biostrings)

pasA1_10k <- readDNAStringSet(filepath = "outputs/fastas/pasA1_ncbi_hits_10kbpFlanks.fasta")

db <- dbConnect(SQLite(), ":memory:")

Seqs2DB(seqs = pasA1_10k[1:10],
        type = "XStringSet",
        dbFile = db, 
        identifier = names(pasA1_10k[1:10]))

synteny <- FindSynteny(db, verbose = T)
pasA1_10k_syn <- AlignSynteny(synteny, db, verbose = T)