library(ape)
library(microseq)

#-----read gff3 using read.gff() from the 'ape' package-----#
pasA1_phage_gff <- read.gff(file = "inputs/phage_genomes/Enterobacteria_phage_YYZ_2008.gff3") %>%
  mutate(type = as.factor(type))

#-----read blast data-----#
phage_blast <- read.csv(file = "inputs/ncbi_phage_blast.csv") %>%
  mutate(qseqid = as.factor(qseqid))

#-----test gff extraction manually-----#
lower <- 13622
upper <- 33705
pasA1_phage_gff_10k <- pasA1_phage_gff %>%
  filter(type == "CDS") %>%
  filter(!(end < lower | start > upper))
  
writeGFF(pasA1_phage_gff_10k, out.file = "outputs/gffs/YYZ_2008_pasA1_10k.gff3")