library(ape)
library(microseq)

#-----reading in phage gff files-----#
file_names <- list.files(path = "inputs/phage_genomes/", pattern = "*.gff3", recursive = T)
pathname <- "inputs/phage_genomes/"
gff_list <- vector(mode = "list", length = length(file_names))
for ( i in 1:length(file_names)){
  print(paste0("Reading in file ", file_names[i]))
  gff_list[[i]] <- ape::read.gff(file = paste0(pathname, file_names[i]))
}
names(gff_list) <- str_extract(file_names, "^[^/.]+")

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

#-----lapply subset_gffs() to extract phage gffs-----#
phage_gff_10k <- lapply(X = split(phage_blast, f = seq(nrow(phage_blast))), 
                        FUN = subset_gff,
                        gffs = gff_list,
                        upstream = 10000,
                        downstream = 10000,
                        a.acc.name = "phage.name")
names(phage_gff_10k) <- phage_blast$names

#-----write out gffs-----#
write.gff.list(gffs = phage_gff_10k, outfilepath = "outputs/gffs/phages/")