import subprocess
import sys

cmd1 = "datasets download genome accession --inputfile assembly_accessions.txt --filename ncbi_pas_blast_hits.zip --include genome,gff3"
subprocess.run(cmd1, shell = True)