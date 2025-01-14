import subprocess
import sys

cmd1 = "/Users/zhudx/ncbi-blast/bin/blastn -db prok_complete_genomes -query inputs/pas_sRNAs.fasta -remote -outfmt \"6 qseqid sseqid sacc sallacc pident length qcovs qstart qend sstart send sstrand evalue bitscore \" -max_target_seqs 5000 -out outputs/prok_complete_genomes_blast.txt"
subprocess.run(cmd1, shell = True)