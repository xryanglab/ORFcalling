prepare_transcripts -g gencode.v19.annotation.gtf -f hg19_genome.fa -o ribocode_annot
RiboCode -a ribocode_annot -c config.txt  -o ribocode