#!/bin/bash

GENOME=$1
THREADS=64
if [[ ! -e "$GENOME" ]]; then
  echo "USAGE: $0 GENOME.fasta"
  echo "(Input file $GENOME not found)"
  exit 1
fi
set -euxo pipefail

mkdir -p genome
echo "# Cleanup"
funannotate-docker clean -i $GENOME  --minlen 1000 -o genome/1.clean.fa
echo "# Sort contigs"
funannotate-docker  sort -i genome/1.clean.fa -b scaffold -o genome/2.sorted.fa
echo "# Mask contigs"
funannotate-docker mask -i genome/2.sorted.fa --cpus $THREADS -o genome/3.mask.fa
echo "# Augustus gene prediction"
funannotate-docker predict -i genome/3.mask.fa -o genome/fun --species "Candida parapsilosis" --strain Y9 --busco_seed_species candida_albicans --cpus $THREADS
echo "# Funannotation"
funannotate-docker annotate -i fun/  --cpus $THREADS
gzip genome/*.fa
