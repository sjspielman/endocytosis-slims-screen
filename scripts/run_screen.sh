#!/bin/bash
set -euo pipefail

cd $(dirname "${BASH_SOURCE[0]}")

# Detect motifs
#echo "Screening proteome..."
#./screen-disordered-proteome.py

# Join with gene names
echo "Joining with gene names..."
./join-gene-names.R

# clean up
#rm temp.csv
