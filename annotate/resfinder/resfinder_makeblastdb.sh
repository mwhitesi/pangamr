#!/usr/bin/env bash
cat *.fsa > resfinder.fasta
makeblastdb -dbtype nucl -in resfinder.fasta -title "Resfinder" -out resfinder 
