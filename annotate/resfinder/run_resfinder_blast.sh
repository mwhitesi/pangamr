#!/usr/bin/env bash
QUERY=$1
blastn -query $QUERY -db resfinder -out resfinder_report.tsv -evalue 0.0001 -outfmt 6 -perc_identity 90 
