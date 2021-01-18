#!/bin/bash

# Download dataset

# Assembly
wget --quiet https://banana-genome-hub.southgreen.fr/filebrowser/download/430 -O ensete_glaucum.assembly.fna

# Structural annotation provide by EVM (GFF3 file)

wget --quiet https://banana-genome-hub.southgreen.fr/filebrowser/download/675 -O EVM.gff

#Functional annotation (tab file)
wget --quiet https://banana-genome-hub.southgreen.fr/filebrowser/download/672 -O ensete_glaucum_product.txt 
wget --quiet https://banana-genome-hub.southgreen.fr/filebrowser/download/674 -O interpro.tar.gz

perl gff2fasta.pl --fasta ensete_glaucum.assembly.fna --gff EVM.gff --prefix ensete_glaucum --verbose

tar -xzf  interpro.tar.gz
perl cnv_interpro.pl --result tsv --prefix ensete_glaucum

perl add_functional_annotation2gff3.pl -product ensete_glaucum_product.txt -go_file ensete_glaucum_go.txt -interpro_file ensete_glaucum_ipr.txt -gff3_file ensete_glaucum.gff3 -prefix ensete_glaucum
