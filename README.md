
## Download dataset

Assembly (FASTA file)

```get
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/430 -O ensete_glaucum.assembly.fna
```

Structural annotation provide by EVM (GFF3 file)

```get
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/670 -O EVM.gff
```

Functional annotation (tab file)

```get
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/672 -O ensete_glaucum_product.txt
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/671 -O ensete_glaucum_evm_ipr.txt
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/673 -O ensete_glaucum_evm_go.txt
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/674 -O interpro.tar.gz
```

## gff2fasta.pl

Convert EvidenceModeler output to fasta for each gene model and rename ID according to the order of appearance in the assembly  

```perl
perl gff2fasta.pl --fasta ensete_glaucum.assembly.fna --gff EVM.gff  --verbose
```
Produce the following file :
 - ensete_glaucum_gene.fna 
 - ensete_glaucum_cds.fna
 - ensete_glaucum_cdna.fna 
 - ensete_glaucum_protein.faa
 - ensete_glaucum.gff3

## Parse interproscan result

```perl
tar -xzf interpro.tar.gz
perl cnv_interpro.pl --result interpro --prefix ensete_glaucum
```

## Run ncbi-blast

 - blast_SwissProt.out
 - blast_TrEMBL.out
 - blast_MUSAC.out


## cnv_blast.pl

Parse multiple Blastp and choose the best product for each gene

```perl
perl cnv_blast.pl --blast blast_SwissProt.out --blast blast_TrEMBL.out --blast blast_MUSAC.out --output ensete_glaucum_product.txt
```

## add_functional_annotation2gff3.pl

```perl
perl add_functional_annotation2gff3.pl -product ensete_glaucum_product.txt -go_file ensete_glaucum_evm_go.txt -interpro_file ensete_glaucum_evm_ipr.txt -gff3_file ensete_glaucum.gff3 -prefix ensete_glaucum
```
