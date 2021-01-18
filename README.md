
## Download dataset

Assembly (FASTA file)

```bash
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/430 -O ensete_glaucum.assembly.fna
```

Structural annotation provide by EVM (GFF3 file)

```bash
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/670 -O EVM.gff
```

Functional annotation (tab file)

```bash
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/672 -O ensete_glaucum_product.txt
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/671 -O ensete_glaucum_evm_ipr.txt
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/673 -O ensete_glaucum_evm_go.txt
wget https://banana-genome-hub.southgreen.fr/filebrowser/download/675 -O interpro.tar.gz
```

## gff2fasta.pl

Convert EvidenceModeler output to fasta for each gene model and rename ID according to the order of appearance in the assembly  

```bash
perl gff2fasta.pl --fasta ensete_glaucum.assembly.fna --gff EVM.gff  --verbose
```
Produce the following file :
 - ensete_glaucum_gene.fna 
 - ensete_glaucum_cds.fna
 - ensete_glaucum_cdna.fna 
 - ensete_glaucum_protein.faa
 - ensete_glaucum.gff3

## Run InterProScan (v 5.41-78.0)

```bash
perl split_fasta.pl --file ensete_glaucum_protein.faa
sed -i "s:*::" fasta*fna
mkdir interpro
for i in $PWD/fasta/*fna; do echo qsub -b y -q normal.q -l mem_free=20G -N iprscan  -V interproscan.sh -i $i -cpu 8 -dp --iprlookup --pathways --goterms -d $PWD/interpro >> iprscan.sh;done
chmod +x iprscan.sh
module load bioinfo/interproscan/5.41-78.0
./iprscan.sh
```

## Parse InterProScan

```bash
tar -xzf interpro.tar.gz
perl cnv_interpro.pl --result interpro --prefix ensete_glaucum
```

Produce the following file :
 - ensete_glaucum_ipr.txt
 - ensete_glaucum_go.txt

## Run ncbi-blast (v2.9)

```bash
blastp -query ensete_glaucum_protein.faa -out blast_SwissProt.out -db <path swissprot db> -evalue 1e-10 -max_target_seqs 5
blastp -query ensete_glaucum_protein.faa -out blast_TrEMBL.out -db <path trembl db>  -evalue 1e-10 -max_target_seqs 5
blastp -query ensete_glaucum_protein.faa -out blast_MUSAC.out -db <path dh pahang db> -evalue 1e-10 -max_target_seqs 5
```

## cnv_blast.pl

Parse multiple Blastp and choose the best product for each gene

```bash
perl cnv_blast.pl --blast blast_SwissProt.out --blast blast_TrEMBL.out --blast blast_MUSAC.out --output ensete_glaucum_product.txt
```

Produce the following file :
 - ensete_glaucum_product.txt

## add_functional_annotation2gff3.pl

```bash
perl add_functional_annotation2gff3.pl -product ensete_glaucum_product.txt -go_file ensete_glaucum_evm_go.txt -interpro_file ensete_glaucum_evm_ipr.txt -gff3_file ensete_glaucum.gff3 -prefix ensete_glaucum
```
