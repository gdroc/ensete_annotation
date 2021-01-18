qsub -b y -q normal.q -N rename perl bin/gff2fasta.pl --fasta /gs7k1/projects/HUBs/banana/Reference_genomes_assemblies/Ensete_glaucum_1.0/ensete_glaucum.assembly.fna --gff /work/droc/ensete/EVM_gff1.all.gff --prefix ensete_glaucum --dir $PWD --verbose


qsub -b y -q normal.q -N test perl bin/cnv_blast.pl --blast /work/droc/ensete/blast_SwissProt.out --blast /work/droc/ensete/blast_TrEMBL.out --blast /work/droc/ensete/blast_DHPahang.out --prefix ensete_glaucum --dir $PWD


