#!/usr/bin/perl
use FindBin;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Getopt::Long;
use Pod::Usage; 
use Cwd;
 
my %conversion = (
    "cluster01" => "Eg08",
    "cluster02" => "Eg01",
    "cluster03" => "Eg02",
    "cluster04" => "Eg07",
    "cluster05" => "Eg03",
    "cluster06" => "Eg09",
    "cluster07" => "Eg05",
    "cluster08" => "Eg04",
    "cluster09" => "Eg06"
);

my $help        = "";
my $gff         = "";
my $prefix      = "ensete_glaucum";
my $fasta       = "";
my $dir         = cwd();
my $verbose     = 0;
my $get_by_name = 0;
my $current_dir = getcwd; 
#Macma304_g029910.1

my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters  
    --gff        gff3 file  [required] 
    --fasta      fasta file [required]
    --prefix     prefix for output file 
    --verbose    
    --help       show this help and exit
/;
GetOptions(
    'gff=s'    => \$gff,
    'fasta=s'  => \$fasta,
    'prefix=s' => \$prefix, 
    'verbose'  => \$verbose,    
    'help|h|?' => \$help
)  or pod2usage(-message => $usage);
if ($help) { pod2usage(-message => $usage); } 
  
if ($gff eq "") {
    print "\nWarn :: --gff is empty. Please specify the GFF file \n\n";
    print $usage;
    exit 0;
}
if ($fasta eq "") {
    print "Warn :: --fasta is empty. Please specify the FASTA file\n";
    print $usage;
    exit 0;
}
 

 
my $logfile = $dir ."/gff2fasta.log";
open(LOG,">$logfile"); 

my $output_gene    = $dir."/".$prefix ."_gene.fna";
my $output_cds     = $dir."/".$prefix ."_cds.fna";
my $output_cdna    = $dir."/".$prefix ."_cdna.fna";
my $output_protein = $dir."/".$prefix ."_protein.faa";
my $output_gff3    = $dir."/".$prefix .".gff3";


# Input file
my $fasta_in = new Bio::SeqIO(
	-file => $fasta,
	-format => 'fasta'
);
my $gff_in = new Bio::Tools::GFF(
	-file => $gff,
	-gff_version => 3
); 
# Output file
###### Output type description ######
# gene 			- the entire gene sequence (including UTRs and introns)
# cds 			- translated sequence (starting with ATG and ending with a stop codon included)
# cdna 			- transcribed sequence (devoid of introns, but containing untranslated exons)
# protein 		- cds translated (includes a * as the stop codon)

my $gene_out = new Bio::SeqIO(
	-file => ">$output_gene",
	-format => 'fasta'
);
my $cds_out = new Bio::SeqIO(
	-file => ">$output_cds",
	-format => 'fasta'
);
my $cdna_out = new Bio::SeqIO(
	-file => ">$output_cdna",
	-format => 'fasta'
);
my $protein_out = new Bio::SeqIO(
	-file => ">$output_protein",
	-format => 'fasta'
);
my $gff_out = new Bio::Tools::GFF(
	-file => ">$output_gff3",
	-gff_version => 3
); 

# Read FASTA
 
my %fasta;
&print("# Read reference fasta file\n");
while (my $seqobj = $fasta_in->next_seq) {
	$fasta{$seqobj->display_id} = $seqobj;
}
$fasta_in->close;


# Read & parse GFF3
my %gene;
my %mrna;
my %exon;
my %cds;
my $cpt; 

&print("# Read gff3 file\n");

while(my $feature = $gff_in->next_feature) {
    my $seq_id = $feature->seq_id();
    $feature->seq_id($conversion{$seq_id});
	if ($feature->primary_tag() eq "gene") {
        my ($gene_id) = $feature->get_tag_values("Name"); 
        push @{$gene{$conversion{$seq_id}}}, $feature; 
    }
	elsif ($feature->primary_tag() eq "mRNA") {
        my ($mrna_id) = $feature->get_tag_values("ID");
        my ($mrna_name_id) = $feature->get_tag_values("Name");
        my ($gene_id) = $feature->get_tag_values("Parent"); 
        $mrna{$gene_id}{$mrna_id} = $feature;
	}
	elsif ($feature->primary_tag() eq "CDS") {
		my ($mrna_id) = $feature->get_tag_values("Parent");  
		push @{$cds{$mrna_id}} , $feature;
		 
	} 
	elsif ($feature->primary_tag() eq "exon") {
		my ($mrna_id) = $feature->get_tag_values("Parent"); 
		push @{$exon{$mrna_id}} , $feature;
	}
} 
$gff_in->close; 
&print("# Write fasta file\n");
my $cpt_gene    = 0;
my $cpt_cds     = 0;
my $start_codon = 0;
my @start;
my $stop_codon  = 0;
my @stop;
my @bad_sequence;



foreach my $seq_id (sort {$a cmp $b} keys %gene) {
    my $cpt = 0;
    my $new_seq_id = $seq_id;
    $new_seq_id =~ s/Eg/chr/;
	my $seqobj = $fasta{$new_seq_id};
	my $strand;
    foreach my $feature (sort {$a->start <=> $b->start} @{$gene{$seq_id}}) {
        $cpt++;
        $cpt_gene++;
        $feature->seq_id($new_seq_id);
        my ($gene_id) = $feature->get_tag_values("ID");
        my $new_gene_id = sprintf( "%s_g%06d", $seq_id, $cpt * 10 );
        my $new_mrna_id = sprintf( "%s_t%06d", $seq_id, $cpt * 10 );
        $feature->remove_tag("ID");
        $feature->remove_tag("Name");
        $feature->add_tag_value("ID",$new_gene_id);
        $feature->add_tag_value("Name",$new_gene_id);
        $gff_out->write_feature($feature);
		my $seqobj_gene = $seqobj->trunc($feature->start,$feature->end);
		my $strand = $feature->strand;
        if ($strand =~ /-/) {
            $seqobj_gene = $seqobj_gene->revcom();
        }	 
		$seqobj_gene->display_id($new_gene_id); 
		$gene_out->write_seq($seqobj_gene);
        foreach my $mrna_id (keys %{$mrna{$gene_id}}) {
            $cpt_cds++;
            my $feature_mrna = $mrna{$gene_id}{$mrna_id};
            $feature_mrna->seq_id($new_seq_id); 
            $feature_mrna->remove_tag("ID");
            $feature_mrna->remove_tag("Name");
            $feature_mrna->remove_tag("Parent");
            $feature_mrna->add_tag_value("ID",$new_mrna_id);
            $feature_mrna->add_tag_value("Name",$new_mrna_id);
            $feature_mrna->add_tag_value("Parent",$new_gene_id);
            $gff_out->write_feature($feature_mrna);
            my $cds_seq;
            my $cdna_seq;
            foreach my $cds (sort {$a->start <=>$b->start} @{$cds{$mrna_id}} ) {
                $cds->seq_id($new_seq_id);
                $cds->remove_tag("ID") if $cds->has_tag("ID");
                $cds->remove_tag("Parent");
                $cds->add_tag_value("Parent",$new_mrna_id);
                $gff_out->write_feature($cds);
                $cds_seq .= $seqobj->subseq($cds->start,$cds->end); 
                $strand = $cds->strand;
                my $subseq = $seqobj->subseq($cds->start,$cds->end);
            }
            foreach my $exon (sort {$a->start <=>$b->start} @{$exon{$mrna_id}} ) {
                $exon->seq_id($new_seq_id);
                $exon->remove_tag("ID") if $exon->has_tag("ID");
                $exon->remove_tag("Parent");
                $exon->add_tag_value("Parent",$new_mrna_id);
                $gff_out->write_feature($exon);	
                $cdna_seq .= $seqobj->subseq($exon->start,$exon->end); 
            }
            my $seqobj_cds = Bio::PrimarySeq->new(
                -seq        => $cds_seq, 
                -display_id => $new_mrna_id
            ); 
            if ($strand =~ /-/) {
                $seqobj_cds = $seqobj_cds->revcom();
            }
            if ($seqobj_cds->seq =~ /^ATG.*/){
                $start_codon++;
                push @start, $mrna_id; 
            } 
            if ($seqobj_cds->seq =~ /.*(TAG|TAA|TGA)$/){
                $stop_codon++;
                push @stop , $mrna_id; 
            }
            
            my $seqobj_protein4 = $seqobj_cds->translate();
            my $seqobj_protein1 = $seqobj_cds->translate(-offset=>1);
            my $seqobj_protein2 = $seqobj_cds->translate(-offset=>2);
            my $seqobj_protein3 = $seqobj_cds->translate(-offset=>3);
            $seqobj_protein4->display_id($prot_id);
            $seqobj_protein1->display_id($prot_id);
            $seqobj_protein2->display_id($prot_id);
            $seqobj_protein3->display_id($prot_id);
            my @seq_prot4 = (split(/\*/,$seqobj_protein4->seq()));
            my @seq_prot1 = (split(/\*/,$seqobj_protein1->seq()));
            my @seq_prot2 = (split(/\*/,$seqobj_protein2->seq()));
            my @seq_prot3 = (split(/\*/,$seqobj_protein3->seq()));
            my $seq_protein;
            if (scalar(@seq_prot4) <= 1) {
                $seq_protein = $seqobj_protein4;
            }
            if (scalar(@seq_prot1) <= 1) {
                $seq_protein = $seqobj_protein1;
            }
            elsif (scalar(@seq_prot2) <=1 ) {
    
                $seq_protein = $seqobj_protein2;
            }
            elsif (scalar(@seq_prot3) <= 1) {
    
                $seq_protein = $seqobj_protein3;
            }
            if ($seq_protein) {
                $cds_out->write_seq($seqobj_cds);
                $protein_out->write_seq($seq_protein);
                my $cdna_seq;
                foreach my $feature (sort{$a->start <=> $b->start} @{$exon{$seq_id}{$mrna_id}}) {
                    $cdna_seq .= $seqobj->subseq($feature->start,$feature->end);
                }
                my $seqobj_cdna = Bio::PrimarySeq->new(
                    -seq        => $cdna_seq,
                    -display_id => $mrna_id
                );
                if ($strand =~ /-/) {
                    $seqobj_cdna = $seqobj_cdna->revcom();
                }
                $cdna_out->write_seq($seqobj_cdna);
            }
        } 
	}  
}
$gene_out->close;
$cds_out->close;
$cdna_out->close;
$protein_out->close;
&print("Output file:\n");
&print("- $output_gene\n");
&print("- $output_cds\n");
&print("- $output_cdna\n");
&print("- $output_protein\n");
&print("- $output_gff3\n\n");
&print("Number of gene(s) : ". $cpt_gene ."\n");
&print("Number of transcript : " . $cpt_cds ." (with ATG ". $start_codon ."; with stop_codon ". $stop_codon .")\n"); 

 

sub print {
    my $cmd = shift;
    print LOG $cmd;
    print $cmd  if $verbose; 
}

 