#!/usr/bin/perl  
use Getopt::Long;
use FindBin;
use Pod::Usage; 
use Bio::SeqIO;
use Bio::Tools::GFF;
  
my $help      = "";
my $file      = ""; 
my $gff3_file = "";
my $go_file   = "";
my $interpro_file = "";
my $product = "";
my $te_file = "";
my %ipr;
my %go;
my %locus_info;
my $prefix   = "ensete_glaucum";
my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters

    -gff3_file   gff3 file to be annotated  
    -prefix      prefix name for fasta output
    -product     locus tag file, required
    -go_file
    -interpro_file
    -help
/;
GetOptions(  
    'gff3_file=s'  => \$gff3_file,   
    'go_file=s'       => \$go_file,
    'interpro_file=s' => \$interpro_file,
    'product=s'          => \$product,
    'prefix=s'     => \$prefix, 
    'help|h|?'     => \$help
)  or pod2usage(-message => $message);
if ($help) { pod2usage(-message => $message); } 

if ($gff3_file eq "") {
    warn "\nWarn :: --file is empty. Please specify the locus tag file\n\n";
    warn $usage;
    exit 0;
}      
my $gff3_out_file = $prefix ."_annoted.gff3"; 
my $gff_in = new Bio::Tools::GFF(
    -file => $gff3_file,
    -gff_version => 3
);
my $gff_out = new Bio::Tools::GFF(
    -file => ">$gff3_out_file",
    -gff_version => 3
);
 

open(IN,$product); 
while (<IN>) {
    chomp;
    my ($primary_id,$product,$ic_code,$completeness,$gene_name,$dbxref) = (split(/\t/,$_)) ; 
    $completeness = "missing_functional_completeness" if $completeness eq "";
    $locus_info{$primary_id} = { 
        product      => $product,
        completeness => $completeness,
        ic_code      => $ic_code,
        gene_name    => $gene_name, 
        dbxref       => $dbxref
    };    
} 
close IN; 
if ($go_file) {
    open(GO,$go_file);
    while (<GO>) {
        chomp;
        my ($primary_id,$go_id) = (split(/\t/,$_));
        push @{$go{$primary_id}} , $go_id;
    }
    close GO;
}
if ($interpro_file) {
    open(IPR,$interpro_file);
    while (<IPR>) {
        chomp;
        my ($primary_id,$ipr_id) = (split(/\t/,$_));
        push @{$ipr{$primary_id}} , $ipr_id;
    }
    close IPR;
}   
while (my $feature = $gff_in->next_feature) { 
    if ($feature->primary_tag() eq "gene") {
        my ($gene_id) = $feature->get_tag_values("ID");
        $gene_id =~ s/_g/_t/;
        $feature->add_tag_value("Note",$locus_info{$gene_id}{product}) if $locus_info{$gene_id}{product};    
        $gff_out->write_feature($feature); 
    }
    elsif ($feature->primary_tag() eq "mRNA") { 
        my ($mrna_id) = $feature->get_tag_values("ID");  
        if (defined $locus_info{$mrna_id}) { 
            $feature->add_tag_value("Note",$locus_info{$mrna_id}{product})  if $locus_info{$mrna_id}{product};
            $feature->add_tag_value("Dbxref",$locus_info{$mrna_id}{dbxref}) if $locus_info{$mrna_id}{dbxref} &&  $locus_info{$mrna_id}{dbxref} ne "N/A";
            $feature->add_tag_value("Ontology_term","CC_functional_completeness:".$locus_info{$mrna_id}{completeness}) if $locus_info{$mrna_id}{completeness}; 
            $feature->add_tag_value("Ontology_term","CC_evidence_code:".$locus_info{$mrna_id}{ic_code}) if $locus_info{$mrna_id}{ic_code};
            if (defined $locus_info{$mrna_id}{gene_name}) {
                my @alias = (split(/,/,$locus_info{$mrna_id}{gene_name}));
                foreach my $alias (@alias) {
                    $feature->add_tag_value("Alias",$alias);
                }
            }

        }

        if (defined $ipr{$mrna_id}){
            foreach my $ipr_id (@{$ipr{$mrna_id}}) {
                $feature->add_tag_value("Dbxref","InterPro:".$ipr_id);   
            }
        }
        if (defined $go{$mrna_id}) {
            foreach my $go_id (@{$go{$mrna_id}}) {
                $feature->add_tag_value("Ontology_term",$go_id);      
            }                
        }
        $gff_out->write_feature($feature); 
    }  
    else {
        $gff_out->write_feature($feature); 
    }  
}
$gff_in->close;
$gff_out->close;