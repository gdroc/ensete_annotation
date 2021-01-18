#!/usr/bin/perl
use FindBin;
use Getopt::Long;
use Pod::Usage; 
use Cwd;
my $result         = "";
my $directory      = cwd();
my $verbose        = 0;
my $current_dir    = getcwd; 
my $correspondance = "";
my $prefix         = "locus"; 
my $help           = "";
my %ipr2go;
my %go2ipr;
my %exist;
my %list_gene; 
my %list;
# List of file to remove
my @remove = (".log.txt",".sequence.txt",".xml.xml",".out.txt",".htmltarball.html.tar.gz",".gff.txt");


my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters  
    --result    Directory containing all the Interpro output files [required]  
    --prefix    Prefix for output 
    --help
/;
GetOptions(
    'result=s'         => \$result, 
    'prefix=s'         => \$prefix,  
    'help|h|?'         => \$help
)  or pod2usage(-message => $usage);
if ($help) { pod2usage(-message => $usage); }  

if ($result eq "") {
    print "\nWarn :: --result is empty\nPlease specify a directory of Interpro Result\n\n";
    print $usage;
    exit 0;
} 

system("mkdir $directory") unless -e $directory;
# Download interpro2go
my $interpro2go = "interpro2go";
my $url = "http://geneontology.org/external2go/interpro2go";
system("wget --quiet $url") unless -e $interpro2go;


# Create list between IPR id and GO id
open(IPR, $interpro2go);
while (<IPR>) {
    chomp;
	# InterPro:IPR000018 P2Y4 purinoceptor > GO:purinergic nucleotide receptor activity, G-protein coupled ; GO:0045028
    my ($ipr_id, $go_id )  = ($_ =~ /^InterPro:(IPR\d+)\s.*\s>\s.*\s;\s+(GO:\d+)/) ; 
    $ipr_id =~ s/InterPro\://;
    push @{$ipr2go{$ipr_id}} , $go_id; 
    push @{$go2ipr{$go_id}} , $ipr_id;  
}
close IPR;

# Scan Directory and push all IPR id for each locus
open(DIR,"ls $result/*.tsv |");
while (my $file = <DIR>) {
    chomp($file); 
    open(IN,$file);  
    while(<IN>) {
        chomp;
        my ($gene_id,$type,$ipr_id,$desc) = (split(/\t/,$_))[0,8,11,12];
        if ($ipr_id) {
            if (defined $exist{$gene_id}{$ipr_id}){
                next;
            }
            else {
                $exist{$gene_id}{$ipr_id} = 1; 
                push @{$list_gene{$gene_id}} , $ipr_id; 
            }
        }
    }
    close IN;
}
close DIR;
 
open(IPR, ">". $directory ."/". $prefix."_ipr.txt");
open(GO, ">". $directory ."/".  $prefix."_go.txt");

foreach my $gene_id (keys %list_gene) { 
    foreach my $ipr_id (@{$list_gene{$gene_id}}){
        print IPR join("\t",$gene_id,$ipr_id),"\n";
        foreach my $go_id (@{$ipr2go{$ipr_id}}) {
            if (defined $exist{$gene_id}{$go_id}){
                next;
            }
            else {
                $exist{$gene_id}{$go_id} = 1;  
                print GO join("\t",$gene_id,$go_id),"\n"; 
            }
        }
    }
}
close GO;
close IPR; 
