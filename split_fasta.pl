#!/usr/bin/perl

use Bio::SeqIO; 
use Getopt::Long;
use FindBin;
use Pod::Usage;
use Cwd;
my $file        = "";
my $prefix      = "fasta";
my $number      = 5000;
my $tab_file    = ""; 
my $directory   = cwd();
my $verbose     = 0;
my $current_dir = getcwd; 

my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters   
    --file     FASTA file, required
    --prefix   Prefix file name for output (default : fasta)
    --dirout   Output directory (default : current_directory)
    --number   Number of sequences to be breaked (default : 5000) 
    --help
/;
GetOptions( 
    'file=s'   => \$file,
    'prefix=s' => \$prefix, 
    'number=i' => \$number, 
    'dirout=s' => \$directory,
    'help|h|?' => \$help
)  or pod2usage(-message => $message);
if ($help) { pod2usage(-message => $message); }
if ($file eq "") {
    print "\nWarn :: --file is empty. Please specify the FASTA file\n\n";
    print $usage;
    exit 0;
}
 
system("mkdir $directory") unless -e $directory;
my $in  = new Bio::SeqIO(
    -file  => $file,
    -format => 'fasta'
);

my $count = 0;
my $fcount = 0;
my $out;
while (my $seq = $in->next_seq) {
    if ($count % $number == 0) {
        $fcount++;
        my $file_out = $directory ."/".$prefix . $fcount.".fna";
        if ($number == 1) {
            $file_out = $directory ."/" .$seq->display_id.".fa";
        }
        $out = new Bio::SeqIO(
            -file => ">$file_out",
            -format=>'fasta'
        );
    }
    $out->write_seq($seq);
    $count++;
}