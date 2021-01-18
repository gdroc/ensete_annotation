#!/usr/bin/perl
use Bio::SeqIO;
use Bio::SearchIO;
use File::Basename;
use FindBin;
use Pod::Usage; 
use Getopt::Long;
use Data::Dumper;
use Cwd;
my $dir         = cwd();
my $verbose     = 0;
my $current_dir = getcwd; 

my $blast;
my ($product,$gene,$source); 
my $output;  

my $help;
my @source;
my %cc_gene;
my %genedb_products; 

my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help

Parameters
    --blast         Blast [required]  
    --output        Output filename 
    --help          Show this help and exit
/;
GetOptions(
     'blast=s@'  => \$blast,      
     'output=s'  => \$output, 
     'help|h|?'  => \$help
)  or pod2usage(-message => $usage);
if ($help) { pod2usage(-message => $usage); }

if ($blast eq "") {
    print "\nWarn :: --blast is empty. Please specify the blastp file \n\n";
    print $usage;
    exit 0;
}

foreach my $blast_result ( @$blast) {
     my ($name,$path,$suffix) = fileparse($blast_result,".out"); 
     ($product,$gene,$source) = &parse_blastp($blast_result,$product,$gene,$name);
     push @source , $source;
}  
open(OUT,">$output"); 
foreach my $locus (sort keys %$product) {
    my $function = "Hypothetical protein";
    my $cds      = "missing_functional_completeness";
    my $evidence_code = "ISS_5"; 
    my @gene;  
    my @match;
    foreach my $source (@source) { 
        push @match, @{$product->{$locus}->{$source}};
		push @gene , @{$gene->{$locus}->{$source}};
    }
    my @sorted =  sort { $a->{evidence_code} cmp $b->{evidence_code} } @match;
     
    my $best_match = $sorted[0];
    $function      = $best_match->{product};
    $cds           = $best_match->{cds_type};
    $evidence_code = $best_match->{evidence_code};
    my $uniprot_name = $best_match->{uniprot_name} ;  
    @gene = @{&sort_uniq(\@gene)};
    my $gene_name = @gene ? join(",",@gene) : "unknown_gene";
    print OUT join("\t",$locus, $function, $evidence_code ,$cds,$gene_name,$uniprot_name),"\n"; 
         
}    
close OUT;  

sub parse_blastp {
    my ($file,$product,$gene_name,$bank) = @_;
    $bank = (split(/\_/,$bank))[1]; 
    my $searchio = new Bio::SearchIO( 
        -format => 'blast',
        -file   => $file
    );   
    my %dbxref; 
    my %product = %$product;
    my %gene_name = %$gene_name; 
    while(my $result = $searchio->next_result()) {     
        my $have_product = join("_","have_product",$bank);
        my $query_length = $result->query_length;
        my $query_name   = $result->query_name;
        my $query_desc   = $result->query_description; 
        my $locus_tag    = $query_name;
        my $rank_hit = 0;
        my $num_hits = $result->num_hits();
        if ($num_hits > 0) {
            while(my $hit = $result->next_hit()) { 
                $rank_hit++;
                last if $rank_hit == 5;
                my $hit_length = $hit->length();
                my $hit_name = $hit->name();
                my ($description,$species,$gene_name,$dbxref,$uniprot_name,$alias);	
                if ($bank =~ "SwissProt" || $bank eq "TrEMBL") { 
                    ($name,$alias) =(split(/\|/,$hit_name))[1,2]; 
                    if ($name eq "") {
                        $name = $hit_name;
                    } 
                    if ($hit->description() =~ /(.*)\sOS=.*\sOX=\d+\sGN=(\S+)\s.*/) { 
                        $description = $1; 
                    	$gene_name = $2;
                    	push @{$gene_name{$locus_tag}{$bank}} , $gene_name if $gene_name;
                    }
                    elsif ($hit->description() =~ /(.*)\sOS=(.*)\sOX=(\d+)\sPE=/) {
                        $description = $1;
                    }
                    $dbxref = join(":", $bank, $name); 
                }			
                else {
                    $name = $hit_name; 
                    $description = $hit->description();
                    $dbxref = join(":", $bank, $name); 
                }	 
                $description = "Hypothetical protein" if $description eq "";    
                my $hsp = $hit->next_hsp();
                my $hsp_query_cum_length = $hsp->length( 'query');
                my $hsp_hit_cum_length = $hsp->length( 'hit');
                my $hsp_cum_length = $hsp->length( 'total');
                my $cum_num_identical = $hsp->num_identical();
                my $strand = $hsp->strand('hit') == 0 ? "+" : "-"; 
                my $query_start = $hsp->start('query');
                my $query_end   = $hsp->end('query');
                my $hit_start   = $hsp->start('subject');
                my $hit_end     = $hsp->end('subject'); 
                my $evalue      = $hit->significance();
                my $qcov        = sprintf('%.2f',($hsp_query_cum_length / $query_length));
                my $scov        = sprintf('%.2f',($hsp_hit_cum_length / $hit_length));
                my $identity    = sprintf('%.2f',($cum_num_identical / $hsp_cum_length));
                my ($product,$cds_type, $evidence_code) = &eval_prediction($qcov,$scov,$identity,$description);  
                push @{$product{$locus_tag}{$bank}} , {
                    product       => $product,
                    cds_type      => $cds_type,
                    qcov 	      => $qcov,
                    scov	      => $scov,
                    identity      => $identity,
                    evidence_code => $evidence_code,
                    hit_name      => $hit_name,
                    description   => $description,
                    uniprot_name  => $dbxref
                };     
            }  
        } 
        else {
            push @{$product{$locus_tag}{$bank}} , { 
                product =>"Hypothetical protein",  
                cds_type => "missing_functional_completeness",
                qcov 	 => "N/A",
                scov	 => "N/A",
                identity => "N/A",
                evidence_code => "ISS_6",
                hit_name      => "N/A",
                description   => "N/A",
                uniprot_name  => "N/A"
            };        
        }
    }	
    return (\%product,\%gene_name,$bank);
}			

sub _max (@) {
    my $i = shift;
    foreach (@_) {
		$i = $_ if $_ > $i;
    }
    return $i;
}


sub _min (@) {
    my $i = shift;
    foreach (@_) {
		$i = $_ if $_ < $i;
    }
    return $i;
}

sub eval_prediction {
    my ( $QCov, $SCov, $identity, $description ) = @_;  
    my $evidence_code = "ISS_5";
    my $product = "Hypothetical protein";
    my $cds_type = "missing_functional_completeness";  
    if( $identity < 0.25 ) {
        $evidence_code = "ISS_5";
        $product = "Hypothetical protein";
    }
    elsif (( $description =~ /hypothetical protein/i ) || ( $description =~ /uncharacterized protein/i ) || ( $description =~ /^unknown protein$/i ) || ( $description =~ /^predicted protein$/i )) {
        if ( $QCov >= 0.8 ) {
            $product = "Conserved hypothetical protein";
            $evidence_code = "ISS_4";
            if ( $SCov >= 0.8 ) {
                $cds_type = "complete";
            }
            else {
                $cds_type = "fragment";
            }
        }
        elsif ( $SCov >= 0.8 ) {
            $product = "Conserved hypothetical protein";
            $evidence_code = "ISS_4";
            $cds_type = "modules";
        }
    }
    elsif (( $description =~ /putative/i ) || ( $description =~ /probable/i ) || ( $description =~ /probably/i ) || ( $description =~ /predicted/i ) || ( $description =~ /^Unknown /i )) {
        if ( $QCov >= 0.8 ) {
            $evidence_code = "ISS_3";
            $product = $description; 
            if ( $product =~ /putative/i ) {
                $product =~ s/putative/Putative/;
            }
            elsif ( $description =~ /probable/i ) {
                $product =~ s/[Pp]robable/Putative/;
            } 
            elsif ( $description =~ /predicted/i ) {
                $product =~ s/[Pp]redicted/Putative/;
            }
            else {$product = $description;}
            if ( $SCov >= 0.8 ) {
                $cds_type = "complete";
            }
            else {
                $cds_type = "fragment";
            }
        }
        elsif ( $SCov >= 0.8 ) {
            $evidence_code = "ISS_3";
            $product = $description;
			if ( $product =~ /putative/i ) {
                $product =~ s/putative/Putative/;
            }
			elsif ( $description =~ /probable/i ) {
                $product =~ s/[Pp]robable/Putative/;
            } 
			elsif ( $description =~ /predicted/i ) {
                $product =~ s/[Pp]redicted/Putative/;
            }
			else {$product = $description;}
            $cds_type = "modules";
        }
    }
    else {
        $evidence_code = "ISS_4";
        $product = "Conserved hypothetical protein";		
        if ( $QCov >= 0.8 ) {
            $evidence_code = "ISS_3";
            $product = "Putative $description";
            if ( $SCov >= 0.8 ) {
                $cds_type = "complete";
            }
            else {
                $cds_type = "fragment";
            }
        }
        elsif ( $SCov >= 0.8 ) {
            $evidence_code = "ISS_3";
            $product = "Putative $description";
            $cds_type = "modules";
        } 
        if (( $identity >= 0.45 ) && ( $QCov >= 0.8 ) &&( $SCov >= 0.8 )) {
            $evidence_code = "ISS_2";
            $product = $description; 
       	}
        if (( $identity >= 0.9 ) && ( $QCov >= 0.9 ) &&( $SCov >= 0.9 )) {
            $product = $description;
            $evidence_code = "ISS_1";
        }
    }
    return ( $product, $cds_type, $evidence_code );
}

 

sub sort_uniq {
    my $array = shift;
    my %seen = ();
    my @array = grep { !$seen{$_}++ } @{$array};
    return \@array;
}

sub encod {
    my  $encod = shift;
	$encod =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!-])/uc sprintf("%%%02x",ord($1))/eg; 
	return $encod;
}
