use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use Parallel::ForkManager;
use Cwd;
require "./komodize.pm";

my $infile;
my $verbose;
my $root_outdir;
my $MAX_PROCESSES;
my $busco_cutoff;
my $busco_database;
my $project_config;
#Getting options from command line

GetOptions ("in|infile=s" => \$infile, # file containing species names or taxonomy IDs 
	    "verbose=i"  => \$verbose,   #as
            "outdir|out=s" => \$root_outdir,
            "cpu=i" => \$MAX_PROCESSES,
            "busco_cutoff=i" => \$busco_cutoff,
	    "busco_database=s" => \$busco_database,
            "project_config=s" => \$project_config)
or die("Error in command line arguments\n");

#read config files
my $komodize_path = dirname (rel2abs($0));
my $parameters = read_config_files ($project_config, $komodize_path);
check_parameters($parameters); #checkin if user setted the parameters right

my $pm = Parallel::ForkManager->new($parameters->{max_processors});

open (In, $parameters->{infile});
my @queries = <In>;

#create directory structure 
create_directory_structure($parameters->{outdir}, $parameters->{extract_mode});
my $log_fh;
open($log_fh, ">>", "$parameters->{outdir}/log.txt") || die("Can't create log.txt");

#open and parse infile
my @species = parse_infile(\@queries);
my @ids = @species;
#print Dumper (@ids);
# Download genome data from NCBI if not already downloaded
#my @ids = get_ids(\@queries, $root_outdir);
#print Dumper ("@ids");

print $log_fh ("Welcome to Komodize Genomes\n");


@ids = get_genomes(\@ids, $parameters->{outdir}, $parameters->{max_processors}, $log_fh, $parameters->{assembly_summary});

#print Dumper (@ids);

my @ids2 = @ids;
@ids = ();
foreach my $id2(@ids2){
  my $id = $id2;
  $id =~ s/\s/\_/g ;
  push @ids, $id ;
}
#print Dumper (@ids);

# Extract protein from genomes
if ($parameters->{extract_mode} eq "protein") {
  PROTEINS:
  foreach my $id (@ids) {
    my $pid = $pm -> start and next PROTEINS;
    if (-e "$parameters->{outdir}/longest_proteins/$id.aa.fasta_longest"){
      print "-> Proteins form $id already extracted\n";
      print $log_fh "-> Proteins form $id already been extracted\n";
    }    
    else{
      print ("-> Extractiong Proteins from $id\n") if ($parameters->{verbose});
      print $log_fh ("-> Extractiong Proteins from $id\n");
      my $log = get_proteins ($parameters->{outdir}, $id);
      if ($log) {
        print $log_fh ("$log\n");
      }
    }  
      $pm -> finish;
  }
    $pm -> wait_all_children;
}
if ($parameters->{extract_mode} eq "orfs"){
# Extract ORFs from genomes
  ORFS:
  foreach my $id(@ids){
    my $pid = $pm ->start and next ORFS;
    if (-e "$parameters->{outdir}/longest_ORFS/$id\_nt.fasta_longest"){
      print "-> ORFS form $id already extracted\n";
      print $log_fh "-> ORFS form $id already been extracted\n";
    }
    else {
      print ("-> Extracting ORFS from $id\n") if ($parameters->{verbose});
      print $log_fh "-> Extracting ORFS from $id \n";
      my $log = get_ORF($parameters->{outdir}, $id);
      if ($log) {
        print $log_fh ("$log\n");
      }
    }
    $pm -> finish;
  }
  $pm ->wait_all_children;
}
#get longest ORF for each locus to summarize per locus
LONGEST_SEQUENCE:
foreach my $id (@ids) {
  my $pid = $pm ->start and next LONGEST_SEQUENCE;
  
  if (-e "$parameters->{outdir}/longest_proteins/$id.aa.fasta_longest"){
    print "-> Genome $id already summarized per locus\n";
    print $log_fh "-> Genome $id already summarized per locus\n";
  }
  else {
    print ("-> Summarizing genome $id per locus\n") if ($parameters->{verbose});
    print $log_fh "-> Summarizing genome $id per locus\n";
    my $log = get_longest_sequence($id, $parameters->{outdir}, $parameters->{extract_mode});
    if ($log) {
      print $log_fh("$log\n");
    }
  }
  $pm -> finish;
}
$pm ->wait_all_children;

if ($parameters->{extract_mode} eq "orfs"){
#translating
  TRANSLATE:
  foreach my $id (@ids) {
    my $pid = $pm ->start and next TRANSLATE;
    my $log = translate($id, $parameters->{outdir}); 
    print ("-> Translating ORFs for genome $id\n") if ($parameters->{verbose});
    print $log_fh "-> Translating ORFs for genome $id\n";
    if ($log) {
      print $log_fh("$log\n");
    }
    $pm -> finish;
  }
  $pm ->wait_all_children;
}
#checking for BUSCO completeness
foreach my $id (@ids) { 
  if (-e "$parameters->{outdir}/busco/run_$id"){
    print "-> Busco completeness already check for genome $id\n";
    print $log_fh "-> Busco completeness already check for genome $id\n";
  }
  else { 
    print "-> Checking genome $id for completeness\n";
    print $log_fh "-> Checking genome $id for completeness\n";
    my $log = run_busco($id, $parameters->{outdir}, $parameters->{busco_database},$parameters->{max_processors},$parameters->{python_path}, $parameters->{busco_path});
    if ($log) {
      print $log_fh ("$log\n");
    }
  }
}

my $completeness;
my $duplicate;
CHECK_COMPLETENESS:
foreach my $id (@ids) {
  my $pid = $pm ->start and next CHECK_COMPLETENESS;
  ($completeness, $duplicate) = check_completeness ($id, $parameters->{outdir});
  print "-> Genome $id completeness is $completeness and has $duplicate% duplicate\n";
  print $log_fh "-> Genome $id completeness is $completeness and has $duplicate% duplicate\n";
  $pm ->finish;
} 
$pm -> wait_all_children;

#annotating genomes using InterproScan
foreach my $id (@ids) {
  print ("-> Annotating genome $id\n") if ($parameters->{verbose});
  print $log_fh "-> Annotating genome $id\n";
  my $log = annotate_genome($id, $parameters->{outdir}, $parameters->{max_processors}, $parameters->{interpro_path});
  if ($log) {
    #print $log_fh("$log\n");
  }
}


#creating KOMODO2-compatible output files
GENE2GO:
foreach my $id (@ids) {
  my $pid = $pm ->start and next GENE2GO;
  print "-> Creat komodo compatible imput GENE to Gene Ontology\n";
  my $log = gene2GO ($id, $parameters->{outdir});
  if ($log) {
    #print $log_fh("$log\n"); 
  $pm ->finish;
  }
  $pm -> wait_all_children;
}

FEATURE2GO:
foreach my $id (@ids) {
  my $pid = $pm ->start and next FEATURE2GO;
  print "-> Creat komodo compatible imput Pfam to Gene Ontology\n";
  my $log = feature2GO ($id, $parameters->{outdir});
  if ($log) {
    #print $log_fh("$log\n");
  $pm ->finish;
  }
}
$pm -> wait_all_children;

PARSER_PFAM:

foreach my $id (@ids) {
  my $pid = $pm ->start and next PARSER_PFAM;
  print "-> Creat komodo compatible imput Pfam\n";
  my $log = parser_pfam ($id, $parameters->{outdir});
  if ($log) {
    #print $log_fh("$log\n");
  $pm ->finish;
  }
}
$pm -> wait_all_children;

#GENE2SUPERFAMILY:
#foreach my $id (@ids) {
#  my $pid = $pm ->start and next GENE2SUPERFAMILY;
#  print "creat komodo compatible imput GENE to Superfamily \n";
#  my $log = gene2superfamily ($id, $parameters->{outdir});
#  if ($log) {
#    print $log_fh("$log\n");
#  $pm ->finish;
#  }
#  $pm -> wait_all_children;
#}
print "Finish! Closing Komodize, see you around\n";
