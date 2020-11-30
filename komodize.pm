use warnings;
use strict;

sub read_config_files { #read config files in the form element = value #comment
  my $project_config = $_[0];
  my $komodize_path = $_[1];
  my %parameters;
  open(my $fh_project_config, "<", "$project_config") || die ("Couldn't open the project configuration file, you may have to chech if the name or permission of the file are correct.\nDetails: $!\n");
  open(my $fh_komodize_config, "<", "$komodize_path/komodize_config") || die ("The configuration file 'potion_config' couldn't be read, please check if the file exists and if its permissions are properly set.\n");
  
  # -=-=-= INPUT FILES =-=-=-
  $parameters{infile} = read_config_file_line ('infile', $fh_project_config);
  
  # -=-=-= OUTPUT DIR =-=-=-
  $parameters{outdir} = read_config_file_line ('outdir', $fh_project_config);

  # -=-=-= PROJECT CONFIGURATION =-=-=-
  $parameters{busco_database} = read_config_file_line ('busco_database', $fh_project_config);
  $parameters{assembly_summary} = read_config_file_line ('assembly_summary', $fh_project_config);
  $parameters{max_processors} = read_config_file_line ('max_processors', $fh_project_config);
  $parameters{verbose} = read_config_file_line ('verbose', $fh_project_config);
  $parameters{extract_mode} = read_config_file_line ('extract_mode', $fh_project_config);
  
  close($fh_project_config);
 
  # -=-=-= PATH TO EXTERNAL PROGRAMS =-=-=-
  $parameters{python_path} = read_config_file_line ('phyton', $fh_komodize_config);
  $parameters{busco_path} =  read_config_file_line ('busco', $fh_komodize_config);
  $parameters{interpro_path} = read_config_file_line ('interpro', $fh_komodize_config);

 close($fh_komodize_config);

 return \%parameters;
}

sub read_config_file_line { # file format element = value
  my ($parameter, $fh_config_file) = @_;

  seek($fh_config_file, 0, 0);              # added to allow any order of parameters in the config files, preventing unfriendly error messages if the user changes the order
  while (my $line = <$fh_config_file>){
    if ($line =~ /^\s*$parameter\s*=/) {    # the string to be searched in the file
      chomp ($line);
      $line =~ s/^\s*$parameter\s*=\s*//;   # removing what comes before the user input
      $line =~ s/#.*$//;                    # removing what comes after the user input (commentaries)
      $line =~ s/\s*$//;                    # removing what comes after the user input (space caracteres)
      if ($line eq 'undef' || $line eq '') { return; }
      else { return $line; }
    }
  }
  return;
}

sub check_parameters { #check for all parameters, 
  my $parameters = shift;

  my $config_path = getcwd();
  $config_path =~ s/\/\w+$/\/config_files/;

  # -=-=-= INPUT FILES =-=-=-
  if (!defined $parameters->{infile}) { die ("No path to the species files was specified in your project's configuration file, please fill the parameter 'infile'.\n"); }
  if (!-s $parameters->{infile}) { die ("The path to your project's species files isn't a valid file, please check if the path in 'infile' is correct: $parameters->{infile}\n"); }
  if (!-r $parameters->{infile}) { die ("You don't have permission to read in your project's species directory, please redefine your permissions.\n"); }

  # -=-=-= OUTPUT DIRECTORY =-=-=-
  if (!defined $parameters->{outdir}) { die ("No path to outdir was specified in project_config, please open this file and fill the parameter 'outdir'.\n"); }
  #if (!-d $parameters->{outdir}) { die ("The path to outdir isn't a valid directory, please check if the path in 'outdir' is correct: $parameters->{outdir}\n"); }
  #if (!-w $parameters->{outdir}) { die ("You don't have permission to write in the Komodize directory, please redefine your permissions for this directory.\n"); }
  
  # -=-=-= OUTPUT DIRECTORY =-=-=-
  if (!defined $parameters->{assembly_summary}) { die ("No path to assembly summary files was specified in your project's configuration file, please fill the parameter 'assembly_summary'.\n"); }
  if (!-s $parameters->{assembly_summary}) { die ("The path to assembly summary files isn't a valid file, please check if the path in 'assembly_summary' is correct: $parameters->{assembly_summary}\n"); }
  if (!-r $parameters->{assembly_summary}) { die ("You don't have permission to read assembly summary file, please redefine your permissions.\n"); }
  
  # -=-=-= PATH TO EXTERNAL PROGRAMS USED BY KOMODIZE =-=-=-
  if (!defined $parameters->{python_path}) { die ("No path to Python was specified in komodize_config at $config_path, please open this file and fill the parameter 'Python3.6'.\n"); }
  if (!-s $parameters->{python_path}) { die ("The executable of Python wasn't found in the specified path, please check if the path is correct: $parameters->{python_path}\n"); }
  if (!-x $parameters->{python_path}) { die ("You don't have permission to execute the Python file specified at komodize_config, please check permissions or replace the file\n"); }

 if (!defined $parameters->{busco_path}) { die ("No path to BUSCO was specified in komodize_config at $config_path, please open this file and fill the parameter 'busco'.\n"); }
  if (!-s $parameters->{busco_path}) { die ("The executable of BUSCO wasn't found in the specified path, please check if the path is correct: $parameters->{busco_path}\n"); }
  if (!-x $parameters->{busco_path}) { die ("You don't have permission to execute the BUSCO file specified at komodize_config, please check permissions or replace the file\n"); }

  if (!defined $parameters->{interpro_path}) { die ("No path to INTERPRO was specified in komodize_config at $config_path, please open this file and fill the parameter 'interpro'.\n"); }
  if (!-s $parameters->{interpro_path}) { die ("The executable of INTERPRO wasn't found in the specified path, please check if the path is correct: $parameters->{interpro_path}\n"); }
  if (!-x $parameters->{interpro_path}) { die ("You don't have permission to execute the INTERPRO file specified at komodize_config, please check permissions or replace the file\n"); }

  # -=-=-= PROJECT CONFIGURATION =-=-=-  
  if (!defined $parameters->{verbose}) { $parameters->{verbose} = 0; } # default value

  if (!defined $parameters->{max_processors}) {die "Max processors not setted. Please set max_processors element in configuration file\n";}

  if (!defined $parameters->{busco_database}) { die ("No path to busco_database was specified in project_config, please open this file and fill the parameter 'busco_database'.\n"); }
  if (!-d $parameters->{busco_database}) { die ("The path to busco_database isn't a valid directory, please check if the path in 'busco_database' is correct: $parameters->{outdir}\n"); }
  #if (!-w $parameters->{potion_dir}) { die ("You don't have permission to write in the Potion directory, please redefine your permissions for this directory.\n"); }
  
  if (!defined $parameters->{extract_mode}) {die "Extract_mode not setted. Please set extract_mode element in configuration file\n";}
  if (($parameters->{extract_mode} !~ /[protein|orfs]/)) {die "You must set the parameter \"extract_mode\" with \"proteins\" or \"orfs\" \n";}
}

sub create_directory_structure {
  my $root_outdir = $_[0];
  my $extract_mode = $_[1];
  if (!-e $root_outdir) { mkdir ($root_outdir) || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n"); 
  }
  else {
    my $i = 0;
    while ($i < 1) {
      print"Directory $root_outdir already exists, do you want to continue? This action can overwrite your files. (yes or not)\n";
      my $check = <STDIN>;
      chomp $check;
      if ($check eq "not" || $check eq "n") { 
        die ("closing komodize_genomes.pm\n");
      }
      elsif ($check eq "yes" || $check eq "y") {
      #  print "You have the power and the consequences are yours\n";
        print "Trying to continue with Komodize\n";
        $i = 1;
      }
      else {
        print "Invalid option, try again\n";
      }
    }
  }
  if (!-e "$root_outdir/genomes") {mkdir ("$root_outdir/genomes") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
  }
  if ($extract_mode eq "orfs"){
    if (!-e "$root_outdir/ORFS") {mkdir ("$root_outdir/ORFS") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
    }
    if (!-e "$root_outdir/longest_ORFS") {mkdir ("$root_outdir/longest_ORFS") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
    }
  }
  if ($extract_mode eq "protein") {
    if (!-e "$root_outdir/proteins") {mkdir ("$root_outdir/proteins") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
    }
  }
  if (!-e "$root_outdir/longest_proteins") {mkdir ("$root_outdir/longest_proteins") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
  }
  if (!-e "$root_outdir/interpro") {mkdir ("$root_outdir/interpro") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
  }
  if (!-e "$root_outdir/komodo2_imputs") {mkdir ("$root_outdir/komodo2_imputs") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
    mkdir ("$root_outdir/komodo2_imputs/feature2GO");
    mkdir ("$root_outdir/komodo2_imputs/gene2GO");
    mkdir ("$root_outdir/komodo2_imputs/pfam");
  }
}


sub parse_infile {

  my @queries;
  my @species = @{$_[0]};
  my $i = 0;

  foreach my $line(@species) {
    chomp $line;
    $queries[$i] = $line;
    # print $line;
    $i++;
  }
  return (@queries);
}

sub get_ids {
  my @species = @{$_[0]};
  use Time::HiRes qw(time);
  my @ids;
  if ( -e "TINY.XML") {die "ABORTING: TINY.XML exists!!!\n"};

  #my $lastupdate="2018/04/06";
  my $lastupdate="1918/04/06";
  foreach my $specie(@species){
    #my $query = '%22mus+musculus%22[orgn]+OR+%22homo+sapiens%22[orgn]'; #%22 é o " (aspas duplas), "+" é o espaço;   #lista tem um limite de 1000 caracteres(aprox?)
    my $ACCESS;
    $specie =~ s/\n//;
    my $query = $specie;
    $query =~ s/\s/+/g;
    $query = "%22$query%22[orgn]";

    #assemble the esearch URL
    my $db		= "genome";
    my $base	= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url		= $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

    #post the esearch URL
    my $output = get($url);
    #parse WebEnv, QueryKey and Count (# records retrieved)

    my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
    my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
    my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
    $efetch_url .= "&query_key=$key";
    $efetch_url .= "&rettype=docsum";
    #print $efetch_url."\n";
    my $genome = get($efetch_url);
    while ($genome =~ /<Item Name="Assembly_Accession" Type="String">(\S+)</g) {
      $ACCESS= $1; 
    }
    if (defined $ACCESS) {
      push @ids, $ACCESS; 
    }
    else {
      print "$specie no contains genome id\n";
    }
  }
  return @ids;
}
sub get {
    my $url=shift @_;
    return `wget -q -O - "$url"`;
}


sub get_genomes { 
  
  my @all_ids = @{$_[0]};
  my $cpu = $_[2];    
  my $root_outdir = $_[1];
  my $log = $_[3];
  my $assembly_summary = $_[4];
  open (Lista, $assembly_summary) or die ("assembly_summary not found\n");
  my %urls;
  my @species_filtered;
  my @all_species_in_assembly_summary;
  while (<Lista>){
    if ($_ =~ m/^#/){next;}
    my @line = split ('\t', $_);
    $urls {$line[7]} = $line[19];
    push @all_species_in_assembly_summary, $line[7];
  }
  my @ids_err;
  open (ID_ERR, ">", "$root_outdir/especies_not_found.txt") or die ("");
  foreach my $id(@all_ids){
    if (grep (/^$id$/, @all_species_in_assembly_summary)){
     # print "I FOUND $id IN ASSEMBLY_SUMMARY_REFSEQ\n";
      push @species_filtered, $id;;
    }
    else {
      print "I DID NOT FIND $id IN ASSEMBLY_SUMMARY_REFSEQ. CHECK IF THE SPECIE NAME IS CORRECT OR IF IT HAS A SUB-SPECIE THIRD NAME. eg *Canis lupus familiaris*\n";
      print $log "I DID NOT FIND $id IN ASSEMBLY_SUMMARY_REFSEQ. CHECK IF THE SPECIE NAME IS CORRECT OR IF IT HAS A SUB-SPECIE THIRD NAME. eg *Canis lupus familiaris*\n";
      push @ids_err, $id;
      print ID_ERR "$id\n";
    }
    
  }
  close ID_ERR;
  #print " I DID NOT FIND THIS SPECIES IN IN ASSEMBLY_SUMMARY_REFSEQ. CHECK IF THE SPECIE NAME IS CORRECT OR IF IT HAS A SUB-SPECIE THIRD NAME. eg *Canis lupus familiaris*\n";
  #foreach my $id(@ids_err){
  #print "$id\n";
  #}

  my @ids = @species_filtered;
  my $pm = Parallel::ForkManager->new($cpu);
  #print Dumper (@ids);
  #my $kk = <STDIN>;
  DOWNLOAD:
    foreach my $id(@ids){
    sleep (2); 
    my $pid = $pm ->start and next DOWNLOAD;
    # my @comp = split ('\s',$id);
    my $id2 = $id;
    $id2 =~ s/\s/\_/g;
    my @vetor = split ('\/', $urls{$id});
    my $var = pop @vetor;
  
    $urls{$id}=~s/\/($var.*?)$/\/$1\/$1_rna.gbff.gz/;
  
    if (-e "$root_outdir/ORFS/$id2.nt.fasta" || -e "$root_outdir/proteins/$id2.aa.fasta"){
      print "-> Genome $id already downloaded\n";
      print $log "-> Genome $id already downloaded\n";
    }
    else {
      print "-> Downloading genome $id\n";
      print $log "-> Downloading genome $id\n";   
      my $dowload = `wget  -q -O $root_outdir/genomes/$id2.gbff.gz  '$urls{$id}'`;
      #print $dowload."\n";
      my $gzip = `gzip -d $root_outdir/genomes/$id2.gbff.gz`;
    }
  $pm -> finish;
  sleep (2);
  }
  $pm -> wait_all_children;
  return (@ids);
}

sub get_proteins {
  
  my $root_outdir = $_[0];
  my $id = $_[1];
  my $path_2_file = "$root_outdir/genomes/$id.gbff";
  my $in = new Bio::SeqIO(-format => 'genbank',
                           -file => $path_2_file);
  my $seq;

  my $name = $id;
  my $i = 1;

  my $outfile_path = "$root_outdir/proteins/$id.aa.fasta";
  #next if (-e $outfile_path);
  open (OUT, ">", "$outfile_path");

  while($seq = $in->next_seq()){
    my $start;
    my $end;
    my $organism;
    my $source;
    my $cds;
    my $flag = 0;
    for my $feat_object ($seq->get_all_SeqFeatures) {
      #print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
        for my $value ($feat_object->get_tag_values($tag)) {
          if ($tag eq "pseudo") { #remove pseudogenes
            $flag = 1;
            next;
          }
        }
        if ( $feat_object->primary_tag eq 'CDS' ) {
          my $cds_object = $feat_object->spliced_seq;
          $cds = $cds_object->seq;
        }
        if ($feat_object->primary_tag eq "source") {
          $source = $id;
          $source =~s/\.gb.*$//g;
          $source =~ s/\s+/_/g;
          for my $tag ($feat_object->get_all_tags) {
            if ($tag eq "organism") {
              for my $value ($feat_object->get_tag_values($tag)) {
                $organism = $value;
              }
            }
          }
          $organism =~ s/\s+/_/g;
          $organism =~ s/\-/_/g;
        }
      }
      if ($feat_object->primary_tag eq "CDS") {
        my $final_line;
        $final_line = ">$source\_$i";
        for my $tag ($feat_object->get_all_tags) {
          if ($tag eq "protein_id") {
            for my $value ($feat_object->get_tag_values($tag)) {
              $final_line = $final_line."||protein_id:$value";
              $i++;
            }
          }
          if ($tag eq "locus_tag") {
            for my $value ($feat_object->get_tag_values($tag)) {
              $final_line = $final_line."||locus_tag:$value";
            }
          }
          if ($tag eq "gene") {
            for my $value ($feat_object->get_tag_values($tag)) {
              $final_line = $final_line."||gene:$value";
            }
          }
          if ($tag eq "translation") {
            for my $value ($feat_object->get_tag_values($tag)) {
              $final_line = $final_line."\n$value\n";
              print OUT $final_line;
              #print ("gene $i parsed\n");
            }
          }
        }
      }
    }
    $i++;
  }
}


sub get_ORF {
  my $root_outdir = $_[0];
  #print $root_outdir."  sub\n";
  my $file = $_[1];
  my $output_dir = "$root_outdir/ORFS/";

  my $len_cutoff = 100;

  my $flags = shift;

  if (! defined $len_cutoff) {
    $len_cutoff = 100;
  }

  if (! defined $flags) {
    $flags = "all";
  }
  my $path_2_file = "$root_outdir/genomes/$file.gbff";
  my $in = new Bio::SeqIO(-format => 'genbank',-file => $path_2_file);
  my $seq;
  my $name = $file;
  my $i = 1;
  open (OUT, ">", "$output_dir/$file.nt.fasta");
    
  while($seq = $in->next_seq()){
    my $sequence = $seq->seq();
    my $start;
    my $end;
    my $organism;
    my $source;
    my $cds;
    my $flag = 0;
    for my $feat_object ($seq->get_all_SeqFeatures) {
#    print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
        for my $value ($feat_object->get_tag_values($tag)) {
          if ($tag eq "pseudo") { #remove pseudogenes
          $flag = 1;
          next;
          }
        }
        if ( $feat_object->primary_tag eq 'CDS' ) {
          my $cds_object = $feat_object->spliced_seq;
          $cds = $cds_object->seq;
        }
        if ($feat_object->primary_tag eq "source") {
          $source = $file;
          $source =~s/\.gbff.*$//g;
          $source =~ s/\s+/_/g;
          for my $tag ($feat_object->get_all_tags) {
            if ($tag eq "organism") {
              for my $value ($feat_object->get_tag_values($tag)) {
                $organism = $value;
              }
            }
          }
          $organism =~ s/\s+/_/g;
          $organism =~ s/\-/_/g;
        }
      }
      if ($feat_object->primary_tag eq "CDS") {
        next if check_cds($cds); #next if cds contains any problem
        print OUT (">$source\_$i");
        for my $tag ($feat_object->get_all_tags) {
          if ($tag eq "protein_id") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||protein_id:$value");
              $i++;
            }
          }
          if ($tag eq "organism") { 
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||organism:$value");
            }
          }
          if ($tag eq "locus_tag") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||locus_tag:$value");
            }
          }
          if ($tag eq "gene") {
            for my $value ($feat_object->get_tag_values($tag)) {
              print OUT ("||gene:$value");
            }
          }
        }
        print OUT ("\n$cds\n");
      }
    }
    $i++;
  }
  
  sub check_cds {
    my $tmp_seq = $_[0];
    my $start_codon = substr($tmp_seq, 0, 3);
    my $stop_codon = substr($tmp_seq, -3);
    if ($start_codon !~ /ATG|GTG/i) {
      print "Not valid start codon\t$tmp_seq\n";
      return (1);
    }
    if ($stop_codon !~ /TAA|TAG|TGA/i) {
      print "Not valid stop codon\t$tmp_seq\n";
      return (1);
    }
    if (length($tmp_seq) % 3 == 0) {
    } else {
      print "length not multiple of three\t$tmp_seq\n";
      return (1);
    }
    if ($tmp_seq =~ /[^ACGT]/i) {
      print "non-standard nucleotides\t$tmp_seq\n";
      return (1);
    }
    if (length($tmp_seq) < $len_cutoff) {
      print "length smaller than cutoff\n";
      return (1);
    }
    return (0);
  }
}

sub get_longest_sequence {
  my $id = $_[0];
  my $root_outdir = $_[1];
  my $extract_mode = $_[2];
  my $seqio_obj;
  #print "LALALALALALLAA $id\n";
  if ($extract_mode eq "orfs"){
    open (OUT, ">", "$root_outdir/longest_ORFS/$id.nt.fasta_longest");

    $seqio_obj = Bio::SeqIO->new(  -format => "fasta",
                                      -file => "$root_outdir/ORFS/$id.nt.fasta"
                                   );
  }
  if ($extract_mode eq "protein"){
    open (OUT, ">", "$root_outdir/longest_proteins/$id.aa.fasta_longest");

    $seqio_obj = Bio::SeqIO->new(  -format => "fasta",
                                      -file => "$root_outdir/proteins/$id.aa.fasta"
                                   );
  }

  my %gene_data; 
  while (my $seq_obj = $seqio_obj->next_seq) {
    my $seq = $seq_obj->seq();
    my $acc = $seq_obj->id();
    my @aux = split(/\|\|/, $acc);
    my $uniq_id;
    my $locus_id;
    my $gene_id;
    my $protein_id;
    my $flag = 0;
    foreach my $element (@aux) {
     if ($element =~ /(^[A-Z][a-z])/) {
        $uniq_id = $element;
        next;
      }
      if ($element =~ /gene:/) {
        $gene_id = $element;
        next;
      }
      elsif ($element =~ /locus_tag:/) {
        $locus_id = $element;
        next;
      }
      elsif ($element =~ /protein_id:/) {
        $protein_id = $element;
        next;
      }
      else {
        $flag = 1;
      }
    }
    if ($flag == 1) {
      next;
      print $acc."\n";
    }
    next if (! defined $protein_id);
    next if ((! defined $gene_id)&&(! defined $locus_id));
    if (defined $gene_id) {
      if (defined $gene_data{$gene_id}{gene_id}{sequence}) {
        my $actual_length = length($gene_data{$gene_id}{gene_id}{sequence});
        my $new_length = length($seq);
        if ($new_length > $actual_length) {
          $gene_data{$gene_id}{gene_id}{sequence} = $seq;
          $gene_data{$gene_id}{gene_id}{protein_id} = $protein_id;
        }
      }
      else {
        $gene_data{$gene_id}{gene_id}{sequence} = $seq;
        $gene_data{$gene_id}{gene_id}{protein_id} = $protein_id;
      }
    }
    if (defined $locus_id) {
      if (defined $gene_data{$locus_id}{locus_id}{sequence}) {
        my $new_length = length($seq);
        my $actual_length = length($gene_data{$locus_id}{locus_id}{sequence});
        if ($new_length > $actual_length) {
          $gene_data{$locus_id}{locus_id}{sequence} = $seq;
          $gene_data{$locus_id}{locus_id}{protein_id} = $protein_id;
        }
      }
      else {
      $gene_data{$locus_id}{locus_id}{sequence} = $seq;
      $gene_data{$locus_id}{locus_id}{protein_id} = $protein_id;
      }
    }
  }
  my %print_flag;
  foreach my $key (keys %gene_data) {
    my $tmp_prot_id;
    if (defined $gene_data{$key}{gene_id}{protein_id}) {
      $tmp_prot_id = $gene_data{$key}{gene_id}{protein_id};
    }
    if (defined $gene_data{$key}{locus_id}{protein_id}) {
      $tmp_prot_id = $gene_data{$key}{locus_id}{protein_id};
    }
    if (defined $print_flag{$tmp_prot_id}) {
      next;
    }
    else {
      if ((defined $gene_data{$key}{gene_id}{sequence})&&(!defined $print_flag{$tmp_prot_id})) {
        print OUT (">$key|$gene_data{$key}{gene_id}{protein_id}\n$gene_data{$key}{gene_id}{sequence}\n");
        $print_flag{$tmp_prot_id} = 1;
        next;
      }
      if ((defined $gene_data{$key}{locus_id}{sequence})&&(!defined $print_flag{$tmp_prot_id})) {
        print OUT (">$key|$gene_data{$key}{locus_id}{protein_id}\n$gene_data{$key}{locus_id}{sequence}\n");
        $print_flag{$tmp_prot_id} = 1;
        next;
      }
    }
  }
}

sub translate {
  my $infile = $_[0];
  my $root_outdir = $_[1];
  chomp $infile;

  my $in  = Bio::SeqIO->new(-file => "$root_outdir/longest_ORFS/$infile.nt.fasta_longest" ,
                       -format => 'Fasta');
  my $out = Bio::SeqIO->new(-file => ">$root_outdir/longest_proteins/$infile.aa.fasta_longest" ,
                       -format => 'Fasta');

  while ( my $seq = $in->next_seq() ) {
    my $prot_obj = $seq->translate(-codontable_id => 1);
    $out->write_seq($prot_obj);
  }
  system ("sed -i s/*\$// $root_outdir/translated/$infile.aa.fasta");
}

sub run_busco{
  my $root_outdir = $_[1];
  if (!-e "$root_outdir/busco") {mkdir ("$root_outdir/busco") || die ("Couldn't create the directory specified as '$root_outdir', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");}
  my $id = $_[0];
  my $busco_database = $_[2];
  my $cpu = $_[3];
  my $python_path = $_[4];
  my $busco_path = $_[5];
  #print $cpu."AQUI\n\n";
  my $input = "$root_outdir/longest_proteins/$id.aa.fasta_longest";
  #print $busco_database."\n";
  system ("$python_path $busco_path -i $input -o $id -l $busco_database -c $cpu -m prot -q -f");
  system ("mv -u run_$id $root_outdir/busco"); 
}

sub check_completeness{
  my $id = $_[0];
  my $root_outdir = $_[1];
  my $completeness;
  my $duplicate;
  open (SUMMARY, "$root_outdir/busco/run_$id/short_summary_$id.txt") || die ("can't open short_summary_$id.txt");
  while (<SUMMARY>){
    my $lines = $_;
    if ($lines =~ m/C:/){
      my @comp = split ('%', $lines);
      #print Dumper (@comp);
      my @completeness_data = split (':', $comp[0]);
      $completeness = $completeness_data[1];
      my @duplicate_data = split (':', $comp[2]);
      $duplicate = $duplicate_data[1];
      #print "Genome $id completeness is $completeness and has $duplicate% duplicate\n";
      return ($completeness, $duplicate);
    }
  }  
}

sub annotate_genome{
  my $id = $_[0];
  my $root_outdir = $_[1];
  my $cpu = $_[2];
  my $interpro_path = $_[3];
  my $input = "$root_outdir/longest_proteins/$id.aa.fasta_longest";
  
  system ("$interpro_path -i $input  -iprlookup -goterms -pa --cpu $cpu -dra -d $root_outdir/interpro -f tsv");

}

sub gene2GO{
  my $id = $_[0];
  my $root_outdir = $_[1];

  #interpro output file, tsv-delimited
  my $infile = "$root_outdir/interpro/$id.aa.fasta_longest.tsv";

  #output directory
  my $outdir_path = "$root_outdir/komodo2_imputs/gene2GO";
  chomp $outdir_path;

  #the bunch of code below is just to generate a suitable output file path
  my @interpro_path = split(/\//, $infile);

  my $file_name = pop @interpro_path;

  my @tmp_name = split(/_/, $file_name);

  my $final_name = join("_", $tmp_name[0], $tmp_name[1]);

  my $outfile_path = join("/", $outdir_path, $final_name);

  $outfile_path = join(".", $outfile_path, "gene2GO", "txt");

  #print ("Generating file $outfile_path\n");

  if (-e $outfile_path) {
    die "Outfile $outfile_path already exists.\n";
    } 
  else {
    open(OUT,">$outfile_path");
  }

  my $feature_name = "GO";

  open(IN, "<$infile") || die($!);

  #my $header = <IN>;

  #my @col_names = split(/\t/, $header);

  print OUT ("Feature\t$feature_name\n");

  #will store gene-centered data
  my %data;

  while (my $line = <IN>) {
    chomp $line;   
    my @aux = split(/\t/, $line);
    my $gene_id = shift @aux;
    foreach my $element (@aux) {
      if ($element =~ /^GO:/) {
#     my $tmp = parse_GO($element);
        if (defined $data{$gene_id}) {
          $data{$gene_id} = join("|", $element, $data{$gene_id});
        } 
        else {
          $data{$gene_id} = $element;
        }
      }
    }
  }
  close IN;

  foreach my $key(keys %data) {
    my $unique = unique_values($data{$key});
    print OUT ("$key\t$unique\n");
  }

  sub unique_values {
    my $tmp = $_[0];
    my %tmp_data;
    my @aux = split(/\|/, $tmp);
    foreach my $element (@aux) {
      if (defined $tmp_data{$element}) {
 
      } 
      else {
        $tmp_data{$element} = 1;
      }
    }
    my @uniquearray = keys %tmp_data;
    my $tmp_out = join(";", @uniquearray);
    return $tmp_out;
  }
  close OUT;

}

sub feature2GO{
  my $id = $_[0];
  my $root_outdir = $_[1];
  my $in = "$root_outdir/interpro/$id.aa.fasta_longest.tsv";
  my @path = split(/\//, $in);
  my $outdir_path = "$root_outdir/komodo2_imputs/feature2GO";
  my $file_name = pop @path;

  my @tmp_name = split(/_/, $file_name);

  my $final_name = join("_", $tmp_name[0], $tmp_name[1]);

  $final_name =~ s/\./_/g;

  my $infile = $in;
  my $outfile_path = join("/", $outdir_path, $final_name);
  my $feature_name = "Pfam";
  chomp $feature_name;

  open(IN, "<$infile") || die($!);
  if (-e $outfile_path) {
    die "Outfile $outfile_path already exists.\n";
    }
  else {
    open(OUT,">$outfile_path");
  }

  #my $header = <IN>;

  #my @col_names = split(/\t/, $header);

  print OUT ("Feature\t$feature_name\n");

  #my %data;
  
  #my $total = 0;

  my $i = 0;

  while (my $line = <IN>) {
    chomp $line;   
    my @aux = split(/\t/, $line);
    my $actual_feature = $aux[3];
    if ((defined $actual_feature)&&($actual_feature eq $feature_name)) {
      my $GOs = get_GO($line);
      if (defined $GOs) { #new feature has GO terms; print them
        $GOs =~ s/\|/;/g;
        print OUT "$feature_name\_$i\t$GOs\n";
        $i++;
      } 
      else {
        next();
      }
    }
  }

  close IN;
  close OUT;
  sub get_GO {
    my $tmp = $_[0];
    my @aux = split(/\t/, $tmp);
    foreach my $element (@aux) {
      if ($element =~ /^GO:/) {
        return $element;
      }
    }
    return undef;
  }
}

sub parser_pfam {

  my $id = $_[0];
  my $root_outdir = $_[1];
  #interpro output file, tsv-delimited
  my $infile = "$root_outdir/interpro/$id.aa.fasta_longest.tsv";

  #which feature you want to summarize (e.g. "Pfam")
  my $feat_class = "Pfam";

  #output directory
  my $outdir_path = "$root_outdir/komodo2_imputs/pfam";
  chomp $outdir_path;

  #the bunch of code below is just to generate a suitable output file path
  my @interpro_path = split(/\//, $infile);

  my $file_name = pop @interpro_path;

  my @tmp_name = split(/_/, $file_name);

  my $final_name = join("_", $tmp_name[0], $tmp_name[1]);

  my $outfile_path = join("/", $outdir_path, $final_name);

  $outfile_path = join(".", $outfile_path, $feat_class, "txt");

  #print ("Generating file $outfile_path\n");

  if (-e $outfile_path) {
    die "Outfile $outfile_path already exists.\n";
  } 
  else {
    open(OUT,">$outfile_path");
  }

  my $feature_name = $feat_class;

  open(IN, "<$infile") || die($!);

  #my $header = <IN>;

  #my @col_names = split(/\t/, $header);

  print OUT ("Feature\t$feature_name\n");

  #my %data;

  #my $total = 0;

  my $i = 0;

  while (my $line = <IN>) {
    chomp $line;   
    my @aux = split(/\t/, $line);
    my $actual_feat_class = $aux[3];
    if ($actual_feat_class eq $feat_class) {
      my $actual_feature = $aux[4];
      if (defined $actual_feature) {
        print OUT ("$i\t$actual_feature\n");
      }
    $i++;
    }
  }

  close IN;
  close OUT;

}



1
