#!/bin/perl
use strict;
use warnings;
use POSIX;
 

#This takes a MAF output from cactus2 and sums up counts of patterns in windows. It requires all samples to have data (i.e. not N, - or lower case). It uses a window based on the outgroup, specified in the input. It doesn't output invariant sites.

my $window_size = 1000000;
my $outgroup = $ARGV[0];

#First, get sample names from the starting phylogeny. Ignore samples starting with Anc
my %samples;
my @samples;
my $sets_processed = 0;
#Bases to keep
my %usable_bases;
$usable_bases{"A"}++;
$usable_bases{"T"}++;
$usable_bases{"C"}++;
$usable_bases{"G"}++;

my %sets;
#Variables refreshed every set;
my %sample_counts;
my $region_start;
my $region_direction;
my $region_chr;
my %data;
my $seq_length;
while (<STDIN>){
  chomp;
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^# hal/){
    my $line = $_;
    $line =~ s/# hal //g;
    my @fields = split(/,/,$line);
    foreach my $field (@fields){
      my @parts = split(/:/,$field);
      my $name = $parts[0];
      $name =~ s/\(//g;
      $name =~ s/\)//g;
      if ($name =~ m/^Anc/){next;}
      $samples{$name}++;
    }
    foreach my $sample (sort keys %samples){
      push(@samples, $sample);
    }
    next;
  }
  #If it starts with s, then load in the sequences.
  if ($_ =~ m/^s/ ){
    my @a = split(/\t/,$_);
    my $sample_chr = $a[1];
    if ($sample_chr =~ m/^Anc/){next;}
    my @sample_chr = split(/\./,$sample_chr);
    my $sample = $sample_chr[0];
    $sample_counts{$sample}++;
    my $chr = $sample_chr;
    $chr =~ s/${sample}.//g;
    my $start = $a[2];
    my $direction = $a[4];
    my $max_length = $a[5];
    if ($direction eq "-"){
      $start = $max_length - $start;
    }
    if ($sample eq $outgroup){
      $region_start = $start;
      $region_chr = $chr;
      $region_direction = $direction;
    }
    #Load in all bases;
    my @bases = split(//,$a[6]);
    foreach my $i (0..$#bases){
      $data{$sample}{$i} = $bases[$i];
    }
    $seq_length = $#bases;
  }
  #If it starts with a, then process sequences and reset hashes.
  if ($_ =~ m/^a/ ){
    $sets_processed++;
    if ($sets_processed % 10000 == 0){print STDERR "Processed $sets_processed alignments...\n";}
    #Check to see if each sample is only represented once.
    foreach my $sample (@samples){
      unless($sample_counts{$sample}){
        goto REFRESH;
      }
      if ($sample_counts{$sample} != 1){
        goto REFRESH;
      }
    }
    my $position = $region_start- 1; #Minus 1 because it increments up before determining window.
    foreach my $i (0..$seq_length){
      #Keep track of the current position in the outgroup reference genome.
      my $ancestral_allele = $data{$outgroup}{$i};
      if ($ancestral_allele ne "-"){
        $position++;
      }
      my $window = floor($position/$window_size);

      
      #Check if each sample has a called genotype, and phase it as ancestral or derived.
      my %total_alleles;
      my %phased_genotypes;
      foreach my $sample (@samples){
        unless ($usable_bases{$data{$sample}{$i}} ){
          goto NEXTBASE;
        }
        $total_alleles{$data{$sample}{$i}}++;
        if ($data{$sample}{$i} eq $ancestral_allele){
          $phased_genotypes{$sample} = "A";
        }else{
          $phased_genotypes{$sample} = "B";
        }
      }
      #Check that there are two alleles in all samples;
      my $total_alleles = keys %total_alleles;
      if ($total_alleles != 2){
        next;
      }
      my $set;
      foreach my $sample (@samples){
        $set.=$phased_genotypes{$sample};
      }
      $sets{$region_chr}{$window}{$set}++;
      NEXTBASE:
    }
    
    REFRESH:
    #Refresh all the global variables that need changing every alignment block.
    undef(%sample_counts);
    undef($region_start);
    undef($region_direction);
    undef($region_chr);
    undef($seq_length);
    undef(%data);
  }
}

#For last alignment in file process it
########
#Check to see if each sample is only represented once.
foreach my $sample (@samples){
  unless($sample_counts{$sample}){
    goto MOVEON;
  }
  if ($sample_counts{$sample} != 1){
    goto MOVEON;
  }
}
my $position = $region_start- 1; #Minus 1 because it increments up before determining window.
foreach my $i (0..$seq_length){
  #Keep track of the current position in the outgroup reference genome.
  my $ancestral_allele = $data{$outgroup}{$i};
  if ($ancestral_allele ne "-"){
    $position++;
  }
  my $window = floor($position/$window_size);
  
  
  #Check if each sample has a called genotype, and phase it as ancestral or derived.
  my %total_alleles;
  my %phased_genotypes;
  foreach my $sample (@samples){
    unless ($usable_bases{$data{$sample}{$i}} ){
      goto SKIPBASE;
    }
    $total_alleles{$data{$sample}{$i}}++;
    if ($data{$sample}{$i} eq $ancestral_allele){
      $phased_genotypes{$sample} = "A";
    }else{
      $phased_genotypes{$sample} = "B";
    }
  }
  #Check that there are two alleles in all samples;
  my $total_alleles = keys %total_alleles;
  if ($total_alleles != 2){
    next;
  }
  my $set;
  foreach my $sample (@samples){
    $set.=$phased_genotypes{$sample};
  }
  $sets{$region_chr}{$window}{$set}++;
  SKIPBASE:
}
MOVEON:
########
print "#samples:\t$samples[0]";
foreach my $sample (1..$#samples){
  print ",$samples[$sample]";
}
print "\n#outgroup:\t$outgroup";
print "\n#windowsize:\t$window_size";
print "\nchr\twindow\tset\tcount";
#Now print out all set counts 
foreach my $chr (sort keys %sets){
  foreach my $window (sort {$a <=> $b} keys %{$sets{$chr}}){
    foreach my $set (sort keys %{$sets{$chr}{$window}}){
      print "\n$chr\t$window\t$set\t$sets{$chr}{$window}{$set}";
    }
  }
}
