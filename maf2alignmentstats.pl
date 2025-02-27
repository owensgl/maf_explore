#!/bin/perl
use strict;
use warnings;
use POSIX;
 
#This takes a MAF. It looks for how many bases are good or have duplicates or are masked or are gaps (In that order). It uses one genome (outgroup) as the reference to keep track of position.
#Check for duplicates
#Checkf for missing copy
#Check for deletions in bases
#Check for masked bases
#Otehrwise its good.

my $window_size = 10000;
my $outgroup = $ARGV[0]; #Genome to use as reference

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
my @region_start;
my @region_direction;
my @region_chr;
my @region_window;
my @region_length;
my %data;
my $seq_length;
my %matches;
my %length;
my %all_regions;
my %good_bases;
my %dup_bases;
my %del_bases;
my %dup_sample;
my %del_sample;
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
    my $current_win = floor($start/$window_size);
    my $ref_length = $a[3];
    my $direction = $a[4];
    if ($sample eq $outgroup){
      push(@region_start, $start);
      push(@region_window,$current_win);
      push(@region_chr, $chr);
      push(@region_direction, $direction);
      push(@region_length,$ref_length);
      $all_regions{$chr}{$current_win}++;
    
    }
  }
  #If it starts with a, then process sequences and reset hashes.
  if ($_ =~ m/^a/ ){
    my %duplicates;
    my %deletions;
    $sets_processed++;
    if ($sets_processed % 10000 == 0){print STDERR "Processed $sets_processed alignments...\n";}
    #Check to see if each sample is only represented once.
    unless($sample_counts{$outgroup}){
      goto SKIP;
    }
    foreach my $sample (@samples){
      unless($sample_counts{$sample}){
        $deletions{$sample}++;
      }elsif ($sample_counts{$sample} != 1){
        $duplicates{$sample}++;
      }
    }
    
    #If no duplications or deletions, record number of good bases;
    my $n_dup = keys %duplicates;
    my $n_del = keys %deletions;
    if (($n_dup == 0) and ($n_del == 0)){
      foreach my $i (0..$#region_start){
        $good_bases{$region_chr[$i]}{$region_window[$i]}+=$region_length[$i];
      }
    }elsif ($n_dup > 0){
      foreach my $i (0..$#region_start){
        $dup_bases{$region_chr[$i]}{$region_window[$i]}+=$region_length[$i];
      }
    }elsif ($n_del > 0){
      foreach my $i (0..$#region_start){
        $del_bases{$region_chr[$i]}{$region_window[$i]}+=$region_length[$i];
      }
    }
    foreach my $sample (sort keys %duplicates){
      foreach my $i (0..$#region_start){
        $dup_sample{$region_chr[$i]}{$region_window[$i]}{$sample}+=$region_length[$i];
      }
    }
    foreach my $sample (sort keys %deletions){
      foreach my $i (0..$#region_start){
        $del_sample{$region_chr[$i]}{$region_window[$i]}{$sample}+=$region_length[$i];
      }
    }
    SKIP:
    #Refresh all the global variables that need changing every alignment block.
    undef(%sample_counts);
    undef(@region_start);
    undef(@region_direction);
    undef(@region_chr);
    undef(@region_window);
    undef(@region_length);
    undef(%data);
  }
}

my %duplicates;
my %deletions;
#Check to see if each sample is only represented once.
unless($sample_counts{$outgroup}){
  goto SKIP2;
}
foreach my $sample (@samples){
  unless($sample_counts{$sample}){
    $deletions{$sample}++;
  }elsif ($sample_counts{$sample} != 1){
    $duplicates{$sample}++;
  }
}
    
#If no duplications or deletions, record number of good bases;
my $n_dup = keys %duplicates;
my $n_del = keys %deletions;
if (($n_dup == 0) and ($n_del == 0)){
  foreach my $i (0..$#region_start){
    $good_bases{$region_chr[$i]}{$region_window[$i]}+=$region_length[$i];
   }
}elsif ($n_dup > 0){
  foreach my $i (0..$#region_start){
    $dup_bases{$region_chr[$i]}{$region_window[$i]}+=$region_length[$i];
  }
}elsif ($n_del > 0){
  foreach my $i (0..$#region_start){
    $del_bases{$region_chr[$i]}{$region_window[$i]}+=$region_length[$i];
  }
}
foreach my $sample (sort keys %duplicates){
  foreach my $i (0..$#region_start){
    $dup_sample{$region_chr[$i]}{$region_window[$i]}{$sample}+=$region_length[$i];
  }
}
foreach my $sample (sort keys %deletions){
  foreach my $i (0..$#region_start){
    $del_sample{$region_chr[$i]}{$region_window[$i]}{$sample}+=$region_length[$i];
  }
}
 


SKIP2:
########
print "#samples:\t$samples[0]";
foreach my $sample (1..$#samples){
  print ",$samples[$sample]";
}
print "\n#outgroup:\t$outgroup";
print "\n#windowsize:\t$window_size";
print "\nchr\twindow\ttype\tsample\tbases";
#Now print out all set counts 
foreach my $chr (sort keys %all_regions){
  foreach my $window (sort {$a <=> $b} keys %{$all_regions{$chr}}){
    my $real_window = $window * $window_size;
    unless ($good_bases{$chr}{$window}){
      $good_bases{$chr}{$window} = 0;
    }
    unless ($dup_bases{$chr}{$window}){
      $dup_bases{$chr}{$window} = 0;
    }
    unless ($del_bases{$chr}{$window}){
      $del_bases{$chr}{$window} = 0;
    }
    print "\n$chr\t$real_window\tsingle_copy\tNA\t$good_bases{$chr}{$window}";
    print "\n$chr\t$real_window\tdup_copy\tNA\t$dup_bases{$chr}{$window}";
    print "\n$chr\t$real_window\tdeleted_copy\tNA\t$del_bases{$chr}{$window}";
    foreach my $sample (@samples){
      unless ($dup_sample{$chr}{$window}{$sample}){
        $dup_sample{$chr}{$window}{$sample} = 0;
      }
      unless ($del_sample{$chr}{$window}{$sample}){
        $del_sample{$chr}{$window}{$sample} = 0;
      }
      print "\n$chr\t$real_window\tdup_copy_sample\t$sample\t$dup_sample{$chr}{$window}{$sample}";
      print "\n$chr\t$real_window\tdeleted_copy_sample\t$sample\t$del_sample{$chr}{$window}{$sample}";
    }
  }
}
