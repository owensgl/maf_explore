#!/bin/perl
use strict;
use warnings;
use POSIX;
#this takes a maf with one sequence per species per alignment and outputs all paired alignment loctions
my %data;
my $target_species = "Sebastes_schlegelii";
my $min_size = 50;
print "sample_1\tchr_1\tstart_1\tlength_1\tsample_2\tchr_2\tstart_2\tlength_2";
while (<STDIN>){
  chomp;
  if ($_ =~ m/^#/){next;}
  if ($_ =~ m/^s/ ){
    my @a = split(/\s+/,$_);
    my $sample_chr = $a[1];
    if ($sample_chr =~ m/^Anc/){next;}
    my @sample_chr = split(/\./,$sample_chr);
    my $sample = $sample_chr[0];
    my $chr = $sample_chr;
    $chr =~ s/${sample}.//g;
    my $start = $a[2];
    my $ref_length = $a[3];
    my $direction = $a[4];
    my $max_length = $a[5];
    if ($direction eq '-'){
      $start = $max_length - $start;
    }
    $data{$sample}{"chr"} = $chr;
    $data{$sample}{"start"} = $start;
    $data{$sample}{"length"} = $ref_length;
  }  if ($_ =~ m/^a/ ){
    unless(scalar(keys(%data)) > 0){next;}
      my $sample1 = $target_species;
      foreach my $sample2 (sort keys %data){
        if ($sample1 eq $sample2){next;}
        if ($data{$sample1}{'length'} < $min_size){next;}
        print "\n$sample1\t$data{$sample1}{'chr'}\t$data{$sample1}{'start'}\t$data{$sample1}{'length'}\t";
        print "$sample2\t$data{$sample2}{'chr'}\t$data{$sample2}{'start'}\t$data{$sample2}{'length'}";
      }
   
    undef(%data);
  }
}


