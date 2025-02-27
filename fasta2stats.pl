#!/bin/perl
use strict;
use warnings;
use POSIX;

#outputs fasta stats
my $filename = $ARGV[0];
#Number of sites without any Ns
#Number of variable sites
#Number of sites with any masking.
my %n;
my %masked;
my %masked_bases;
$masked_bases{"a"}++;
$masked_bases{"c"}++;
$masked_bases{"t"}++;
$masked_bases{"g"}++;
my %base;
$base{"A"}++;
$base{"T"}++;
$base{"C"}++;
$base{"G"}++;
my $length;
my %data;
my %counts;
my $samples;
open FILE, $filename;
while(<FILE>){
  chomp;
  if ($_ =~ m/^>/){
    $samples++;
    next;
  }
  my @a = split(//,$_);
  foreach my $i (0..$#a){
    $length=$#a+1;
    if ($a[$i] eq "N"){
      $n{$i}++;
      next;
    }
    if ($a[$i] eq "-"){next;}
    if ($masked_bases{$a[$i]}){
      $masked{$i}++;
    }
    if($base{uc($a[$i])}){
      print STDERR "$a[$i]\n";
      $data{$i}{uc($a[$i])}++;
    }
    $counts{$i}++
  }
}
my $N_bases;
my $total_masked_bases;
my $variable_bases;
my $full_data;
foreach my $i (0..$length){
  if ($n{$i}){
    $N_bases++;
  }
  if ($masked{$i}){
    $total_masked_bases++;
  }
  if($data{$i}){
    my $alleles = keys %{$data{$i}};
    if ($alleles >1){
      $variable_bases++;
    }
  }
  unless($counts{$i}){
    $counts{$i} =0;
  }
  if ($counts{$i} == $samples){
    $full_data++;
  }
}
close FILE;
print "$filename\t$full_data\t$N_bases\t$total_masked_bases\t$variable_bases\t$length";
