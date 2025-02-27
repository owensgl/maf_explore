#!/bin/perl
use strict;
use warnings;
use POSIX;


#Count the number of Ns in windows for a reference genome in fasta

my $window_size = 10000;
my %data;
my %total_bases;
my $current_chr;
my $current_pos;
my $current_window;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^>/){
    my $chr = $_;
    $chr =~ s/>//;
    $current_chr =$chr;
    $current_pos = 0;
  }else{
    my @a = split(//,$_);
    foreach my $i (0..$#a){
      $current_pos++;
      $current_window = floor($current_pos/$window_size);
      if ($a[$i] eq "N"){
        $data{$current_chr}{$current_window}++;
      }
      $total_bases{$current_chr}{$current_window}++;
    }
  }
}
print "chr\twindow\tbases\tNs";
foreach my $chr (sort keys %total_bases){
  foreach my $window (sort keys %{$total_bases{$chr}}){
    unless($data{$chr}{$window}){
      $data{$chr}{$window} = 0;
    }
    print "\n$chr\t$window\t$total_bases{$chr}{$window}\t$data{$chr}{$window}";
  }
}
