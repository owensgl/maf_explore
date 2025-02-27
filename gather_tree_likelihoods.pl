#!/bin/perl
use strict;
use warnings;

#Take a list of tree likelihoods and process it into one file. 
#Pipe in an ls of trees, eg PGA_scaffold_31__9.1407933.n.9999.trees
print "chr\tstart\tdirection\ttree_ID\trel_lik\tlik\tmax_lik";
while(<STDIN>){
  chomp;
  my @a = split(/\./,$_);
  my $chr = $a[0];
  my $start = $a[1];
  my $direction = $a[2];
  my $length = $a[3];
  open FILE, $_;
  my $max_lik = -99999999;
  my %liks;
  while(my $row = <FILE>){
    my @b = split(/ /,$row);
    my $tree_ID = $b[2];
    my $lik = $b[3];
    $lik =~ s/lh=//g;
    if ($lik > $max_lik){
      $max_lik = $lik;
    }
    $liks{$tree_ID} = $lik;
  }
  foreach my $tree_ID (sort keys %liks){
    my $rel_lik = $max_lik - $liks{$tree_ID};
    print "\n$chr\t$start\t$direction\t$tree_ID\t$rel_lik\t$liks{$tree_ID}\t$max_lik";
  }
  close FILE;
}


