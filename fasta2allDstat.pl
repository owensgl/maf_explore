#!/bin/perl
use strict;
use warnings;

#This takes a fasta file and calculates all possible D statistics with one set outgroup
my $outgroup = $ARGV[0];
my $current_base;
my $current_name;
my %data;
my %site_counts;
my $total_samples;
my %possible_bases;
$possible_bases{"A"}++;
$possible_bases{"T"}++;
$possible_bases{"G"}++;
$possible_bases{"C"}++;
print STDERR "Loading in fasta sequence...";
while(<STDIN>){
  chomp;
  if ($_ =~ m/>/){
    #It's a header line, so read in the name
    my $name = $_;
    $name =~ s/>//g;
    $current_name = $name;
    $current_base = 0;
    $total_samples++;
  }else{
    #It's a sequence row so read it in. 
    my @a = split(//,$_);
    foreach my $i (0..$#a){
      $data{$current_name}{$current_base} = $a[$i];
      $current_base++;
    }
  }
}
#Check the possible number of combinations
my $total_comparisons = ($total_samples-1)*($total_samples-1)*($total_samples-1);

print STDERR "\n$total_comparisons possible sets";
my $tested_sets = 0;
print STDERR "\nAdding D stats together";
#Now do all possible D statistics
foreach my $sample1 (sort keys %data){
  foreach my $sample2 (sort keys %data){
    foreach my $sample3 (sort keys %data){
      $tested_sets++;
      if ($tested_sets % 1000 == 0){print STDERR "\n$tested_sets";}
      if ($sample1 eq $sample2){next;}
      if ($sample1 eq $sample3){next;}
      if ($sample3 eq $sample2){next;}
      if ($sample1 eq $outgroup){next;}
      if ($sample2 eq $outgroup){next;}
      if ($sample3 eq $outgroup){next;}
      my @samples = ($sample1, $sample2, $sample3, $outgroup);
      foreach my $i (0..($current_base-1)){
        my %tmp_bases;
        my %tmp_oriented;
        foreach my $n (0..3){
          $tmp_bases{$n} = uc($data{$samples[$n]}{$i});
          unless($possible_bases{$tmp_bases{$n}}){goto NEXTLINE;}
        }
        foreach my $n (0..2){
          if ($tmp_bases{$n} eq $tmp_bases{3}){
            $tmp_oriented{$n} = "A";
          }else{
            $tmp_oriented{$n} = "B";
          }
        }
        #Check if abba or baba exists
        my $state = "$tmp_oriented{0}$tmp_oriented{1}$tmp_oriented{2}";
        
        if ($state eq "ABB"){
          #Check if derived alleles match
          if ($tmp_bases{1} eq $tmp_bases{2}){
            $site_counts{$sample1}{$sample2}{$sample3}{"ABBA"}++;
          }
        }

        if ($state eq "BAB"){
          #Check if derived alleles match
          if ($tmp_bases{0} eq $tmp_bases{2}){
            $site_counts{$sample1}{$sample2}{$sample3}{"BABA"}++;
          }
        }

        NEXTLINE:
      }
    }
  }
}

foreach my $sample1 (sort keys %data){
  foreach my $sample2 (sort keys %data){
    foreach my $sample3 (sort keys %data){
      if ($sample1 eq $sample2){next;}
      if ($sample1 eq $sample3){next;}
      if ($sample3 eq $sample2){next;}
      if ($sample1 eq $outgroup){next;}
      if ($sample2 eq $outgroup){next;}
      if ($sample3 eq $outgroup){next;}
      unless($site_counts{$sample1}{$sample2}{$sample3}{"ABBA"}){
        $site_counts{$sample1}{$sample2}{$sample3}{"ABBA"} = 0;
      }
      unless($site_counts{$sample1}{$sample2}{$sample3}{"BABA"}){
        $site_counts{$sample1}{$sample2}{$sample3}{"BABA"} = 0;
      }
      print "\n$sample1\t$sample2\t$sample3\t$outgroup\t";
      print "$site_counts{$sample1}{$sample2}{$sample3}{'ABBA'}\t";
      print "$site_counts{$sample1}{$sample2}{$sample3}{'BABA'}\t";
    }
  }
}

