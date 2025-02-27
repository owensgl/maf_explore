#!/bin/perl
use strict;
use warnings;
use POSIX;

#Converts maf to fasta format, one file per alignment block

my $reference_species = "Sebastes_aleutianus"; # Used to label file names.
my %sequences;
my $ref_name;
while(<STDIN>){
  if ($_ =~ /^#/){next;}
  if ($_ =~ /^s/){
    #Load in sequences.
    my @a = split(/\s+/,$_);
    my @scaffold_info = split(/\./,$a[1]);
    my $start = $a[2];
    my $length = $a[3];
    my $dir = $a[4];
    my $max_length = $a[5];
    my $species = $scaffold_info[0];
    if ($species eq $reference_species){
      my $scaffold = $scaffold_info[1];
      if ($scaffold_info[2]){
        foreach my $i (2..$#scaffold_info){
          $scaffold .= $scaffold_info[$i];
        }
      }
      my $pos = $start;
#      if ($dir eq "-"){
#        $pos = $max_length - $start;
#      }
      my $dir_code = "p";
      if ($dir eq "-"){
        $dir_code = "n";
      }
      $ref_name = "$scaffold.$pos.$dir_code.$length";
    }
    my $seq = $a[6];
    $sequences{$species} = $seq;
  }elsif ($_ =~ /^a/){
    unless ($ref_name){next;}
print STDERR "Printing: $ref_name\n";
    open (my $file_name, '>', "$ref_name.fa");
    my $count;
    foreach my $species (sort keys %sequences){
      if ($count){
        print $file_name "\n>$species\n$sequences{$species}";
      }else{
        print $file_name ">$species\n$sequences{$species}";
        $count++;
      }
    }
    close $file_name;
    undef(%sequences);
    undef($ref_name);
  }
}
unless($ref_name){
  exit;
}
open (my $file_name, '>', "$ref_name.fa");
my $count;
foreach my $species (sort keys %sequences){
  if ($count){
    print $file_name "\n$species\n$sequences{$species}";
  }else{
    print $file_name "$species\n$sequences{$species}";
    $count++;
  }
}
close $file_name;
