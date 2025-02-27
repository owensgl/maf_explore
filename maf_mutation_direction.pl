#!/bin/perl
use strict;
use warnings;
use POSIX;

#Command for processing filtered mafs

#../../bin/mafOrder rockfish_highqual_8genomes.filt.unique.maf species_8.txt /dev/stdout | ../../bin/mafTools/bin/mafStrander --maf /dev/stdin --seq Sebastes_aleutianus --strand + | ../../bin/phast/bin/maf_parse -o MAF -g ../hetsites/ref_genomes/tmp.hetlist.txt  /dev/stdin

my $hetfile = $ARGV[0];

print STDERR "Loading site lists...\n";
open(HETSITES, "gunzip -c $hetfile |");
my %ref;
my %alt;
my %before;
my %after;
while(<HETSITES>){
  chomp;
  my @a = split(/\t/,$_);
  my $chr = $a[2];
  my $pos = $a[3];
  my $ref = $a[4];
  my $alt = $a[5];
  my @context = split(//,$a[6]);
  $ref{$chr}{$pos} = $ref;
  $alt{$chr}{$pos} = $alt;
  $before{$chr}{$pos} = uc($context[0]);
  $after{$chr}{$pos} = uc($context[2]);
} 
close HETSITES;

print STDERR "Finished loading sites...\n";
print STDERR "Processing MAF...\n";
my $matching_site;
my %current_counts;
my $current_chr;
my @current_pos;
my %sites;
my $counter = 0;
while (<STDIN>){
  chomp;
  my $line = $_;
  if ($line =~ m/^a/){
    #Print results from previous sets;
    if (%current_counts){
      foreach my $chr (keys %current_counts){
        foreach my $pos (keys %{$current_counts{$chr}}){
          unless($current_counts{$chr}{$pos}{$ref{$chr}{$pos}}){
            $current_counts{$chr}{$pos}{$ref{$chr}{$pos}} = 0;
          }
          unless($current_counts{$chr}{$pos}{$alt{$chr}{$pos}}){
            $current_counts{$chr}{$pos}{$alt{$chr}{$pos}} = 0;
          }
          #check if both alleles are present in relatives
          if (($current_counts{$chr}{$pos}{$ref{$chr}{$pos}} > 0) and ($current_counts{$chr}{$pos}{$alt{$chr}{$pos}} > 0)){
            #print "\n$chr\t$pos\ttranspolymorphism\tNA";
            $sites{"TR"}++;
          }
          #If ref is in relatives
          elsif ($current_counts{$chr}{$pos}{$ref{$chr}{$pos}} > 0){
            #print "\n$chr\t$pos\t$ref{$chr}{$pos}:$alt{$chr}{$pos}\t$current_counts{$chr}{$pos}{$ref{$chr}{$pos}}";
            $sites{"$before{$chr}{$pos}$ref{$chr}{$pos}$after{$chr}{$pos}:$before{$chr}{$pos}$alt{$chr}{$pos}$after{$chr}{$pos}"}++;

          }
          #If alt is in relatives
          elsif ($current_counts{$chr}{$pos}{$alt{$chr}{$pos}} > 0){
            #print "\n$chr\t$pos\t$alt{$chr}{$pos}:$ref{$chr}{$pos}\t$current_counts{$chr}{$pos}{$alt{$chr}{$pos}}";
            $sites{"$before{$chr}{$pos}$alt{$chr}{$pos}$after{$chr}{$pos}:$before{$chr}{$pos}$ref{$chr}{$pos}$after{$chr}{$pos}"}++;

          }else{
            #print "\n$chr\t$pos\tneither\tNA";
            $sites{"NA"}++;

          }
          $counter++;
          if ($counter % 10000 == 0){
            print STDERR "$hetfile $counter sites processed..\n";
          }
        }
      }
    }
    undef($matching_site);
    undef(%current_counts);
    undef($current_chr);
    undef(@current_pos);
  }
  if ($line =~ m/^s/){
    unless($matching_site){
      my @a = split(' ',$line);
      my $chr = $a[1];
      unless($ref{$chr}){next;}
      my $start = $a[2];
      my @bases = split(//,$a[6]);
      foreach my $i (0..$#bases){
        my $pos = $start + $i + 1;
        if ($ref{$chr}{$pos}){
          $matching_site++;
          $current_chr = $chr;
          $current_pos[$i] = $pos;
        }
      }
      next;
    }else{
      my @a = split(' ',$line);
      my @bases = split(//,$a[6]);
      foreach my $i (0..$#bases){
        unless($current_pos[$i]){next;}
        if ($ref{$current_chr}{$current_pos[$i]}){
          if ($bases[$i] eq "-"){next;}
          if ($bases[$i] eq "N"){next;}
          $bases[$i] = uc($bases[$i]);
          $current_counts{$current_chr}{$current_pos[$i]}{$bases[$i]}++;
          # print "\n$current_chr\t$current_pos[$i]\t$ref{$current_chr}{$current_pos[$i]}:$alt{$current_chr}{$current_pos[$i]}\t";
          # print "$bases[$i]";
        }
      }
    }

  }
}
print "site_type\tn";
foreach my $type (sort keys %sites){
  print "\n$type\t$sites{$type}";
}
