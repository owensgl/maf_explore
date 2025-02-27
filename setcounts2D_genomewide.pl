#!/bin/perl
use strict;
use warnings;
use POSIX;

#This script takes the output of maf2setcounts.pl and does all D stat tests including jackknifing with set block size.
my @samples;
my $outgroup;
my %search_strings;
my %abba_data;
my %baba_data;


while(<STDIN>){
  chomp;
  if ($_ =~ m/^#samples/){
    my @tmp = split(/\t/,$_);
    my @tmp_samples = split(/,/,$tmp[1]);
    @samples = @tmp_samples;
    $outgroup = $samples[0];

    #Make the ABBA and BABA search terms 
    foreach my $i (1..$#samples){
      foreach my $j (1..$#samples){
        if ($i == $j){next;}
        foreach my $k (1..$#samples){
          if ($i == $k){next;}
          if ($j == $k){next;}
          #Match ABBA
          my $ABBA_string = "A";
          foreach my $x (1..$#samples){
            if ($x == $i){
              $ABBA_string .= "A";
            }elsif($x ==$j){
              $ABBA_string .= "B";
            }elsif($x == $k){
              $ABBA_string .= "B";
            }else{
              $ABBA_string .= "\\w";
            }
          }
          my $BABA_string = "A";
          foreach my $x (1..$#samples){
            if ($x == $i){
              $BABA_string .= "B";
            }elsif($x ==$j){
              $BABA_string .= "A";
            }elsif($x == $k){
              $BABA_string .= "B";
            }else{
              $BABA_string .= "\\w";
            }
          }
          $search_strings{$i}{$j}{$k}{"ABBA"}= $ABBA_string;
          $search_strings{$i}{$j}{$k}{"BABA"}= $BABA_string;
        }
      }
    }

  }
  if ($_ =~ m/^#/){next;}
  if ($_ =~ m/^chr\t/){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $window = $a[1];
  my $set = $a[2];
  my $count = $a[3];
  #Go through each possible comparison of samples and count ABBA and baba
  foreach my $i (sort keys %search_strings){
    foreach my $j (sort keys %{$search_strings{$i}}){
      foreach my $k (sort keys %{$search_strings{$i}{$j}}){
        my $search_string_abba = $search_strings{$i}{$j}{$k}{'ABBA'};
        my $search_string_baba = $search_strings{$i}{$j}{$k}{'BABA'};
        if ($set =~ /$search_string_abba/){
          $abba_data{$i}{$j}{$k}{$chr}{$window} += $count;
        }elsif ($set =~ /$search_string_baba/){
          $baba_data{$i}{$j}{$k}{$chr}{$window} += $count;
        }
      }
    }
  }
}
print "s1\ts2\ts3\ts4\tchr\twindow\tabba\tbaba\tspecies_set";
#Print out total D stat and do the jackknife boostrap.
foreach my $i (sort keys %search_strings){
  foreach my $j (sort keys %{$search_strings{$i}}){
    foreach my $k (sort keys %{$search_strings{$i}{$j}}){
      my $abba = 0;
      my $baba = 0;
      my $jackknife_D;
      my @abba;
      my @baba;
      my $species_alphabetical;
      my %tmp_list;
      $tmp_list{$samples[$i]}++;
      $tmp_list{$samples[$j]}++;
      $tmp_list{$samples[$k]}++;
      my @ordered_list = sort keys %tmp_list;
      $species_alphabetical = "$ordered_list[0]:$ordered_list[1]:$ordered_list[2]";
      #Make sure the keys for abba and babs are filled in with zeros
      foreach my $chr (sort keys %{$abba_data{$i}{$j}{$k}}){
        foreach my $window (sort keys %{$abba_data{$i}{$j}{$k}{$chr}}){
          unless($abba_data{$i}{$j}{$k}{$chr}{$window}){
            $abba_data{$i}{$j}{$k}{$chr}{$window} = 0;
          }
          unless ($baba_data{$i}{$j}{$k}{$chr}{$window}){
            $baba_data{$i}{$j}{$k}{$chr}{$window} = 0;
          }
        }
      }
      foreach my $chr (sort keys %{$baba_data{$i}{$j}{$k}}){
        foreach my $window (sort keys %{$baba_data{$i}{$j}{$k}{$chr}}){
          unless ($baba_data{$i}{$j}{$k}{$chr}{$window}){
            $baba_data{$i}{$j}{$k}{$chr}{$window} = 0;
          }
          unless($abba_data{$i}{$j}{$k}{$chr}{$window}){
            $abba_data{$i}{$j}{$k}{$chr}{$window} = 0;
          }
        }
      }
      foreach my $chr (sort keys %{$abba_data{$i}{$j}{$k}}){
        foreach my $window (sort keys %{$abba_data{$i}{$j}{$k}{$chr}}){
          print "\n$samples[$i]\t$samples[$j]\t$samples[$k]\t$outgroup\t$chr\t$window\t$abba_data{$i}{$j}{$k}{$chr}{$window}\t$baba_data{$i}{$j}{$k}{$chr}{$window}\t$species_alphabetical";
        }
      }
    }
  }
}

