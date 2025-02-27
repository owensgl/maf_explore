#!/bin/perl
use strict;
use warnings;
use POSIX;
use Statistics::Descriptive;
use Math::CDF;

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
  }elsif ($_ =~ m/^#outgroup/){
    my @tmp = split(/\t/,$_);
    $outgroup = $tmp[1];
    #Make the ABBA and BABA search terms 
    foreach my $i (0..$#samples){
      if ($samples[$i] eq $outgroup){next;}
      foreach my $j (0..$#samples){
        if ($samples[$j] eq $outgroup){next;}
        if ($i == $j){next;}
        foreach my $k (0..$#samples){
          if ($samples[$k] eq $outgroup){next;}
          if ($i == $k){next;}
          if ($j == $k){next;}
          #Match ABBA
          my $ABBA_string;
          foreach my $x (0..$#samples){
            if ($x == $i){
              $ABBA_string .= "A";
            }elsif($x ==$j){
              $ABBA_string .= "B";
            }elsif($x == $k){
              $ABBA_string .= "B";
            }elsif ($samples[$x] eq $outgroup){
              $ABBA_string .= "A";
            }else{
              $ABBA_string .= "\\w";
            }
          }
          my $BABA_string;
          foreach my $x (0..$#samples){
            if ($x == $i){
              $BABA_string .= "B";
            }elsif($x ==$j){
              $BABA_string .= "A";
            }elsif($x == $k){
              $BABA_string .= "B";
            }elsif($samples[$x] eq $outgroup){
              $BABA_string .= "A";
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
print "s1\ts2\ts3\ts4\tabba\tbaba\td\tjackknife_d\tstd_err\tzscore\tspecies_set";
#Print out total D stat and do the jackknife boostrap.
foreach my $i (sort keys %search_strings){
  foreach my $j (sort keys %{$search_strings{$i}}){
    foreach my $k (sort keys %{$search_strings{$i}{$j}}){
      my $abba = 0;
      my $baba = 0;
      my $jackknife_D;
      my @abba;
      my @baba;
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
          $abba += $abba_data{$i}{$j}{$k}{$chr}{$window};
          push(@abba,$abba_data{$i}{$j}{$k}{$chr}{$window});
        }
      }
      foreach my $chr (sort keys %{$baba_data{$i}{$j}{$k}}){
        foreach my $window (sort keys %{$baba_data{$i}{$j}{$k}{$chr}}){
          $baba += $baba_data{$i}{$j}{$k}{$chr}{$window};
          push(@baba,$baba_data{$i}{$j}{$k}{$chr}{$window});
        }
      }
      my $d = ($abba - $baba) / ($abba + $baba);
      my @d_jackknives;
      foreach my $x (0..$#abba){
        my $abba_tmp;
        my $baba_tmp;
        foreach my $y (0..$#abba){
          if ($x == $y){next;}
          $abba_tmp+=$abba[$y];
          $baba_tmp+=$baba[$y];
        }
        my $d_tmp = ($abba_tmp - $baba_tmp) / ($abba_tmp + $baba_tmp);
        push(@d_jackknives, $d_tmp);
      }
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(\@d_jackknives);
      my $mean_d = $stat->mean();
      my $var_d  = $stat->variance();
      my $stderr_d = sqrt($var_d * ($#d_jackknives+1));
      my $z = abs($mean_d)/$stderr_d;
      
      my $species_alphabetical;
      my %tmp_list;
      $tmp_list{$samples[$i]}++;
      $tmp_list{$samples[$j]}++;
      $tmp_list{$samples[$k]}++;
      my @ordered_list = sort keys %tmp_list;
      $species_alphabetical = "$ordered_list[0]:$ordered_list[1]:$ordered_list[2]";
        

      print "\n$samples[$i]\t$samples[$j]\t$samples[$k]\t$outgroup\t$abba\t$baba\t$d\t$mean_d\t$stderr_d\t$z\t$species_alphabetical";
    }
  }
}



