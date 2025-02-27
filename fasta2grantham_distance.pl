#!/bin/perl
use warnings;
use strict;

my $aaDistanceMatrix = "/home/owens/working/Sebastes/bin/AAdistMatrix_Grantham.cnv";
my $fastaFile = $ARGV[0];
my $reference_sample = $ARGV[1]; #Species to compare all others to.

# Read aa distances into hash by combining all possible aa pairs
my $aa1;
my $aa2;
my $tempAA;
my $numberOfAAs;
my @aaList;
my $distance;
my %aaPairwiseDistance = ();
my $problemFound = 0;
my $inputLine;
open (AAMATRIX, '<', "$aaDistanceMatrix") or die "$!";
$inputLine = <AAMATRIX>; 		# read and process header line
chomp $inputLine;
while ($inputLine =~m/\t/) {
	($tempAA, $inputLine) = split(/\t/, $inputLine, 2);
	push(@aaList, $tempAA);
}
push(@aaList, $inputLine);

foreach $tempAA (@aaList) {
}
while ($inputLine = <AAMATRIX>) {
	chomp $inputLine;
	$a = 0;
	($aa1, $inputLine) = split(/\t/, $inputLine, 2);
	while ($inputLine =~m/\t/) {
		$a = $a + 1;
		($distance, $inputLine) = split(/\t/, $inputLine, 2);
		$aa2 = $aa1 . $aaList[$a];
		$aaPairwiseDistance{$aa2} = $distance;
	}
	$a = $a + 1;
	$aa2 = $aa1 . $aaList[$a];
	$aaPairwiseDistance{$aa2} = $inputLine;
}
close AAMATRIX;

open (FASTA, '<', "$fastaFile") or die "Could not open input fasta file. $!";
my $ID;
my $sequence;
my $a;
my %SequenceArray = ();
while ($inputLine = <FASTA>) {
	chomp $inputLine;
	$inputLine =~ s/\r//g;
	if ($inputLine =~ m/^\s*$/) {next;}	#Ignore empty lines
	if ($inputLine =~ m/>/) {
		$sequence = "";
		($a, $ID) = split(/>/, $inputLine);
	}
	else {
		$sequence = $sequence . $inputLine;
		$SequenceArray{$ID} = $sequence;
	}
}
my @alleleList = sort (keys %SequenceArray);

# Calculate and output pairwise average distance between alleles
my $lengthAlleleArray = $#alleleList;
my $overallSequenceLength = length($SequenceArray{$alleleList[0]});
my $outputFilePairwise = $fastaFile . "_PairwiseDistanceList.txt";
my $sequenceLength;
my $distanceSum;
my $i;
my $j;
my $aaIndex;
my $problematicAA;

open (PAIRS, '>', "$outputFilePairwise") or die "$!\n";
print PAIRS "Species_1\tSpecies_2\tposition\tDistance\n";			# print header line

foreach my $sample (sort keys %SequenceArray){
  if ($sample eq $reference_sample){next;}
  foreach my $aaIndex (0..length($SequenceArray{$sample})){
print "$aaIndex\t$reference_sample\t$sample\n";
    $tempAA = substr($SequenceArray{$sample}, $aaIndex, 1) . substr($SequenceArray{$reference_sample}, $aaIndex, 1);
    if ($tempAA =~ m/[^ARNDCQEGHILKMFPSTWYV]+/) {
      $problemFound = 1;
      $problematicAA = $tempAA;
      $sequenceLength = $sequenceLength - 1;
    } else {
      $distance= $aaPairwiseDistance{$tempAA};
      my $position = $aaIndex+1;
      print PAIRS "\n$sample\t$reference_sample\t$position\t$distance";
    }
  }
}
