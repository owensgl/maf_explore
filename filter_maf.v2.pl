#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#This script takes a maf file and filters blocks based on requiring that each sample be present in one or two copies, from the two haplotypes. You give it a list of sample names. For samples with two haplotypes, there should be two lines in the list file (one for each haplotype) It looks at each block and asks how many copies of each haplotype are present. If neither haplotype is present then it skips. If either haplotype has more than one copy, then it skips. Otherwise it will output the block as is. The log file gives you some information about why it skipped each block (i.e. who was missing or had too many copies) and it also gives some final counts on how many blocks were retained.

#In relaxed mode it keeps a block if one haplotype is present but the other is missing. If relaxed mode is off, then it requires each haplotype to be present (and single copy).
my $relaxed_mode = 1;

# Usage check
if (@ARGV < 3) {
    die "Usage: $0  <maf_file> <species.haplotype_list_file> <log_file>\n";
}

my $maf_file = $ARGV[0];
my $list_file = $ARGV[1];
my $log_file = $ARGV[2];

# Open log file for writing
open(my $log_fh, '>', $log_file) or die "Could not open log file '$log_file': $!\n";

# Read the list of desired species.haplotype strings
my @desired_sp_hap = ();
my %species_haps = ();  # For relaxed mode: track haplotypes by species

open(my $list_fh, '<', $list_file) or die "Could not open list file '$list_file': $!\n";
while (my $line = <$list_fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;  # Skip empty lines
    push @desired_sp_hap, $line;

    # For relaxed mode: build mapping of species to their haplotypes
    if ($relaxed_mode) {
        my ($species, $haplotype) = split(/\./, $line, 2);
        push @{$species_haps{$species}}, $haplotype if defined $species && defined $haplotype;
    }
}
close($list_fh);

if (!@desired_sp_hap) {
    die "No species.haplotype strings found in list file.\n";
}

# Print mode and filtering criteria
my $mode_str = $relaxed_mode ? "RELAXED MODE" : "STRICT MODE";
print $log_fh "$mode_str - Filtering for blocks with ";
if ($relaxed_mode) {
    print $log_fh "at least one haplotype per species from: ";
} else {
    print $log_fh "exactly one occurrence of each: ";
}
print $log_fh join(", ", @desired_sp_hap), "\n";

# Open the MAF file (handles both regular and gzipped files)
my $maf_fh;
if ($maf_file =~ /\.gz$/) {
    # Open gzipped file
    open($maf_fh, '-|', "gzip -dc $maf_file") or die "Could not open gzipped MAF file '$maf_file': $!\n";
} else {
    # Open regular file
    open($maf_fh, '<', $maf_file) or die "Could not open MAF file '$maf_file': $!\n";
}

# Variables to track current block
my @current_block = ();
my %sp_hap_counts = ();
my %species_counts = ();  # For relaxed mode
my $in_block = 0;
my $block_number = 0;
my $blocks_kept = 0;
my $blocks_filtered = 0;


print "##maf   version=1\n\n";
while (my $line = <$maf_fh>) {
    # Add line to current block
    push @current_block, $line;

    # Start of a new alignment block
    if ($line eq "a\n") {
        $in_block = 1;
        %sp_hap_counts = map { $_ => 0 } @desired_sp_hap;
        %species_counts = () if $relaxed_mode;
        $block_number++;
        next;
    }

    # Process sequence lines to count species.haplotypes
    if ($in_block && $line =~ /^s\s+/) {
        chomp $line;
        my @fields = split(/\s+/, $line);
        my $seq_name = $fields[1];
        my @parts = split(/\./, $seq_name);

        if (@parts >= 3) {
            my $species = $parts[0];
            my $haplotype = $parts[1];
            my $sp_hap = "$species.$haplotype";

            # Check if this is one of our desired species.haplotypes
            if (exists $sp_hap_counts{$sp_hap}) {
                $sp_hap_counts{$sp_hap}++;

                # For relaxed mode: also track by species
                if ($relaxed_mode) {
                    $species_counts{$species}{$haplotype} = 1;
                }
            }
        }
    }

    # End of alignment block (blank line or EOF)
    if ($in_block && ($line eq "\n" || eof($maf_fh))) {
        $in_block = 0;

        # Determine if we keep this block
        my $keep_block = 1;
        my @reasons = ();

        if ($relaxed_mode) {
            # RELAXED MODE: Check if at least one haplotype is present for each species
            foreach my $species (keys %species_haps) {
                # Count how many of this species' haplotypes are present
                my $present_haps = 0;
                foreach my $hap (@{$species_haps{$species}}) {
                    if (exists $species_counts{$species}{$hap}) {
                        $present_haps++;
                    }
                }

                # Get total haplotypes for this species
                my $total_haps = scalar @{$species_haps{$species}};

                # In relaxed mode, we need at least one haplotype per species
                if ($present_haps == 0) {
                    push @reasons, "species $species has no haplotypes present (needed at least 1)";
                    $keep_block = 0;
                }
                # Also check for excess copies (more than 1 of same haplotype)
                foreach my $sp_hap (@desired_sp_hap) {
                    if ($sp_hap =~ /^$species\./ && $sp_hap_counts{$sp_hap} > 1) {
                        push @reasons, "has $sp_hap_counts{$sp_hap} copies of $sp_hap (wanted 0 or 1)";
                        $keep_block = 0;
                    }
                }
            }
        } else {
            # STRICT MODE: Check for exactly one of each desired species.haplotype
            foreach my $sp_hap (@desired_sp_hap) {
                if ($sp_hap_counts{$sp_hap} == 0) {
                    push @reasons, "missing $sp_hap";
                    $keep_block = 0;
                } elsif ($sp_hap_counts{$sp_hap} > 1) {
                    push @reasons, "has $sp_hap_counts{$sp_hap} copies of $sp_hap (wanted 1)";
                    $keep_block = 0;
                }
            }
        }

        # Output block or rejection reason
        if ($keep_block) {
            print @current_block;
            $blocks_kept++;
        } else {
            print $log_fh "Block $block_number rejected: ", join(", ", @reasons), "\n";
            $blocks_filtered++;
        }

        # Reset for next block
        @current_block = ();
        next;
    }
}
close($maf_fh);

# Print summary statistics to STDERR
my $total_blocks = $blocks_kept + $blocks_filtered;
my $kept_percent = $total_blocks ? sprintf("%.2f%%", ($blocks_kept / $total_blocks) * 100) : "0.00%";
my $filtered_percent = $total_blocks ? sprintf("%.2f%%", ($blocks_filtered / $total_blocks) * 100) : "0.00%";

print STDERR "Summary: $blocks_kept blocks kept ($kept_percent), $blocks_filtered blocks filtered out ($filtered_percent), from $total_blocks total blocks\n";

close($log_fh);
