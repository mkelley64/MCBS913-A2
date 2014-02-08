#
#   Mark Kelley
#   Assignment 2
#   MCBS913 - Spring 2014
#   02/06/10
#

package FastaORFUtils;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(translateDNAToProtein reverseComplement translateReadingFrame);

our $VERSION = "1.00";

#++++++++++++++++++++ DNA Lookup Table ++++++++++++++++++++++++++
my %genetic_code = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # Leucine
'TTG' => 'L', # Leucine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '_', # Stop
'TAG' => '_', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '_', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # Leucine
'CTC' => 'L', # Leucine
'CTG' => 'L', # Leucine
'CTT' => 'L', # Leucine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # Glutamine
'CAG' => 'Q', # Glutamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # Isoleucine
'ATC' => 'I', # Isoleucine
'ATT' => 'I', # Isoleucine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # Glutamic Acid
'GAG' => 'E', # Glutamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G'  # Glycine
);


#+++++++++++++++++++++++++++++++++++++++++++++++
#           translateToProtein
#
sub translateDNAToProtein
{
    my($sequence) = @_;
    
    my $protein = "";

    for (my $i=0; $i<length($sequence)-2; $i+=3) {
        my $codon = substr($sequence,$i,3);
        
        $codon = uc($codon);
        
        if (exists $genetic_code{$codon}) {
            $protein .= $genetic_code{$codon};
        } else {
            $protein .= "*";
        }   
    }
    
    return $protein;
}

sub reverseComplement
{
    my($sequence) = @_;

    my $revCom = reverse $sequence;
    
    $revCom =~ tr/ACGTacgt/TGCAtgca/;
    
    return $revCom;
}

sub translateReadingFrame
{
    my($sequence, $start, $end) = @_;
    
    unless ($end) {
        $end = length($sequence);
    }
    
    my $frameSequence = substr($sequence, $start, $end - $start);

    return $frameSequence;
}

# dummy subroutine
#my $str = shift;

#    return $str;

1;
