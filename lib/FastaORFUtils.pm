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
our @EXPORT_OK = qw( translateDNAToProtein
                     reverseComplement
                     translateReadingFrame
                     getHeaderParts
                     getLongestORF);

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
#           translateDNAToProtein
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

#+++++++++++++++++++++++++++++++++++++++++++++++
#           reverseComplement
#
sub reverseComplement
{
    my($sequence) = @_;

    my $revCom = reverse $sequence;
    
    $revCom =~ tr/ACGTacgt/TGCAtgca/;
    
    return $revCom;
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#           translateReadingFrame
#
sub translateReadingFrame
{
    my($sequence, $start, $end) = @_;
    
    unless ($end) {
        $end = length($sequence);
    }
    
    # return empty string if sequence is too short for end index
    return "" if ($end <= $start);
    
    my $frameSequence = substr($sequence, $start, $end - $start);

    return $frameSequence;
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#           getHeaderParts
#
sub getHeaderParts
{
    my($header) = @_;
    my @headerParts = split(" ", $header);
    my $seqId = shift(@headerParts);
    my $body = "";
    
    if ($seqId && length($seqId) > 1) {
        $seqId = substr($seqId, 1);
        $body = substr($header, length($seqId)+1);
        $body =~ s/^\s+|\s+$//g  #trim white space
    
    } else {
        $seqId = "<unknown>";  
    }
    
    return ($seqId, $body);
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#           getLongestORF
#
sub getLongestORF
{
    my($seq) = @_;
    
    my $posIndex = 0;
    my $posIndexLongest = 0;
    my $sequenceLongest = "";
    
    # match codons until a stop codon, not greedy
    while ($seq =~ m/((\w\w\w)*?)(TAA|TAG|TGA)/ig) {
        # print "Found '$1' at position $posIndex\n";
        if (length($1) > length($sequenceLongest)) {
            $posIndexLongest = $posIndex;
            $sequenceLongest = $1;
        }
        
        $posIndex = pos($seq);
    }
    
    # check leftover string
    my $leftover = substr($seq, $posIndex);
    $leftover =~ m/((\w\w\w)*)/i;
    
    if (defined($1) && length($1) > length($sequenceLongest)) {
        $posIndexLongest = $posIndex;
        $sequenceLongest = $1;
    }
    
    # print "Leftover: " . substr($seq, $posIndex) . " at position $posIndex\n";
    
    return ($sequenceLongest, $posIndexLongest);
}

# dummy subroutine
#my $str = shift;

#    return $str;

1;
