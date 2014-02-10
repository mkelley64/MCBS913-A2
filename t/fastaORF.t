#!/usr/bin/perl -w
#
#   Mark Kelley
#   Assignment 2
#   MCBS913 - Spring 2014
#   02/06/10
#
use 5.10.0;
use warnings;
use strict;

use Test::More qw( no_plan );

use lib '../lib';

use FastaORFUtils qw( translateDNAToProtein
                      reverseComplement
                      translateReadingFrame
                      getLongestORF
                      getHeaderParts );


# Test translateDNAToProtein

is(translateDNAToProtein('TCA'), 'S', 'Works for single replacement');
is(translateDNAToProtein('TCATTCGAT'), 'SFD', 'Works for multiple replacement');
is(translateDNAToProtein('TAATAGTGA'), '___', 'Handles STOP codons');
is(translateDNAToProtein('CTN'), '*', 'Handles single unknown nucleotide');
is(translateDNAToProtein('CTNNNCNNN'), '***', 'Handles multiple unknown nucleotides');
is(translateDNAToProtein('tcatTcgAT'), 'SFD', 'Handles lower case nucleotides');
is(translateDNAToProtein('AT'), '', 'Handles nucleotide strings with length < 3');


# Test reverseComplement

is(reverseComplement('CCGGAAAAAAAATTTATATAT'), 'ATATATAAATTTTTTTTCCGG', 'Works for multiple nucleotides all upper case');
is(reverseComplement('CCGGaaaaaaaatttatatat'), 'atatataaattttttttCCGG', 'Works for multiple nucleotides mixed case');
is(reverseComplement('CCGGAAAAAANNTTTATATAn'), 'nTATATAAANNTTTTTTCCGG', 'Handles unknown nucleotides');


# Test translateReadingFrame

is(translateReadingFrame('CCGGAAAAAAAATTTATATAT', 0), 'CCGGAAAAAAAATTTATATAT', 'Works for reading frame 0');
is(translateReadingFrame('CCGGAAAAAAAATTTATATAT', 1), 'CGGAAAAAAAATTTATATAT', 'Works for reading frame 1');
is(translateReadingFrame('CCGGAAAAAAAATTTATATAT', 2), 'GGAAAAAAAATTTATATAT', 'Works for reading frame 2');
is(translateReadingFrame('CAT', 3), '', 'Handles nucleotide string length == substring length');
is(translateReadingFrame('CAT', 4), '', 'Handles nucleotide string length < substring length');


# Test getLongestORF

my($seq, $len) = getLongestORF('TGA');
is($seq, '', 'Handles only stop codon: Sequence');
is($len, 0, 'Handles only stop codon: Length');

($seq, $len) = getLongestORF('CATCAT');
is($seq, 'CATCAT', 'Handles no stop codons: Sequence');
is($len, 0, 'Handles no stop codons: Length');

($seq, $len) = getLongestORF('CATCATTAGCAT');
is($seq, 'CATCAT', 'Handles one stop codon: Sequence');
is($len, 0, 'Handles one stop codon: Length');

($seq, $len) = getLongestORF('GGGCATCATTAGCATTAAGGG');
is($seq, 'GGGCATCAT', 'Handles two stop codons, first frame longest: Sequence');
is($len, 0, 'Handles two stop codons, first frame longest: Length');

($seq, $len) = getLongestORF('GGGTAGGGGCATCATTAAGGG');
is($seq, 'GGGCATCAT', 'Handles two stop codons, 2nd frame longest: Sequence');
is($len, 6, 'Handles two stop codons, 2nd frame longest: Length');

($seq, $len) = getLongestORF('TAGCATTAAGGGCATCAT');
is($seq, 'GGGCATCAT', 'Handles two stop codons, last frame longest: Sequence');
is($len, 9, 'Handles two stop codons, last frame longest: Length');

($seq, $len) = getLongestORF('TAGCATTAAGGGCATCATG');
is($seq, 'GGGCATCAT', 'Handles two stop codons, last frame longest: Sequence');
is($len, 9, 'Handles two stop codons, last frame longest: Length');

($seq, $len) = getLongestORF('TAGCATTAAGGGCATCATGG');
is($seq, 'GGGCATCAT', 'Handles two stop codons, last frame longest, one extra char: Sequence');
is($len, 9, 'Handles two stop codons, last frame longest, two extra chars: Length');


# Test getHeaderParts

my($seqId, $seqBody) = getHeaderParts('>simple_1');
is($seqId, 'simple_1', 'Works for one word header: ID');
is($seqBody, '', 'Works for one word header: Body');

($seqId, $seqBody) = getHeaderParts('>simple_1 more stuff');
is($seqId, 'simple_1', 'Works for header with body: ID');
is($seqBody, 'more stuff', 'Works for header with body: Body');

($seqId, $seqBody) = getHeaderParts('>');
is($seqId, '<unknown>', 'Handles missing sequence id: ID');
is($seqBody, '', 'Handles missing sequence id: Body');

