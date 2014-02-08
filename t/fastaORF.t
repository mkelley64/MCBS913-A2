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

use FastaORFUtils qw(translateDNAToProtein reverseComplement translateReadingFrame);


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