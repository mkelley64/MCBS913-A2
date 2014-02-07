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

use FastaORFUtils qw(translateDNAToProtein);


# Test translateDNAToProtein

is(translateDNAToProtein('TCA'), 'S', 'Works for single replacement');

is(translateDNAToProtein('TCATTCGAT'), 'SFD', 'Works for multiple replacement');

is(translateDNAToProtein('CTN'), '*', 'Handles single unknown nucleotide');

is(translateDNAToProtein('CTNNNCNNN'), '***', 'Handles multiple unknown nucleotides');