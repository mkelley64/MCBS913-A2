#!/usr/bin/perl -w
#
# fastaORF fastafile
#
#   Mark Kelley
#   Assignment 2
#   MCBS913 - Spring 2014
#   02/06/10
#
my $usageMsg = q(   Usage: fastaORF fastafile

        Extract each sequence from a fastafile into a single string.
        Create reverse complement of sequence
        Find longest ORF

        Output for each sequence is <one line per frame>:
          
        sequenceId  frame  longestOrfLength  startPosition
          
        Output sent to standard output. );
          
use 5.10.0;
use warnings;
use strict;

use lib 'lib';
use FastaORFUtils;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
&checkUsage();


#++++++++++++ Main Routine +++++++++++++++++++++++++++

my $seqFile = $ARGV[0];

open (IN, $seqFile)  or die "Unable to open: ".$seqFile ;

# first line better be a sequence header
my $header = <IN>;

if (substr( $header, 0, 1 ) ne '>') {
    print "********* ERROR ********* Expecting header, got:\n $header";
    print " is this a fasta file??? ";
    &checkUsage();
    exit;
}

# for each sequence
while ($header) {
    
    # get sequenceId
    my $seqId = FastaORFUtils::getSequenceId($header);
      
    my $seq = ""; 
    my $inLine = <IN>;

    # read in all input lines of bases to create sequence string
    while ($inLine && substr($inLine, 0, 1 ) ne '>') {
        chomp($inLine);     # remove line feed
        $seq = $seq . $inLine;
        $inLine = <IN>;
    }
    
    my %readingFrames;
    
    # generate forward reading frames
    for (my $index=0; $index<3; $index++) {
        $readingFrames{$index} = FastaORFUtils::translateReadingFrame($seq, $index);
    }
    
    # get reverse complement
    my $revCom = FastaORFUtils::reverseComplement($seq);
    
    # generate reverse reading frames
    for (my $index=0; $index<3; $index++) {
        $readingFrames{"r$index"} = FastaORFUtils::translateReadingFrame($revCom, $index);
    }
    
    my %longestORFs;
    my $longestORFLength = 0;
    
    foreach my $key (keys %readingFrames) {
        my($orf, $len) = FastaORFUtils::getLongestORF($readingFrames{$key});
        $longestORFs{$key} = $orf;
    }
    
    
    # temp print results
    for (my $index=0; $index<3; $index++) {
        my $key = $index;
        my $val = $longestORFs{$key};
        say "$seqId  $key   " .length($val) . "  $val";
    }
    
    for (my $rindex=0; $rindex<3; $rindex++) {
        my $key = "r$rindex";
        my $val = $longestORFs{$key};
        say "$seqId  $key  " .length($val) . "  $val";
    }
    
    print "\n";
    # end temp
    
    # call translation subroutine
    # my $proteinSequence = FastaORFUtils::translateDNAToProtein($seq);

    # output results
    #print $header;
    # say $proteinSequence;

    #--------------------------------------------------------
    $header = $inLine;    # last line read is either next header or null
}


#+++++++++++++++++++++++++++++++++++++++++++++++
#                checkUsage
#
sub checkUsage()
{
    if ( @ARGV == 0 || $ARGV[0] eq "-h" ) {
        print STDERR "$usageMsg\n";
        exit;
    }
}






