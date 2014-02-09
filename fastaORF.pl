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
    
    # create db of ORFs
    my %longestORFdb;
    my $longestORFLength = 0;
    
    foreach my $key (keys %readingFrames) {
        my($orf, $start) = FastaORFUtils::getLongestORF($readingFrames{$key});
        my $orfLength = length($orf);
        
        if ($orfLength > $longestORFLength) {
            $longestORFLength = $orfLength;
        }

        $longestORFdb{$key} = { "orf" => $orf,
                                "start" => $start};
    }
    
    # print results to console
    print $header;
    
    for (my $index=0; $index<3; $index++) {
        my $frameId = $index;
        printORFdata ($seqId, $frameId, $longestORFdb{$frameId}, $longestORFLength);
    }
    
    for (my $rindex=0; $rindex<3; $rindex++) {
        my $frameId = "r$rindex";
        printORFdata ($seqId, $frameId, $longestORFdb{$frameId}, $longestORFLength);
    }
    
    print "\n";
    
    # write "best" protein sequence to file
    foreach my $key (keys %longestORFdb) {
        my $orf = $longestORFdb{$key}->{'orf'};
        
        if (length($orf) == $longestORFLength) {
            my $outputFile = $seqId . ".out";
            
            unless (open(OUTPUT, ">$outputFile")) {
                print "Can not write to $outputFile";
            }
            
            print OUTPUT ">$seqId\_$key\n";
            print OUTPUT FastaORFUtils::translateDNAToProtein($orf);
            
            close(OUTPUT);
            
            next;
        } 
    }
    
    #--------------------------------------------------------
    $header = $inLine;    # last line read is either next header or null
}

#+++++++++++++++++++++++++++++++++++++++++++++++
#                printORFdata
#
sub printORFdata
{
    my $sequenceId = $_[0];
    my $frameId = $_[1];
    my $orfData = $_[2];
    my $maxLength = $_[3];
    
    my $orfLength = length($orfData->{'orf'});
    my $maxFlag = ($orfLength == $maxLength) ? "*": " ";
    
    say "$maxFlag$sequenceId\t$frameId\t$orfLength\t$orfData->{'start'}";
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






