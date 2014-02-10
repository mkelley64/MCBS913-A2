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

# db for protein output file
my @proteinData;

# for each sequence
while ($header) {
    
    # get sequenceId, sequenceBody
    my($seqId, $seqBody) = FastaORFUtils::getHeaderParts($header);
      
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
    
    # create db of longest ORF per frame
    my %longestORFdb;
    my $longestORFLength = 0;
    
    foreach my $frameId (keys %readingFrames) {
        my($orf, $start) = FastaORFUtils::getLongestORF($readingFrames{$frameId});
        my $orfLength = length($orf);
        
        if ($orfLength > $longestORFLength) {
            $longestORFLength = $orfLength;
        }

        $longestORFdb{$frameId} = { "orf"   => $orf,
                                    "start" => $start };
    }
    
    # create array of longest ORF per sequence
    foreach my $frameId (keys %longestORFdb) {
        my $orf = $longestORFdb{$frameId}->{"orf"};
       
        if (length($orf) == $longestORFLength) {
            my $data = { "seqId"    => $seqId,
                         "seqBody"  => $seqBody,
                         "frameId"  => $frameId,
                         "orf"      => $orf };
            
            push(@proteinData, $data);
        }
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
    
    
    
    #--------------------------------------------------------
    $header = $inLine;    # last line read is either next header or null
}


# write "best" protein sequences to file    
my $outputFile = $seqFile . ".out";
    
unless (open(OUTPUT, ">$outputFile")) {
    print "Can not write to $outputFile";
    exit;
}

foreach (@proteinData)
{
      my $protein = $_;
      print OUTPUT ">$protein->{'seqId'}\_$protein->{'frameId'} $protein->{'seqBody'}\n";
      print OUTPUT FastaORFUtils::translateDNAToProtein($protein->{'orf'});
      print OUTPUT "\n\n";
}
    
close(OUTPUT);

exit;

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






