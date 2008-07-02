#!/usr/bin/perl -w
use Getopt::Std;
use Getopt::Long;
use strict;

use constant USAGE =><<END;
SYNOPSIS:

 clean_bed.pl <input_bed_file> [OPTIONS]

DESCRIPTION:

Takes a BED file and merges overlapping segments of each track. Optionally
merges segments of a given track where separated by gaps.

OPTIONS:

       --help
             Prints usage message.
       --output <output_bed_file>
             Outputs the results to the given file
       --gaps <int>
             Merges segments separated by <int> bp gaps

EXAMPLES:

 clean_bed.pl input.bed
 clean_bed.pl input.bed -o output.bed
 clean_bed.pl input.bed -o output.bed -g 100
 
AUTHOR:

Brady Olsen

END

sub dieUsage() {
    print USAGE;
    exit();
}

sub dieHelp() {
    print "Usage: clean_bed.pl <input_bed_file> [OPTIONS]\n";
    print "Parameters:\n";
    print "--help\n";
    print "    Prints usage message.\n";
    print "--output <output_bed_file>\n";
    print "    Outputs the results to the given file\n";
    print "--gaps <int>\n";
    print "    Merges segments separated by <int> bp gaps\n";
    exit();
}

# get list of chromosomes contained in the given bed lines
sub getChromosomes() {
    my ($linesRef) = @_;
    my @lines = @$linesRef;

    my @chrs = ();
    foreach my $line (@lines) {
        my @data = split(/\s/, $line);
        my $lineChr = $data[0];
        my $shouldAdd = 1;
        foreach my $chr (@chrs) {
            if ($chr eq $lineChr) { # already in list
                $shouldAdd = 0;
            }
        }
        if ($shouldAdd) {
            push(@chrs, $lineChr);
        }
    }
    return @chrs;
}

# gets the lines with the given chromosome
sub getChromosomeLines() {
    my ($linesRef, $chr) = @_;
    my @lines = @$linesRef;
    
    my @chrLines = ();
    foreach my $line (@lines) {
        my @data = split(/\s/, $line);
        my $lineChr = $data[0];
        if ($chr eq $lineChr) {
            push(@chrLines, $line);
        }
    }
    return @chrLines;
}

# sorts the given line by start position
sub getSortedLines() {
    my ($linesRef) = @_;
    my @lines = @$linesRef;
    
    my %hash = ();
    foreach my $line (@lines) {
        my @data = split(/\s/, $line);
        my $start = $data[1];
        $hash{$start} = $line;
    }
    
    my @sortedLines = ();
    foreach my $key (sort(keys %hash)) {
        push(@sortedLines, $hash{$key});
    }
    return @sortedLines;
}

sub getMinStart() {
    my ($linesRef) = @_;
    my @lines = @$linesRef;
    
    my $minLine = "";
    my $minStart = 0;
    foreach my $line (@lines) {
        my @data = split(/\s/, $line);
        my $start = $data[1];
        if ($start < $minStart) {
            $minStart = $start;
            $minLine = $line;
        }
    }
    return $minLine;
}

sub getWeightedScore() {
    my ($lineA, $lineB) = @_;
    
    my @dataA = split(/\s/, $lineA);
    my $lenA = $dataA[2] - $dataA[1];
    my $scoreA = $dataA[4];
    
    my @dataB= split(/\s/, $lineB);
    my $lenB = $dataB[2] - $dataB[1];
    my $scoreB = $dataB[4];

    return ($scoreA * $lenA + $scoreB * $lenB) / ($lenA + $lenB);
}

# returns an array of the lines with overlaps combined into one line
sub getContinuousLines() {
    my ($sortedLinesRef) = @_;
    my @sortedLines = @$sortedLinesRef;
    
    my @mergedLines = ();
    for (my $a = 0; $a < @sortedLines; $a++) {              
        my $lineA = $sortedLines[$a];
        my @dataA = split(/\s/, $lineA);
        my $startA = $dataA[1];
        my $endA = $dataA[2];
        
        # look for a merge
        my $mergeFound = 0;
        foreach (my $b = $a + 1; $b < @sortedLines; $b++) {
            my $lineB = $sortedLines[$b];
            my @dataB = split(/\s/, $lineB);
            my $startB = $dataB[1];
            my $endB = $dataB[2];
            if (($startA >= $startB && $startA <= $endB) ||
                ($endA >= $startB && $endA <= $endB) ||
                ($startA <= $startB && $endA >= $endB) ||
                ($startB <= $startA && $endB >= $endA)) {
                # overlap found
                my $mergeStart = ($startA < $startB) ? $startA : $startB;
                my $mergeEnd = ($endA > $endB) ? $endA : $endB;
                my $mergeScore = &getWeightedScore($lineA, $lineB);
                
                # copy merged data into line B
                $dataB[1] = $mergeStart;
                $dataB[6] = $mergeStart; # repeat for color
                $dataB[2] = $mergeEnd;
                $dataB[7] = $mergeEnd; # repeat for color
                $dataB[4] = $mergeScore;
                $sortedLines[$b] = join("\t", @dataB) . "\n";
                $mergeFound = 1;
                last;
            }
        }
        if (!$mergeFound) {
            push(@mergedLines, $lineA);
        }
    }
    return @mergedLines;
}

sub getGaplessLines() {
    my ($sortedLinesRef, $gapThreshold) = @_;
    my @sortedLines = @$sortedLinesRef;
    
    my @mergedLines = ();
    for (my $a = 0; $a < @sortedLines; $a++) {              
        my $lineA = $sortedLines[$a];
        my @dataA = split(/\s/, $lineA);
        my $startA = $dataA[1];
        my $endA = $dataA[2];
        
        # look for a gap smaller than threshold
        my $mergeFound = 0;
        foreach (my $b = $a + 1; $b < @sortedLines; $b++) {
            my $lineB = $sortedLines[$b];
            my @dataB = split(/\s/, $lineB);
            my $startB = $dataB[1];
            my $endB = $dataB[2];
            
            if (($startA < $startB && $endA < $startB) &&
                ($startB - $endA <= $gapThreshold)) {
                # gap found
                my $mergeStart = $startA;
                my $mergeEnd = $endB;
                my $mergeScore = &getWeightedScore($lineA, $lineB);
                
                # copy merged data into line B
                $dataB[1] = $mergeStart;
                $dataB[6] = $mergeStart; # repeat for color
                $dataB[2] = $mergeEnd;
                $dataB[7] = $mergeEnd; # repeat for color
                $dataB[4] = $mergeScore;
                $sortedLines[$b] = join("\t", @dataB) . "\n";
                $mergeFound = 1;
                last;
            }
        }
        if (!$mergeFound) {
            push(@mergedLines, $lineA);
        }
    }
    return @mergedLines;
}

sub getTrack() {
    my ($line) = @_;
    my @data = split(/\s/, $line);
    return $data[3];
}

sub main() {
    # get core parameters
    &dieUsage() if (@ARGV == 0);
    my $inputFilename = $ARGV[0];
    if (!(-e $inputFilename)) {
        print "File not found: $inputFilename\n\n";
        &dieHelp();
    }
    
    # get optional parameters
    my $outputFilename = $inputFilename;
    my $help = 0;
    my $maxGap = 0;
    GetOptions(
        "help" => \$help,
        "gap=s" => \$maxGap,
        "output=s" => \$outputFilename) or &dieHelp();
    &dieUsage() if $help;

    open (FILE_IN, "<$inputFilename");
    my @linesIn = <FILE_IN>;
    close(FILE_IN);
    
    open(FILE_OUT, ">$outputFilename");
    for (my $i = 0; $i < @linesIn; $i++) {
        if ($linesIn[$i] =~ /^chr/) { # new track            
            # read down to find last item in this track
            my $track = &getTrack($linesIn[$i]);
            my $bottom = $i;
            for (my $j = $i + 1; $j < @linesIn; $j++) {
                if ($linesIn[$j] =~ /^chr/) {
                    if (&getTrack($linesIn[$j]) eq $track) {
                        $bottom = $j;
                    } else {
                        last;
                    }
                }
            }
            
            # read in lines to bottom of track
            my @lineBuffer = ();
            for (my $j = $i; $j <= $bottom; $j++) {
                if ($linesIn[$j] =~ /^chr/) { # data line 
                    push(@lineBuffer, $linesIn[$j]);
                } else {
                    print FILE_OUT "$linesIn[$j]";
                }
            }
        
            # process lines for this track found in lineBuffer
            my @chrs = &getChromosomes(\@lineBuffer);
            foreach my $chr (@chrs) {
                my @chrLines = &getChromosomeLines(\@lineBuffer, $chr);
                my @sortedLines = &getSortedLines(\@chrLines);
                my @continuousLines = &getContinuousLines(\@sortedLines);
                my @gaplessLines = &getGaplessLines(\@continuousLines, $maxGap);
                foreach my $line (@gaplessLines) {
                    print FILE_OUT "$line";
                }
            }
            
            # processed through to bottom, continue loop from here
            $i = $bottom;
        } else { # comment
            print FILE_OUT "$linesIn[$i]";
        }
    }
    close(FILE_OUT);
}

&main();