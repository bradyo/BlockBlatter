#!/usr/bin/perl -w

# =============================================================================
# DOCUMENTATION
# =============================================================================
# synopsis:
#   blats sequences in a fasta file (source) against a folder of genome fasta
#   files and generates a bed file from the hits
# input:
#   source file (string) - the name of the source fasta file containing block
#       sequences
#   genome directory (string) - the directory containing genome fasta files,
#       the files in the directory must contain headers as "chrN" where N is
#       the corresponding chromosome
#   output file (string) - the filename to write the output BED data to
#   config (string) - an identifier matching a config section in the
#       config.ini file. Default is "windows"
# output:
#   a bed file representing the genomic hits of the source file sequences.
#   Resulting psl files are written to the psl directory and resulting block
#   sequences are written to the blocks directory

use strict;
use PslItem;
use Bio::Seq;
use Bio::SeqIO;

# default system parameters
use constant DEFAULT_FA_TO_TWO_BIT_PROGRAM => 'blat/faToTwoBit.exe';
use constant DEFAULT_GF_SERVER_PROGRAM => 'blat/gfServer.exe';
use constant DEFAULT_GF_SERVER_NAME => 'localhost';
use constant DEFAULT_GF_SERVER_PORT => 50000;
use constant DEFAULT_GF_CLIENT_PROGRAM => 'blat/gfClient.exe';

# default paramters
use constant CONFIG => 'config/config.ini';
use constant CONFIG_MIN_PERCENT_LENGTH => 'config/minPercentLength.ini';
use constant CONFIG_COLORS => 'config/colors.ini';
use constant DEFAULT_MIN_PERCENT_LENGTH => 10;
use constant MIN_UNIQUE_PERCENT_LENGTH => 18;
use constant DEFAULT_COLOR => '0,0,0';
use constant DEFAULT_MAX_MERGE_GAP => 100;

# load default configuration
my $faToTwoBitProgram = DEFAULT_FA_TO_TWO_BIT_PROGRAM;
my $gfServerProgram = DEFAULT_GF_SERVER_PROGRAM;
my $gfServerName = DEFAULT_GF_SERVER_NAME;
my $gfServerPort = DEFAULT_GF_SERVER_PORT;
my $gfClientProgram = DEFAULT_GF_CLIENT_PROGRAM;
my $minPercentLength = DEFAULT_MIN_PERCENT_LENGTH;
my $maxMergeGap = DEFAULT_MAX_MERGE_GAP;

# hashes for configuration settings
my %trackColors = ();
my %minPercentLengths = ();

&main(@ARGV);

# =============================================================================
# MAIN
# =============================================================================

sub main() {
    my ($blocksFile, $genomeDir, $outputFile, $config) = @_;
       
    # override default config with those given in config/config.ini    
    &loadConfig($config) if (defined($config));
    
    # load hashes from config files
    my %trackColors = &getConfigHash('config/colors.ini');
    my %minPercentLengths = &getConfigHash('config/minPercentLength.ini');
    
    print("blocksFile = $blocksFile\n");
    print("genomeDir = $genomeDir\n");
    print("outputFile = $outputFile\n");
    print("\n");
    
    print("faToTwoBitProgram = $faToTwoBitProgram\n");
    print("gfServerProgram = $gfServerProgram\n");
    print("gfServerName = $gfServerName\n");
    print("gfServerPort = $gfServerPort\n");
    print("gfClientProgram = $gfClientProgram\n");
    print("minPercentLength = $minPercentLength\n");
    print("maxMergeGap = $maxMergeGap\n");
    print("\n");
    
    print("trackColors:\n");
    foreach my $key (keys(%trackColors)) {
        print("$key => $trackColors{$key}\n");
    }
    print("\n");
    
    print("minPercentLengths:\n");
    foreach my $key (keys(%minPercentLengths)) {
        print("$key => $minPercentLengths{$key}\n");
    }
    print("\n");
       
    # start blat server for genome
    print "building 2bit file:\n";
    print "$genomeDir/*.fa -> genome.2bit\n";
    my $genomeFile = &create2Bit($genomeDir);
    print "done\n\n";
    
    print "starting blat server using $genomeFile\n";
    &startBlatServer($genomeFile);
    print "done\n\n";
    
    # blat the blocks, write output to "psl" folder
    print "blatting sequences in $blocksFile\n";
    &gfPolyClient($genomeDir, $blocksFile);
    print "done\n\n";
    
    # delete genome 2bit file
    unlink($genomeFile);
    
    # make a bed file from output psl folder
    print "making $outputFile BED file from psl directory\n";
    &makeBed("psl", $outputFile);
    print "done\n\n";
    
    # remove overlaps and small gaps in output bed file
    print "removing overlaps and gaps from $outputFile\n";
    &cleanBed($outputFile, $outputFile);
    print "done\n\n";
}

# =============================================================================
# PRIMARY SUBROUTINES
# =============================================================================

sub create2Bit() {
    my ($genomeDir) = @_;
    my $genomeFile = "$genomeDir/genome.2bit";
    # list fasta files
    my @files = &getFaFiles($genomeDir);

    my $cmd = "\"$faToTwoBitProgram\" @files $genomeFile –ignoreDups";
    print "executing: $cmd\n";
    system($cmd);
    return $genomeFile;
}

sub getFaFiles() {
    my ($dir) = @_;
    opendir(DIR, $dir);
    my @files = grep(/\.fa$/, readdir(DIR));
    my @fullFiles = ();
    foreach my $file (@files) {
        push(@fullFiles, $dir."/".$file);
    }
    closedir(DIR);
    return @fullFiles;
}

sub startBlatServer() {
    my ($genomeFile) = @_;
    my $cmd = "\"$gfServerProgram\" start $gfServerName $gfServerPort "
        . "–stepSize=5 $genomeFile";
    print "executing: $cmd\n";
    system($cmd);
}

sub gfPolyClient() {
    my ($databaseDir, $fastaFilename) = @_;

	# load sequences from fasta query
	my @blockSeqs = ();
	my $in = Bio::SeqIO->new(-file => "$fastaFilename", -format => 'Fasta');
	while (my $seq = $in->next_seq()) {
		push(@blockSeqs, $seq);
	}
	
	# create directories if they dont already exist
    mkdir('blocks') if (!(-e 'blocks'));
    opendir(DIR, 'blocks');	
	while (my $filename = readdir(DIR)) {
		unlink("blocks/$filename") if $filename =~ m/\.fa$/;
    }
    
	mkdir('psl') if (!(-e 'psl'));
    opendir(DIR, 'psl');	
	while (my $filename = readdir(DIR)) {
		unlink("psl/$filename") if $filename =~ m/\.psl$/;
    }
	
	foreach my $blockSeq (@blockSeqs) {
		# extract block name from sequence
		my $blockName = $blockSeq->id;
		print "blating: $blockName. ";
		
		# write a fasta file for this block to blat with
		my $blockFilename = "blocks/$blockName.fa";
		my $seqOut = Bio::SeqIO->new('-file' => ">$blockFilename",
			'-format' => 'Fasta');
		$seqOut->write_seq($blockSeq);
		
		# blat sequences in fasta query
		my $pslFilename = "psl/$blockName.psl";
        my $cmd = "\"$gfClientProgram\" $gfServerName $gfServerPort "
            . "$databaseDir $blockFilename $pslFilename "
            . "-minScore=0 -minIdentity=0";
        print "executing: $cmd\n";
		system($cmd);
	}
}

sub makeBed() {
	my ($pslFolder, $outputFilename) = @_;

    open(OUT, ">$outputFilename");
        
	print OUT "track name=\"Gap\" visibility=1\n";
	print OUT "track name=\"Blocks\" visibility=3 itemRgb=\"On\"\n";
	
	opendir(DIR, $pslFolder);	
	while (my $filename = readdir(DIR)) {
		# skip files not ending in psl
		next if ($filename !~ m/\w+.psl/);
		
		# handle each psl file
		my ($trackName) = ($filename =~ m/(\w+)\.\w+/);
		my $pslFilename = "$pslFolder/$filename";
		print "processing: $pslFilename\n";
		
        	# filter each line by percent identity > 90
            my @fragments = ();
            print OUT "#$trackName\n";
            open(PSL_FILE_IN, "<$pslFilename");
            while (my $line = <PSL_FILE_IN>) {
                next if (substr($line, 0, 1) !~ m/\d/); # not data line
                
                # make a psl object from the line of data
                my $pslItem = new PslItem($line);
                my $chr = $pslItem->{tName};
                my $strand = $pslItem->{strand};
                my $color = &getColor($trackName);
                my $percentIdentity = $pslItem->getPercentIdentity();
                if ($percentIdentity > 90 && $chr !~ m/(random|Un)/) {
                    # output continuous fragments greater than 10% query size
                    my %fragments = $pslItem->getContinuousFragments();
                    while (my ($start, $end) = each(%fragments)) {
                        my $percentLength = ($end - $start) / $pslItem->{qSize}
                            * 100;
                        my $minPercentLength = &getMinPercentLength($trackName);
                        if ($percentLength > $minPercentLength) {				
                            print OUT "$chr\t$start\t$end\t$trackName\t"
                                . "$percentIdentity\t$strand\t$start\t"
                                . "$end\t$color\n";	
                        }
                    }
                }
            }
            close(PSL_FILE_IN);
	}
	closedir(DIR);
}

sub cleanBed() {
    my ($inputFilename, $outputFilename, $maxGap) = @_;
    
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
        } else { 
            print FILE_OUT "$linesIn[$i]";
        }
    }
    close(FILE_OUT);
}

# =============================================================================
# SECONDARY SUBROUTINES
# =============================================================================

sub loadConfig() {
    my ($config) = @_;
    
    # look for config section in config file
    my $filename = CONFIG;
    open(IN, "<$filename");
    while (my $line = <IN>) {
        if ($line =~ /\[$config config\]/) {
            $line = <IN>;
            while (defined($line) && $line !~ /\[.*\]/) {
                if ($line =~ /faToTwoBitProgram/) {
                    ($faToTwoBitProgram) = ($line =~ /=\s*(.+)/);
                }
                elsif ($line =~ /gfServerProgram/) {
                    ($gfServerProgram) = ($line =~ /=\s*(.+)/);
                }
                elsif ($line =~ /gfClientProgram/) {
                    ($gfClientProgram) = ($line =~ /=\s*(.+)/);
                }
                elsif ($line =~ /gfServerName/) {
                    ($gfServerName) = ($line =~ /=\s*(.+)/);
                }
                elsif ($line =~ /gfServerPort/) {
                    ($gfServerPort) = ($line =~ /=\s*(.+)/);
                }
                $line = <IN>;
            }
        }
        
        if ($line =~ /\[filters\]/) {
            $line = <IN>;
            while (defined($line) && $line !~ /\[.*\]/) {
                if ($line =~ /minPercentLength/) {
                    ($minPercentLength) = ($line =~ /=\s*(.+)/);
                }
                elsif ($line =~ /maxMergeGap/) {
                    ($maxMergeGap ) = ($line =~ /=\s*(.+)/);
                }
                $line = <IN>;
            }
        }    
    }
    close(IN);
}

sub getConfigHash() {
    my ($filename) = @_;
    my %hash = ();    
       
    open(IN, "<$filename") or print("could not open config file: $filename\n");
    my @lines = ();
    while (my $line = <IN>) {
        $line =~ /^(.+)\s*=\s*(.+)$/;
        my $key = $1;
        my $value = $2;
        $hash{$key} = $value;
    }
    close(IN);
    return %hash;
}

sub getColor() {
	my ($trackName) = @_;
	my $color = DEFAULT_COLOR;
	if (exists $trackColors{$trackName}) {
		$color = $trackColors{$trackName}
	}
	return $color;
}

sub getMinPercentLength() {
	my ($trackName) = @_;
	my $minLength = DEFAULT_MIN_PERCENT_LENGTH;
	if (exists $minPercentLengths{$trackName}) {
		$minLength = $minPercentLengths{$trackName};
	} elsif ($trackName =~ m/unique/) {
		$minLength = MIN_UNIQUE_PERCENT_LENGTH;
	}
	return $minLength;
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




