package PslItem;
use strict;

sub new {
	my ($class, $line) = @_;
	my $self = {
        match => 0,
        repMatch => 0,
        misMatch => 0,
        ns => 0,
        qNumInsert => 0,
        qGapBases => 0,
        tNumInsert => 0,
        tGapBases => 0,
        strand => "",
        qName => "",
        qSize => 0,
        qStart => 0,
        qEnd => 0,
        tName => "",
        tStart => 0,
        tEnd => 0,
        blockCount => 0,
        blockSizes => "",
        qStarts => "",
        tStarts => ""
    };
	bless $self, $class;
    
    if ($line) {
        $self->loadFromString($line);
    }
	return $self;
}

sub round() {
    my ($numberStr) = @_;
    my $roundedStr = $numberStr;
    if ($numberStr =~ m/\./) {
        my ($majorStr, $minorStr) = ($numberStr =~ m/(\d*)\.(\d+)/);
        if (&isRoundUp($minorStr)) {
            $roundedStr = $majorStr + 1;
        } else {
            $roundedStr = $majorStr;
        }
    }
     
    return $roundedStr;
}

sub isRoundUp() {
    my ($minorStr) = @_;
    my $roundUp = 0;
    my $firstDigit = substr($minorStr, 0, 1);
    if ($firstDigit eq "") {
        $roundUp = 1; 
    } elsif ($firstDigit > 5) {
        $roundUp = 1;
    } elsif ($firstDigit < 5) {
        $roundUp = 0;
    } elsif ($firstDigit == 5) {
        if (length($minorStr) == 1) {
            $roundUp = 1;
        } else {
            $roundUp = &isRoundUp(substr($minorStr, 1));
        }
    }
    return $roundUp;
}

sub loadFromString() {
    my ($self, $str) = @_;

    my @data = split(/[\ \t\n]+/, $str);
    $self->{match}    = $data[0];
    $self->{misMatch} = $data[1];
    $self->{repMatch} = $data[2];
    $self->{ns} = $data[3]; # n's
    
    $self->{qNumInsert} = $data[4]; # gap count
    $self->{qGapBases} = $data[5];
    $self->{tNumInsert} = $data[6]; # gap count
    $self->{tGapBases} = $data[7];
    
    $self->{strand} = $data[8];
    
    $self->{qName} = $data[9];    
    $self->{qSize}  = $data[10];
    $self->{qStart} = $data[11];
    $self->{qEnd}   = $data[12];
    
    $self->{tName} = $data[13];
    $self->{tSize} = $data[14];
    $self->{tStart} = $data[15];
    $self->{tEnd}   = $data[16];
    
    $self->{blockCount} = $data[17];
    $self->{blockSizes} = $data[18];
    $self->{qStarts} = $data[19];
    $self->{tStarts} = $data[20];
}

sub print() {
    my ($self) = @_;
    print "tName=$self->{tName} ";
    print "strand=$self->{strand} ";
    print "match=$self->{match} ";
    print "\n";
}

# returns the percent identity of a psl data structure
sub getPercentIdentity() {
    my ($self) = @_;
    my $isMrna = 1;
    
    my $qAlignSize = $self->{qEnd} - $self->{qStart};
    my $tAlignSize = $self->{tEnd} - $self->{tStart};
    my $min = ($qAlignSize < $tAlignSize) ? $qAlignSize : $tAlignSize;
    if ($min <= 0) {
        return 0;
    }

    my $sizeDiff = $qAlignSize - $tAlignSize;
    if ($sizeDiff < 0) {
        $sizeDiff = ($isMrna) ? 0 : -$sizeDiff;
    }

    my $insertFactor = $self->{qNumInsert};
    if (!$isMrna) {
        $insertFactor += $self->{tNumInsert};
    }

    my $milliBad = 0;
    my $total = $self->{match} + $self->{repMatch} + $self->{misMatch};
    if ($total != 0) {
        $milliBad = 1000 * ($self->{misMatch} + $insertFactor +
            &round(3 * log(1 + $sizeDiff))) / $total;
    }
    
    return 100.0 - $milliBad * 0.1; # percent identity
}

# returns the score of a psl data structure
sub getScore() {
    my ($self) = @_;
    return $self->{match} + ($self->{repMatch}>>1) - $self->{misMatch} -
        $self->{qNumInsert} - $self->{tNumInsert};
}

sub getContinuousFragments() {
    my ($self) = @_;
    
    my $minFragSize = 0.1 * $self->{qSize};
    my $maxGapSize = 0.01 * $self->{qSize};
    
    my @blocks = split(",", $self->{tStarts});
    my @blockSizes = split(",", $self->{blockSizes});

    # group pieces into fragments separated by gaps < 1% qSize
    my @frags = ($blocks[0]);
    my @fragSizes = ();
    for (my $i = 1; $i < @blocks; $i++) {
        my $gap = $blocks[$i] - ($blocks[$i-1] + $blockSizes[$i-1]);
        if ($gap > $maxGapSize) {
            # end the fragment
            my $fragSize = ($blocks[$i-1] + $blockSizes[$i-1]) - $frags[-1];
            push(@fragSizes, $fragSize);
            
            # start new fragment
            push(@frags, $blocks[$i]);
        }
    }
    my $fragSize = ($blocks[-1] + $blockSizes[-1]) - $frags[-1];
    push(@fragSizes, $fragSize);
    
    # return only fragments that exceed threshold size
    my %fragments = ();
    for (my $i = 0; $i < @frags; $i++) {
        if ($fragSizes[$i] > $minFragSize) {
            $fragments{$frags[$i]} = $frags[$i] + $fragSizes[$i];
        }
    }
    
    return %fragments; 
}


1;
