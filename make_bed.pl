#!/usr/bin/perl -w
use strict;
use PslItem;

my %trackColors = (
	"Default" => "0,0,0",
	"Block1" => "255,66,0",
	"Block10" => "255,146,88",
	"Block11" => "11,55,11",
	"Block12" => "77,225,255",
	"Block13" => "77,255,77",
	"Block14" => "77,44,00",
	"Block15" => "88,55,110",
	"Block16" => "11,77,11",
	"Block17" => "22,22,121",
	"Block18" => "255,255,66",
	"Block19" => "0,121,0",
	"Block2" => "55,196,55",
	"Block20" => "255,33,169",
	"Block21" => "225,225,00",
	"Block22" => "110,66,33",
	"Block23" => "66,22,88",
	"Block24" => "55,110,55",
	"Block25" => "121,66,33",
	"Block26" => "255,99,169",
	"Block27" => "77,169,255",
	"Block28" => "99,110,144",
	"Block29" => "255,144,99",
	"Block3" => "22,110,255",
	"Block30" => "255,196,88",
	"Block31" => "121,255,121",
	"Block32n" => "144,110,169",
	"Block33" => "110,225,255",
	"Block34" => "144,144,00",
	"Block35" => "99,225,255",
	"Block36" => "77,110,77",
	"Block37" => "55,167,55",
	"Block38" => "110,225,169",
	"Block39n" => "215,215,99",
	"Block4" => "255,164,0",
	"Block40" => "255,225,110",
	"Block41" => "33,33,33",
	"Block42" => "255,169,77",
	"Block43" => "225,110,55",
	"Block44" => "99,99,66",
	"Block45" => "255,00,144",
	"Block46" => "99,99,66",
	"Block47" => "00,00,169",
	"Block48" => "169,00,255",
	"Block49" => "99,99,99",
	"Block5" => "0,169,255",
	"Block6" => "0,88,0",
	"Block7" => "99,77,169",
	"Block8" => "55,55,225",
	"Block9" => "33,144,144",
	"PseudoX" => "99,99,99",
	"TAR" => "169,00,66",
	"Telomere" => "66,169,99",
);

my %minPercentLengths = (
	"Default" => 10,
	"Block14" => 15,
	"Block26" => 18,
	"Block32n" => 20,
	"gaattg_repeats" => 18,
	"overlapping" => 18,
	"PseudoX" => 18
);

sub getColor() {
	my ($trackName) = @_;
	my $color = $trackColors{"Default"};
	if (exists $trackColors{$trackName}) {
		$color = $trackColors{$trackName}
	};
	return $color;
}

sub getMinPercentLength() {
	my ($trackName) = @_;
	my $minLength = $minPercentLengths{"Default"};
	if (exists $minPercentLengths{$trackName}) {
		$minLength = $minPercentLengths{$trackName};
	} elsif ($trackName =~ m/unique/) {
		$minLength = 18;
	}
	return $minLength;
}

sub outputBed() {
	my ($trackName, $pslFilename) = @_;
	
	# filter each line by percent identity > 90
	my @fragments = ();
	print "#$trackName\n";
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
			# output all continuous fragments greater than 10% query size
			my %fragments = $pslItem->getContinuousFragments();
			while (my ($start, $end) = each(%fragments)) {
				my $percentLength = ($end - $start) / $pslItem->{qSize} * 100;
				my $minPercentLength = &getMinPercentLength($trackName);
				if ($percentLength > $minPercentLength) {				
	     			print "$chr\t$start\t$end\t$trackName\t";
					print "$percentIdentity\t$strand\t$start\t$end\t$color\n";	
				}
			}
		}
	}
	close(PSL_FILE_IN);
}

sub main() {
	my ($pslFolder, $bedFilename) = @_;

	print "track name=\"Gap\" visibility=1\n";
	print "track name=\"Blocks\" visibility=3 itemRgb=\"On\"\n";
	
	opendir(DIR, $pslFolder);	
	while (my $filename = readdir(DIR)) {
		# skip files not ending in psl
		next if ($filename !~ m/\w+.psl/);
		
		# handle each psl file
		my ($trackName) = ($filename =~ m/(\w+)\.\w+/);
		my $pslFilename = "$pslFolder/$filename";
		print STDERR "processing: $pslFilename\n";
		&outputBed($trackName, $pslFilename);
	}
	closedir(DIR);
}

&main(@ARGV);

