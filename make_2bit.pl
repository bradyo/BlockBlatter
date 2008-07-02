
my $genomePath = shift;
my $genomeFilename = shift;

if (!defined($genomeFilename)) {
    $genomeFilename = "genome.2bit";
}

my $faToTwoBitProgram = "C:/blat/faToTwoBit.exe";

my $filename = &create2Bit($genomePath);
print "\ncreated $filename\n";

sub create2Bit() {
    my ($genomeDir) = @_;
    my $genomeFile = "$genomeDir/genome.2bit";
    #my $genomeFile = "genome.2bit";
    # list fasta files
    my @files = &getFaFiles($genomeDir);

    my $cmd = "\"$faToTwoBitProgram\" @files $genomeFile –ignoreDups";
    print "executing: $cmd\n";
    #system($cmd);
    return $genomeFile;
}

sub getFaFiles() {
    my ($dir) = @_;
    opendir(DIR, $dir);
    my @files = grep(/\.fa$/, readdir(DIR));
    my @fullFiles = ();
    foreach my $file (@files) {
        $dir =~ s/\\/\//g;
        push(@fullFiles, $dir."/".$file);
    }
    return @fullFiles;
}