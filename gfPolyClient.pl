#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

# run: perl gfPolyClient.pl server port databaseDir blocks.fa
sub main() {
	my ($server, $port, $databaseDir, $fastaFilename) = @_;

	# load sequences from fasta query
	my @blockSeqs = ();
	my $in = Bio::SeqIO->new(-file => "$fastaFilename", -format => 'Fasta');
	while (my $seq = $in->next_seq()) {
		push(@blockSeqs, $seq);
	}
	
	# create directories if they dont already exist
	mkdir('blocks') if (!(-e 'blocks'));
	mkdir('psl') if (!(-e 'psl'));
	
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
		system("gfClient $server $port $databaseDir $blockFilename " .
			   "$pslFilename -minScore=0 -minIdentity=0");
	}
}

&main(@ARGV);



