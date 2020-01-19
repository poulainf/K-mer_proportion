#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use Smart::Comments;
use autodie;
use warnings FATAL => 'uninitialized';

# Input files
my $seq_file = shift;
my $n = shift;
my ($outname) = $seq_file =~ m/(.*)_cds/g;

### Seq analysis
my %seq_for= read_fasta( $seq_file );

### Split values
my @keys = keys %seq_for;
my $n_seq = ( @keys / $n) + 1;

### Print output fasta file
SPLIT_FILES:
for my $i(1..$n){
	my $outfile = "$outname-$i.cds.fasta";
	open my $out, '>', $outfile;

	for	my $y(1..$n_seq){ 
		my $id = shift @keys;
		say {$out} ">", "$id";
		
		for  my $in_id (keys %{$seq_for{$id}}){
			print {$out} "$seq_for{$id}{$in_id}";
		}
		
		print {$out} "\n";
		
	}
	next SPLIT_FILES;
}

sub read_fasta { ### Elapsed time |===[%]

	my $infile = shift ;	

	### Reading input file: $infile

	open my $in, '<', $infile;

	my $seq_id;
	my $seq;
	my $seq_virus;

	LINE:
				
	while (my $line = <$in> ){
		chomp $line;
		
		# at each '>' char...
		if ($line =~ m/^>[\w\d]*\|([\d\w\.]+\ )/xms) {
			
		 # add current seq to hash  (if any)
		 if($seq) {
			$seq_for{$seq_virus}{$seq_id} = $seq;
			$seq = q{};
			}
		
			# extract new seq_id
			$seq_id = $1;
			
			($seq_virus) = $seq_id =~ m/([\w\d]+\.*\d*)_\w+_.*/xms;
		
			next LINE;
			
		}
		# elongate current seq (seq can be broken sevetal lines)
		$seq .= $line;
	}

	#add last seq to hash (if any)
	$seq_for{$seq_virus}{$seq_id} = $seq if $seq;

	close $in;
	return %seq_for;
}
