#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use Algorithm::Numerical::Shuffle qw/shuffle/;
use List::Util qw(sum);
use warnings FATAL => 'uninitialized';


unless (@ARGV==1) {
        die << "EOT";
Usage: $0 <infile.fasta> <infile.gb>

This tool estimate the trikmer enrichment ratio of CCN NCC NNCCNN NCT NCG NCA NGG GGN NNGGNN NNTGNN NNAGNN NNCGNN TCN NTC NNTCNN NTT NTA NTG NGA GAN NNGANN NNTANN NNAANN NNCANN.
For genic sequences. 
        1-The software format fasta sequences/
        2-The software collect for each ID the associated virus Specie, Family, Genus, Group, from file.gb
        3-The software randomize sequence and count nucleotide motif (1000 times)
        4-The software count observed nucleotide motif in original sequence
        6-The software calculate the enrichment ratio between obs/exp
        7-The software print the result with Virus informations

Example: $0 test.fasta test.gb
EOT
}

### Input files
my $seq_file = shift;
my ($outname) = $seq_file =~ m/(.*)\.fasta/g;


### Seq analysis
my %seq_for= read_fasta( $seq_file );


### Start printing
my $outfile1 = join "" , $outname,'_Kmer-counting1000_gene.tsv';
open my $out, '>',  $outfile1;
my @kmers = split / /, q(CCN NCC NNCCNN NCT NCG NCA NGG GGN NNGGNN NNTGNN NNAGNN NNCGNN TCN NTC NNTCNN NTT NTA NTG NGA GAN NNGANN NNTANN NNAANN NNCANN);
say {$out} join  "\t", "#ID","Virus",@kmers, "Gene", "Protein", "Location","Group", "Genus", "Family", "Specie";

# Build Kmer pathern

my @kmer_pathern = split / /, q(CC. .CC ..CC.. .CT .CG .CA .GG GG. ..GG.. ..TG.. ..AG.. ..CG.. TC. .TC ..TC.. .TT .TA .TG .GA GA. ..GA.. ..TA.. ..AA.. ..CA..);


### Start Virus genome analysis

for my $id (keys %seq_for){		### Elapsed time |===[%]
	my $seq = $seq_for{$id}{"seq"};
	my $virus = $seq_for{$id}{"virus"};
		
	print {$out} "$id\t$virus\t";
	my $seq_length = length $seq;
		
	for my $kmer (@kmer_pathern){
		my $l = length $kmer;
		my $compt = 0;
		
		for my $y (1..1000){	
			my $random_seq = join '', shuffle split //, $seq;
			for (my $i = 0 ; $i < $seq_length; $i +=3 ){
				my $codon = substr ($random_seq, $i, $l);
				$compt++ if ($codon =~ m/$kmer/xms);
			}
		}
		
		my $pExp = $compt/1000;
		$compt = 0;
		
		for (my $i = 0 ; $i < $seq_length; $i +=3 ){
			my $codon = substr ($seq, $i, $l);
			$compt++ if ($codon =~ m/$kmer/xms);
			
		}
		
		my $result;
		if ($pExp) {$result = sprintf "%5.3f", $compt/$pExp }  else {$result = 0};
		print {$out} "$result\t";	
	}
	
	if ($seq_for{$id}{"seq_info"} =~ m/\[gene=(\w*)\]/xms){ print {$out} "$1\t"} else {print {$out} "\t"};
	if ($seq_for{$id}{"seq_info"} =~ m/\[protein=([\w\d\ \.]+)\]\ /xms){ print {$out} "$1\t"} else {print {$out} "\t"};
	if ($seq_for{$id}{"seq_info"} =~ m/\[location=([\w\d\ \.\(\)\<\>]+)\]\ /xms){ print {$out} "$1\t"} else {print {$out} "\t"};
	say {$out} join "\t","yoyo";
}

close $out;


# Functions
sub read_fasta {

	my $infile = shift ;	

	### Reading input file: $infile
	### Elapsed time |===[%]

	open my $in, '<', $infile;

	my $seq_id;
	my $seq;
	my $virus;
	my $seq_info;



	LINE:
		
		
	while (my $line = <$in> ){
		chomp $line;
		
		# at each '>' char...
		if ($line =~ m/>lcl\|(.*)_cds_([\w\d\.]+)\ (.*)/xms) {
			
		 # add current seq to hash  (if any)
		 if($seq) {
			$seq_for{$seq_id} = {
				seq => $seq,
				virus => $virus,
				seq_info => $seq_info,
			};
			
			$seq = q{};
			}
		
			# extract new seq_id
			$seq_id = $2;
			$virus = $1;
			$seq_info = $3;
			next LINE;
			
		};
		# elongate current seq (seq can be broken sevetal lines)
		$seq .= $line;
	};

	#add last seq to hash (if any)
	$seq_for{$seq_id} = {
			seq => $seq,
			virus => $virus,
			seq_info => $seq_info,
	};

	close $in;
	
	# %seq_for

	return %seq_for;
}
