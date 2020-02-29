#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use List::AllUtils 'mesh';
use Algorithm::Numerical::Shuffle qw/shuffle/;
use List::Util qw(sum);
use warnings FATAL => 'uninitialized';


unless (@ARGV==2) {
	die << "EOT";
Usage: $0 <infile.fasta> <infile.gb>

This tool estimate the trikmer enrichment ratio of CCN NCC NNCCNN NCT NCG NCA NGG GGN NNGGNN NNTGNN NNAGNN NNCGNN TCN NTC NNTCNN NTT NTA NTG NGA GAN NNGANN NNTANN NNAANN NNCANN NGT NGC NAC NAT NAA NAG .
	1-The software format fasta sequences/
	2-The software collect for each ID the associated virus Specie, Family, Genus, Group, from file.gb
	3-The software randomize sequence and count nucleotide motif (100 times)
	4-The software count observed nucleotide motif in original sequence
	6-The software calculate the enrichment ratio between obs/exp
	7-The software print the result with Virus informations

Example: $0 test.fasta test.gb
EOT
}

### Input files
my $seq_file = shift;
my $info_file = shift;
my ($outname) = $seq_file =~ m/(.*)\.fasta/g;


### Seq analysis
my %seq_for= read_fasta( $seq_file );
my %info_for = taxon_id( $info_file );

### Start printing 
my $outfile1 = join "" , $outname,'_Kmer-counting100_global.tsv';
open my $out, '>',  $outfile1;
my @kmers = split / /, q(CCN NCC NNCCNN NCT NCG NCA NGG GGN NNGGNN NNTGNN NNAGNN NNCGNN TCN NTC NNTCNN NTT NTA NTG NGA GAN NNGANN NNTANN NNAANN NNCANN NGT NGC NAC NAT NAA NAG );
say {$out} join  "\t", "#ID", @kmers, "Group", "Genus", "Family", "Specie";

# Build Kmer pathern

my @kmer_pathern = split / /, q(CC. .CC ..CC.. .CT .CG .CA .GG GG. ..GG.. ..TG.. ..AG.. ..CG.. TC. .TC ..TC.. .TT .TA .TG .GA GA. ..GA.. ..TA.. ..AA.. ..CA.. .GT .GC .AC .AT .AA .AG );
							 


### Start Virus genome analysis

for my $id (keys %seq_for){ ### Elapsed time |===[%]
		
	print {$out} "$id\t";
	my $seq = $seq_for{$id};
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
		
		my $result = sprintf "%5.3f", $compt/$pExp;
		print {$out} "$result\t";	
	}
	
	say {$out} join "\t",($info_for{$id}{"group"}, $info_for{$id}{"famille"}, $info_for{$id}{"genus"}, $info_for{$id}{"specie"});
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

	LINE:
		
		while (my $line = <$in> ){
		chomp $line;
		
		# at each '>' char...
		if ($line =~ m/>(\w*.\d)/xms) {
			
		 # add current seq to hash  (if any)
		 if($seq) {
			$seq_for{$seq_id} = $seq;
			$seq = q{};
			}
		
			# extract new seq_id
			$seq_id = $1;
			next LINE;
			
		}
		# elongate current seq (seq can be broken sevetal lines)
		$seq .= $line;
	}

	#add last seq to hash (if any)
	$seq_for{$seq_id} = $seq if $seq;

	close $in;

	return %seq_for;
}

sub taxon_id {		
	my $infile = shift ;	

	### Reading input file: $infile
	### Elapsed time |===[%]

	open my $in, '<', $infile;

	
	my $ID;
			my $group;
			my $famille;
			my $specie;
			my $genus;
			my $orga_line = q{};
			my %orga_liner;
	LONE:		
	while (my $line = <$in> ){
		chomp $line;

		if ($line =~ m/VERSION\ *(\w*.\d)\ *.*/xms) { 
			$ID = $1; 
			next LONE;
		};
		
		if ($line =~ m/^\ \ ORGANISM/xms){$orga_line = $line ; next LONE }; 		
		if ($line =~ m/^REFERENCE/xms){
			$orga_line =~ s/\ \ +/\ /g;
			$orga_liner{$ID} = $orga_line;
			if ($orga_line =~ m/ORGANISM\ ([\w\ \d-]*).*\ Viruses/xms){$specie = $1};
			if ($orga_line =~ m/Viruses;([\w\ ]*)/xms) {$group = $1};
			if ($orga_line =~ m/positive/xms){$group = join " ",$group,"+"};
			if ($orga_line =~ m/negative/xms){$group = join " ",$group,"-"};
			if ($orga_line =~ m/Viruses.*([A-Z]\w*ae)/xms ) {$famille=$1};
			if ($orga_line =~ m/Viruses;.*;([\ *\w*]*)\./xms ){$genus = $1};
			if ($group =~ m/ssRNA\ positive/xms) {$group = "ssRNA viruses +"};
			if ($group =~ m/unclassified\ bacterial\ viruses/xms) {$group = "Phages"};
			if ($group =~ m/unclassified\ virophages/xms) {$group = "Phages"};
			if ($group =~ m/unclassified\ phages/xms) {$group = "Phages"};
			if ($group =~ m/unclassified\ archaeal\ viruses/xms) {$group = "Archaeal viruses"};
			if ($group =~ m/Retro/xms) {$group = "ssRNA-RT"};
			if ($group =~ m/Virus\ families\ not\ assigned\ to\ an\ order/xms) {$group = "Unassigned Viruses"};
			if ($group =~ m/unclassified\ viruses/xms) {$group = "Unassigned Viruses"};
			if ($group =~ m/unassigned\ viruses/xms) {$group = "Unassigned Viruses"};
			if ($group =~ m/unclassified\ RNA\ viruses/xms) {$group = "Unassigned Viruses"};
			if ($group =~ m/environmental\ samples/xms) {$group = "Environmental samples"};
		    if (!($group)) {$group = "Unassigned Viruses"};
			$orga_line = q{};
				
			next LONE;
			}; 
			
		if ($orga_line) { $orga_line .= $line; next LONE };
		 
		if ($line =~ m/^\/\//xms){
			$info_for{$ID} = { 
					group => $group,
					famille => $famille,
					genus => $genus,
					specie=>$specie};
			
			$orga_line = q{};
			
			next LONE;
			
			};
			
		next LONE;
	};
	
	$info_for{$ID} = { 
					group => $group,
					famille => $famille,
					genus => $genus,
					specie=>$specie};
		
	close $in;


open my $out, '>',  "ORGANISM_lines.txt";

	for	my $id (keys %orga_liner){
		my $org_line = $orga_liner{$id};
		say {$out} join "\t", $id, $org_line;
	}
	
	close $out;
	
	return %info_for;
}
