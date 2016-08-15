#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;

system("clear");


my $outDir = "03-extract_exons";
my $fileName = "$outDir/transcript.pos";
my $fileOut = "$outDir/panel.bed";


# read in input file
open (IN, "< $fileName") or die $!;
my @inputFile;
while (<IN>) {
	chomp;
	my @tmp = split(/\t/, $_);
	push(@inputFile, \@tmp);
}
close (IN);


# for each gene produce  hash conatining gene as key and list of cds as value
my %chr = ();
my @allCds;
my %genes = ();
foreach my $line (@inputFile){

	# extract key
	my @line = @$line;
	my $curGene = $line[0];
		
	# remember chr for later use
	my $chr = $line[1];
	${$chr{$curGene}} = \$chr;
	
	# form array of exon positions
	my @exonPositions;
	for (my $counter = 4; $counter < scalar(@line); $counter++){
		my $exon = $line[$counter];
		push @exonPositions, $exon;
	}

	# save positions to hash
	push @{$genes{$curGene}}, \@exonPositions;
}


# check for longest transcript
my %longTrans = ();
for my $key (keys %genes){
	
	# make counter and compare variable
	my $length = 0;
	my $compare = 0;
	my @longest = 0;
	my @listOfList =  @{$genes{$key}};

	foreach my $list (@listOfList){

		# get length of list
		$length = @$list;
		if ($length > $compare){
			$compare = $length;
			@longest = @$list;
		}
	}
	
	# save longest transcript of each gene
	push (@{$longTrans{$key}}, @longest);
}


## save gene in bed file
## open output file
open (OUT, "> $fileOut");

foreach my $curGene (keys %longTrans){
	
	# get chromosom
	my $chrRef = ${$chr{$curGene}};
	my $chr = $$chrRef;
	
	#get array of exons
	my @exons = @{$longTrans{$curGene}};
	
	# save each exon as bed entry and add exon number
	my $counter = 1;
	my $numberExons = @exons;
	foreach my $exon (@exons){
	
		my ($start, $end) = split(/-/, $exon);
	
		say OUT "chr$chr\t$start\t$end\t$curGene\t$counter/$numberExons";
	
		# increment counter
		$counter++;
	}
}



close (OUT);














































