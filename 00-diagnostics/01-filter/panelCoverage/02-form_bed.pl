#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;

#system("clear");


my $outDir = $ARGV[0];
my $fileName = "$outDir/$ARGV[1]-transcripts.pos";

# create directory for bed files
mkdir("$outDir/bed");
#my $fileOut = "$outDir/Panel.bed";

# define a flanking region of each gene to include splice site variations
my $flanking = 25;


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


foreach my $curGene (keys %longTrans){
	
	# create variable for line collection for later print
	my @outLines;
		
	# get chromosom
	my $chrRef = ${$chr{$curGene}};
	my $chr = $$chrRef;
	
	#get array of exons
	my @exons = @{$longTrans{$curGene}};
	chomp @exons;	
	# save each exon as bed entry and add exon number
	my $counter = 1;
	my $numberExons = @exons;
	
	foreach my $exon (@exons){
		chomp $exon;
		my ($start, $end) = split(/-/, $exon);
		
		
		#add flanking to exons to include splice sites
		# exclude genes with no known exon
		next if ($start eq "" || $end eq "");
		$start = $start - $flanking;
		$end = $end + $flanking;
		
		# save in file		
		my $line = "chr$chr\t$start\t$end\t$curGene\t$counter/$numberExons";
		push (@outLines, $line);


		# increment counter
		$counter++;
	}
	
	# if lines for outfile available print in file
	
	
	
	if (@outLines){
		# open out file
		my $fileOut = "$outDir/bed/$curGene.bed";
		open (OUT, "> $fileOut");
		say OUT join ("\n", @outLines);
		close (OUT);
	}
}


















































