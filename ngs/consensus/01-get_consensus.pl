#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils  'first_index';
use Scalar::Util qw(looks_like_number);


=head1 NAME

=head1 SYNOPSIS

Options:

	-help		brief help message
	-man		full documentation
	-vcfIn		the vcf input file containing different caller reads
	-out		the name of the output vcf file (default = out.vcf)
	-consensusName	name of the consensus call in the vcf file (default = consensus)
	-minimumHits	integer defining how many call need at least detect a variant (default = 2)

=head1 DESCRIPTION

=cut


######## define variabales ########
my $help;
my $man;
my $vcfFile;
my $vcfOut = "out.vcf";
my $consensusName = "consensus";
my $minimumHits = 2;



# parse command line options
my $result = GetOptions (	
		"help"			=> \$help,
		"man"			=> \$man,
		"vcfIn=s" 		=> \$vcfFile,
		"out=s"			=> \$vcfOut,
		"consensusName=s"	=> \$consensusName,
		"minHits=i"		=> \$minimumHits,
		);

pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);

######## make entry that program is running
say "Start running consesus calling.";

#### open vcf in and out file ####
die if ($vcfFile eq "");
open (VCFIN, "< $vcfFile") or die "Can't open $vcfFile.";
open (my $outFile, "> $vcfOut") or die "Can't create $vcfOut.";

# define variables
my %callerMap;

#### read in file and make entry depending on number of caller that found the specific variant
while (<VCFIN>) {

	# remove white space at end of line and empty lines if any
	chomp;
	next if ($_ eq "");
	
	#########################################
	######## manipulate and save header lines
	
	
	# if header found check if it is header for columns 
	# if so modify else print it out
	if (/^#/) {
		
		###############
		#### check for the final header line and exclude all caller lines and replace them with consensus
		if (/^#CHROM/) {
		
			# make entry for new added consensus variable
			say $outFile '##INFO=<ID=CONS,Number=2,Type=Integer,Description="Indicating the number of calls supporting this variant">';	
			say $outFile '##INFO=<ID=ORIG,Number=1,Type=Character,Description="The variantcaller the consensus line originates from."';
		
			# split header to exclude different caller
			my @line = split (/\t/, $_);
			my @caller = splice (@line, 9);
			push (@line , $consensusName);
			say $outFile join ("\t", @line);
						
			# save caller as key in hash for late identification which caller was used
			my $counter = 0;
			foreach (@caller){
			$callerMap{ $counter } = $caller[$counter];
			$counter++;
			}
			next;
		}

		############
		#### print out all header
		say $outFile $_;
		next;
	}
	
	
	
	######################
	######## variant lines
	 
	######## check each called variant and extract best quality call, use those values 
	######## and add in FORMAT how many out of all caller found that variant
	
	# split line to extract colums for editing
	my @line = split (/\t/, $_);
	
	# save first 7 parts to keep in out variable since they won't get manipulated
	my @lineOut = splice (@line,0, 7);

	
	# divide rest into INFO, FORMAT and caller part for distinct processing
	my $info = shift @line; 
	my $format = shift @line;
	my @caller = @line;


	## search for position of quality entry (GQ = genotype quality in phred scale) in format array
	my @format= split (":", $format);
	my $indexGQ = first_index { /GQ/ } @format;

	############
	#### compare different GQ scores and take the entry with the highes score
	my @compare;
	my $failCase;
	# split entries to get GQ values and save them in array for comparison		
	foreach my $key (keys %callerMap) {
		my $callerName = $callerMap{$key};
		my @callerValues = split (":", $caller[$key]);
		$compare [$key] = $callerValues[$indexGQ];
		
		## check if more then one entry available. -> Valid entries available even if GQ score 
		# not available. needed  for later use
		if ($#callerValues > 1) {
			$failCase = $key;
		}
		
	}
	
	#### compare score values in array, and extract index of largest (index = key in hash)
	## save the number of caller calling this variant
	## in case no GQ value available take last entry with real entries
	my $counter = 0;
	my $indexBest = "";
	my $comp = -1;
	my $numberCaller = 0;
	foreach my $curCallerScore (@compare){
		if (looks_like_number($curCallerScore)) {
			## get the number of caller calling this variant
			$numberCaller++;

			if ($curCallerScore > $comp) {
				$comp = $curCallerScore;
				$indexBest = $counter;
			}
		}
		$counter++;
	}

	# save out call with largest score. If indexBest not awailable (no values for GQ availabel)
	# take last entry with valid entries (saved in failCase)
	if ($indexBest eq ""){
		$indexBest = $failCase;
	}
	
	my $outCaller = $caller[$indexBest];
	
	
	
	
	
	########################################################
	######## gather information to save in vcf file ########
	########################################################	
	
	# check that call is supported by sufficient caller then print it out
	
	if ($numberCaller >=	 $minimumHits){
		## add to info field the caller the current line originate from
		# save where the called variant values origin from 
		my $infoOut = $info . ";ORIG=$callerMap{$indexBest}"; 
	
		# save the number of caller that found that variant and the total number of caller used
		$infoOut .= ";CONS=$numberCaller/" . ($#compare + 1);	

		#### save current consensus line in file
		say $outFile join ("\t", @lineOut) . "\t" . $infoOut . "\t" . $format . "\t" . $outCaller;
	}
}

close (VCFIN);
close ($outFile);





#### say that program is finished
say "Consensus calling successuflly finished."









































