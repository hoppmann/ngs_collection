#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use File::Basename;


######## define external programs and resources
my $samtools = "/data/programs/bin/ngs/samtools/1.3.1/bin/samtools";
my $hg19 = "/data/public_resources/bundle/2.8/ucsc.hg19.fasta";

# get commandline arguments
my $bedFile = $ARGV[0];
my $bamFile = $ARGV[1];
my $outDir = $ARGV[2];

$bedFile or die "bed file missing.\n";
$bamFile or die "bamfile dir missing.\n";
$outDir or die "out dir missing.\n";


# define variables
my $minCov = 30;


# get name for outFile
my $fileName = basename($bedFile);
$fileName =~ s/\..*$//;
my $outFile = "$outDir/$fileName.csv";


# check if outfile exist else create it
if ( ! -d $outDir ) {
	mkdir ($outDir);
}




# get a list of all usable bam files
die "Bamfile $bamFile not available.\n" if ( ! -f $bamFile);


###################################################
#### check that gene has entry, else scip gene ####
###################################################

open (IN, "<", $bedFile);
my @bed = <IN>;
close (IN);

my $counter = 0;
foreach my $curLine (@bed) {

	chomp $curLine;
	next if ( $curLine eq '' || $curLine =~ /^\s*$/);
	$counter++;

}

die "Skipping $bedFile: emtpy file!\n" if ($counter == 0);







# define hash to save coverage information
# key = bamList;
# value = array( all lines in ped file as array )
# all lines = array ( entries in that line )
my %infoGather;

# check if bamfile exists
( -e $bamFile ) or die "Bam file $bamFile not found.";


######################
######## check if bam file has index, else create index
if ( ! -f "$bamFile.bai" ){
	# make status update
	say "#### indexing file $bamFile ####";

	#prepare index command
	my $indexCmd = "$samtools \\\n";
	$indexCmd .= "index \\\n";
	$indexCmd .= "$bamFile";

	# execute index command
	system ($indexCmd);
}

#########
######## extract coverage according to bedfile locations
########

# prepare command
my $covExtractCmd = "$samtools \\\n";
$covExtractCmd .= "bedcov \\\n";
$covExtractCmd .= "--reference $hg19 \\\n";
$covExtractCmd .= "$bedFile \\\n";
$covExtractCmd .= "$bamFile";


####### calculate mean coverage
foreach (`$covExtractCmd`) {
	chomp $_;
	my @splitLine = split ("\t", $_);
	my $mean = $splitLine[5] / ($splitLine[2] - $splitLine[1]);

	# save mean array with rest of line info
	push (@splitLine, $mean);
	push (@{$infoGather{$bamFile}}, \@splitLine);

}


###############################
######## create output ########
###############################


#write initial header
my $upperLines = "usedBam\t#goodExons\n";
my $lowerLines;
my %exonInfo = ();
my $sumCov = 0;
my $sumAll = 0;


# extract info and print interesting parts out
foreach my $curBamFile (keys %infoGather) {


	# init variables
	my $goodExon = 0;
	my $badExon = 0;
	my @curBamArrayOfLines = @{$infoGather{$curBamFile}};


	# make first header for lower line
	$lowerLines .= $curBamFile . "\n";
	my $goodLines = 0;
	my $badLines = 0;


	# for each bam file extract for each exon if good or bad covered
	foreach (@curBamArrayOfLines){
		my @curLine = @$_;

		# check if exon is sufficient covered
		my $meanCov = ($curLine[$#curLine]);

		if ($meanCov >= 30) {

			$goodExon++;

			# save good exons for later print
			$goodLines .= $curLine[4] . "\t" . sprintf ("%.2f", $meanCov) . "\n";

		} else {

			$badExon++;

			#save bad exons for later print
			$badLines .= $curLine[4] . "\t" . sprintf ("%.2f", $meanCov) . "\n";

			push (@{$exonInfo{$curLine[4]}}, basename($curBamFile));
		}
	}





	# percentual overall coverage
	my $numAllExon = $#curBamArrayOfLines+1;
	$sumCov += $goodExon;
	$sumAll += $numAllExon;

	# prepare entries for upper lines
	$upperLines .= $curBamFile . "\t" . "$goodExon/" . ($numAllExon) . "\n" ;

	# prepare entries for lowe lines
	$lowerLines .= "Good covered exons:\nExon\tmeanCov\n" . $goodLines . "\n"  . "Bad covered exons:\nExon\tmeanCov\n" . $badLines . "\n\n";



}


	my $thirdLines = "Exon\t#BadBam\n";
	my $totalBam = 1;
	for my $key (keys %exonInfo){
		my $badBam = @{$exonInfo{$key}};
		$thirdLines .= $key . "\t" . $badBam . "/" . $totalBam . "\t" . join ("\t", @{$exonInfo{$key}}) . "\n";
	}

	## check that not divided by 0
	$sumAll = 1 if ($sumAll == 0);

	my $coverLine = "Overall coverage = " . sprintf ("%.2f", ($sumCov / $sumAll) * 100 ) . "%\n";
	#### open file to write in
	open (my $fhOut, ">", $outFile) or die "Can't create resultfile $outFile\n";
	say $fhOut $coverLine;
	say $fhOut $upperLines;
	say $fhOut $thirdLines;
	say $fhOut $lowerLines;
	close $fhOut;
