#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use DBI;
use File::Basename;
use File::Path;



=head1 NAME

=head1 SYNOPSIS

Options:

	-help		brief help message
	-man		full documentation
	-database	name of the gemini database to use (mandatory option)
	-panel		name of the panel file to use (tab seperated list of gene and corresponding transcript)
	-outDir		name of the output directory (default = "out")
	-lookup		if choosen additionaly all chosen genes are extrated disregard any inheritance
	-torrent	if chosen excludes all frameshift_variants
	-NoAutoOut	do not choose out dir name from input DB ()
	

=head1 DESCRIPTION

=cut

# parse command line options
my $help;
my $man;
my $panel;
my $dbName;
my $outDir;
my $lookup;
my $torrent;
my $autoOut;

my $result = GetOptions (	
				"help"		=> \$help,
				"man"		=> \$man,
				"panel=s"	=> \$panel,
				"database=s"	=> \$dbName,
				"outDir=s"	=> \$outDir,
				"lookup"	=> \$lookup,
				"torrent"	=> \$torrent,
				"noAutoOut"	=> \$autoOut,
			);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);


# define standard variales
my $gemini = "/data/ngs/bin/gemini/anaconda/bin/gemini";
my $columns = "/h/hoppmann/scripts/ngs/filtering/00-columns.sh";


####################
######## check panel file and database


# check if database is given and if exists
$dbName or die "No database given. Use option -database\n";
die "Database $dbName not found.\n" if (! -f $dbName);



##########################################################################
######## read in column file, and remove bash variable assignment ########
##########################################################################

open (COLUMNS, "<", $columns);
my $cols = <COLUMNS>;
close (COLUMNS);

# remove bash variable assignment
chomp $cols;
$cols =~ s/columns=//g;
$cols=~ s/"//g;

#$cols = ("gene");

#####################################################
######## read in panel file and save in hash ########
#####################################################

# read in panel file
my $panelName;
my %panelInfo;
if ($panel){
	open (FILEIN, "<", $panel);
	my @curFile = <FILEIN>;
	close(FILEIN);

	## save read in panel in hash key = gene value = transcript
	foreach (@curFile){
		# remove white space and split, then save in hash
		chomp $_;
		next if ($_ eq "");
		my @splitLine = split ("\t", $_);
		my $geneName = $splitLine[0];
		my $transcript = $splitLine[1];
		
		# remove version number from transcript
		if ( $transcript){
			$transcript =~ s/\..*$//;
		
			# save info in hash for later use
			$panelInfo{$geneName} = $transcript;
		} else {
			$panelInfo{$geneName} = "";
		}
	}
	
	
	$panelName = basename($panel);
	$panelName =~ s/\..*$//;
	
	say ("Panel file \"$panelName\" successfully read in.");
}



#####################################################################
######## extract patient ID from database and prepare out dir #######
#####################################################################

my $cmd = "$gemini query -q \"select * from samples \" $dbName";
my @patLine = split ("\t", `$cmd`);
my $pat = $patLine[2];

######## create out dir
# check if out dir is given
$outDir = "out" if (! $outDir);

# if autoOut option chosen save pat ID as outDir
$outDir = "$outDir/$pat" if ( ! $autoOut);

# create out dirs
mkpath $outDir if ( ! -d $outDir);

#####################################################################
######## prepare neccessary information for filter production #######
#####################################################################
# out file names
my $outDomHigh = "$outDir/$pat-DOM-HIGH.xls";
my $outDomMed = "$outDir/$pat-DOM-MED.xls";
my $outRecHigh = "$outDir/$pat-REC-HIGH.xls";
my $outRecMed = "$outDir/$pat-REC-MED.xls";

# arrays containing filter specifications
my @domHigh = ("DOM", "HIGH");
my @domMed = ("DOM", "MED");
my @recHigh = ("REC", "HIGH");
my @recMed = ("REC", "MED");
my @lookupHigh = ("LUP", "HIGH");
my @lookupMed = ("LUP", "MED");
my @lookupAll = ("LUP", "ALL");
my @compHetHigh = ("COMP", "HIGH");
my @compHetMed = ("COMP", "MED");

# save in hash
my @difFilters;
if ($panel){
	@difFilters = (\@domHigh, \@domMed, \@recHigh, \@recMed, \@compHetHigh, \@compHetMed, \@lookupHigh, \@lookupMed, \@lookupAll);
} else {
	@difFilters = (\@domHigh, \@domMed, \@recHigh, \@recMed, \@compHetHigh, \@compHetMed);
}


#@difFilters = (\@domHigh);







#########################################
######## compute comp het filter ########
#########################################



foreach my $curFilter (@difFilters){

	#init variables
	my $preSql;
	my @sqlQuery;
	my $postSql;
	
	# prepare general things
	$preSql = "$gemini query -q \" \\\n";
	
	################
	######## SQL query
	
	push (@sqlQuery, "select $cols from variants where \\\n");
	push (@sqlQuery, "(aaf_adj_exac_nfe <= 0.01 or aaf_adj_exac_nfe is null) \\\n");
	
	# exclude all frameshift variants if torrent runs are fitlered
	push (@sqlQuery, "and impact != 'frameshift_variant' \\\n") if ($torrent);
	
	
	
	######################
	######## check if impact severity is HIGH or MED
	if (@$curFilter[1] eq "HIGH"){
		push (@sqlQuery, "and (impact_severity == 'HIGH' or impact_severity is null) \\\n");

	} elsif (@$curFilter[1] eq "MED") {

		push (@sqlQuery, "and (impact_severity == 'MED' or impact_severity == 'HIGH' or impact_severity is null) \\\n");

	} elsif (@$curFilter[1] eq "ALL") {
	} 



	##############################
	####### include panel genes and transcripts (if transcript available)
	
	if ($panel){
		my $setOr = 0;
		foreach my $curGene (keys %panelInfo){
			# every time besides first time add the or expression
			if (!$setOr) {
				push (@sqlQuery, "and (");
				$setOr = 1;
			} else {
				push (@sqlQuery, "or ");
			}
		
			if ($panelInfo{$curGene}) {
				push (@sqlQuery, "(gene = '$curGene' and transcript = '$panelInfo{$curGene}') \\\n");
			} else {
				push (@sqlQuery, "(gene = '$curGene') \\\n");
			}

		}
		
		# close bracket to conclude block of panel genes
		push (@sqlQuery, ") \\\n");
	}


	
	
	#################
	######## close sql query
	push (@sqlQuery, "\" \\\n");



	###################
	######## non SQL filters
	# check if dominant or recessive inheritance or just plain lookup
	if (@$curFilter[0] eq "REC"){
	
		$postSql .= "--gt-filter \"gt_types.$pat == HOM_ALT\" \\\n";
	
	} elsif (@$curFilter[0] eq "DOM" || @$curFilter[0] eq "COMP" ){
	
		$postSql .= "--gt-filter \"gt_types.$pat == HET\" \ \\\n";
	
	} elsif (@$curFilter[0] eq "LUP" ) {
	
	}


	###################
	######## general stuff
	$postSql .= "--header \\\n";
	
	# if panel file given add panel name to output
	my $outFileName;
	if ($panelName) {
		$outFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]";
	} else {
		$outFileName = "$pat-@$curFilter[0]-@$curFilter[1]";
	}
	$postSql .= "$dbName > $outDir/$outFileName.xls";



	######################
	######## combine different query parts and execute query
	my $filterCmd = $preSql;
	foreach (@sqlQuery){
		$filterCmd .= $_;
	} 
	$filterCmd .= $postSql;
	

	# print out info	
	say "Filtering $outFileName";

#say $filterCmd;
	# execute filter command
	system ($filterCmd);


}


######################################
######## comp Het post filter ########
######################################

foreach my $curFilter (@difFilters){
	
	# if filter not specified for comp het scip to next entry
	next if (@$curFilter[0] !~ "COMP");
	
	
	# get input fiel name
	# if panel file given add panel name to output
	my $inFileName;
	if ($panelName) {
		$inFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]";
	} else {
		$inFileName = "$pat-@$curFilter[0]-@$curFilter[1]";
	}
	my $inFile = "$outDir/$inFileName.out";


	
	# check if inFile realy exist
	( -e $inFile) or die "File $inFile for comp het extraction not found.";
	
	
	# open file to extract all genes
	open (INFILE, "<", $inFile);
	my @fileIn = <INFILE>;
	close (INFILE);
	
	# create out file
	open (my $fileOut, "> $inFile") or die "Can't open file for writing.";
	
	# extract gene names
	my @genes;
	foreach my $curLine (@fileIn) {
		chomp $curLine;
		if ( $curLine !~ m/^gene/) {
			my @splitLine = split ("\t", $curLine);
			push (@genes, $splitLine[0]);
		} else {
			# print out header
			say $fileOut $curLine;
		}
	}
	
	
	# extract all uniq genes	   
	my @uniqueGenes = do { my %seen; grep { !$seen{$_}++ } @genes };

	# for each unique gene check if there are 2 or more entries. if so save as comp het
	foreach my $curGene (@uniqueGenes) {
		my @hits = grep {/^$curGene/} @fileIn;

		# check if more then 1 hits found. if so save.
		if ($#hits + 1 > 1){
			say $fileOut join ("\n", @hits)
		}
	}
	
	# close new comp Het file
	close ($fileOut);

}


# print out info that filtering finished succesfully
say "Filtering successfully  finished!";




































