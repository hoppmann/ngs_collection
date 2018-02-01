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
	-outDir		name of the output directory (default = "04-filter")
	-lookup		TODO: if choosen only the gene of choice will be looked uped (to be implemented)
	-patID		the patient ID to query for, if not given the patID is extracted automatically
	-trio		if chosen an trio analysis will be performed
	-parentalIDs	the IDs of the parents to be used for filtering (requires exact 2 entries space seperated)
	-printCmd	if chosen instead of executing the filter, the filter command is printed to screen


=head1 DESCRIPTION

=cut

# parse command line options
my $help;
my $man;
my $panel;
my $dbName;
my $outFolder;
my $lookup;
my $torrent;
my $patID;
my $trio;
my @parentalIDs;
my $cmdOnly;

my $result = GetOptions (
				"help"			=> \$help,
				"man"			=> \$man,
				"panel=s"		=> \$panel,
				"database=s"		=> \$dbName,
				"outDir=s"		=> \$outFolder,
				"lookup=s"		=> \$lookup,
				"patID=s"		=> \$patID,
				"trio"			=> \$trio,
				"parentalIDs=s{2}"	=> \@parentalIDs,
				"printCmd"		=> \$cmdOnly,
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

#####################################################################
######## extract patient ID from database and prepare out dir #######
#####################################################################


my @patID;
if (! $patID){

	my $cmd = "$gemini query -q \"select * from samples \" $dbName";
	my @tmpArray = split ("\n", `$cmd`);
	foreach my $tmpLine (@tmpArray){
		my @splitLine = split ("\t", $tmpLine);
		push (@patID, $splitLine[2]);

	}
} else {
	push (@patID, $patID);
}


######## create out dir
# check if out dir is given
$outFolder = "04-filter" if (! $outFolder);


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



#
# say join ("\t", keys %panelInfo);
# die;



#############################################################
######## for each patient ID run all filtering steps ########
#############################################################
foreach my $pat (@patID){



	##########################################################################
	######## read in column file, and remove bash variable assignment ########
	##########################################################################

	######## get all ids of interest (index, relative 1, relative 2 )
	my @allIDs = ($patID, @parentalIDs);

	## read in column file
	open (COLUMNS, "<", $columns);
	my $cols = <COLUMNS>;
	close (COLUMNS);

	# remove bash variable assignment
	chomp $cols;
	$cols =~ s/columns=//g;
	$cols =~ s/"//g;


	if (!$trio){
		$cols =~ s/\(gt_depths\)\.\(\*\)/gt_depths\.$pat/g;
		$cols =~ s/\(gts\)\.\(\*\)/gts.$pat/g;
		$cols =~ s/\(gt_types\)\.\(\*\)/gt_types.$pat/g;
	} elsif ($trio) {


		#### remove all instances of gt_* and add new ones per ID
		$cols =~ s/, \(gt_depths\)\.\(\*\)//g;
		$cols =~ s/, \(gts\)\.\(\*\)//g;
		$cols =~ s/, \(gt_types\)\.\(\*\)//g;

		## for each gt_* run loop over all IDs in oder to cluster them
		foreach my $curID (@allIDs){
			$cols .= ", gt_depths.$curID";
		}

		foreach my $curID (@allIDs){
			$cols .= ", gts.$curID";
		}

		foreach my $curID (@allIDs){
			$cols .= ", gt_types.$curID";
		}

	}



	# if noAutoOut option chosen save pat ID as outDir

	my $outDir = "$outFolder/$pat";

	# create out dirs
	mkpath $outDir if ( ! -d $outDir);






	#####################################################################
	######## prepare neccessary information for filter production #######
	#####################################################################
	# out file names; depeding on the filter situation
#	my $outDomHigh;
#	my $outDomMed;
#	my $outRecHigh;
#	my $outRecMed;
#	if (! $trio ) {
#		$outDomHigh = "$outDir/$pat-DOM-HIGH.out";
#		$outDomMed = "$outDir/$pat-DOM-MED.out";
#		$outRecHigh = "$outDir/$pat-REC-HIGH.out";
#		$outRecMed = "$outDir/$pat-REC-MED.out";
#	} else {
#
#		$outDomHigh = "$outDir/$pat-DOM-HIGH-trio.out";
#		$outDomMed = "$outDir/$pat-DOM-MED-trio.out";
#		$outRecHigh = "$outDir/$pat-REC-HIGH-trio.out";
#		$outRecMed = "$outDir/$pat-REC-MED-trio.out";
#
#	}
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
#		push (@sqlQuery, "(aaf_adj_exac_nfe <= 0.01 or aaf_adj_exac_nfe is null or aaf_adj_exac_nfe >= 0.99) \\\n");

		push (@sqlQuery, "(aaf_gnomad_afr < 0.01 or aaf_gnomad_afr > 0.99 or aaf_gnomad_afr is null) \\\n");
		push(@sqlQuery, "and (aaf_gnomad_eas < 0.01 or aaf_gnomad_eas > 0.99 or aaf_gnomad_eas is null) \\\n");
		push(@sqlQuery, "and (aaf_gnomad_nfe < 0.01 or aaf_gnomad_nfe > 0.99 or aaf_gnomad_nfe is null) \\\n");
		push(@sqlQuery, "and (aaf_gnomad_sas < 0.01 or aaf_gnomad_sas > 0.99 or aaf_gnomad_sas is null) \\\n");





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
					push (@sqlQuery, "and ( \\\n");
					$setOr = 1;
				} else {
					push (@sqlQuery, "or ");
				}


		######## Transcripts are excluded due to the fact, that VEP only emmits the transcript with the most severe impact.
#				if ($panelInfo{$curGene}) {
#					push (@sqlQuery, "(gene = '$curGene' and transcript = '$panelInfo{$curGene}') \\\n");
#				} else {
					push (@sqlQuery, "(gene = '$curGene') \\\n");
#				}

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
		# prepare command for single patien case


		#### Legend of gt_type numbers
		# 0 = HOM_REF
		# 2 = ./. which means uncalled, happens due to merge steps. Can mostly be interpreted as HOM_REF
		# 1 = HET
		# 3 = HOM_ALT



		if (! $trio){
			if (@$curFilter[0] eq "REC"){

				$postSql .= "--gt-filter \"gt_types.$pat == HOM_ALT\" \\\n";

			} elsif (@$curFilter[0] eq "DOM" || @$curFilter[0] eq "COMP" ){

				$postSql .= "--gt-filter \"gt_types.$pat == HET\" \ \\\n";

			} elsif (@$curFilter[0] eq "LUP" ) {

			}

		} elsif ($trio){

			# recessive case: parents are HET and patient HOM_ALT
			if (@$curFilter[0] eq "REC"){

				$postSql .= "--gt-filter \"gt_types.$pat == HOM_ALT and gt_types.$parentalIDs[0] == HET and gt_types.$parentalIDs[1] == HET \" \\\n";


			# dominant case: parents are either HOM_REF or ./. and patient HET
			# will also be used for COMP_HET filtering in later case
			} elsif (@$curFilter[0] eq "DOM" || @$curFilter[0] eq "COMP" ){



				$postSql .= "--gt-filter \"gt_types.$pat == 1 and (gt_types.$parentalIDs[0] == 0 or gt_types.$parentalIDs[0] == 2) and (gt_types.$parentalIDs[1] == 0 or gt_types.$parentalIDs[1] == 2)\" \ \\\n";

			# in case of lookup extract all variants without any restrictions
			} elsif (@$curFilter[0] eq "LUP" ) {

			}
		}




		###################
		######## general stuff
		$postSql .= "--header \\\n";

		# if panel file given add panel name to output
		# mark out file if trio data was analysed
		my $outFileName;
		if ($panelName) {
			if ($trio){
				$outFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]-trio";
			}else {
				$outFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]";
			}
		} else {
			if ($trio){
				$outFileName = "$pat-@$curFilter[0]-@$curFilter[1]-trio";
			} else {
				$outFileName = "$pat-@$curFilter[0]-@$curFilter[1]";
			}
		}


		$postSql .= "$dbName > $outDir/$outFileName.out";



		######################
		######## combine different query parts and execute query
		my $filterCmd = $preSql;
		foreach (@sqlQuery){
			$filterCmd .= $_;
		}
		$filterCmd .= $postSql;


		# print out info
		say "Filtering $outFileName";

#	 	say $filterCmd . "\n";
		# execute filter command
		if ($cmdOnly){
			say $filterCmd;
		} else {
			system ($filterCmd);
		}

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
			if ($trio){
				$inFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]-trio";
			} else {
				$inFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]";
			}
		} else {
			if ($trio) {
				$inFileName = "$pat-@$curFilter[0]-@$curFilter[1]-trio";
			} else {
				$inFileName = "$pat-@$curFilter[0]-@$curFilter[1]";
			}
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

}
# print out info that filtering finished succesfully
say "Filtering successfully finished!";
