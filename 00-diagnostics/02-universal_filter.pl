#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use DBI;
use File::Basename;
use File::Path;
use experimental 'smartmatch';


=head1 NAME

=head1 SYNOPSIS

Options:


	-help		brief help message

	-man		full documentation

	-database	name of the gemini database to use (mandatory option)

	-panel		name of the panel file to use (tab seperated list of gene and corresponding transcript)

	-outDir		name of the output directory (default = "04-filter")

	-noLUP		if chosen skips the lookup step in panel filterin

	-lookup	if choosen only the gene of choice will be looked uped, multiple genes lookups by comma seperated list

	-screen		if chosen the lookup outout is printed on screen instead of a file

	-all		if chosen all genotype information of all patients are shown in an lookup, not only the chosen patient

	-patID		the patient ID to query for, if not given the patID is extracted automatically

	-trio		if chosen an trio analysis will be performed

	-parentalIDs	the IDs of the parents to be used for filtering (requires exact 2 entries space seperated)

	-cmdOnly	if chosen instead of executing the filter, the filter command is printed to screen


=head1 DESCRIPTION

=cut

# parse command line options
my $help;
my $man;
my $panel;
my $dbName;
my $outFolder;
my $lup;
my $lookup;
my $all;
my $screen;
my $torrent;
my $patID;
my $trio;
my @parentalIDs;
my $cmdOnly;
my $maf=0.01;

my $result = GetOptions (
		"help"				=> \$help,
		"man"				=> \$man,
		"panel=s"			=> \$panel,
		"database=s"		=> \$dbName,
		"outDir=s"			=> \$outFolder,
		"lookup=s"			=> \$lookup,
		"screen"			=> \$screen,
		"all"				=> \$all,
		"patID=s"			=> \$patID,
		"trio"				=> \$trio,
		"parentalIDs=s{2}"	=> \@parentalIDs,
		"cmdOnly"			=> \$cmdOnly,
		"noLUP"				=> \$lup,
);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);


# define standard variales
my $gemini = "/data/programs/bin/ngs/gemini/anaconda/bin/gemini";
my $columns = "/data/programs/scripts/hoppmann/00-diagnostics/pipeline/00-columns.sh";


####################
######## check panel file and database


# check if database is given and if exists
$dbName or die "No database given. Use option -database\n";
die "Database $dbName not found.\n" if (! -f $dbName);







#######################################################################
######## extract patient ID (from database) and prepare out dir #######
#######################################################################


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
$outFolder = "09-filter" if (! $outFolder);










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
		# remove all non-alphanumerical letters to avoid bug
		$geneName =~ s/\W//g;

		## split mode of inheritance in case of more then one MOI
		my @mois;
		if ($splitLine[1] && $splitLine[1] ne "" && $splitLine[1]!~ /^\s$/) {
			@mois = split (",", $splitLine[1]);
			@mois = grep(s/\s*$//g, @mois);
		}
		if (@mois){
			# save info in hash for later use
			$panelInfo{$geneName} = \@mois;
		} else {
			$panelInfo{$geneName} = "";
		}
	}


	$panelName = basename($panel);
	$panelName =~ s/\..*$//;

	say ("Panel file \"$panelName\" successfully read in.");
}








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
	my @cols = <COLUMNS>;
	close (COLUMNS);

	# remove bash variable assignment
	chomp @cols;


	## exclude comment lines starting with "#"
	my @temp;
	foreach my $curCol (@cols) {
		if ($curCol eq "" || $curCol =~ /^\s*#/ || $curCol =~ /^$/) {
			next;
		} else {
			push @temp, $curCol;
		}
	}
	my $cols = join ", ", @temp;




	####################
	######## adapt genotyp output
	if ($trio) {

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




	} elsif ($all) {

		#### use wildcards as found in file

	} else {
		$cols =~ s/\(gt_depths\)\.\(\*\)/gt_depths\.$pat/g;
		$cols =~ s/\(gts\)\.\(\*\)/gts.$pat/g;
		$cols =~ s/\(gt_types\)\.\(\*\)/gt_types.$pat/g;
	}


	# if noAutoOut option chosen save pat ID as outDir

	my $outDir = "$outFolder/$pat";

	# create out dirs
	mkpath $outDir if ( ! -d $outDir);


















	#####################################################################
	######## prepare neccessary information for filter production #######
	#####################################################################

	# arrays containing filter specifications
	my @domHigh = ("DOM", "HIGH");
	my @domMed = ("DOM", "MED");
	my @recHigh = ("REC", "HIGH");
	my @recMed = ("REC", "MED");
	my @lookupHigh = ("LUP", "HIGH");
	my @lookupMed = ("LUP", "MED");
	my @lookupAll = ("LUP", "ALL");
	my @compHetHigh = ("COMP", "HIGH");
	my @compHetMed = ("COMP", "MED_HIGH");
	my @xlHigh = ("XL", "HIGH");
	my @xlMed = ("XL", "MED");

	# save in hash
	my @difFilters;
	if ($panel){
		if ($lup){
			@difFilters = (\@domHigh, \@domMed, \@recHigh, \@recMed, \@xlHigh, \@xlMed, \@compHetHigh, \@compHetMed, \@lookupHigh, \@lookupMed, \@lookupAll);
		} else {
			@difFilters = (\@domHigh, \@domMed, \@recHigh, \@recMed, \@xlHigh, \@xlMed, \@compHetHigh, \@compHetMed);
		}
	} else {
		@difFilters = (\@domHigh, \@domMed, \@recHigh, \@recMed, \@compHetHigh, \@compHetMed);
	}

	# @difFilters = (\@compHetMed, \@compHetHigh);





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

		push (@sqlQuery, "(aaf_gnomad_afr < " . $maf . " or aaf_gnomad_afr > " . (1-$maf) . " or aaf_gnomad_afr is null) \\\n");
		push(@sqlQuery, "and (aaf_gnomad_eas < $maf or aaf_gnomad_eas > " . (1-$maf) . " or aaf_gnomad_eas is null) \\\n");
		push(@sqlQuery, "and (aaf_gnomad_nfe < $maf or aaf_gnomad_nfe > " . (1-$maf) . " or aaf_gnomad_nfe is null) \\\n");
		push(@sqlQuery, "and (aaf_gnomad_sas < $maf or aaf_gnomad_sas > " . (1-$maf) . " or aaf_gnomad_sas is null) \\\n");















		######################
		######## check for chosen impact_severity
		if (@$curFilter[1] eq "HIGH"){
			push (@sqlQuery, "and (impact_severity == 'HIGH' or impact_severity is null) \\\n");

		} elsif (@$curFilter[1] eq "MED") {

			push (@sqlQuery, "and (impact_severity == 'MED' or impact_severity is null) \\\n");

		} elsif (@$curFilter[1] eq "MED_HIGH") {

			## for the case of COMP_HET all HIGH and MED variants are needed for usefull filtering
			push (@sqlQuery, "and (impact_severity == 'HIGH' or impact_severity == 'MED' or impact_severity is null) \\\n");

		} elsif (@$curFilter[1] eq "ALL") {

		}





















		##############################
		####### include panel genes and transcripts (if transcript available)

		if ($panel){
			my $setOr = 0;


			foreach my $curGene (keys %panelInfo){

				#####################
				#### add genes that have concurrend MOI

				## check for dominat and COMP_HET genes
				if (@$curFilter[0] eq "DOM" ) {


					#### following lines should probably replace other if clauses due to experimental design
					# if ((grep /^AD$/, @{$panelInfo{$curGene}}) || ! exists ${$panelInfo{$curGene}}[0] ){
					# 	my $pan = $panelInfo{$curGene};
					# 	say $curGene;
					# }



					if ("AD" ~~ @{$panelInfo{$curGene}} || $panelInfo{$curGene} eq ""){
						# every time besides first time add the or expression
						if (!$setOr) {
							push (@sqlQuery, "and ( \\\n(gene = '$curGene') \\\n");
							$setOr = 1;
						}  else {
							push (@sqlQuery, "or (gene = '$curGene') \\\n");
						}
					}




					## check for recessive genes
				} elsif (@$curFilter[0] eq "REC" || @$curFilter[0] eq "COMP") {
					if ("AR" ~~ @{$panelInfo{$curGene}} || $panelInfo{$curGene} eq ""){
						# every time besides first time add the or expression
						if (!$setOr) {
							push (@sqlQuery, "and ( \\\n(gene = '$curGene') \\\n");
							$setOr = 1;
						}  else {
							push (@sqlQuery, "or (gene = '$curGene') \\\n");
						}
					}



					## check if lookup if so add all genes
				} elsif (@$curFilter[0] eq "LUP") {
					# every time besides first time add the or expression
					if (!$setOr) {
						push (@sqlQuery, "and ( \\\n(gene = '$curGene') \\\n");
						$setOr = 1;
					}  else {
						push (@sqlQuery, "or (gene = '$curGene') \\\n");
					}




					## check if x-Linked case is chosen
				} elsif (@$curFilter[0] eq "XL") {
					if ("XL" ~~ @{$panelInfo{$curGene}} || "XLD" ~~ @{$panelInfo{$curGene}} || "XLR" ~~ @{$panelInfo{$curGene}}){
						# every time besides first time add the or expression
						if (!$setOr) {
							push (@sqlQuery, "and ( \\\n(gene = '$curGene') \\\n");
							$setOr = 1;
						}  else {
							push (@sqlQuery, "or (gene = '$curGene') \\\n");
						}
					}
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
		# prepare command for single patient case


		#### Legend of gt_type numbers
		# 0 = HOM_REF
		# 2 = ./. which means uncalled, happens due to merge steps. Can mostly be interpreted as HOM_REF
		# 1 = HET
		# 3 = HOM_ALT

		#####################
		######## no trio analysis
		if (! $trio){

			## recessive case: HOM_ALT variant
			if (@$curFilter[0] eq "REC"){

				$postSql .= "--gt-filter \"gt_types.$pat == HOM_ALT\" \\\n";

				## dom or comp_het case: HET variant
			} elsif (@$curFilter[0] eq "DOM" || @$curFilter[0] eq "COMP" ){

				$postSql .= "--gt-filter \"gt_types.$pat == HET\" \ \\\n";

				## Lookup case: all variants that aren't ./.
			} elsif (@$curFilter[0] eq "LUP" || @$curFilter[0] eq "XL") {
				$postSql .= "--gt-filter \"gt_types.$pat != 2\" \ \\\n";
			}


			####################
			######## trio analysis

		} elsif ($trio){


			if (@$curFilter[0] eq "REC"){

				# recessive case: parents are HET and patient HOM_ALT
				$postSql .= "--gt-filter \"gt_types.$pat == HOM_ALT and gt_types.$parentalIDs[0] == HET and gt_types.$parentalIDs[1] == HET \" \\\n";



			} elsif (@$curFilter[0] eq "DOM" ){

				# dominant case: parents are either HOM_REF or ./. and patient HET
				$postSql .= "--gt-filter \"gt_types.$pat == 1 and (gt_types.$parentalIDs[0] == 0 or gt_types.$parentalIDs[0] == 2) and (gt_types.$parentalIDs[1] == 0 or gt_types.$parentalIDs[1] == 2)\" \ \\\n";

			} elsif (@$curFilter[0] eq "COMP") {

				# COMP_HET case: patient is HET and parents mustn't be HOM_ALT
				$postSql .= "--gt-filter \"gt_types.$pat == 1 and (gt_types.$parentalIDs[0] != 3) and (gt_types.$parentalIDs[1] != 3)\" \ \\\n";




			} elsif (@$curFilter[0] eq "LUP" ) {

				# in case of lookup extract all variants without any restrictions
				$postSql .= "--gt-filter \"(gt_types.(*).(!= 2).(any))\" \ \\\n";

			}
		}















		###################
		######## general stuff / naming


		$postSql .= "--header \\\n";

		# if panel file given add panel name to output
		# mark out file if trio data was analysed
		my $outFileName;
		if ($panelName) {
			if ($trio){
				# if trio and panel
				$outFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]-trio";
			} else {
				# if not trio but panel
				$outFileName = "$pat-$panelName-@$curFilter[0]-@$curFilter[1]";
			}
		} else {
			if ($trio){
				# if trio and panel
				$outFileName = "$pat-@$curFilter[0]-@$curFilter[1]-trio";
			} else {
				# if not trio but panel
				$outFileName = "$pat-@$curFilter[0]-@$curFilter[1]";
			}
		}


		$postSql .= "$dbName ";
		if (! $screen ){
			$postSql .= "> $outDir/$outFileName.out";
		}


		######################
		######## combine different query parts and execute query
		my $filterCmd = $preSql;
		foreach (@sqlQuery){
			$filterCmd .= $_;
		}
		$filterCmd .= $postSql;


		# print out info
		if (!$lookup){
			say "Filtering $outFileName";
		}
		#	 	say $filterCmd . "\n";
		# execute filter command
		if (! $cmdOnly && ! $lookup) {
			# in case a panel is filtered check if a gene name is specified, else gemini has an sql error
			if ($panel) {
				if ($filterCmd =~ "gene ="){
					system ($filterCmd);
				} else {
					say "Filter skipped since there are no genes in this panel.";
				}
			} else {
				system ($filterCmd);
			}
		} elsif ($cmdOnly && !$lookup){


			say "$filterCmd\n\n\n"
		}

	}

















	#################################
	######## performe lookup ########
	#################################


	if ($lookup){
		#### get genes for lookup
		my @lookup = split(/,/,$lookup);


		#### print out Status
		say "Running lookup for " . join (", ",@lookup);


		#### make lookup query
		my $cmd = "$gemini query -q \\\n";
		$cmd .= "\" select $cols from variants where \\\n";

		#### add genes for lookup
		my $setOr = 0;
		my $geneList;
		foreach my $curGene (@lookup){
			$geneList .= "_$curGene";
			if (!$setOr) {
				$cmd .= "(gene == '$curGene') \\\n";
				$setOr = 1;
			}  else {
				$cmd .= "or \\\n(gene == '$curGene') \\\n";
			}
		}
		$cmd .= "\" \\\n";

		#### exclude variants where patID has 2 (./.) as genotypes
		if (!$all){
			$cmd .= "--gt-filter \"gt_types.$pat != 2\" \ \\\n";
		}


		#### add header
		$cmd .= "--header \\\n";

		#### add database
		$cmd .= "$dbName \\\n";

		#### prepare out file
		my $outFile = "$outDir/$pat\_LUP$geneList.out";

		say $outFile;
		if (! $screen){
			$cmd .= "> $outFile";
		}


		if ($lookup && ! $cmdOnly){
			system $cmd;
		} else {
			say $cmd;
		}

	}


















	######################################
	######## comp Het post filter ########
	######################################
	if (!$lookup ){
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

			#### create out files
			my $fileOutName = $inFile;
			$fileOutName =~ s/-raw//;
			open (my $fileOutNonDeNovo, "> $fileOutName") or die "Can't open file $fileOutName for writing.";

			#### create out file for cases with de novo muations
			my $fileOutDeNovo;
			if ($trio){
				my $fileOutNameDeNovo = $inFile;
				$fileOutNameDeNovo =~ s/.out/-DeNovo.out/;
				open ($fileOutDeNovo, "> $fileOutNameDeNovo") or die "Can't open file $fileOutNameDeNovo for writing.";
			}




			########################
			######## sort each line in hash to corresponding gene entry and prepare header
			my %allGenes;
			my @gtIndex;
			foreach my $curLine (@fileIn) {
				chomp $curLine;

				# if not header split line and save in hash based on gene names
				if ( $curLine !~ m/^gene/) {
					my @splitLine = split ("\t", $curLine);
					my $curGene = $splitLine[0];
					push (@{$allGenes{$curGene}}, $curLine);
					# push (@genes, $splitLine[0]);

				} else {

					#### if header, print out header
					say $fileOutNonDeNovo $curLine;
					say $fileOutDeNovo $curLine if ($trio);

					############
					#### extract index of patient and parental gt_types in case of trio filtering
					if ($trio){
						my @header = split ("\t", $curLine);

						# get index of patient gt_types
						my $curGt = "gt_types." . $patID;
						push (@gtIndex, grep { $header[$_] =~ $curGt} 0...$#header);

						# get index of parental gt_types
						foreach (@parentalIDs){
							my $curGt = "gt_types." . $_;
							push (@gtIndex, grep { $header[$_] =~ $curGt} 0...$#header);
						}
					}
				}
			}









			####################
			######## check for each gene if multiple fitting variants are available

			if ($trio){

				foreach my $curGene (keys %allGenes){

					## retrieve gene array for readability
					my @curGeneVarList = @{$allGenes{$curGene}};

					## define variable used to remember indeces
					my @matchDeNovoLines;
					my @matchNonDeNovoLines;



					## only use genes with more then one entry.
					## possible comp_het genes
					my $variantsInGene = @curGeneVarList;



					#### exclude all genes with "None" -> they mess up the analysis since in real there are diffenrent genes included here. This nomenclature is introduced by the way GATK writes indels... Needs to be fixed anyways
					next if ( $variantsInGene <= 1 || $curGene =~ /None/);


					#######################################
					######## run COMP_HET in case of a trio

					## compare each line if inheritance pattern matches
					for (my $counter = 0; $counter < $variantsInGene; $counter++) {

						#### get gt_types for current line to compare with remaining lins of that gene
						my $firstLine = $curGeneVarList[$counter];
						my @splitLine = split "\t", $firstLine;
						my @gtTypeCurLine;
						foreach (@gtIndex){
							push @gtTypeCurLine, $splitLine[$_];
						}



						#### iterate over remaining lines and check if inheritance matches restrictions
						for (my $counter2 = $counter + 1; $counter2 < $variantsInGene; $counter2++) {


							#### get gt_types for next line
							my $secondLine = $curGeneVarList[$counter2];
							my @splitLine = split "\t", $secondLine;
							my @gtTypeNextLine;
							foreach (@gtIndex){
								push @gtTypeNextLine, $splitLine[$_];
							}


							########################################################
							######## check that inheritance is not violated ########
							########################################################

							#### case 1 -> one partent is Het in both variants
							#### case 2 -> both parents are WT (0 and 2)
							#### case 3 -> parent 1 and 2 are het in same variant
							#### case 4 -> no de novo mutaiton is allowed
							#### case 4 is omited for DeNovo-file

							#### write gtTypes explicit for ease of reading
							my $gtPar1F = $gtTypeCurLine[1];
							my $gtPar1S = $gtTypeNextLine[1];
							my $gtPar2F = $gtTypeCurLine[2];
							my $gtPar2S = $gtTypeNextLine[2];




							#####################################
							######## check cases for de_novo file
							#	###########
							if ( $gtPar1F == 1 && $gtPar1S == 1 ) {
							} elsif ( $gtPar2F == 1 && $gtPar2S == 1 ) {
								#################
								#### check case 2
							} elsif ((( $gtPar1F == 0 || $gtPar1F == 2 ) && ( $gtPar1S == 0 || $gtPar1S == 2 ))  &&   (( $gtPar2F == 0 || $gtPar2F == 2 ) && ( $gtPar2S == 0 || $gtPar2S == 2 ))  ) {
								#################
								#### check case 3
							} elsif ( ($gtPar1F == 1 &&  $gtPar2F == 1) || ($gtPar1S == 1 &&  $gtPar2S == 1) ) {
								#### if no exclusion case fits save index to write out in basic file
							} else {
								push (@matchDeNovoLines, $counter);
								push (@matchDeNovoLines, $counter2);
							}




							#########################################
							######## check cases for non-de_novo file
							#	###########
							#	#### check case1
							#	## par1 is HET in both positions
							if ( $gtPar1F == 1 && $gtPar1S == 1 ) {
								# say "Par 1 HET in both $gtPar1F $gtPar1S";

							} elsif ( $gtPar2F == 1 && $gtPar2S == 1 ) {
								#################
								#### check case 2
							} elsif ( (( $gtPar1F == 0 || $gtPar1F == 2 ) && ( $gtPar1S == 0 || $gtPar1S == 2 ))  &&   (( $gtPar2F == 0 || $gtPar2F == 2 ) && ( $gtPar2S == 0 || $gtPar2S == 2 )) ) {
								#################
								#### check case 3
							} elsif ( ($gtPar1F == 1 &&  $gtPar2F == 1) || ($gtPar1S == 1 &&  $gtPar2S == 1) ) {
								###########
								#### check case 4
							} elsif ( (( $gtPar1F == 0 || $gtPar1F == 2 ) && ( $gtPar1S == 0 || $gtPar1S == 2 )) || (( $gtPar2F == 0 || $gtPar2F == 2 ) && ( $gtPar2S == 0 || $gtPar2S == 2 )) ) {
								#### if no exclusion case fits save index to write out in basic file
							} else {
								# say "This should be correct: $gtPar1F $gtPar1S || $gtPar2F $gtPar2S ";
								push (@matchNonDeNovoLines, $counter);
								push (@matchNonDeNovoLines, $counter2);
							}



						} ## end of second line iteration
					} ## end of first line iteration





					######## extract unique indeces and write out line
					my @uniqueLines = do { my %seen; grep { !$seen{$_}++ } @matchNonDeNovoLines };
					foreach my $index (@uniqueLines){
						say $fileOutNonDeNovo $curGeneVarList[$index];
					}


					######## extract unique indeces and write out line
					@uniqueLines = do { my %seen; grep { !$seen{$_}++ } @matchDeNovoLines };
					foreach my $index (@uniqueLines){
						say $fileOutDeNovo $curGeneVarList[$index];
					}


				} ## end iteration over each gene

				# close comp Het out file
				close ($fileOutDeNovo);
				close ($fileOutNonDeNovo);

			}	#### end of trio COMP_HET filter




			if (! $trio ){

				foreach my $curGene (keys %allGenes){


					## retrieve gene array for readability
					my @curGeneVarList = @{$allGenes{$curGene}};

					## define variable used to remember indeces
					my @matchDeNovoLines;
					my @matchNonDeNovoLines;



					## only use genes with more then one entry.
					## possible comp_het genes
					my $variantsInGene = @curGeneVarList;



					#### exclude all genes with "None" -> they mess up the analysis since in real there are diffenrent genes included here. This nomenclature is introduced by the way GATK writes indels... Needs to be fixed anyways
					next if ( $variantsInGene <= 1 || $curGene =~ /None/);

					#### print out variants for genes
					say $fileOutNonDeNovo join "\n", @curGeneVarList;

				} #### close for each gene loop


				# close comp Het out file
				close ($fileOutDeNovo) if ($trio);
				close ($fileOutNonDeNovo);
			}





		} ## end iteration over each filter
	} ## end if not lookup
} ## end iterate over each patient



















# print out info that filtering finished succesfully
say "Filtering successfully finished!";
