#!/usr/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use Parallel::ForkManager;
use File::Basename;

=head1 NAME

=head1 SYNOPSIS

Options:

	-help		brief help message
	-man		full documentation
	-input		file containing names of xsq files
	-cpu		number of available CPUs (default: autodetection)
	-bam		list with bamfiles to use (only to use if (1) - (3) is skipped)
	-logfile	set name of logfiel
	-start		first task to start with (see below)
	-end		last task to execute (see below)
	-Individuals	choose this option if combined vcf file should als be extracted into individuals besides families.
	-family		if chosen supplement file (tab deliminated) containing family id and individual ids (#famID  indID) all lines starting with "#" will be excluded
	-cmdonly	choose this option to run pipeline without executing commands
	-out		name of out folder (default = "out")
	-tools		gives out the directories of the tools and refrence files
	-hard		if chosen hard filtering is used instead of VQSR
	
	Parts of the Pipeline
		SOLID specific:
			(1) 	XSQ-Conversion
			(2) 	Mapping
			(3) 	sam-bam conversion and merging;
		GATK pipeline:
			(4) 	mark duplicates
			(5) 	indel realignment
			(6) 	base quality score recalibration BQSR
			(7) 	create matrices
			(8) 	variant calling and VQSR
			(9) 	extraction of individuals and or families
			(10)	VEP-annotation for later use with gemini
			(11)	Load families in gemini DB


Examples:

	######## Ped file
	file containing the following columns tab deliminated in shown order

	#famID  indID   patID   matID   m/f     pheno
	oran    group1  group3  -9      2       2
	oran    group2  grpup3  -9      1       2
	oran    group3  group5  group6  1       1
	oran    group5  -9      -9      1       1
	oran    group6  -9      -9      2       1

	fam 	=> family ID, needed to group by families and for identification for Gemini load
	indID	=> ID for each individual
	patID	=> ID of the father
	matID	=> ID of the mother
	m/f	=> identifyier if male (1) or female (2)
	pheno	=> identifier if case (2) or control (1)

	missing is "-9" or "null"


=head1 DESCRIPTION

=cut

# parse command line options
my $help;
my $man;
my $input;
my $cpu;
my $log = "logfile.log";
my $start = 1;
my $end = 20;
my $cmdOnly;
my $extractInd;
my $pedIN;
my $bamIN;
my $outFolder = "out";
my $tools;
my $hard;

# save commandline options for storage in Log file
my @chosenOptions = @ARGV;

my $result = GetOptions (	"help"			=> \$help,
				"man"			=> \$man,
				"input=s"		=> \$input,
				"cpu=i"			=> \$cpu,
				"logfile=s"		=> \$log,
				"start=i"		=> \$start,
				"end=i"			=> \$end,
				"cmdOnly"		=> \$cmdOnly,
				"Individuals!"		=> \$extractInd,
				"family=s"		=> \$pedIN,
				"bam=s"			=> \$bamIN,
				"out=s"			=> \$outFolder,
				"tools!"		=> \$tools,
				"hard!"			=> \$hard,
			);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);


###############################
######## define pathes ########

#programs
my $ngsPlumbingDir = "/data/ngs/bin/ngs_plumbing-0.13.1/xsq-convert";
my $shrimp2Dir = "/data/ngs/bin/SHRiMP_2_2_3";
my $samtools = "/data/ngs/bin/samtools-1.1/bin/samtools";
my $picardDir = "/data/ngs/bin/picard-tools-1.119";
my $bedtoolsDir = "/data/ngs/bin/bedtools2-2.20.1/bin";
my $vcftools = "/data/ngs/bin/vcftools_0.1.12b/bin/vcftools";
my $VEP = "/data/ngs/bin/variant_effect_predictor/variant_effect_predictor.pl";
my $gemini = "/data/ngs/bin/gemini/bin/gemini";
my $vtMaster = "/data/ngs/bin/vt-master/vt";


# References
my $hg19Dir = "/data/ngs/programs/bundle_2.8/ucsc.hg19.fasta";
my $indelMills = "/data/ngs/programs/bundle_2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf";
my $indel1kgp = "/data/ngs/programs/bundle_2.8/1000G_phase1.indels.hg19.vcf";
my $snps1kgp = "/data/ngs/programs/bundle_2.8/1000G_phase1.snps.high_confidence.hg19.vcf";
my $dbsnp = "/data/ngs/programs/bundle_2.8/dbsnp_138.hg19.vcf";
my $exonFile = "/data/ngs/exome_targets/SureSelect_V4_S03723314_Regions.bed";
my $hapmap = "/data/ngs/programs/bundle_2.8/hapmap_3.3.hg19.vcf";
my $omni = "/data/ngs/programs/bundle_2.8/1000G_omni2.5.hg19.vcf";



#GATK
my $gatk_dir = "/data/ngs/bin/GATK/3.3-0/";
#my $gatk_dir = "/data/ngs/programs/GenomeAnalysisTK-2.8-1-g932cd3a/";
my $gatk = $gatk_dir . "GenomeAnalysisTK.jar";
my $gatk_key = $gatk_dir . "yong.li_uniklinik-freiburg.de.key";

#other
my $tmpDir = "/scratch/global/tmp";

#### open Log-file and ERROR-file
open (my $LOG, "> $log") or die "\nCan't open $log.\n";
open (my $ERROR, "> ERROR.log") or die "\nCan't open file ERROR.log.";

#### define global variables
#create directories for output
my %dirs;

######## Write chosen options in log
say_out ($LOG, "Options chosen:");
my $chosenOptions =  join (" ", @chosenOptions);
say_out ($LOG, "$chosenOptions \n");

######## if tools chosen give out tool pathes and then die;
if ($tools){
	say "######## TOOLS ######## \n$ngsPlumbingDir\n$shrimp2Dir\n$samtools\n$picardDir\n$picardDir\n$bedtoolsDir\n$vcftools\n$VEP\n$gemini \n\n";
	say "######## References ######## \n$hg19Dir\n$indelMills\n$indel1kgp\n$dbsnp\n$exonFile \n";
	die "\n";
}


#########################################
######## get available resources ########
#########################################

#### Add starting date and time
my $date = localtime ();
say_out ($LOG, "Pipeline started: $date.\n");
######## check available CPUs if not specified
say"";
if (!$cpu){
	say_out($LOG, "Using automatically defined number of CPUs.");
	$cpu = `cat /proc/cpuinfo |grep processor |wc -l`;
	chomp $cpu;
	say_out ($LOG, "Available CPUs:\t$cpu");
} else {
	say_out($LOG, "Using user defined number of CPUs.\nAvailableCPUs:\t$cpu");
}

######### extract available memory

open (IN, "< /proc/meminfo");
my $totalMem = <IN>;
close IN;
chomp $totalMem;
my @totalMem = split (/ /, $totalMem);
$totalMem = int (($totalMem[$#totalMem-1]/1000000)-10);
say_out  ($LOG, "Available RAM:\t".$totalMem."GB\n");


######## make general out folder
$dirs{'out'} = "out";
if (!-d "$dirs{'out'}"){
	mkdir $dirs{'out'};
}

############################################################################
######## check for/read in existance of files depending to be given ########
############################################################################


###############
######## check and read in ped file
#define global variable
my %family;
my %ped;

if ($pedIN) {
		
	#open and read inped-like-file
	open (INFILE, "< $pedIN") or die_command ("Can't open $pedIN. \n$!");
	my @fileIn = <INFILE>;
	close INFILE;
	chomp @fileIn;

	#declare variables
	%family = ();
	
	# assign individual to famID using hash of arrays			
	foreach my $line (@fileIn) {
		#skip if line is header or empty
		if ($line =~ /^#/ or $line =~ /^\W$/) {
			next;
		} else {
			my @line = split (/\t/, $line);
			push (@{$family{$line[0]}}, $line[1]); 
			push (@{$ped{$line[1]}},  $line[2], $line[3], $line[4], $line[5]);
		}
	}
}


#################################
######## read in list of Bamfiles if option is chosen

#init variable for global use	
my %bamFilesFromFile = ();


if ($bamIN) {
	# open list containing bam files if file exists
	open (BAMIN, "< $bamIN") or die_command ("Can not open $bamIN.");
	while (<BAMIN>) {
	
	#skip if line is empty or header
		if ($_ =~ /^#/ or $_ =~ /^$/) {
			next;
		} else {
			
			chomp $_;
			my $file = $_;

			#check if file exists, else exit and give out Error Message
			if (!-e $_) {
				die_command ("$_ does not exist.");
			}

			#extract patients name				
			my $_ = basename($_);
			$_ =~ s/.bam$//;
#			push (@patients, $_);
				
			#save file destination for later use
			$bamFilesFromFile{$_} = $file;
		}
	}
	close BAMIN;
}

########################################
######## convert XSQ to csfasta ########
########################################

#folder for xsq-convert step	
####create output directory for xsq-convert step
make_dir ("xsq-convert", "01-xsq-convert", "out");

if ($start > 1 or $end < 1) {
	say_out ($LOG, "***** (1) Skipping xsq conversion");
} else {
	
	#check if bam option is used
	(!$bamIN) or die_command ("Use option '-bam' only if (1) - (3) is skipped.");
	
	#open file with input xsq if starting with (1)
	if (!$input){
		if ($start < 4) {
			die_command("Input obligatory");
		}
	} elsif (!-f $input){
		die_command("Inputfile not found.");
	}


	open (INPUT, "< $input") or die "Can't open $input\n";
	my @inputFile = <INPUT>;
	close INPUT;
	chomp @inputFile;
	
																	#TODO include cmd forming in execution loop
	#### prepare command for xsq-conversion
	my @cmd;
	my @currentFile;
	my $count = 1;
	foreach my $xsq (@inputFile){
		my $cmd_xsqConvert;
		$cmd_xsqConvert = "$ngsPlumbingDir $xsq --format FASTA -s color --dir $dirs{'xsq-convert'} ";
		push (@cmd, $cmd_xsqConvert);
		push (@currentFile, $xsq);
		$count++;
	}

	#make log entry
	my $print = "######## Converting files ########";
	my $string;
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out($LOG, "\n$string\n$print\n$string\nStart running xsq-convert in up to $cpu threads.");

	#get Version number for xsq-convert
	my $version = `$ngsPlumbingDir -V 2>&1`;
	chomp $version;
	say_out($LOG, "Version: $version.\n");

	#### set maximum number of threads ####

	my $pm = Parallel::ForkManager->new($cpu);


	#fork and execute xsq-convert
	for (my $count = 0; $count <= $#cmd; $count++) {

		#### make log entry
		say_file ($LOG, "$cmd[$count]");
	
		#generate pid for forking and execute
		my $pid = $pm->start and next;
		
		$cmdOnly or system($cmd[$count]);
		say_out ($LOG, "Done with ".basename($currentFile[$count]).".");

		$pm->finish;
	}

	# waiting for all children to finish
	$pm->wait_all_children;

	#make grap-able log entry 
	say_out ($LOG, "=> Conversion done.\n");
}

#################################
######## MAPPING SHRIMP2 ########
#################################

#################
####create output directory for mapping step
make_dir ("mapping", "02-mapping", "out");

#skip if chosen
if ($start > 2 or $end < 2) {
	say_out ($LOG, "***** (2) Skipping mapping");
} else {
	
	#check if bam option is used
	(!$bamIN) or die_command ("Use option '-bam' only if (1) - (3) is skipped.");
	
	#get directories created by xsq-convert
	#open directory
	opendir my $dh, $dirs{'xsq-convert'} or die "$0: opendir: $!\n";

	#extract subdirectories and die if non exist
	my @dirs = grep {-d "$dirs{'xsq-convert'}/$_" && ! /^\.{1,2}$/} readdir($dh);
	@dirs or die_command ("No directories containing fasta exist. \n");

	################
	#### get maximum number of threads for mapping

	#### maximum by memory

	my $cpuPerThread = 4;

	#48GB for hg19 in ram; each thread needs hg19 seperately -> possible threads = total ram / 48GB
	my $workingMem = int($totalMem/48);
	if ($workingMem == 0){
		$workingMem = 1;
	}

	#### maximum by cpu
	my $cpuThreads = int($cpu/$cpuPerThread);
	if ($cpuThreads == 0) {
		$cpuThreads = 1;
	}

	#get lowest number of possible threads (either defined by ram or by number of cpu/4 (4 thread per sequence))
	my $maxThreads;
	if ($workingMem < $cpuThreads) {
		$maxThreads = $workingMem;
	} else {
		$maxThreads = $cpuThreads;
	}

	#recompute possible number of threads
	$cpuPerThread = int ($cpu / $maxThreads);

	#get SHRIMP version
	my $version = ` $shrimp2Dir/bin/gmapper-cs --help 2>&1| grep SHRiMP`;
	my @version = split(/\n/, $version);

	#Write in Log file
	my $print = "######## Mapping files using SHRiMP ########";
	my $string= ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\nVersion: $version[0]\nRunning $maxThreads threads in parallel.\nUsing $cpuPerThread threads per sequence.");

	#set maximum number of pids for SHRiMP
	my $pm = Parallel::ForkManager->new($maxThreads);
	
	my @lanes;
	foreach my $currentDir (@dirs){
	
		#if not jet exists create directory
		if  (!-d "$dirs{'mapping'}/$currentDir"){
			mkdir "$dirs{'mapping'}/$currentDir";
		}
	
	
		say_out ($LOG, "\n****Running on files in directory $currentDir\n");
	
		#save directoryname for later use
		push (@lanes, $_);
	
		#extract different csfasta files
		my @files = glob ("$dirs{'xsq-convert'}/$currentDir/*F3.csfasta");

		foreach (@files){
			select(undef, undef, undef, 0.1);

			#prepare command
			my $currentF3 = $_;
			my $currentF5 = $_;		
			$currentF5 =~ s/F3.csfasta$/F5-DNA.csfasta/;
		
			#extract current patient
			my $currentPatient = basename($_);
			$currentPatient =~ s/_\d*_F3.csfasta$//;
			say_out ($LOG, "Running on patient $currentPatient");
		
			#build command
			my $cmd = "$shrimp2Dir/bin/gmapper-cs ";
			$cmd .= "-1 $currentF3 \\\n";
			$cmd .= "-2 $currentF5 \\\n";
			$cmd .= "-N $cpuPerThread \\\n";
			$cmd .= "-p opp-in \\\n";
			$cmd .= " $hg19Dir \\\n";
			$cmd .= ">$dirs{'mapping'}/$currentDir/$currentPatient.sam \\\n";
			$cmd .= "2>$dirs{'mapping'}/$currentDir/$currentPatient.log ";
		
			#write command in log-file
			say_file ($LOG, $cmd, "\n");
		
			#start forking	
			my $pid = $pm->start and next;
		
			#execute command
			$cmdOnly or system($cmd);
		
			#end forking
			$pm->finish;
		}
	}

	# waiting for children to finish
	$pm->wait_all_children;

	#make LOG entry
	select(undef, undef, undef, 0.1);
	say_out ($LOG, "\n=> Mapping done.\n");
}

########################################################
######## convert sam to bam and merge bam Files ########
########################################################
####create output directory for mapping step
make_dir ("merge", "03-merge", "out");


if ($start > 3 or $end < 3) {
	say_out ($LOG, "***** (3) Skipping sam-bam conversion and merging");
} else {

	#check if bam option is used
	(!$bamIN) or die_command ("Use option '-bam' only if (1) - (3) is skipped.");
	
	#### get patients
	my @patients;
	
	#get directories created by xsq-convert
	#open directory
	opendir my $dh, $dirs{'mapping'} or die "$0: opendir: $!\n";

	#extract subdirectories
	my @dirs = grep {-d "$dirs{'mapping'}/$_" && ! /^\.{1,2}$/} readdir($dh);
	foreach my $currentDir (@dirs){

		#extract different csfasta files
		my @files = glob ("$dirs{'mapping'}/$currentDir/*.sam");
		foreach (@files){
		
		#extract current patient
		my $currentPatient = basename($_);
		$currentPatient =~ s/.sam$//;

			#save unique patientnames for later use
			my $controll = 0;
			foreach (@patients){
				if ($currentPatient eq $_){
					$controll = 1;
				}
			}
	
			#if patient is unique save for later use
			if ($controll == 0){
				push (@patients, $currentPatient);
			}
		}
	}
		
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if sam files exist in $dirs{'mapping'}. Possibly the mapping didn't work correct.");
	}

	#calculate numbers of thread for 4 cpus per thread
	my $cpuPerThread = 4;
	my $maxThreads = int($cpu / $cpuPerThread);

	#get samtools version
	my $version = `$samtools --help 2>&1| grep Version`;
	chomp $version;

	#Write in Log file
	my $print = "######## Converting sam to bam and merging files ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\n$version\nRunning $maxThreads threads in parallel.\nUsing $cpuPerThread threads per sequence.");

	say_out ($LOG, "#### Converting files\n");

	# init variables
	my $failure = 0;

	#######################
	#### convert sam to bam

	#set maximum number of threads to number of patients
	my $pm = Parallel::ForkManager->new($maxThreads);

	#get directories created by xsq-convert
	#open directory
	opendir $dh, $dirs{'mapping'} or die "$0: opendir: $!\n";

	#extract subdirectories
	@dirs = grep {-d "$dirs{'mapping'}/$_" && ! /^\.{1,2}$/} readdir($dh);

	#create Directory per Lane
	foreach (@dirs) {
		my $currentDir = $_;
		if  (!-d "$dirs{'merge'}/$currentDir"){
			mkdir "$dirs{'merge'}/$currentDir";
		}
	}



	#open directory with mapping results
	opendir $dh, $dirs{'mapping'} or die "$0: opendir: $!\n";
	@dirs = grep {-d "$dirs{'mapping'}" && ! /^\.{1,2}$/} readdir($dh);
	my @lanes;
	foreach (@dirs) {
		#save directoryname for later use
		push (@lanes, $_);
	}
	#sort lanes-array for subsequent checking for unique patients
	@lanes = sort @lanes;

	#open different subdirs to extract samfiles
	my @files = glob ("$dirs{'mapping'}/$lanes[0]/*.sam");

	#prepare command for samtool convert
	my @currentPath;
	foreach my $currentPatient (@patients){
	
		#make Log entries
		say_out ($LOG, "\nConverting $currentPatient.");
	
	
		my $cmd;
		foreach my $currentLane (@lanes) {
	
			say "Running on $currentLane";
			#build convert sam command
			$cmd = "$samtools \\\n";
			$cmd .= "view -b \\\n";
			$cmd .= "-@ $cpuPerThread \\\n";
			$cmd .= "$dirs{'mapping'}/$currentLane/$currentPatient.sam \\\n";
			$cmd .= "> $dirs{'merge'}/$currentLane/$currentPatient.bam ";
	
			#build sort sam command	
			my $cmdSort = "$samtools \\\n";
			$cmdSort .= "sort -@ $cpuPerThread \\\n";
			$cmdSort .= "$dirs{'merge'}/$currentLane/$currentPatient.bam \\\n";
			$cmdSort .= "$dirs{'merge'}/$currentLane/$currentPatient.sorted ";
		
			#make log entry
			say_file($LOG, "$cmd \n");
			say_file ($LOG, "$cmdSort \n");

			#execute convertion command
			if (-f "$dirs{'mapping'}/$currentLane/$currentPatient.sam") {
				my $pid = $pm -> start and next;
				$cmdOnly or system ($cmd);
				$cmdOnly or system ($cmdSort);
				$pm -> finish;
			} else {
				#make header if first error
				if ($failure == 0){
					say_file ($ERROR, "\n#### Converting sam to bam.\n");
					$failure = 1;
				}
				say_out ($ERROR, "ERROR: File $dirs{'mapping'}/$currentLane/$currentPatient.sam not found. Skipping!");
				say_file ($LOG, "ERROR: File $dirs{'mapping'}/$currentLane/$currentPatient.sam not found. Skipping!");
			}
	
		}
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#make log entry
	say_out ($LOG, "\n=> Done Converting.\n");

	###############################
	#### merging bam files into one

	#init variables
	$failure = 0;

	#set maximum number of threads
	$pm = Parallel::ForkManager->new($maxThreads);

	#write in log file
	say_out ($LOG, "#### Merging bam files and indexing\n");

	foreach my $currentPatient (@patients) {

		#build merge command
		my $check = 0;
		my $cmd;
		$cmd = "$samtools \\\n";
		$cmd .= "merge \\\n";
		$cmd .= "-f -@ $cpuPerThread \\\n";
		$cmd .= "$dirs{'merge'}/$currentPatient.bam \\\n";
		#add sam file from different lanes
		foreach my $currentLane (@lanes) {
			$cmd .= "$dirs{'merge'}/$currentLane/$currentPatient.sorted.bam \\\n";

			#check if file exists. If not skip and make ERROR entry;
			if (!-f "$dirs{'merge'}/$currentLane/$currentPatient.sorted.bam") {
				$check++;
			}
		}
	
		#build index bam command
		my $cmdIdx = "$samtools \\\n";
		$cmdIdx .= "index \\\n";
		$cmdIdx .= "$dirs{'merge'}/$currentPatient.bam ";
	
		if ($check == 0){
			say_out ($LOG, "Merging and indexing: $currentPatient\n");
			say_file ($LOG, $cmd);
			say_file ($LOG, "$cmdIdx \n");
			my $pid = $pm -> start and next;
			$cmdOnly or system ($cmd);
			$cmdOnly or system ($cmdIdx);
			$pm -> finish;
		} else {
			#make ERROR header if first error
			if ($failure == 0){
				say_file ($ERROR, "\n#### Merging bam files\n");
				$failure = 1;
			}
			say_file ($ERROR, "ERROR: Not for all lanes samfiles found. Skipping merging $currentPatient.");
			say_out ($LOG, "ERROR: Not for all lanes samfiles found. Skipping merging $currentPatient.");
		}
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "\n=> Done merging.\n");


	#exclude non converted patients from further processing
	my @temp;
	foreach (@patients) {
		my $currentPatient = $_;
		if (-f "$dirs{'merge'}/$currentPatient.bam") {
			push (@temp, $currentPatient);
		}
	}
}

###############################
######## GATK Pipeline ########
###############################

# set java ram Xmx
my $ram = 10;
my $xmx = "Xmx".$ram."g";


####################################################################
######## Picard mark PCR duplicates and add readgroup names ########
####################################################################
####create output directory for mapping step
make_dir ("duplicate", "04-mark-duplicates", "out");

#init variable for global use	
#my %bamFilesFromFile = ();

#skip if wanted
if ($start > 4 or $end < 4) {
	say_out ($LOG, "***** (4) Skipping mark duplicates and add readgroup names");
} else {

	#set variables
	my @patients;
#	my %bamFilesFromFile = ();
	
	#get bam files if bam option is chosen
	
	if ($bamIN) {
		
		@patients = keys %bamFilesFromFile;
		
	} else {
	
		###############
		####extract number and name of patients

		#extract bam files
		my @files = glob ("$dirs{'merge'}/*.bam");
		foreach (@files){

		#extract current patient
		my $currentPatient = basename($_);
		$currentPatient =~ s/.bam$//;

			#save unique patientnames for later use
			my $controll = 0;
			foreach (@patients){
				if ($currentPatient eq $_){
					$controll = 1;
				}
			}

			#if patient is unique save for later use
			if ($controll == 0){
				push (@patients, $currentPatient);
			}
		}
	
		#if still no patients available stop processing
		if (!@patients){
			die_command ("No patients available. \nCheck if out files exist in $dirs{'merging'}. Possibly the mapping didn't work correct.");
		}
	}

	###########
	####extract maximum number of possible threads
	
	#calculate numbers of thread defined by given ram (10GB per job)
	my $maxThreads = int ($totalMem / $ram);
	if ($maxThreads > $cpu) {
		$maxThreads = $cpu;
	}

	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);

	#get samtools version
	my $version = `java -$xmx -jar $picardDir/MarkDuplicates.jar --help 2>&1 |grep Version`;
	chomp $version;


	#Write in Log file
	my $print = "######## Picard mark duplicates ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\n$version\nRunning up to $maxThreads threads in parallel.");

	######## Run mark PCR duplicate

	foreach my $currentPatient (@patients){
		
		# create folder for intermediate files
		make_dir ($currentPatient, $currentPatient, "duplicate");

		#create mark duplicate command
		my $cmd = "java -$xmx -Djava.io.tmpdir=$tmpDir \\\n";
		$cmd .= "-jar $picardDir/MarkDuplicates.jar \\\n";
		if ($bamIN) {
			$cmd .= "INPUT=$bamFilesFromFile{$currentPatient} \\\n"
		} else {
			$cmd .= "INPUT= $dirs{'merge'}/$currentPatient.bam \\\n";
		}
		$cmd .= "OUTPUT=$dirs{'duplicate'}/$currentPatient/$currentPatient.dup.bam \\\n";
		$cmd .= "METRICS_FILE=$dirs{'duplicate'}/$currentPatient/$currentPatient.metrics.txt CREATE_INDEX=true \\\n";
		$cmd .= "VALIDATION_STRINGENCY=LENIENT \\\n";
		$cmd .= "2> $dirs{'duplicate'}/$currentPatient/$currentPatient.dup.log ";
	
		#create add readgroups name
		my $cmdRG = "java -$xmx -Djava.io.tmpdir=$tmpDir \\\n";
		$cmdRG .= "-jar $picardDir/AddOrReplaceReadGroups.jar \\\n";
		$cmdRG .= "INPUT=$dirs{'duplicate'}/$currentPatient/$currentPatient.dup.bam \\\n";
		$cmdRG .= "OUTPUT=$dirs{'duplicate'}/$currentPatient.bam \\\n";
		$cmdRG .= "RGLB=$currentPatient RGPL=SOLID RGPU=$currentPatient RGSM=$currentPatient \\\n";
		$cmdRG .= "CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \\\n";
		$cmdRG .= "2> $dirs{'duplicate'}/$currentPatient/$currentPatient.add.log";
	
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "*Mark duplicate\n$cmd\n\n*Adding Readgroups\n$cmdRG");
	
		#run command in parallel mode
	
		my $pid = $pm -> start and next;
			$cmdOnly or system ($cmd);
			$cmdOnly or system ($cmdRG);
			$pm -> finish;
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "\n=> Done checking for duplicates.\n");
}

#########################################
######## GATK: indel realignment ########
#########################################

####create output directory for mapping step
make_dir ("realign", "05-realign", "out");

#skip if wanted
if ($start > 5 or $end < 5) {
	say_out ($LOG, "***** (5) Skipping indel realignment");
} else {

	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'duplicate'}/*.bam");
	foreach (@files){

	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.bam$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if out files exist in $dirs{'duplicate'}. Possibly the mapping didn't work correct.");
	}
	
	###########
	####extract maximum number of possible threads
	
	#calculate numbers of thread for 4 cpus per thread (Exponential decay per thread = most decay till 4/5) or #cpu/#patients
	my $cpuByPatients = int($cpu / ($#patients+1));
	my $cpuByTesting = 4;
	my $cpuPerThread;
	if ($cpuByPatients > $cpuByTesting) {
		$cpuPerThread = $cpuByPatients;
	} else {
		$cpuPerThread = $cpuByTesting;
	}

	#calculate numbers of thread defined by ram (10GB per thread) and cpu (4 per thread)
	my $maxThreads;
	my $threadsByRam = int ($totalMem / $ram);
	my $threadsByCpu = int($cpu / $cpuPerThread);
	if ($threadsByRam > $threadsByCpu) {
		$maxThreads = $threadsByCpu;
	} else {
		$maxThreads = $threadsByRam;
	}

	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);




	# Write version in Log file
	my $version = `java -$xmx -jar $gatk --version`;
	chomp $version;

	my $print = "######## Realigning indels ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}

	say_out ($LOG, "\n$string\n$print\n$string\nGATK version: $version\nRunning $maxThreads threads in parallel.\nUsing $cpuPerThread cpu per thread.");


	foreach my $currentPatient (@patients) {
	
		# create folder for intermediate files
		mkdir "$dirs{'realign'}/$currentPatient";
	
		# Step 1 create table of indels => RealignerTargetCreator
		my $cmd = "java -$xmx -jar $gatk -K $gatk_key -et NO_ET \\\n";
		$cmd .= "-T RealignerTargetCreator \\\n";
		$cmd .= "-R $hg19Dir \\\n";
		$cmd .= "-known $indelMills \\\n-known $indel1kgp \\\n";
		$cmd .= "--filter_mismatching_base_and_quals \\\n";
		$cmd .= "--filter_bases_not_stored \\\n";
		$cmd .= "--filter_reads_with_N_cigar \\\n";
		$cmd .= "-nt $cpuPerThread \\\n";
		$cmd .= "-I $dirs{'duplicate'}/$currentPatient.bam \\\n";
		$cmd .="-o $dirs{'realign'}/$currentPatient/$currentPatient.bam.list \\\n";
		$cmd .= "-log $dirs{'realign'}/$currentPatient/$currentPatient.targetcreator.log ";

		#Step2 realign reads in these targets => IndelReader
		my $cmdRealign = "java -$xmx -jar $gatk -K $gatk_key -et NO_ET \\\n";
		$cmdRealign .= "-T IndelRealigner \\\n";
		$cmdRealign .= "--filter_mismatching_base_and_quals \\\n";
		$cmdRealign .= "--filter_bases_not_stored \\\n";
		$cmdRealign .= "--filter_reads_with_N_cigar \\\n";
		$cmdRealign .= "-I $dirs{'duplicate'}/$currentPatient.bam \\\n-R $hg19Dir \\\n";
		$cmdRealign .= "-targetIntervals $dirs{'realign'}/$currentPatient/$currentPatient.bam.list \\\n";
		$cmdRealign .= "-o $dirs{'realign'}/$currentPatient.bam \\\n";
		$cmdRealign .= "-log $dirs{'realign'}/$currentPatient/$currentPatient.realign.log ";
	
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmd\n");
		say_file ($LOG,"$cmdRealign\n");
		my $pid = $pm -> start and next;
			$cmdOnly or system ($cmd);
			$cmdOnly or system ($cmdRealign);
	 		$pm -> finish;
	
	}


	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "\n=> Indel Realignment done.\n");
}

########################################
######## GATK: BaseRecalibrator ########
########################################
####create output directory for mapping step
make_dir ("BQSR", "06-BQSR", "out");

#skip if wanted
if ($start > 6 or $end < 6) {
	say_out ($LOG, "***** (6) Skipping base quality score recalibration BQSR");
} else {
	
	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'realign'}/*.bam");
	foreach (@files){

	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.bam$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if out files exist in $dirs{'realign'}. Possibly the mapping didn't work correct.");
	}
	
	########################
	#### calculate and set maximum number of available threads
	
	#calculate numbers of thread for 4 cpus per thread (Exponential decay per thread = most decay till 4/5) or #cpu/#patients
	my $cpuByPatients = int($cpu / ($#patients+1));
	my $cpuByTesting = 4;
	my $cpuPerThread;
	if ($cpuByPatients > $cpuByTesting) {
		$cpuPerThread = $cpuByPatients;
	} else {
		$cpuPerThread = $cpuByTesting;
	}
	
	#calculate numbers of thread defined by ram (10GB per thread) and cpu (4 per thread)
	my $maxThreads;
	my $threadsByRam = int ($totalMem / $ram);
	my $threadsByCpu = int($cpu / $cpuPerThread);
	if ($threadsByRam > $threadsByCpu) {
		$maxThreads = $threadsByCpu;
	} else {
		$maxThreads = $threadsByRam;
	}
	
	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);





	# Write version in Log file
	my $version = `java -Xmx4g -jar $gatk --version`;
	chomp $version;

	#write in log file
	my $print = "######## Base recalibration ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}

	say_out ($LOG, "\n$string\n$print\n$string\nGATK version: $version\nRunning $maxThreads threads in parallel.\nUsing $cpuPerThread threads per sequence.");

	foreach my $currentPatient (@patients) {

		#generate subdirectories
		mkdir "$dirs{'BQSR'}/$currentPatient";
	
		####generate commands
		#BaseReaclibrator 1: count covariates
																		#TODO add nct and check best amound out
		my $cmdCovar = "java -$xmx -jar $gatk -l INFO \\\n";
		$cmdCovar .= "-T BaseRecalibrator -K $gatk_key -et NO_ET \\\n";
		$cmdCovar .= "-R $hg19Dir \\\n";
		$cmdCovar .= "-knownSites $dbsnp \\\n-knownSites $indelMills \\\n-knownSites $indel1kgp \\\n";
		$cmdCovar .= "-solid_nocall_strategy PURGE_READ -sMode SET_Q_ZERO_BASE_N \\\n";
		$cmdCovar .= "-I $dirs{'realign'}/$currentPatient.bam \\\n";
		$cmdCovar .= "-nct $cpuPerThread \\\n";
		$cmdCovar .= "-o $dirs{'BQSR'}/$currentPatient/$currentPatient.recal.grp \\\n";
		$cmdCovar .= "-log $dirs{'BQSR'}/$currentPatient/$currentPatient.recal.log";
	
		#BaseRecalibrator 2: analyze covariates remaining after recalibration
	
		my $cmdPost = "java -$xmx -jar $gatk -l INFO \\\n";
		$cmdPost .= "-T BaseRecalibrator -K $gatk_key -et NO_ET \\\n";
		$cmdPost .= "-R $hg19Dir \\\n";
		$cmdPost .= "-knownSites $dbsnp \\\n-knownSites $indelMills \\\n-knownSites $indel1kgp \\\n";
		$cmdPost .= "-solid_nocall_strategy PURGE_READ -sMode SET_Q_ZERO_BASE_N \\\n";
		$cmdPost .= "-BQSR $dirs{'BQSR'}/$currentPatient/$currentPatient.recal.grp \\\n";
		$cmdPost .= "-I $dirs{'realign'}/$currentPatient.bam \\\n";
		$cmdPost .= "-nct $cpuPerThread \\\n";
		$cmdPost .= "-o $dirs{'BQSR'}/$currentPatient/$currentPatient.post-recal.grp \\\n";
		$cmdPost .= "-log $dirs{'BQSR'}/$currentPatient/$currentPatient.post-recal.log";
	
		#BaseRecalibrator 3: generate before/after plots
	
		my $cmdPlot = "java -$xmx -jar $gatk -l INFO \\\n";
		$cmdPlot .= "-T AnalyzeCovariates -K $gatk_key -et NO_ET \\\n";
		$cmdPlot .= "-R $hg19Dir \\\n";
		$cmdPlot .= "-before $dirs{'BQSR'}/$currentPatient/$currentPatient.recal.grp \\\n";
		$cmdPlot .= "-after $dirs{'BQSR'}/$currentPatient/$currentPatient.post-recal.grp \\\n";
		$cmdPlot .= "-plots $dirs{'BQSR'}/$currentPatient/$currentPatient.recal.plots.pdf \\\n";
		$cmdPlot .= "-log $dirs{'BQSR'}/$currentPatient/$currentPatient.plots.log";

		#BaseRecalibrator 4: print reads
	
		my $cmdPrint = "java -$xmx -jar $gatk -l INFO \\\n";
		$cmdPrint .= "-T PrintReads -K $gatk_key -et NO_ET \\\n";
		$cmdPrint .= "-BQSR $dirs{'BQSR'}/$currentPatient/$currentPatient.recal.grp \\\n";
		$cmdPrint .= "-R $hg19Dir \\\n";
		$cmdPrint .= "-nct $cpuPerThread \\\n";
		$cmdPrint .= "-I $dirs{'realign'}/$currentPatient.bam \\\n";
		$cmdPrint .= "-o $dirs{'BQSR'}/$currentPatient.bam \\\n";
		$cmdPrint .= "-log $dirs{'BQSR'}/$currentPatient/$currentPatient.print.log";
	
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmdCovar\n");
		say_file ($LOG,"$cmdPost\n");
		say_file ($LOG,"$cmdPlot\n");
		say_file ($LOG,"$cmdPrint\n");
	
		#execute command
		my $pid = $pm -> start and next;
			$cmdOnly or system ($cmdCovar);
			$cmdOnly or system ($cmdPost);
			$cmdOnly or system ($cmdPlot);
			$cmdOnly or system ($cmdPrint);
	 	$pm -> finish;
	
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "\n=> Base Recalibration done.\n");
}

#########################################
######## Create diverse matrices ########
#########################################


####create output directory for mapping step
make_dir ("stats", "07-stats", "out");

#skip if wanted
if ($start > 7 or $end < 7) {
	say_out ($LOG, "***** (7) Skipping create matrices");
} else {
	
	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'BQSR'}/*.bam");
	foreach (@files){
	
	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.bam$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if input files exist in $dirs{'BQSR'}.");
	}

	########################
	#### calculate and set maximum number of available threads
	
	#calculate numbers of thread defined by given ram (10GB per job)
	my $maxThreads = int ($totalMem / $ram);
	if ($maxThreads > $cpu) {
		$maxThreads = $cpu;
	}
	
	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);

	# Write version in Log file
	my $version = `java -Xmx4g -jar $gatk --version`;
	chomp $version;

	#write in log file
	my $print = "######## Create Matrices ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\nGATK version: $version\nRunning up to $maxThreads threads in parallel.\n");

	############################################
	######## Building and executing calling command 
		
	foreach my $currentPatient (@patients) {
		
		#generate subdirectories
		mkdir "$dirs{'stats'}/$currentPatient";
		
		#prepare command "Mean Depth"
		my $cmdMeanDepth = "$samtools depth -b $exonFile $dirs{'BQSR'}/$currentPatient.bam ";
		$cmdMeanDepth .= "| awk '{sum+=\$3;cnt++} END {print \"$dirs{'BQSR'}/$currentPatient.bam \\t\"sum/cnt}' ";
		$cmdMeanDepth .= "> $dirs{'stats'}/$currentPatient/$currentPatient.meanDepth.txt";
		
		#prepare command CollectAlignmentSummaryMetrics (Picard)
		my $cmdAsMetrics = "java -$xmx -Djava.io.tmpdir=$tmpDir \\\n";
		$cmdAsMetrics .= "-jar $picardDir/CollectAlignmentSummaryMetrics.jar \\\n";
		$cmdAsMetrics .= "REFERENCE_SEQUENCE=$hg19Dir \\\n";
		$cmdAsMetrics .= "VALIDATION_STRINGENCY=LENIENT \\\n";
		$cmdAsMetrics .= "INPUT=$dirs{'BQSR'}/$currentPatient.bam \\\n";
		$cmdAsMetrics .= "OUTPUT=$dirs{'stats'}/$currentPatient/$currentPatient.ASMetric.txt \\\n";
		$cmdAsMetrics .= "2> $dirs{'stats'}/$currentPatient/$currentPatient.ASMetric.log ";


		#prepare command CollectGcBiasMetrics (Picard)
		my $cmdGc = "java -$xmx -Djava.io.tmpdir=$tmpDir \\\n";
		$cmdGc .= "-jar $picardDir/CollectGcBiasMetrics.jar \\\n";
		$cmdGc .= "REFERENCE_SEQUENCE=$hg19Dir \\\n";
		$cmdGc .= "VALIDATION_STRINGENCY=LENIENT \\\n";
		$cmdGc .= "CHART_OUTPUT=$dirs{'stats'}/$currentPatient/$currentPatient.BCBMetrics.pdf \\\n";
		$cmdGc .= "INPUT=$dirs{'BQSR'}/$currentPatient.bam \\\n";
		$cmdGc .= "OUTPUT=$dirs{'stats'}/$currentPatient/$currentPatient.BCBMetrics.txt \\\n";
		$cmdGc .= "2> $dirs{'stats'}/$currentPatient/$currentPatient.GCBMetrics.log";


		#prepare command CollectInsertSiteMetrics (Picard)
		my $cmdIns = "java -$xmx -Djava.io.tmpdir=$tmpDir \\\n";
		$cmdIns .= "-jar $picardDir/CollectInsertSizeMetrics.jar \\\n";	
		$cmdIns .= "VALIDATION_STRINGENCY=LENIENT \\\n";
		$cmdIns .= "INPUT=$dirs{'BQSR'}/$currentPatient.bam \\\n";
		$cmdIns .= "OUTPUT=$dirs{'stats'}/$currentPatient/$currentPatient.ISMetrics.txt \\\n";
		$cmdIns .= "HISTOGRAM_FILE=$dirs{'stats'}/$currentPatient/$currentPatient.ISMetrics.pdf \\\n";
		$cmdIns .= "2> $dirs{'stats'}/$currentPatient/$currentPatient.ISMetrics.log ";

		#prepare command Percent of Coverage
		
		my $cmdCoverage = "$bedtoolsDir/coverageBed -abam $dirs{'BQSR'}/$currentPatient.bam \\\n";
		$cmdCoverage .= "-b $exonFile -hist \\\n";
		$cmdCoverage .= "> $dirs{'stats'}/$currentPatient/$currentPatient.coverageBed.hist.txt ";

		
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmdMeanDepth\n");
		say_file ($LOG,"$cmdAsMetrics\n");
		say_file ($LOG,"$cmdGc\n");
		say_file ($LOG,"$cmdIns\n");
		say_file ($LOG, "$cmdCoverage\n");
	
		#execute command
		my $pid = $pm -> start and next;
			$cmdOnly or system ($cmdMeanDepth);
			$cmdOnly or system ($cmdAsMetrics);
			$cmdOnly or system ($cmdGc);
			$cmdOnly or system ($cmdIns);
			$cmdOnly or system ($cmdCoverage);

		 	#remove coverageBed.hist.txt
		 	$cmdOnly or unlink ("$dirs{'stats'}/$currentPatient/$currentPatient.coverageBed.hist.txt") or say_out ($LOG, "WARNING: Couldn't remove file: \n$dirs{'stats'}/$currentPatient/$currentPatient.coverageBed.hist.txt");

	 		$pm -> finish;
	 		
	 
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "\n=> Collecting metrices done.\n");
}

########################################################
######## Variant calling using Haplotype caller ########
########################################################

####create output directory for calling step
make_dir ("call", "08-variant-calling", "out");

#skip if wanted
if ($start > 8 or $end < 8) {
	say_out ($LOG, "***** (8) Skipping Variant calling");
} else {

	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'BQSR'}/*.bam");
	foreach (@files){
	
	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.bam$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if input files exist in $dirs{'BQSR'}.");
	}
	
	
	########################
	#### calculate and set maximum number of available threads
	
	#calculate numbers of thread for 4 cpus per thread (Exponential decay per thread = most decay till 4/5) or #cpu/#patients
	my $cpuByPatients = int($cpu / ($#patients+1));
	my $cpuByTesting = 4;
	my $cpuPerThread;
	if ($cpuByPatients > $cpuByTesting) {
		$cpuPerThread = $cpuByPatients;
	} else {
		$cpuPerThread = $cpuByTesting;
	}
	
	#calculate numbers of thread defined by ram (10GB per thread) and cpu (4 per thread)
	my $maxThreads;
	my $threadsByRam = int ($totalMem / $ram);
	my $threadsByCpu = int($cpu / $cpuPerThread);
	if ($threadsByRam > $threadsByCpu) {
		$maxThreads = $threadsByCpu;
	} else {
		$maxThreads = $threadsByRam;
	}
	
	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);


	#########################
	##### Write version in Log file
	my $version = `java -Xmx4g -jar $gatk --version`;
	chomp $version;

	#write in log file
	my $print = "######## Variant calling and filtering ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\nGATK version: $version\nRunning up to $maxThreads threads in parallel.\nUsing $cpuPerThread threads per sequence.");

	############################################
	######## Building and executing calling command 

	foreach my $currentPatient (@patients){
	
		#generate subdirectories
		mkdir "$dirs{'call'}/$currentPatient";
	
		#prepare command for SNP calling
		my $cmdSnp = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
		$cmdSnp .= "-T HaplotypeCaller  \\\n";
		$cmdSnp .= "-R $hg19Dir \\\n";
		$cmdSnp .= "-D $dbsnp \\\n";
		$cmdSnp .= "-L $exonFile \\\n";
		$cmdSnp .= "-ERC GVCF \\\n";
		$cmdSnp .= "-variant_index_type LINEAR \\\n";
		$cmdSnp .= "-variant_index_parameter 128000 \\\n";
		$cmdSnp .= "-nct $cpuPerThread \\\n";
		$cmdSnp .= "-I $dirs{'BQSR'}/$currentPatient.bam \\\n";
		$cmdSnp .= "-o $dirs{'call'}/$currentPatient/$currentPatient-raw.vcf \\\n";
		$cmdSnp .= "-log $dirs{'call'}/$currentPatient/$currentPatient-raw_calling.log";	
	
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmdSnp\n");
	
		#execute command
		my $pid = $pm -> start and next;
		$cmdOnly or system ($cmdSnp);
		$pm -> finish;
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	##########################
	######## Running GenotypeGVCF
	
	#create folder for combined variants
	mkdir "$dirs{'call'}/combined";
	
	#prepare command 	
	my $cmdGenotypeGvcf = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
	$cmdGenotypeGvcf .= "-T GenotypeGVCFs \\\n";
	$cmdGenotypeGvcf .= "-R $hg19Dir \\\n";
	$cmdGenotypeGvcf .= "-D $dbsnp \\\n";
	$cmdGenotypeGvcf .= "-stand_call_conf 30 \\\n";
	$cmdGenotypeGvcf .= "-stand_emit_conf 30 \\\n";
	$cmdGenotypeGvcf .= "-nt $cpu \\\n";
	
	foreach my $currentPatient (@patients) {
		$cmdGenotypeGvcf .= "-V  $dirs{'call'}/$currentPatient/$currentPatient-raw.vcf \\\n";
	}
	
	$cmdGenotypeGvcf .= "-o $dirs{'call'}/combined/combined_raw.vcf \\\n";
	$cmdGenotypeGvcf .= "-log $dirs{'call'}/GenotypeGVCF.log";
	
	#make log entry
	say_out ($LOG, "Running GenotypeGVCF.\n");
	say_file ($LOG, "$cmdGenotypeGvcf\n");
	
	#execute command
	$cmdOnly or system ($cmdGenotypeGvcf);
	
	#log entrs DONE
	say_out ($LOG, "=> GenotypeGVCF done.\n");
	
	#########################
	######## variant filtering	
	
	#prepare command
	my $cmdFilter = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
	$cmdFilter .= "-T VariantFiltration \\\n";
	$cmdFilter .= "-R $hg19Dir \\\n";
	$cmdFilter .= "--filterExpression \"QD < 2.0\" \\\n";
	$cmdFilter .= "--filterName \"QDFilter\" \\\n";
	$cmdFilter .= "--filterExpression \"MQ < 40.0\" \\\n";
	$cmdFilter .= "--filterName \"MQFilter\" \\\n";
	$cmdFilter .= "--filterExpression \"FS > 60.0\" \\\n";
	$cmdFilter .= "--filterName \"FSFilter\" \\\n";
	$cmdFilter .= "--filterExpression \"HaplotypeScore > 13.0\" \\\n";
	$cmdFilter .= "--filterName \"HaplotypeScoreFilter\" \\\n";
	$cmdFilter .= "--filterExpression \"MQRankSum < -12.5\" \\\n";
	$cmdFilter .= "--filterName \"MQRankSumFilter\" \\\n";
	$cmdFilter .= "--filterExpression \"ReadPosRankSum < -8.0\" \\\n";
	$cmdFilter .= "--filterName \"ReadPosRankSumFilter\" \\\n";
	$cmdFilter .= "-V $dirs{'call'}/combined/combined_raw.vcf \\\n";
	$cmdFilter .= "-o $dirs{'call'}/combined/combined_filtered.vcf \\\n";
	$cmdFilter .= "-log $dirs{'call'}/combined/VariantFiltration.log";
	
	
	#make log entry
	say_out ($LOG, "Running VariantFiltration.\n");
	say_file ($LOG, "$cmdFilter\n");
	
	#execute command
	if ($hard){	
		$cmdOnly or system ($cmdFilter);
		
		#log entrs DONE
		say_out ($LOG, "=> VariantFiltration done.\n");
	}
	
		
	###############################
	######## Exclude filtered reads
	
	# get version number and write in LOG
	$version = `$vcftools | grep "VCF"`;
	chomp $version;
	
	#prepare command
	my $cmdExclude = "$vcftools \\\n";
	$cmdExclude .= "--remove-filtered-all \\\n";
	$cmdExclude .= "--recode \\\n";
	$cmdExclude .= "--vcf $dirs{'call'}/combined/combined_filtered.vcf \\\n";
	$cmdExclude .= "--out $dirs{'call'}/combined/combined_filtered_passed ";

	#execute command
	if ($hard){
		#make log entry
		say_out ($LOG, "Running Excluding filtered reads (vcftools).\nVersion: $version.\n");
		say_file ($LOG, "$cmdExclude\n");
		
		
		$cmdOnly or system ($cmdExclude);
		
		#log entry DONE
		say_out ($LOG, "=> Excluding filtered variants done.\n");
	}
	
	
	###################################
	######## VQSR (VariantRecalibrator)

	#create folder for combined variants
	mkdir "$dirs{'call'}/VQSR";
		
	######## prepare commands
	#SNP
	my $cmdRecalSNP = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
	$cmdRecalSNP .= "-T VariantRecalibrator \\\n";
	$cmdRecalSNP .= "-nt $cpu \\\n";
	$cmdRecalSNP .= "-R $hg19Dir \\\n";
	$cmdRecalSNP .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \\\n";
	$cmdRecalSNP .= "-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \\\n";
	$cmdRecalSNP .= "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $snps1kgp \\\n";
	$cmdRecalSNP .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \\\n";
	$cmdRecalSNP .= "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\\n";
	$cmdRecalSNP .= "-mode SNP \\\n";
	$cmdRecalSNP .= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n";
	$cmdRecalSNP .= "-rscriptFile $dirs{'call'}/VQSR/recalibrate_SNP_plots.R \\\n";
	$cmdRecalSNP .= "-recalFile $dirs{'call'}/VQSR/recalibrate_SNP.recal \\\n";
	$cmdRecalSNP .= "-tranchesFile $dirs{'call'}/VQSR/recalibrate_SNP.tranches \\\n";
	$cmdRecalSNP .= "-input $dirs{'call'}/combined/combined_raw.vcf";
	
	
	#### apply recal
	my $cmdApplySnp = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
	$cmdApplySnp .= "-T ApplyRecalibration \\\n";
	$cmdApplySnp .= "-R $hg19Dir \\\n";
	$cmdApplySnp .= "-mode SNP \\\n";
	$cmdApplySnp .= "--ts_filter_level 99.0 \\\n";
	$cmdApplySnp .= "-recalFile $dirs{'call'}/VQSR/recalibrate_SNP.recal \\\n";
	$cmdApplySnp .= "-tranchesFile $dirs{'call'}/VQSR/recalibrate_SNP.tranches \\\n";
	$cmdApplySnp .= "-input $dirs{'call'}/combined/combined_raw.vcf \\\n";
	$cmdApplySnp .= "-o $dirs{'call'}/VQSR/recalibrated_snps_raw_indels.vcf ";
	
	#make log entry
	say_out ($LOG, "Running VQSR on SNPs.\n");
	say_file ($LOG, "VQSR (SNPs)\n\n$cmdRecalSNP \n\nApply (SNPs)\n\n$cmdApplySnp \n");
	
	#execute command
	 if (!$hard) {
		$cmdOnly or system ($cmdRecalSNP);
		$cmdOnly or system ($cmdApplySnp);
	}
	
	#INDEL
	
	my $cmdRecalIndel = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
	$cmdRecalIndel .= "-T VariantRecalibrator \\\n";
	$cmdRecalIndel .= "-nt $cpu \\\n";
	$cmdRecalIndel .= "-R $hg19Dir \\\n";
	$cmdRecalIndel .= "-resource:mills,known=true,training=true,truth=true,prior=12.0 /data/ngs/programs/bundle_2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf \\\n";
	$cmdRecalIndel .= "-an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum \\\n";
	$cmdRecalIndel .= "-mode INDEL \\\n";
	$cmdRecalIndel .= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n";
	$cmdRecalIndel .= "--maxGaussians 4 \\\n";
	$cmdRecalIndel .= "-recalFile $dirs{'call'}/VQSR/recalibrate_INDEL.recal \\\n";
	$cmdRecalIndel .= "-tranchesFile $dirs{'call'}/VQSR/recalibrate_INDEL.tranches \\\n";
	$cmdRecalIndel .= "-rscriptFile $dirs{'call'}/VQSR/recalibrate_INDEL_plots.R \\\n";
	$cmdRecalIndel .= "-input $dirs{'call'}/VQSR/recalibrated_snps_raw_indels.vcf";

	my $cmdApplyIndel = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
	$cmdApplyIndel .= "-T ApplyRecalibration \\\n";
	$cmdApplyIndel .= "-R $hg19Dir \\\n";
	$cmdApplyIndel .= "-mode INDEL \\\n";
	$cmdApplyIndel .= "--ts_filter_level 99.0 \\\n";
	$cmdApplyIndel .= "-recalFile $dirs{'call'}/VQSR/recalibrate_INDEL.recal \\\n";
	$cmdApplyIndel .= "-tranchesFile $dirs{'call'}/VQSR/recalibrate_INDEL.tranches \\\n";
	$cmdApplyIndel .= "-input $dirs{'call'}/VQSR/recalibrated_snps_raw_indels.vcf \\\n";
	$cmdApplyIndel .= "-o $dirs{'call'}/VQSR/combined_recalibrated.vcf ";
	
	
	
	#execute command
	
	if (!$hard){
	
		#make log entry
		say_out ($LOG, "Running VQSR on INDEL.\n");
		say_file ($LOG, "VQSR (SNPs)\n\n$cmdRecalIndel \n\nApply (SNPs)\n\n$cmdApplyIndel \n");
	
		$cmdOnly or system ($cmdRecalIndel);
		$cmdOnly or system ($cmdApplyIndel);
	
		#log entry DONE
		say_out ($LOG, "=> Excluding VQSR done.\n");
	}
	
	

}


##################################################
######## extract individuals and families ########
########        evaluate variants         ########
##################################################


####create output directory for mapping step
make_dir ("extract", "09-Extract-idividuals", "out");

#skip if wanted
if ($start > 9 or $end < 9) {
	say_out ($LOG, "***** (9) Skipping Extract Individuals");
} else {

	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'BQSR'}/*.bam");
	foreach (@files){
	
	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.bam$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if ($extractInd){
		if (!@patients){
			die_command ("No patients available. \nCheck if input files exist in $dirs{'BQSR'}.");
		}
	}
	
	#########################
	##### Write version in Log file
	
	# get version number and write in LOG
	my $version = `$vcftools | grep "VCF"`;
	chomp $version;

	#write in log file
	my $print = "######## Extract individuals and families ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\nVcftools version: \n$version\n");
	
	
	############################
	######## extract individuals

	# if chosen extract individuals to family vcfs
	if ($pedIN) {
#		
 		#prepare cmd
		foreach my $key (keys %family) {
			my $cmdExtractFam = "$vcftools \\\n";
			$cmdExtractFam .= "--recode \\\n";
			foreach my $currentPatient (@{$family{$key}}) {
				$cmdExtractFam .= "--indv $currentPatient \\\n";
			say $currentPatient;				
			}
			if ($hard){
				$cmdExtractFam .= "--vcf $dirs{'call'}/combined/combined_filtered_passed.recode.vcf \\\n";
			} else {
				$cmdExtractFam .= "--vcf $dirs{'call'}/VQSR/combined_recalibrated.vcf \\\n";			
			}
			$cmdExtractFam .= "--out $dirs{'extract'}/$key";

			#make log entry
			say_out ($LOG, "Running on $key\n");
			say_file ($LOG, "$cmdExtractFam\n");

			#execute command
			$cmdOnly or system ($cmdExtractFam);
			$cmdOnly or system ("mv $dirs{'extract'}/$key.recode.vcf $dirs{'extract'}/$key.vcf")
		}
	}


	#### extract all individuals from vcf if not other chosen		
	if ($extractInd){
		#make log entry
		say_out ($LOG, "Running extract individuals (vcftools).\nVersion: $version.\n");
		
		foreach my $currentPatient (@patients){
			# prepare command
			my $cmdExtract = "$vcftools \\\n";
			$cmdExtract .= "--recode \\\n";
			$cmdExtract .= "--indv $currentPatient \\\n";
			if ($hard) {
				$cmdExtract .= "--vcf $dirs{'call'}/combined/combined_filtered_passed.recode.vcf \\\n";
			} else {
				$cmdExtract .= "--vcf $dirs{'call'}/VQSR/combined_recalibrated.vcf \\\n";
			}
			$cmdExtract .= "--out $dirs{'extract'}/$currentPatient";
		
			#make log entry
			say_out ($LOG, "Running on $currentPatient\n");
			say_file ($LOG, "$cmdExtract\n");
		
			#execute command
			$cmdOnly or system ($cmdExtract);
			$cmdOnly or system ("mv $dirs{'extract'}/$currentPatient.recode.vcf $dirs{'extract'}/$currentPatient.vcf");
		}

		#log entries DONE
		say_out ($LOG, "=> Extracting individuals done.\n");
	}
	
	
	
	#######################
	#### variant evaluation



	########################
	#### calculate and set maximum number of available threads
	
	#calculate numbers of thread for 4 cpus per thread (Exponential decay per thread = most decay till 4/5) or #cpu/#patients
	my $cpuByPatients = int($cpu / ($#patients+1));
	my $cpuByTesting = 4;
	my $cpuPerThread;
	if ($cpuByPatients > $cpuByTesting) {
		$cpuPerThread = $cpuByPatients;
	} else {
		$cpuPerThread = $cpuByTesting;
	}
	
	#calculate numbers of thread defined by ram (10GB per thread) and cpu (4 per thread)
	my $maxThreads;
	my $threadsByRam = int ($totalMem / $ram);
	my $threadsByCpu = int($cpu / $cpuPerThread);
	if ($threadsByRam > $threadsByCpu) {
		$maxThreads = $threadsByCpu;
	} else {
		$maxThreads = $threadsByRam;
	}
	
	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);



	#get all vcf files and loop variant eval over them
	@files = glob ("$dirs{'extract'}/*.vcf");
	
	#create file for eval results
	my $outDir = "$dirs{'extract'}/stats";
	if  (!-d "$outDir"){
		mkdir "$outDir";
	}
	
	foreach (@files) {
		
		#create name for out file
		my $currentPatient = basename($_);
		$currentPatient =~ s/.vcf$//;

		#create cmd for variant eval
		my $cmdEval = "java -$xmx -jar $gatk -l INFO -K $gatk_key -et NO_ET \\\n";
		$cmdEval .= "-R $hg19Dir \\\n";
		$cmdEval .= "-D $dbsnp \\\n";
		$cmdEval .= "-T VariantEval \\\n";
		$cmdEval .= "-nt $cpuPerThread \\\n";
		$cmdEval .= "--eval:set1 $_ \\\n";
		$cmdEval .= "-o $outDir/$currentPatient.eval.xls \\\n";
		$cmdEval .= "-ST Sample -noST -noEV -EV CountVariants \\\n";
		$cmdEval .= "-EV TiTvVariantEvaluator -EV CompOverlap ";
		
		
		
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmdEval\n");
	
		#execute command
		my $pid = $pm -> start and next;
		$cmdOnly or system ($cmdEval);
		$pm -> finish;
	
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	say_out ($LOG, "\nVariant Evaluation finished!\n");		
}

################################
######## VEP Annotation ########
################################

####create output directory for mapping step
make_dir ("VEP", "10-VEP", "out");

#skip if wanted
if ($start > 10 or $end < 10) {
	say_out ($LOG, "***** (10) Skipping VEP annotation");
} else {

	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'extract'}/*.vcf");
	foreach (@files){
	
	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.vcf$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if input files exist in $dirs{'extract'}.");
	}
	
	########################
	#### calculate and set maximum number of available threads
	
	#calculate numbers of thread for 4 cpus per thread (Exponential decay per thread = most decay till 4/5) or #cpu/#patients
	my $cpuByPatients = int($cpu / ($#patients+1));
	my $cpuByTesting = 6;
	my $cpuPerThread;
	if ($cpuByPatients > $cpuByTesting) {
		$cpuPerThread = $cpuByPatients;
	} else {
		$cpuPerThread = $cpuByTesting;
	}
	my $maxThreads = int($cpu / $cpuPerThread);
	
	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);



	#########################
	##### Write version in Log file
	my $version = `$VEP -help | grep version`;
	chomp $version;
	#write in log file
	my $print = "######## VEP annotation ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\nVEP version: $version\nRunning up to $maxThreads threads in parallel.\nUsing $cpuPerThread threads per sequence.");


	#########################
	#### VEP annotation & preparation

	foreach my $currentPatient (@patients){
		
		#prepare cmd sed & vt -> needed prior to VEP; demanded by Gemini (0.13 and later)
		my $cmdPrep = "zless $dirs{'extract'}/$currentPatient.vcf \\\n";
		$cmdPrep .="| sed 's/ID=AD,Number=./ID=AD,Number=R/' \\\n";
		$cmdPrep .="| $vtMaster decompose -s - \\\n";
		$cmdPrep .= "| $vtMaster normalize -r $hg19Dir -o $dirs{'VEP'}/$currentPatient-sed-norm-decomp.vcf -";
		
		#prepare cmd VEP annotation
		my $cmdVEP = "$VEP \\\n";
		$cmdVEP .= "-i $dirs{'VEP'}/$currentPatient-sed-norm-decomp.vcf \\\n";
		$cmdVEP .= "-o $dirs{'VEP'}/$currentPatient-VEP.vcf \\\n";
		$cmdVEP .= "-fork $cpuByTesting \\\n";
		$cmdVEP .= "--cache \\\n";
		$cmdVEP .= "--species homo_sapiens \ \\\n";
		$cmdVEP .= "--cache_version 77 \\\n";
		$cmdVEP .= "--assembly GRCh37 \\\n";
		$cmdVEP .= "--force_overwrite \\\n";
		$cmdVEP .= "--sift b \\\n";
		$cmdVEP .= "--polyphen b \\\n";
		$cmdVEP .= "--symbol \\\n";
		$cmdVEP .= "--numbers \\\n";
		$cmdVEP .= "--biotype \\\n";
		$cmdVEP .= "--total_length \\\n";
		$cmdVEP .= "--vcf \\\n";
		$cmdVEP .= "--offline \\\n";
		$cmdVEP .= "--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \n";
		
		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmdPrep\n");
		say_file ($LOG, "$cmdVEP\n");
	
		#execute command
		my $pid = $pm -> start and next;
			$cmdOnly or system ($cmdPrep);
			$cmdOnly or system ($cmdVEP);
			#remove temporary vcf files from sed norm decomp
			say_file ($LOG, "Removing $dirs{'VEP'}/$currentPatient-sed-norm-decom.vcf.");
			unlink ("$dirs{'VEP'}/$currentPatient-sed-norm-decomp.vcf");
		$pm -> finish;
	
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "=> VEP-annotation done.\n");
}



################################
######## Load in Gemini ########
################################

####create output directory for mapping step
make_dir ("databases", "11-databases", "out");

#skip if wanted
if ($start > 11 or $end < 11) {
	say_out ($LOG, "***** (11) Skipping loading in gemini");
} else {

	###############
	####extract number and name of patients
	my @patients;
	
	#extract bam files
	my @files = glob ("$dirs{'VEP'}/*.vcf");
	foreach (@files){
	
	#extract current patient
	my $currentPatient = basename($_);
	$currentPatient =~ s/.vcf$//;

		#save unique patientnames for later use
		my $controll = 0;
		foreach (@patients){
			if ($currentPatient eq $_){
				$controll = 1;
			}
		}

		#if patient is unique save for later use
		if ($controll == 0){
			push (@patients, $currentPatient);
		}
	}
	
	#if still no patients available stop processing
	if (!@patients){
		die_command ("No patients available. \nCheck if input files exist in $dirs{'VEP'}.");
	}
	
	########################
	#### calculate and set maximum number of available threads
	
	#calculate numbers of thread for 6 cpus per thread (Exponential decay per thread = most decay till 6) or #cpu/#patients
	my $cpuByPatients = int($cpu / ($#patients+1));
	my $cpuByTesting = 6;
	my $cpuPerThread;
	if ($cpuByPatients > $cpuByTesting) {
		$cpuPerThread = $cpuByPatients;
	} else {
		$cpuPerThread = $cpuByTesting;
	}
	my $maxThreads = int($cpu / $cpuPerThread);
	
	#set maximum number of threads
	my $pm = Parallel::ForkManager->new($maxThreads);


	#########################
	##### Write version in Log file
	my $version = `$gemini --version 2>&1`;
	chomp $version;
	
	#write in log file
	my $print = "######## Load in Gemini databases ########";
	my $string = ();
	for (my $i=0; $i<length $print; $i++){
		$string .= "#";
	}
	say_out ($LOG, "\n$string\n$print\n$string\nGemini version: $version\nRunning up to $maxThreads threads in parallel.\nUsing $cpuPerThread threads per sequence.");

	#########################
	#### load in gemini database

	foreach my $currentPatient (@patients) {

		my $currentOut = $currentPatient;
		$currentOut =~ s/-VEP$//;
		say $currentOut;

		my $cmdLoad = "$gemini load \\\n";
		$cmdLoad .= "-t VEP \\\n";
		$cmdLoad .= "--cores $cpuPerThread \\\n";
		$cmdLoad .= "-v $dirs{'VEP'}/$currentPatient.vcf \\\n";
		$cmdLoad .= "-p $pedIN \\\n";
		$cmdLoad .= "$dirs{'databases'}/$currentPatient.db";

		#make log entry
		say_out ($LOG, "\nRunning on $currentPatient\n");
		say_file ($LOG, "$cmdLoad\n");
	
		#execute command
		my $pid = $pm -> start and next;
		$cmdOnly or system ($cmdLoad);
		$pm -> finish;
	
	}

	#wait for all children to finish
	$pm -> wait_all_children;

	#log entrs DONE
	say_out ($LOG, "=> Load in Gemini database done done.\n");

}

#### Add ending date and time
$date = localtime ();
say_out ($LOG, "Pipeline finished: $date.");








#############################
######## Subroutines ########
#############################



################
#### die command

sub die_command {
	my $message = $_[0];
	$message .= "\n";
	
	#give out message and write in LOG-file
	say_out ($LOG, "\nFAILURE:");
	say_out ($LOG, $message);
	say_out ($LOG, "#### Execution halted! ####");
	die "\n";
}

############################
#### print on screen and log

sub say_out {
	my $file = shift @_;
	say "@_";
	say $file "@_";
}


##################
#### print in file

sub say_file {
	my $file = shift @_;
	say $file "@_";
}


###########################
####create output directory

sub make_dir {
	my $key = shift;
	my $newPath = shift;
	my $mainPath = shift;
	
	#check if dir exist, if not create
	$dirs{$key} = "$dirs{$mainPath}/$newPath";
	if  (!-d "$dirs{$key}"){
		mkdir "$dirs{$key}";
	}
}

