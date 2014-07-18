#! /usr/bin/perl

#######################################################
# CFIA de novo microbial genome assembly pipeline
#
# Authored by Adam Koziol and Michael Knowles
#
# NOTE: As of now, this script will need certain programs installed in order to run properly
#	1) Velvet and velvet optimizer. Velvet will have to be compiled manually in order to enable multithreading and kmer sizes up to 115
#	#2) prokka and all of its dependencies
#	#3) ffp
#	#4) abacas
#	#5) BLAST
#	6) Trimmomatic
#
# We're also assuming that you're running this in a Linux environment. If not, then we have no idea if it will work.

####################################################### 

# Much of this code is influenced by Torsten Seemann's prokka 1.5.2
# Also, thanks are given to Nick Petronella at Health Canada for code examples and guidance

use Cwd;
#use warnings;
use strict;
use File::Copy;
use File::Path qw(make_path remove_tree);
use Time::Piece;
use Time::Seconds;
use Time::HiRes qw/ time sleep /;

# This start time will be used in calculating the total time of the run
my $start_time = time;
# This script assumes that the raw sequence files are Illumina outputted paired end, .fastq, and still zipped. Additional functionality to include different naming conventions and formats may eventually be included
running_time("Welcome to the CFIA microbial genome assembly pipeline.");

# Determine the number of processors present
my @cpus = CPUs();
running_time("There are @cpus CPU(s) in your system.");

# Runs the folderer subroutine to extract the .fastq files and place them in appropriate directories
# Checks to see if there are .gz archives present in the current directory
my $path = getcwd;
chdir ("$path/sequences");
my $folderer_check = glob "*.gz";

# If there are .gz archives in the sequences folder, the folderer subroutine is called
if ($folderer_check){
	running_time("Your paired-end .fastq files will now be extracted.");
	folderer();
}

=com
# Check to see if there are sequence-containing folders present in the $path - if they are full, then the program knows that folderer has already been completed
# Due to the nature of the subroutine, the check for trimmed sequences can already be stored as a variable
my ($fastq_check, $trimmed_check, $folders) = fastq_check();
my @fastq_check = @$fastq_check;
my @folders = @$folders;
if (@fastq_check) {
	running_time("Fastq files are extracted and in folders.")
}

# Because the fastq_check subroutine returns array references, the scalar $trimmed_check must be converted to an array
# Run checking subroutine to determine if trimming has already been perfomed
my @trimmed_check = @$trimmed_check;

# Run trimmomatic to process the .fastq files to remove adapter contamination and low quality sequence
if (@trimmed_check) {
	running_time("Quality trimming your sequences.");
	trimmomatic(@trimmed_check);
} else {
	running_time("Trimming is complete.");
}

# Perform checks to see if velvetoptimizer has previously run and has either assembled the contigs completely, or was interrupted. In the case of an interruption
# folders with partial data are deleted, and velvetoptimizer is run again.  
velvet_run(@fastq_check);

# Perform checks on quality metrics to determine whether untrimmed or trimmed assemblies are best. Creates a summary spreadsheet
assembly_report(\@fastq_check, \@folders);

=cut

# Return the run time
my $end_time = time;
my $total_time = 1000 * ($end_time - $start_time);
my $milliseconds = sprintf("%.0f", $total_time);
print "The total run time was $milliseconds milliseconds.\n";


##########################################################
# This subroutine allows the printing of the local time each time it is invoked - it allows for the user to see how much time has passed
sub running_time 
{
my $time =  localtime;
my $hms = "[" . $time->hms . "] @_\n";
print $hms;
}

##########################################################
# This subroutine counts the number of CPUs present
sub CPUs
{
	my @cpus = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
	@cpus = 1;
#	system("awk '/^processor/ { N++} END { print N }' /proc/cpuinfo > cpus.temp");
#	open(INPUT, "<", "cpus.temp");
#	my @cpus = <INPUT>;
#	chomp @cpus;
#	unlink "cpus.temp";
	return @cpus;
}


##########################################################
# This subroutine extracts .fastq files, places the decompressed files in folders, and deletes the .gz archives
sub folderer
{
my @files = glob("*.gz");
foreach my $file (@files){
	(my $folder = $file) =~ s/_.*//g;
	(my $unzip = $file) =~ s/.gz//g;
	(my $filename = $file) =~ s/S.+_L001_|.gz//g;
	running_time("Now extracting $filename");
	mkdir $folder;
	move($file,"$folder/$file");
	system("gzip -d $folder/$file");
	if (-f $unzip){unlink "$folder/$file"};
	move("$folder/$unzip","$folder/$filename"); 
	}
}

=com
##########################################################
# This subroutine checks to see if folders in $path have .fastq files
sub fastq_check
{
	# Get the path for file/folder manipulations
	my $path = getcwd;
	# Initialize
	my ($mpath, @fastq_1, @fastq_2, @fastq_folders, @trimming_incomplete, @folders) = "";
	
	# Open all folders in $path 
	opendir(DIR, $path) or die "can't find $path: $!";
	while (defined(my $file = readdir(DIR))) {
		# Ignore special files
		next if $file =~ /^\.\.?$/;
		# Ignore folders without sequence data
		next if $file =~ /Best_Assemblies/;
		if (-d $path."/".$file){
			$mpath = $path;
			# Clear the arrays each time through the loop
			(@fastq_1, @fastq_2) = "";
			# The space was throwing off velvet
			$mpath =~ s/ /\\ /g; 
			chdir $path."/".$file;
			# Grab any untrimmed fastq files and put into appropriate arrays
			@fastq_1 = glob("*R1_001.fastq");
			@fastq_2 = glob("*R2_001.fastq");
			# If both arrays have values, then add the $path and filename to an array
			if (@fastq_1 and @fastq_2) {
				push(@fastq_folders, $path . "/" . $file);
				push(@folders, $file);	
			}
			if (@fastq_1 and @fastq_2 and not -e("." . $file)) {
				push(@trimming_incomplete, $path . "/" . $file);
			}
		}
	}
	chdir $path;
	return (\@fastq_folders, \@trimming_incomplete, \@folders);
}

##########################################################
# This subroutine checks to see if folders in $path have trimmed .fastq files
sub trimmomatic
{
	my @files = @_;
	my $path = getcwd;
	my (@fastq_1, @fastq_2, $file_name) = '';
	foreach (@files) {
		chdir ("$_");
		@fastq_1 = glob("*_R1_001.fastq");
		($file_name = $fastq_1[0]) =~ s/_R1_001.fastq//g;
		@fastq_2 = glob("*_R2_001.fastq");
		system ("trimmomatic-0.30.jar PE -threads 24 -phred33 -trimlog $file_name.log $fastq_1[0] $fastq_2[0] $file_name.paired1.fastq $file_name.unpaired1.fastq $file_name.paired2.fastq $file_name.unpaired2.fastq ILLUMINACLIP:/home/blais/Bioinformatics/Trimmomatic-0.30/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");	
		# Creates a hidden file once the script has finished. This will ensure that if the script is prematurely terminated, that partial files will not be considered "acceptable", and the folder is passed over
		open(OUTPUT, ">", ".$file_name");
		close OUTPUT;
		chdir("$path");
	} 	
}

##########################################################
sub velvet_check 
{
	my($file_locations, $trim_level) = @_;
	my @files = @$file_locations;
	# Get the path for file/folder manipulations
	my $path = getcwd;
	
	# Initialize
	my (@untrimmed_unassembled, @trimmed_unassembled, @auto_data, @partial_auto_data, @log, $tail, $check, $junk) = '';
	
	# Run through the loop
	foreach my $file (@files) {
		# Errors keep popping up due to empty elements in the arrays. By checking that each $file has a length, only elements containing data will be used
		
		if (length $file) {
			chdir ($file);
			# This will remove the Log_file.txt that clog up the main directory after several failed runs
			my @log_file = glob("*_Logfile.txt");
				foreach (@log_file) {
					unlink $_;
				}
			# If there's any auto_data folders that have not been renamed or deleted, they are removed, and the location of this incompletely
			# velvetoptimizer-processed folder are added to the @unassembled_auto_data array
			@partial_auto_data = glob("auto_data*");
			if (@partial_auto_data) {
				if ($trim_level eq "untrimmed") {
					push(@untrimmed_unassembled, $file);
				} elsif ($trim_level eq "trimmed") {
					push(@trimmed_unassembled, $file);
				}
				foreach (@partial_auto_data) {
					chdir($file);
					remove_tree($file . "/" . $_) or warn "Could not delete $file/$_\n";
				}
			}
			
			# Checks to see if velvetoptimizer has run correctly on the untrimmed sequences- there is a check for the presence of the renamed auto_data folder
			# If the untrimmed_auto_data folder isn't present, or if there are multiple copies, then the program assumes there was a problem and removes the
			# folders and velvetoptimizer will run
			@auto_data = glob("$trim_level" . "_auto_data*");
			if (((not @auto_data) or (length @auto_data > 1))) {
				if ($trim_level eq "untrimmed") {
					push(@untrimmed_unassembled, $file);
				} elsif ($trim_level eq "trimmed") {
					push(@trimmed_unassembled, $file);
				}
				if (@auto_data) {
					foreach (@auto_data) {
						chdir ($file);
						remove_tree($file . "/" . $_) or warn "Could not delete $auto_data[0]\n";
					}
				}
			}
			# This checks for the presence of "Log" in the auto_data folders. If the Log file doesn't have the proper format, which can happen if the run
			# was interrupted, then the auto_data folders are deleted and velvetoptimizer will run
			if (length @auto_data == 1) {
				chdir($file . "/" . $auto_data[0]);
				if (-e("Log")) {
					@log = '';
					open(INPUT, "<", "Log");
					@log = (<INPUT>);
					close INPUT;
					my $tail = pop @log;
					if ($tail) {
						($check, $junk) = split(" ", $tail);
					}
					# Log files from prematurely terminated runs seem to have empty lines at the bottom. If the last element of an array created from
					# the log is examined, and it doesn't have a length, then the program assumes that there was a problem with the run, and adds
					# the folder details to the appropriate array
					if (length $check and ($check eq "Final")) {
					} elsif (not (length $check) or not ($check eq "Final")) {
						if ($trim_level eq "untrimmed") {
							push(@untrimmed_unassembled, $file);
							chdir($file);
							remove_tree($file . "/" . $auto_data[0]) or warn "Could not delete $auto_data[0]\n";
						} elsif ($trim_level eq "trimmed") {					
							push(@trimmed_unassembled, $file);
							chdir($file);
							remove_tree($file . "/" . $auto_data[0]) or warn "Could not delete $auto_data[0]\n";
						}	
					}
				}
			}
		}
	}
	chdir($path);
	if ($trim_level eq "untrimmed") {
		return (\@untrimmed_unassembled);
	} elsif ($trim_level eq "trimmed") {		
		return (\@trimmed_unassembled);
	}
}	

##########################################################
sub velvet 
{
	my ($file, $trim_level, $cpus) = @_;
	my @cpus = @$cpus;
	my $path = getcwd;
	my (@fastq_1, @fastq_2, @auto_data, $rename, $check, $junk) = '';
	# Confirm that there are values in $file
	if (length $file) {
		chdir($file);
		# depending on $trim_level, the required files are different
		if ($trim_level eq "untrimmed") {
			@fastq_1 = glob("*R1_001.fastq");
			@fastq_2 = glob("*R2_001.fastq");
		} elsif ($trim_level eq "trimmed") {
			@fastq_1 = glob("*.paired1.fastq");
			@fastq_2 = glob("*.paired2.fastq");
		} else {
			running_time("Unknown trimming requirements");
		}
		# System call to velvetoptimizer - currently using @cpus to determine how many cores to use
		#running_time("Performing velvet assembly of $file sequences.");
		system("VelvetOptimiser.pl -s 45 -e 111 -f '-shortPaired -fastq $file/$fastq_1[0] $file/$fastq_2[0]' -t @cpus -c 'n50*Lcon/tbp+log(Lbp)'");
		
		# I struggled to find a way to incorporate the velvet_check subroutine here to keep from reusing code, but eventually, I decided to copy
		# the code over here and change it around a bit in order to get all the checks working properly	
		@auto_data = glob ("auto_data_*");
		if (length @auto_data == 1) {
			my @log_file = glob("*_Logfile.txt");
			my $auto_data_folder = "$file" . "/" . "$auto_data[0]";
			my $file_name = "Logfile.txt";
			#print "Logfile $log_file[0]\n";
			print "Path of Logfile: $file/$log_file[0]\n";
			print "Path of moved Logfile: $auto_data_folder/$file_name\n";
			move("$file/$log_file[0]", "$auto_data_folder/$file_name") or warn "$!\n";
			chdir($file . "/" . $auto_data[0]);
			# A number of the files present in the folder are not required for our purposes, so they will be removed in order to save on file space
			if (-e("Graph")) {unlink "Graph";}
			if (-e("Graph2")) {unlink "Graph2";}
			if (-e("PreGraph")) {unlink "PreGraph";}
			if (-e("Sequences")) {unlink "Sequences";}
			if (-e("stats.txt")) {unlink "stats.txt";}
			if (-e("Log")) {
				my @log = '';
				open(INPUT, "<", "Log");
				@log = (<INPUT>);
				close INPUT;
				my $tail = pop @log;
				if ($tail) {
					($check, $junk) = split(" ", $tail);
				}
				# Log files from prematurely terminated runs seem to have empty lines at the bottom. If the last element of an array created from
				# the log is examined, and it doesn't have a length, then the program assumes that there was a problem with the run, and adds
				# the folder details to the appropriate array
				if (length $check and ($check eq "Final")) {
					chdir($file);
					($rename = $auto_data[0]) =~ s/auto_data/$trim_level\_auto_data/g;
					move("$file/$auto_data[0]", "$file/$rename");
				} else {
					chdir($file);
					remove_tree($file . "/" . $auto_data[0]) or warn "Could not delete $auto_data[0]\n";
				}
			} else {
				if (@auto_data) {
					chdir($file);
					foreach (@auto_data) {
						remove_tree($file . "/" . $_) or warn "Could not delete $_\n";
					}
				}
			}	
		}
	}	
chdir($path);
}

##########################################################
sub velvet_run 
{
	# Run velvet optimizer to process the trimmed and untrimmed files
	my @fastq_check = @_;
	#my @fastq_check = @$fastq_check;
	my @check = ("untrimmed", "trimmed");
	my (@untrimmed_velvet_check, $untrimmed_velvet_check, @trimmed_velvet_check, $trimmed_velvet_check) = '';
	# For both trimmed and untrimmed sequences
	foreach my $trim (@check) {
		# The velvet_check array is populated using the checks in velvet_check. @fastq_check check is used as the input as it should have all the folders being
		# processed in these analyses - so I don't need to create a new array
		# I found it the best way to deal with the two separate trim levels was to split everything up based on trim - that way each trim would have its own
		# variables. This probably could have been avoided, but because this finally works after a few days of searching for bugs, there is absolutely no
		# way I am changing this now. 
		if ($trim eq "untrimmed") {
			$untrimmed_velvet_check= velvet_check(\@fastq_check, $trim);
			@untrimmed_velvet_check = @$untrimmed_velvet_check;
		} elsif ($trim eq "trimmed") {
			$trimmed_velvet_check= velvet_check(\@fastq_check, $trim);
			@trimmed_velvet_check = @$trimmed_velvet_check;
		}
		# There were a few annoying issues I had a hard time resolving here - @(un)trimmed_velvet_check always seemed to have an undefined value to matter whether there were
		# any real values or not. Normally, the if (@(un)trimmed_velvet_check) loop would resolve this, but I found the length $_was necessary.
		if ($trim eq "untrimmed") {
			# This loop will call the velvet subroutine only if array is populated, otherwise, it assumes that the assembly is already complete
			if (@untrimmed_velvet_check) {
				foreach (@untrimmed_velvet_check) {
					if (length $_) {
						running_time("Performing velvet assembly of $_ $trim sequences.");
						velvet($_, $trim, \@cpus);
					} else {
						running_time("$trim sequences assembled.");		
					}
				}	
			} else {
				running_time("$trim sequences assembled.");
			}
		} elsif ($trim eq "trimmed") {
			if (@trimmed_velvet_check) {
				foreach (@trimmed_velvet_check) {
					if (length $_) {
						running_time("Performing velvet assembly of $_ $trim sequences.");
						velvet($_, $trim, \@cpus);
					} else {
					running_time("$trim sequences assembled.");		
					}
				} 
			} else {
				running_time("$trim sequences assembled.");
			}
		}
	}	
}

##########################################################
sub assembly_report 
{
	# initialize variables
	my ($fastq_check, $folders) = @_;
	my @fastq_check = @$fastq_check;
	my @folders = @$folders;
	my (@auto_data_folder, @untrimmed_auto_data, @trimmed_auto_data, @best_assemblies, @best_trim, $best_trim) = '';
	my ($assembly_score, $best_assembly_score, $best_assembly, $assembly) = 0;
	my (@number_of_contigs, @n50, @longest_contigs, @total_bases);
	my @check = ("untrimmed", "trimmed");
	# get the path for file/folder manipulations
	my $path = getcwd;
	my $count = 0;
	
	# Open the assembly report spreadsheet and write the headings
	open (OUTPUT, ">Assembly_report.csv");
	print OUTPUT "Strain\tTrim\tNumber of Contigs\tLongest Contig\tN50\tTotal Bases\n";
	
	# Because I want to compare the assembly scores for trim levels of the same strain, the original loop is the strains
	# rather than the $trim as per usual 
	foreach my $folder (@folders) {
		# The best assembly score needs to be reset inbetween the different strains
		$best_assembly_score = 0;
		foreach my $trim (@check) {
		chdir ($path . "/" . $folder);
		my @auto_data_folder = glob("$trim\_auto_data*");
		chdir ($auto_data_folder[0]);
		
		# This log file has all the necessary metrics for the assembly report, including the Assembly score.
		# The assembly score is calculated as n50*Lcon/tbp+log(Lbp) - The n50 times the number of long contigs (contigs >1000 bp) divided
		# by the total bases in all contigs plus the log of the number of bases in long contigs. Since this number was used to determine
		# the optimal coverage cutoff in velvetoptimizer, it should be sufficient to determine which assembly is better
		open(LOG, "<", "Logfile.txt");
			while (<LOG>) {
				# Since velvetoptimizer iterates through a few times with different parameters, there are multiple instances
				# of N50, total bases, etc. I want only the text after "Final optimised assembly details:" in the file
				if (/details:/) { 
					# This second while (<LOG>) allows me to get only the text after the above regex is satisfied 
					while (<LOG>) {
						if ($_ =~ m/Assembly\sscore:\s(\d+\W\d+)/) {
							# I'm not sure if this statement is necessary - it was in the example code I used to write this block
							last if (/^$/);
							# The $1 variable is the pattern match from the above regex BUT only between the parentheses.
							# the Assembly\sscore:\s is ignored, while (\d+\W\d+) is set to $1
							$assembly_score = $1;
							chomp $assembly_score;
							}
						}
					}
				} 
			close LOG;
			# Checks to see if $assembly_score is larger than $best_assembly_score for that strain, if so, then $best_assembly_score
			# is then set to $assembly_score. The trim level is also recorded as $best_trim			
			if ($assembly_score > $best_assembly_score) {
				$best_assembly_score = $assembly_score;
				# This is set up in such as way that iterating through this array will allow to chdir into the appropriate
				# folders very easily
				$best_assembly =  ("$path" . "/" . "$folder" . "/" . "$auto_data_folder[0]");
				$best_trim = $trim
			}
		
		} 
		# Arrays of the best assemblies and best trim levels are created, such that each array is filled in the same order
		# eg strain1 best assembly and best trim in $best_assemblies[0] and $best_trim[0], respectively
		push (@best_assemblies, $best_assembly);
		push (@best_trim, $best_trim);
	} 
	
	running_time("Best assemblies are: ");
	
	# Each best assembly Log file is open and the pertinent data are extracted using similar methods as above
	foreach $assembly (@best_assemblies) {
		chdir ($assembly);
		open(LOG, "<", "Logfile.txt");
			while (<LOG>) {
				if (/details:/) { 
					while (<LOG>) {
						if ($_ =~ m/Total\snumber\sof\scontigs:\s(\d+)/) {
							chomp $1;
							push (@number_of_contigs, $1);
						} elsif ($_ =~ m/n50:\s(\d+)/) {
							chomp $1;
							push (@n50, $1);
						} elsif ($_ =~ m/length\sof\slongest\scontig:\s(\d+)/) {
							chomp $1;
							push (@longest_contigs, $1);
						} elsif ($_ =~ m/Total\sbases\sin\scontigs:\s(\d+)/) {
							chomp $1;
							push (@total_bases, $1);
						}	
					}
				}
			}
		# The best assemblies are moved to the Best_Assemblies folder (created here) in order to facilitate
		# future analyses and data manipulations
		system ("mkdir -p $path/Best_Assemblies");
		copy ("$assembly/contigs.fa", "$path/Best_Assemblies/$folders[$count].fa");
		
		# The assembly report is filled out
		print OUTPUT "$folders[$count]\t$best_trim[$count]\t$number_of_contigs[$count]\t$longest_contigs[$count]\t$n50[$count]\t$total_bases[$count]\n";
		
		# Summary is printed to the terminal
		print "\t   $best_trim[$count] $folders[$count] with $number_of_contigs[$count] contigs, an N50 of $n50[$count], and an assembly size of $total_bases[$count] bases\n";
		$count++;
	} 
	close OUTPUT;
	chdir ("$path");
}			

##########################################################
exit; 