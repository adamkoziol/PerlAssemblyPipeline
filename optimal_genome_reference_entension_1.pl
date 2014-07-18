#! usr/bin/perl

#use warnings;
use strict;
use Cwd;
use Spreadsheet::ParseExcel;Untitled event
use Spreadsheet::WriteExcel;
use File::Copy;


# initialize variables
my (@blast_values, @unique_names, @cumulative_bit_score, $best_reference, $best_bit_score, $sorted, $subject, $file, @subject_id, $subject_id, @bit_score, $bit_score, @sorted_values, @sorted_names, @input) = '';
my %seen = ();
my ($count, $cumulative_bit_score) = 0;
my $path = getcwd;

# removing the .csv file at this point was necessary, as this script opens the .csv file to append (>>) later on. This might be avoided if the file was created in this spot instead - For instance, I could create the file, and write the headers, which should be "Strain", "Best Reference", "Cumulative Bit Score"
unlink ("best_references.csv");

# this assumes that the query files have the ".fa"
my @genome_files = glob ("*.fa");

foreach $file (@genome_files) {
	# resetting all variables to "defined, but empty". Shouldn't be necessary - test without this next time this script is run
	(@input, @blast_values, @sorted_values, $subject_id, $bit_score, @subject_id, @bit_score, @unique_names, @sorted_names, $sorted, $count, $cumulative_bit_score, $best_bit_score, $best_reference) = '';
	# re-initializing the %seen hash. This shouldn't be necessary -  see above
	%seen = ();

	# system call to run blastn for each individual query file agains the database of reference genomes. A few things are necessary for this to work: 1) reference genomes - these can be obtained using the entrez.pl script, 2) a blast database of the reference genomes (e.g. formatdb ...)
	# note that this is using a custom output based on format 6 - this only outputs the gi of the reference genomes (I think) and the cumulative bitscore
	system ("blastn -query $file -db C_jejuni_reference_genomes.fas -outfmt '6 sseqid bitscore' -out $file.csv");

	open (INPUT, "<$file.csv");
	@input = (<INPUT>);
	close INPUT;
	foreach (@input) {
		push (@blast_values, (split(/\n/, $_))); # this line takes the @input array and splits it on the newline character. I can't remember exactly why I had to add this, but I think it had something to do with the format of the data
	}

	# this sorts the data alphabetically based on the sseqid of the reference genome
	@sorted_values = sort(@blast_values);

	foreach (@blast_values) {
		($subject_id, $bit_score) = split(/\t/, $_); # splits the data at the tab between sseqid and bitscore
		push(@subject_id, $subject_id);
		push(@bit_score, $bit_score);
	}

	foreach (@bit_score) {
		$_ =~ s/ //g; # removing spaces in the bitscore
	}

	# as there are many hits returned from the blast analysis, this command cycles through the %seen hash to eliminate duplicate names from the array
	@unique_names = grep { ! $seen{$_}++} @subject_id;

	# the unique names are then sorted alphabetically
	@sorted_names = sort(@unique_names);

	foreach $sorted (@sorted_names) {
		($count, $cumulative_bit_score) = 0;
		foreach $subject (@subject_id) {
			# this essentially goes down the bitscore column in the .csv file. If the sseqids between $sorted and $subject match, then set $bit_score to the bitscore in the current row ($count). Add $bit_score to $cumulative_bit_score, and increment the $count. Otherwise, just increment the $count
			if ($sorted eq $subject) {
				$bit_score = $bit_score[$count];
				$cumulative_bit_score = $bit_score + $cumulative_bit_score;
				$count++;
			} else {
				$count++;
			}
		}
		# the reference genome that has the highest cumulative bitscore is the best reference genome. Therefore, if the $cumulative_bit_score for the current $file is better than the previous $best_bit_score, then $best_bit_score is swapped for $cumulative_bit_score, and $best_reference is updated to be the current $sorted
		if ($cumulative_bit_score > $best_bit_score) {
			$best_bit_score = $cumulative_bit_score;
			$best_reference = $sorted;
		}

		print "$file $best_reference $best_bit_score\n";
	}
	print "\n";
	open (OUTPUT, ">>best_references.csv");
	print OUTPUT "$file\t $best_reference\t $best_bit_score\n";
	close OUTPUT;
	}
