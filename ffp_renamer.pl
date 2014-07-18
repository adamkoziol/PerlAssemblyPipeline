#! /usr/bin/perl
use Cwd;
use warnings;
use strict;

my (@files, $strain, $short_name, @tree, %hash, $values, $keys) = '';

#this sprintf command forces the $count variable to be displayed as "000" instead of "0"
my $count = sprintf("%03d", 0);
my $path = getcwd;

chdir "$path/00_FFP_renamer";

@files = glob("*.fasta"); #gets all the ".fasta" files into @files - if there's any ".fa" or ".fas" or other file extensions, this won't work

# using the three-part file opening scheme
open (OUTPUT, ">", "18_renamed_species.txt"); 

foreach (@files) {
	($strain = $_) =~ s/.fasta//g; #removes the trailing ".fasta" from the filenames
	($short_name = $strain) =~ s/$strain/_$count\_/; # replaces the filename with "_$count_" or "_000_ up to _number of elements in @files_"	
	$count++;
	print OUTPUT $short_name, "\n";
	$hash{$strain} = $short_name # this hash associates the $strain name variable with the $short_name variable. This will be useful later on
}

# this ffp command is a pipeline of the various ffp program available. Right now, it's using an l-mer size of 18; modify as required
system ("ffpry -l 18 *.fasta | ffpcol | ffprwn |ffpjsd -p 18_renamed_species.txt | ffptree > 18_tree");

open (INPUT, "<", "18_tree");

@tree = (<INPUT>);
close INPUT;

open (OUTPUT, ">", "18_renamed_tree");
foreach (@tree) {
	while (($keys, $values) = each(%hash)) { # this gets the keys and values from %hash set to $keys and $values, respectively
		$_ =~ s/$values/$keys/g; # substitutes the $value for the $key from the %hash if the $value is found in the tree file
		}
	print OUTPUT $_;
}
