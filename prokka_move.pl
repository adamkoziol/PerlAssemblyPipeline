#! /usr/bin/perl

use warnings;
use strict;
use Cwd;
use File::Copy;

my (@speciesTab, @speciesFolder) = '';

# get the path for file/folder manipulations
my $path = getcwd;
#o open all folders in 
opendir(DIR, $path) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	next if $file =~ /^\.\.?$/;
	if (-d $path."/".$file){
		push (@speciesFolder, $path."/".$file);
		push (@speciesTab, $file);
		}
}

closedir(DIR);

foreach my $species (@speciesTab) {
		
	print $species;
	system ("mkdir -p $path/00_Prokka_predictions");	
	chdir ("$path/$species");
	system ("cp *.ffn $path/00_Prokka_predictions");
	chdir ("$path");
}
