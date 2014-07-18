#!usr/bin/perl

use warnings;
use strict;
use Cwd;
use File::Copy;


my (@speciesFolder, @speciesTab, @folder_contents) = '';

# get the path for file/folder manipulations
my $path = getcwd;
#o open all folders in 

chdir ("$path/E_coli_trimmed");
opendir(DIR, $path) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	next if $file =~ /^\.\.?$/;
	next if $file =~ /Velvet_Assembled_contigs/;
	next if $file =~ /Velvet_Log_Files/;
	if (-d $path."/".$file){
		push (@speciesFolder, $path."/".$file);
		push (@speciesTab, $file);
		}
}
closedir(DIR);
#print @speciesTab;


foreach (@speciesTab) {
	chdir ("$path/E_coli_trimmed/$_");
	#rename "contigs.fa", "$_" . ".fa";
	system ("mkdir -p $path/E_coli_trimmed/Velvet_Assembled_contigs");
	system ("mkdir -p $path/E_coli_trimmed/Velvet_Log_files");	
	@folder_contents = glob "*auto_data*";
	foreach my $autodata (@folder_contents) {
		chdir ("$path/$_/$autodata");	
		copy "contigs.fa", "$path/E_coli_trimmed/Velvet_Assembled_contigs";
		copy "Log", "$path/E_coli_trimmed/Velvet_Log_files";
		chdir ("$path/E_coli_trimmed/Velvet_Assembled_contigs/");
		rename "contigs.fa", "$_" . ".fa";
		chdir ("$path/E_coli_trimmed/Velvet_Log_files");
		rename "Log", "$_" . "_Log" . ".txt";
		chdir "$path";
		}
	}
