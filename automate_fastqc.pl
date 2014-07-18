#! usr/bin/perl

use warnings;
use strict;
use Cwd;
use File::Copy;

# initialize variables
my (@Results, @total, $total, $strainName, $strainName_no_, $protein_count, @protein_count, $worksheet, @speciesFolder, @speciesTab, @strainFolder, @strainName, @log) = '';

my ($fastq_1, $fastq_2) = '';

# get the path for file/folder manipulations
my $path = getcwd;
#o open all folders in 
opendir(DIR, $path) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	next if $file =~ /^\.\.?$/;
	next if $file =~ /Backups/;
	next if $file =~ /test/;
	if (-d $path."/".$file){
		push (@speciesFolder, $path."/".$file);
		push (@speciesTab, $file);
		}
}

foreach (@speciesTab) {
	#print $_, "\n";
	system ("mkdir -p $path/Trimmed/$_");
	chdir ("$path/$_");
	$fastq_1 = "$_" . "_R1_001.fastq";
	$fastq_2 = "$_" . "_R2_001.fastq";
	print $fastq_1, "\n";
	system ("trimmomatic-0.30.jar PE -threads 24 -phred33 -trimlog $_.log $fastq_1 $fastq_2 $_.paired1.fastq $_.unpaired1.fastq $_.paired2.fastq $_.unpaired2.fastq ILLUMINACLIP:/home/blais/Bioinformatics/Trimmomatic-0.30/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36");	
	#system ("fastq_quality_filter -t 30 -i $fastq_1 -o $path/Trimmed/$_/$fastq_1"); 
}
