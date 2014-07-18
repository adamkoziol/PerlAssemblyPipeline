#! /usr/bin/perl
use Cwd;
use File::Copy;
use warnings;
use strict;

my $PipePath="/home/blais/Bioinformatics/PIPELINE_Nick_HC";

# get the path for file/folder manipulations
my $path = getcwd;
#o open all folders in 
my ($mpath, @fastq_1, @fastq_2) = "";
opendir(DIR, $path) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	next if $file =~ /^\.\.?$/;
	next if $file =~ /Backups/;
	next if $file =~ /00_Unprocessed/;
	next if $file =~ /Trimmed/;
	next if $file =~ /E_coli_trimmed/;
	next if $file =~ /E_coli_untrimmed/;
	if (-d $path."/".$file){
		#Currently with 24 cores
		$mpath = $path;
		(@fastq_1, @fastq_2) = "";
		$mpath =~ s/ /\\ /g; # The space was throwing off velvet
		chdir $path."/".$file;
		#system ("pwd");
		@fastq_1 = glob("*R1_001.fastq");
		@fastq_2 = glob("*R2_001.fastq");
#		system("cd $path/$file");
		#print $fastq_1[0],"\n";
		system("mkdir -p $path/00_Unprocessed/$file; cd $path/00_Unprocessed/$file; perl $PipePath/VelvetOptimiser-2.2.5/VelvetOptimiser.pl -s 65 -e 111 -f '-shortPaired -fastq $mpath/$file/$fastq_1[0] $mpath/$file/$fastq_2[0]' -t 24 -c 'n50*Lcon/tbp+log(Lbp)'");
		
	}
}


