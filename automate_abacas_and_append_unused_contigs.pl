#! usr/bin/perl

use warnings;
use strict;
use Cwd;
use File::Copy;

my @input = '';
my $seqio = '';
my $path = getcwd;
my $line = '';
my @sequence = '';
my @C_Jejuni1 = ("00.2425");

chdir ("$path/C_jejuni_contigs_and_reference");
my @C_Jejuni = glob "*";

chdir ("$path/C_coli_contigs_and_reference");
my @C_coli = glob "*";

chdir ("$path/C_fetus_contigs_and_reference");
my @C_fetus = glob "*";

chdir ("$path/Arcobacter_contigs_and_reference");
my @arco = glob "*";

my @folder_contents = '';

=com
####################C.jejuni
foreach (@C_Jejuni) {
	(@input, @sequence, $seqio) = '';	
	chdir ("$path/C_jejuni_contigs_and_reference/$_");
	system ("abacas -r CJ_NCTC_11168.fasta -q contigs.fa -p nucmer -m -b -N -a -o $_");
	system ("cat $_.NoNs.fasta $_.contigsInbin.fas > $_.concatenated.fas"); 	
	
	open (INPUT, "<$_.concatenated.fas") or die "That file doesn't exist!";
	@input = (<INPUT>);
	open (OUTPUT, ">$_.appended.contigs.fa") or die "I couldn't write the output!"; 
	print OUTPUT $input[0];
	foreach $line (@input) {
		next if $line =~ />/;		
		chomp $line;
		push (@sequence, $line);
	}

	$seqio = join ('', @sequence);
	
	for (my $pos = 0 ; $pos < length($seqio) ; $pos += 60 ) {
		print OUTPUT (substr($seqio, $pos, 60)), "\n";
		}
	copy "$_.appended.contigs.fa", "$path/appended_contigs/C_jejuni/$_.ordered.fasta";	
	}
=cut	
#####################C.coli
foreach (@C_coli) {
	(@input, @sequence, $seqio) = '';	
	chdir ("$path/C_coli_contigs_and_reference/$_");
	system ("abacas -r *.fasta -q contigs.fa -p nucmer -m -b -N -a -o $_");
	system ("cat $_.NoNs.fasta $_.contigsInbin.fas > $_.concatenated.fas"); 	
	
	open (INPUT, "<$_.concatenated.fas") or die "That file doesn't exist!";
	@input = (<INPUT>);
	open (OUTPUT, ">$_.appended.contigs.fa") or die "I couldn't write the output!"; 
	print OUTPUT $input[0];
	foreach $line (@input) {
		next if $line =~ />/;		
		chomp $line;
		push (@sequence, $line);
	}

	$seqio = join ('', @sequence);
	
	for (my $pos = 0 ; $pos < length($seqio) ; $pos += 60 ) {
		print OUTPUT (substr($seqio, $pos, 60)), "\n";
		}
	copy "$_.appended.contigs.fa", "$path/appended_contigs/C_coli/$_.ordered.fasta"	
	}
=com
#####################C.fetus
foreach (@C_fetus) {
	(@input, @sequence, $seqio) = '';	
	chdir ("$path/C_fetus_contigs_and_reference/$_");
	system ("abacas -r ../../424819875_C_fetus_venerealis.fasta -q contigs.fa -p nucmer -m -b -N -a -o $_");
	system ("cat $_.NoNs.fasta $_.contigsInbin.fas > $_.concatenated.fas"); 	
	
	open (INPUT, "<$_.concatenated.fas") or die "That file doesn't exist!";
	@input = (<INPUT>);
	open (OUTPUT, ">$_.appended.contigs.fa") or die "I couldn't write the output!"; 
	print OUTPUT $input[0];
	foreach $line (@input) {
		next if $line =~ />/;		
		chomp $line;
		push (@sequence, $line);
	}

	$seqio = join ('', @sequence);
	
	for (my $pos = 0 ; $pos < length($seqio) ; $pos += 60 ) {
		print OUTPUT (substr($seqio, $pos, 60)), "\n";
		}
	copy "$_.appended.contigs.fa", "$path/appended_contigs/C_fetus/$_.ordered.fasta"	
	}
		
#####################Arcobacter
foreach (@arco) {
	(@input, @sequence, $seqio) = '';	
	chdir ("$path/Arcobacter_contigs_and_reference/$_");
	system ("abacas -r ../../Arcobacter_RM4018.fasta -q contigs.fa -p nucmer -m -b -N -a -o $_");
	system ("cat $_.NoNs.fasta $_.contigsInbin.fas > $_.concatenated.fas"); 	
	
	open (INPUT, "<$_.concatenated.fas") or die "That file doesn't exist!";
	@input = (<INPUT>);
	open (OUTPUT, ">$_.appended.contigs.fa") or die "I couldn't write the output!"; 
	print OUTPUT $input[0];
	foreach $line (@input) {
		next if $line =~ />/;		
		chomp $line;
		push (@sequence, $line);
	}

	$seqio = join ('', @sequence);
	
	for (my $pos = 0 ; $pos < length($seqio) ; $pos += 60 ) {
		print OUTPUT (substr($seqio, $pos, 60)), "\n";
		}
	copy "$_.appended.contigs.fa", "$path/appended_contigs/Arcobacter/$_.ordered.fasta";	
	close OUTPUT;	
	}
=cut	
	



		
