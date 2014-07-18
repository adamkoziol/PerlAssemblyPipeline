#! usr/bin/perl

#use warnings;
use strict;
use Cwd;
use Spreadsheet::ParseExcel;
use Spreadsheet::WriteExcel;
use File::Copy;


# initialize variables
my (@Results, $total, $strainName, $strainName_no_, $protein_count, @contigs, @protein_count, $worksheet, @speciesFolder, @speciesTab, @strains, $strain_name, @log, @gene_count, $gene_count, @contig_count, $contig_count, @path, $strainPath, $gene_name, $tail, @total, @n50, @max) = '';

# get the path for file/folder manipulations
my $path = getcwd;
my $count = 0;

open (OUTPUT, ">Assembly_report.csv");
print OUTPUT "Strain\t Number of Contigs\t Longest Contig\t N50\t Total Bases\t Number of Predicted Genes\n";


#chdir ("$path/Velvet_Log_files");
#@log = glob("*.txt");

chdir ("$path/Velvet_Assembled_contigs");
@contigs = glob("*.fa");

=commy @predicted_genes = glob("*.ffn");

foreach (@predicted_genes) {
	($strain_name = $_) =~ s/.ordered.ffn//g;
	push(@strains, $strain_name);
	#print $strain_name,"\n";
}

foreach my $genes (@predicted_genes) {
	$gene_count = 0;		
	open (GENES, "<", $genes) or die "No file!";
	#print "These are the contig_files ", @contig_files, "\n\n";
	# $gene_count is re-initialized to 0 for each loop iteration		
	# as the contigs.fa file is read by perl ...		
	while (<GENES>) {
		# if a ">" is encountered ...			
		if (/>/)  {
		# then increment the count			
			$gene_count++;
			}
		}
	close GENES;		
	print $gene_count,"\n";		
	push (@gene_count, $gene_count);		
}

=cut

foreach my $contigs (@contigs) {
	#print $contigs, "\n";	
	($strain_name = $contigs) =~ s/.fa//g;
	push(@strains, $strain_name);	
	$contig_count = 0;		
	open (CONTIGS, "<", $contigs) or die "No file!";
	#print "These are the contig_files ", @contig_files, "\n\n";
	# $gene_count is re-initialized to 0 for each loop iteration		
	# as the contigs.fa file is read by perl ...		
	while (<CONTIGS>) {
		# if a ">" is encountered ...			
		if (/>/)  {
		# then increment the count			
			$contig_count++;
			}
		}
	close CONTIGS;		
	#print $contig_count,"\n";		
	push (@contig_count, $contig_count);		
}							

foreach (@strains) {
	chdir ("$path/Velvet_Log_files");	
	print $_, "\n";	
	open(INPUT, "<$_" . "_Log.txt"); 
	@log = (<INPUT>);
	$tail = pop @log;
	$tail =~ s/,//g;
	my ($junk1, $junk2, $junk4, $junk5, $junk6, $junk7, $junk8, $junk9, $n50, $junk10, $max, 			$junk11, $total) = split(" ", $tail);
	push (@total, $total);
	push (@n50, $n50);
	push (@max, $max);	
	close INPUT;
	print OUTPUT "$_\t $contig_count[$count]\t $max[$count]\t $n50[$count]\t $total[$count]\t $gene_count[$count]\n";
	$count++;
	
}			
			
	

=com

# invoke the parse::excel module and tell it to open the categories.xls spreadsheet or quit
my $parser = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse('categories.xls');
if (!defined $workbook) {
	die $parser->error(), ".\n";
}

# this loop seems to essentially allow access to every cell in a worksheet while it is running
# this seems inefficient - I think a better way would be to proceed in a similar fashion to the way I populated the spreadsheet later
# using a hash of arrays and a counter. That way only one loop would need to be written, as long as the counter incremented properly
# I may revisit this
for my $worksheet ($workbook->worksheets()) {
	
	# define the $row_min and $row_max variables that determine where the first row with data ($row_min) and the last row with data ($row_max)
	# prevents having to provide row values for where to start/stop reading
	# I don't think I actually needed to do this, as the ranges still had to be defined everytime, and it didn't make anything easier
	my ($row_min, $row_max) = $worksheet->row_range();
	
	# read the data from the spreadsheet into the appropriate arrays
	# this is done in separate, but similar loops, so only the first instance is commented
	# motilityHigh array gets populated
	# the range to read is between rows 5 and 14	
	for my $row (5 .. $row_max) {
		for ($countCol = 1; $countCol <=20; $countCol+=2){			
			for my $col ($countCol) {
				# values from the cells are added to the scalar $cell			
				my $cell = $worksheet->get_cell($row, $col);
				# as long as data exists in the cell, continue adding to $cell 			
				next unless $cell;
				# add each value of $cell to the appropriate array			
				push (@Results, $cell->value());
			}
		}
	}
}


# for some reason, though I have no idea why, when I didn't shift the arrays, they would yield an empty cell when printed in spreadsheets in
# the loops below
shift @Results;

# as there were some inconsistencies in the naming schemes between the various inputs from which this program reads, this loop changes hypens and
# underscores to periods
foreach (@Results) {
	$_ =~ s/-/./g;
	$_ =~ s/_/./g;
}

# this array of suffixes is going to be used to name the tabs of the new spreadsheet
# my @array_of_suffixes = ['mH', 'mL', 'bH', 'bL', 'aH', 'aL', 'iH', 'iL', 'H005H', 'H005L', 'H01H', 'H01L', 'H1H', 'H1L', 'H2H', 'H2L', 'H4H', 'H4L', 'H8H', 'H8L'];

# this hash of arrays links the key to the appropriate array.
# using a hash allows for the association of the array with a text shortcut for naming tabs 
my %hash_of_arrays = (
	'Results' => [@Results], 
);

# the original code for the hash of arrays 
=com
foreach my $scalar (sort keys %hash_of_arrays) {
	print "The name of the array is $scalar \n";
	foreach (@{$hash_of_arrays{$scalar}}) {
		print "\t$_\n";
	}
}


# creates a new excel spreadsheet with the name assembly_report.xls
my $workbook_out = Spreadsheet::WriteExcel->new("ProkkaAssemblyReport.xls");

# set up the excel file
# creates a variable $scalar that is the keys from the % hash of arrays and acts in a similar fashion to $_ 
# and therefore allows the looping through of every "component" of the array
# so for every key in the hash, the keys are sorted and ...
my $path = getcwd;
foreach my $i (@speciesTab){
	my @strainName = '';
	opendir(DIR, "$path/$i") or die "can't find $path/$i: $!";
	# these arrays are always re-initialized to empty through each iteration of the loop to ensure that only the data from
	# that particular interation is contained within the arrays	
	($total,  $protein_count) = '';
	# adds worksheets to the workbook named based on $scalar	
	$worksheet = $workbook_out->add_worksheet($i);	
	
	# $strainName is essentially set to be all the data from all the arrays in the hash of arrays, so all the strain names found in the original
	# spreadsheet categories.xls, but these data are associated with the appropriate key 
#	print scalar @strainName, "\n";
	while (defined(my $strainName = readdir(DIR))) {
			next if $strainName =~ /^\.\.?$/;
			if (-d $path ."/". $i."/".$strainName){
#			push @strainTab, $i."/".$strainName;

			$strainName =~ s/.ordered.prokka//g;
			#print $strainName. "\n";
			push @strainName, "$strainName";			
			
			
			#$count = @strainName;
			#print $count, "\n";
			
			
			# the array header exists because it seemed easier to use the write_row function to write an array to the spreadsheet
			# using spreadsheet:writeexcel rather than to try and figure out if I could use that function to enter rows manually 		
			#my @header = ("\t","\t","\t", "Velvet");#, "\t","\t","\t", "\t", "CLC Genomics Workbench");
			
			# this variable exists because of the inconsistencies in the naming. It depends on when the analyses were performed
			# if the naming conventions have already been streamlined, then this may cause problems finding some sequences
	#		$strainName_no_ = $strainName;
	#		$strainName_no_ =~ s/\./_/g;
			# the description of the data being entered into the spreadsheet
			my @categories = ("strain", "total bases", "prokka proteins"); #"\t", "n50", "longest 			contig", "number of contigs", "total bases", "prokka proteins");
			# write the data to the spreadsheet
			# this requires a reference to an array (\@header)		
			# $worksheet->write_row(0,0, \@header);
			$worksheet->write_row(0,0, \@categories);
=com

			# enter the folder containing the ouputs from the assembly of the appropriate strain
			# ENSURE THAT YOUR NAMING SCHEME MATCHES THIS!!!		
			chdir ("$path/$strainName");
			# the output from velvetg is simply called "Log". 
			# the log file is opened - maybe I should have included an "or die" option?
			open (INPUT, "<Log");
			# the contents of log are read into an array		
			@log = (<INPUT>);
			# as the log file contains multiple iterations through the program, we only want the final (best) iteration.
			# these data are found at the end (tail) of the file (array) 
			$tail = pop @log;
			# the various statistics regarding the assembly are separated by commas		
			$tail =~ s/,//g;
			# the important data are read into variables, while the rest are discarded
			my ($junk1, $junk2, $junk4, $junk5, $junk6, $junk7, $junk8, $junk9, $n50, $junk10, $max, $junk11, $total) = split(" ", $tail);
			# appropriate data are pushed to the appropriate arrays
			push (@total, $total);
			push (@n50, $n50);
			push (@max, $max);
			
	
			# count contigs function
			# the contigs are found in the contigs.fa file
			@contig_files = glob "*.fa";
			$gene_count = 0;
			foreach my $contigs (@contig_files) {
				open (INPUT, "<", $contigs) or die "No file!";
				#print "These are the contig_files ", @contig_files, "\n\n";
				# $gene_count is re-initialized to 0 for each loop iteration		
				# as the contigs.fa file is read by perl ...		
				while (<INPUT>) {
					# if a ">" is encountered ...			
					if (/>/)  {
					# then increment the count			
					$gene_count++;
					}
				}
			}
			# add the final tally of $gene_count to individual component of the contig_count array		
			#print $gene_count, "\n";		
			push (@contig_count, $gene_count);
=cut
			# base counter
			# the outputs are in a different directory named "$strainName.prokka"
			#chdir ("$path/$i/$strainName" . ".ordered.prokka");
			#open (INPUT, "<",getcwd."/"."$strainName".".ordered" .".fna") or die $!;
			#my $total = 0;
			#readline(INPUT);
			#while (<INPUT>) {
			#	$total += length($_);
			#}
			#push (@total, $total);			 
	
			# protein count
			# similar to the contig count function
			# the outputs are in a different directory named "$strainName.prokka"
			
			# protein files have a .faa extension
			#open (INPUT, "<",getcwd."/"."$strainName".".ordered" .".faa");
			#my $protein_count = 0;
			#while (<INPUT>) {
			#	if (/>/)  {
			#		$protein_count++;
			#	}			
			#}
			#push (@protein_count, $protein_count);
			#if (($protein_count > 1000) && ($protein_count < 2500)) {
			#	my $strainPath = getcwd."/"."$strainName".".ordered" .".faa";
			#	push (@path, $strainPath);
#				print $protein_count;
#			}
#			print @path, "\n";
			#write the outputs
			#$worksheet->write_col(2,1, \@n50);
			#$worksheet->write_col(2,2, \@max);

			#$worksheet->write_col(2,3, \@contig_count);
			
=com	
		}
		
	}

	#$worksheet->write_col(1,2, \@protein_count);
	$worksheet->write_col(1,1, \@total);
	$worksheet->write_col(0,0, \@strainName);
	$worksheet->write_col(1,3, \@path);
	#not sure why I need this but it clears the arrays of values
	splice (@total);
	splice (@protein_count);
	splice (@path);
}

