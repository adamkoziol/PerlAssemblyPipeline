#!/usr/bin/perl
use Cwd;
use File::Copy;

$path = getcwd();
@files = glob("*.gz");
foreach $file (@files){
	($folder = $file) =~ s/_.*//g;
	($unzip = $file) =~ s/.gz//g;
	($filename = $file) =~ s/S.+_L001_|.gz//g;
	print $filename."\n";
	mkdir $folder;
	move($file,"$folder/$file");
	system("gzip -d $folder/$file");
	if (-f $unzip){unlink "$folder/$file"};
	move("$folder/$unzip","$folder/$filename") 
}
