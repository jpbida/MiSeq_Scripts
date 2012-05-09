#!/usr/bin/perl
use File::Basename;

$filename = @ARGV[0];
$basename = basename($filename,".cpp");
$count=1;


open(FILE,$filename);
open(OUTFILE,">".$basename.$count.".cpp");
while($line = <FILE>) {
  if($line =~ /\/\/ part/){
    close(OUTFILE);
    $count++;    
    open(OUTFILE,">".$basename.$count.".cpp");
  }
  else {
    print OUTFILE $line;
  }
}
close(OUTFILE);
