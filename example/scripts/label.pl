####################################################################################################
### label.pl
###
### 2012 Jp Bida
###
### Given a file with a single sequence on each line, this script compares each sequence to a wt
### sequence and creates a labeled based on the positon changed
###
### inputs: seq.lst wt_seq
### see label.sh for an example
####################################################################################################
#!/usr/bin/perl

my $input = $ARGV[0];
my $wt = $ARGV[1];

print "$input\n";
open(FH,"<$input");

my @lines=<FH>;

foreach my $line (@lines){
chomp($line);
my $name="";
for $i (0 .. (length($wt)-1)) {
   $source_base = substr($wt,$i,1);
   $str_base    = substr($line,$i,1);

  if ($source_base ne $str_base) {
$name.="$source_base$i$str_base:";  
}
}
if($name eq ""){$name="wt"}
print ">$name\n$line\n";

}
