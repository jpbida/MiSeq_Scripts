#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Bio::Seq::Quality;
my $tail2="AAAGAAACAACAACAACAAC";
### Create an index of sequence ids from fasta file ###
my $lib=$ARGV[0];
my $io = Bio::AlignIO->new(-file   => "$lib",
                           -format => "clustalw" );

while(my $aln=$io->next_aln){
foreach $seq ($aln->each_seq) {
   ### Get ID Position from ideal ###
if($first==0){$seq1=$seq->seq;
$line =~ m{($tail2)}g;
my $pos = pos($line) - length $1;
print "Pos: $pos\n";
}
         if($first==1){$seq2=$seq->seq;}
         if($first==2){$seq3=$seq->seq;}
}
}

my $lib_reads=Bio::SeqIO->new( -file=> $lib, -format=>'fasta');

### Get the ID create a hash of the reverse compliment ###

my $read1=$ARGV[1];
my $read2=$ARGV[2];
my $id_reads =Bio::SeqIO->new( -file   => $read1,
           -format => 'fastq',
                 );


sub getMember{
my ($seq,$begin,$end)=@_;
### find tail2 ###
my $hash;
my $id="NONE";
my $seqout="NONE";
my $all=$seq;
### find T7 ###
if($seq=~m/$begin([AGCTagct]*)([AGCTagct]{8})$end/g){
$id=$2;
$seqout="$begin$1$2$end";
}
### if t7 and tail2 don't exist ###

### Take reverse complement ##

### find tail2 ###

### find T7 ###

### if t7 and tail2 don't exist ###
#print fail #
$hash->{id}=$id;
$hash->{seq}=$seqout;
$hash->{all}=$all;
return $hash;
}

##sort by distance
my %d;
@d{@inputs} = map { abs } adistr("pattern", @inputs);
my @d = sort { $d{$a} <=> $d{$b} } @inputs;
