#!/opt/local/bin/perl
use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq::Quality;
my $tail2="AAAGAAACAACAACAACAAC";
my $t7="TTCTAATACGACTCACTATAGG";
### Create an index of sequence ids from fasta file ###
my $lib=$ARGV[0];
my $io = Bio::AlignIO->new(-file   => "$lib",
                           -format => "clustalw" );

while(my $aln=$io->next_aln){
my $first=0;
my $seq1;
my $seq2;
my $seq3;
my $pos1=0;
my $pos2=0;
my $id0="";
my $id1="N";
my $id2="N";
my $id3="N";
foreach my $seq ($aln->each_seq) {
    $id3="N";
   ### Get ID Position from ideal ###
if($first==0){
$seq1=$seq->seq;
$id0=$seq->id;
$seq1 =~ m{($tail2)}g;
my $m1=$1;
$pos2 = pos($seq1) - length $m1;
$seq1=$seq->seq;
if($seq1 =~ m{($t7)}g){
my $m2=$1;
$pos1 = pos($seq1);
print "Pos: $pos1-$pos2\n";
}
}
         if($first==1){$seq2=$seq->seq;
        $id1=substr($seq2,$pos1,($pos2-$pos1));
         
         }
         if($first==2){
             $seq3=$seq->seq;
        $id3=substr($seq3,$pos2-8,15);
        $id2=substr($seq3,$pos1,($pos2-$pos1));
open(FH,">RNA/seq_$id3");
print FH ">wt1\n$id1\n>wt2\n$id2\n";}
close(FH);
$first++;
}
}
