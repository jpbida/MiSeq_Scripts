#!/opt/local/bin/perl
use strict;
use Bio::SeqIO;
use Bio::Seq::Quality;
use Bio::AlignIO;
#TGGA
#ACTC
#ATCT
my $tail2="AAAGAAACAACAACAACAAC";
### Create an index of sequence ids from fasta file ###
my $read1=$ARGV[0];
my $read2=$ARGV[1];
my $id_reads =Bio::SeqIO->new( -file   => $read1,
           -format => 'fastq',
                 );
my $end_reads =Bio::SeqIO->new( -file   => $read2,
           -format => 'fastq',
                 );


my $i=0;
my $w=0;
while((my $aln1=$id_reads->next_seq) && (my $aln2=$end_reads->next_seq)){
  #  $aln1->force_flush('1');
#my $rc=$aln1->revcom;
#my $rsq=$rc->seq;
my $rsq="AAA";
my $ID="N";
my $mid="N";
#$rsq=~m/([AGCT]{8})$tail2([AGCT]{4})/g;
my $ID="AAA";
my $mid="AAA";
#my $ID=$1;
#my $mid=$2;
my $clus=$aln1->id;
#$clus=~m/([\d]+\:[\d]+)$/;
my $clus_id=$1;
### Get Tail2 Id ###
if($ID ne "" && $mid ne ""){
my $seq2=$aln2->seq;
my $clus2=$aln2->id;
#$clus2=~m/([\d]+\:[\d]+)$/;
my $clus_id2=$1;
#if($clus_id2 ne $clus_id){print "Cluster IDs out of order\n";}
#unless(-d $mid){mkdir $mid or die;}
#open(FH,">>$mid/seq_$ID.fasta");
#print FH "\n>frag\n$seq2\n";
#close(FH);
$i++;
}
if($w % 1000==0){print "f1:$i:$w\n";}
$w++;
}

