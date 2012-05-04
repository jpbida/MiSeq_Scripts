#!/opt/local/bin/perl
use Bio::SeqIO;
#use Bio::Factory::EMBOSS;
#use Bio::SeqIO;
#use Bio::AlignIO;
#use Bio::Tools::pSW;
#use Bio::PrimarySeq;
#use Bio::Tools::Run::Alignment::Clustalw;
#use Bio::Tools::Run::Alignment::TCoffee;
#use Bio::Tools::Run::StandAloneBlast;
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

### Read in a FASTA file ###
$lib=Bio::SeqIO->new(-file=>"/Users/jpbida/labwork/oligos/library_final.txt",'-format'=>'Fasta');




### Sequnece can be accessed from the sequence object ###
#while ( my $seq = $lib->next_seq() ) {
        #print $seq->seq()."\nRevCom";
        #my $rc=$seq->revcom();
        #print $rc->seq()."\n";
         #   }
         my @libseqs;
while(my $seq=$lib->next_seq()){
    my $hash;
        my $rc=$seq->revcom();
$id1=getMember($seq->seq(),"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");
$id2=getMember($rc->seq(),"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");
if($id1->{id}=~m/NONE/g){
if($id2->{id}=~m/NONE/g){
}else{
$hash=$id2;
}
}else{
$hash=$id1;
}

## Extract seqeunces from the reads ###
if($hash->{seq}=~m/TTCTAATACGA/g){
push(@libseqs,$hash);
}
}

### Foreach well in the 96 well plate ####
my @cols=("A" .. "H");
my @rows=("01","02","03","04","05","06","07","08","09","10","11","12");
print "Rows: $#rows\n";
print "Cols: $#cols\n";
my @readseqs;
foreach my $row (@rows){
    foreach my $col (@cols){
        my $name="$col$row";
        print "$name\n";
        my $seq1="";
        my $seq2="";
        my $seq1rc="";
        my $seq2rc="";
$reads1=Bio::SeqIO->new(-file=>"./all1.seqs",'-format'=>'Fasta');
$reads2=Bio::SeqIO->new(-file=>"./all2.seqs",'-format'=>'Fasta');
while ( my $seq = $reads1->next_seq() ) {
my $id=$seq->id();
if($id=~m/$name/i){
$seq1=$seq->seq();
$seq1p=$seq->revcom();
$seq1rc=$seq1p->seq();
break;
}
}
while ( my $seq = $reads2->next_seq() ) {
my $id=$seq->id();
if($id=~m/$name/i){
$seq2p=$seq->revcom();
$seq2=$seq2p->seq();
$seq2rc=$seq->seq();
break;
}
}
### If any of the possible sequences has an EteRNA ID add the sequence ID to the alingment file ####
### Get Tail2 or T7 ###
my $id1=getMember($seq1,"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");
my $id2=getMember($seq1rc,"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");
my $id3=getMember($seq2,"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");
my $id4=getMember($seq2rc,"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");
my $hash;
$hash->{well}=$name;
$hash->{id}="NONE";
$hash->{seq1}="NONE";
$hash->{seq2}="NONE";

if($id1->{id}=~m/NONE/g){
if($id2->{id}=~m/NONE/g){
if($id3->{id}=~m/NONE/g){
if($id4->{id}=~m/NONE/g){


}else{
$hash->{id}=$id4->{id};
$hash->{seq1}=$id2->{all};
$hash->{seq2}=$id4->{all};
print $id2->{id}.":".$id4->{id}."\n";


}
}else{
$hash->{id}=$id3->{id};
$hash->{seq1}=$id1->{all};
$hash->{seq2}=$id3->{all};
print $id3->{id}.":".$id1->{id}."\n";


}


}else{
$hash->{id}=$id2->{id};
$hash->{seq1}=$id2->{all};
$hash->{seq2}=$id4->{all};
print $id2->{id}.":".$id4->{id}."\n";


}
}else{
$hash->{id}=$id1->{id};
$hash->{seq1}=$id1->{all};
$hash->{seq2}=$id3->{all};
print $id1->{id}.":".$id3->{id}."\n";

}

if($hash->{seq1}=~m/TTCTAATACGA/g || $hash->{seq2}=~m/TTCTAATACGA/g){
push(@readseqs,$hash);
}

}
}

my $w=0;
foreach $mem(@readseqs){
    foreach $lib(@libseqs){
if($lib->{id} eq $mem->{id}){
open(FH,">seqs_$w.fasta");
print FH ">lib".$mem->{well}."\n$lib->{seq}\n>mycro1\n$mem->{seq1}\n>mycro2\n$mem->{seq2}\n";
close(FH);
break;
system("clustalw2 seqs_$w.fasta");
$w++;
}
    }
}
