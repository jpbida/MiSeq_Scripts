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
return $hash;
}

### Read in a FASTA file ###
$lib=Bio::SeqIO->new(-file=>"/Users/jpbida/labwork/oligos/library_final.txt",'-format'=>'Fasta');
$reads=Bio::SeqIO->new(-file=>"./all.seqs",'-format'=>'Fasta');

while ( my $seq = $lib->next_seq() ) {

#        print $seq->seq()."\n";
#        getMember($seq->seq(),);
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

my @readseqs;
while ( my $seq = $reads->next_seq() ) {
#        print $seq->seq()."\n";
#        getMember($seq->seq(),);
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
push(@readseqs,$hash);
}
}

my $w=0;
foreach $mem(@readseqs){

    foreach $lib(@libseqs){
if($lib->{id} eq $mem->{id}){
open(FH,">seqs_$w.fasta");
$w++;
print FH ">lib\n$lib->{seq}\n>mycro\n$mem->{seq}\n";
close(FH);
break;
}
    }
}

### Foreach sequence identify the sequence ID and find the matches 
### in the library.  Perform both reverse complement and regular.

my @lib=<FH>;
my @d;
my $st=0;
foreach my $line(@lib){
if($line=~m/\>/g){
if($st!=0){
push(@d,$str);
$str="";
}else{
$st=1;
$str="";
}

}

}

my @lines=<SQ>;
my $seq1="";
my $start=0;
foreach my $line (@lines){
chomp($line);
if($line=~m/\>/){
    if($start==0){
$start=1;
$seq1="";
}else{
##my $seq2=getMember($seq1,"TTCTAATACGACTCACTATAGG","AAAGAAACAACAACAACAAC");

if($seq2!=-1){

my %d;
@d{@inputs} = map { abs } adistr("pattern", @inputs);
my @d = sort { $d{$a} <=> $d{$b} } @inputs;

# this should produce the same alignment as 'water'
#my $factory = new Bio::Tools::Run::Alignment::Clustalw('ktuple' => 2, 
#                           'matrix' => 'BLOSUM');
$seq3=$seq2;
#my $aln = $factory->pairwise_alignment($seq2,$seq3);
#my $d = $stat->D_JukesCantor($aln)
#print "Alignment: $d\n";


}
$seq1="";
}

}else{
$seq1.=$line;
}

}
