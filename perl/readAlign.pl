#!/opt/local/bin/perl
use Bio::AlignIO;
my $file=$ARGV[0];
my $io = Bio::AlignIO->new(-file   => "$file",
                           -format => "clustalw" );


while(my $aln=$io->next_aln){
print $file." ".$aln->length." ".$aln->percentage_identity." ";
my $insertion=0;
my $deletion=0;
my $first=0;
my $seq1="";
my $seq2="";
my $seq3="";
foreach $seq ($aln->each_seq) {
         if($first==0){$seq1=$seq->seq;}
         if($first==1){$seq2=$seq->seq;}
         if($first==2){$seq3=$seq->seq;}
}

         my $sq=$seq->seq;
         print $sq." ";

                  if($first==0){
                      $first=1;
                    #### Remove all missing ####


                  $insertion++ while($sq =~m/\-/g);}else{
                  $deletion++ while($sq=~m/\-/g);
                  }
                  }
print $insertion." ".$deletion."\n"; 
}
