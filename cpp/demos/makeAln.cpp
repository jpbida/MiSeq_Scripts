// FRAGMENT(includes)
#define SEQAN_PROFILE // enable time measurements
#include <seqan/file.h>
#include <iostream>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/basic.h>
using namespace seqan;
typedef StringSet<CharString> THaystacks;
// FRAGMENT(open_file)
int main (int argc, char const * argv[])
{
	SEQAN_PROTIMESTART(loadTime);
//Read in command line arguments pointing to the paired end reads and the expected library
//Example
//      ./makeAln tail2reads 3primereads library.fasta

	MultiSeqFile multiSeqFile1;
	MultiSeqFile multiSeqFile2;
	MultiSeqFile multiSeqFile3;
	if (argc < 3 || (!open(multiSeqFile1.concat, argv[1], OPEN_RDONLY) ||
        !open(multiSeqFile2.concat,argv[2],OPEN_RDONLY) || !open(multiSeqFile3.concat,argv[3],OPEN_RDONLY)))
		return 1;


// FRAGMENT(guess_and_split)
	AutoSeqFormat format1;
	guessFormat(multiSeqFile1.concat, format1);
	split(multiSeqFile1, format1);
	
	AutoSeqFormat format2;
	guessFormat(multiSeqFile2.concat, format2);
	split(multiSeqFile2, format2);
	
    AutoSeqFormat format3;
	guessFormat(multiSeqFile3.concat, format3);
	split(multiSeqFile3, format3);

// FRAGMENT(reserve)
	unsigned seqCount1 = length(multiSeqFile1);
    unsigned seqCount2 = length(multiSeqFile2);
    unsigned seqCount3 = length(multiSeqFile3);

// FRAGMENT(read_sequences)
	String<char> seq1;
	String<char> seq2;
	String<char> seq3;
	CharString qual1;
	CharString qual2;
	CharString id1;
	CharString id2;
    THaystacks haystacks;

// Index library sequences for alignment 
std::cout << "Indexing Sequences..";
for(unsigned j=0; j< seqCount3; j++){
		assignSeq(seq3, multiSeqFile3[j], format3);    // read sequence
        appendValue(haystacks, seq3);
}
Index<THaystacks> index(haystacks);
Finder<Index<THaystacks> > finder(haystacks);
std::cout << "completed" << std::endl;

for (unsigned i = 0; i < seqCount1; ++i)
	{
		assignSeq(seq1, multiSeqFile1[i], format1);    // read sequence
		assignQual(qual1, multiSeqFile1[i], format1);  // read ascii quality values
		assignSeqId(id1, multiSeqFile1[i], format1);   // read sequence id
		
        assignSeq(seq2, multiSeqFile2[i], format2);    // read sequence
		assignQual(qual2, multiSeqFile2[i], format2);  // read ascii quality values
		assignSeqId(id2, multiSeqFile2[i], format2);   // read sequence id
reverseComplement(seq1);
//Get ID and MID
//		std::cout << seq1 << std::endl;
String<char> ndl = "[AGCT]{8}AAAGAAACAACAACAACAAC[AGCT]{4}";
Finder<String<char> > finder1(seq1);
Pattern<String<char>, WildShiftAnd> pattern1(ndl);

if(find(finder1, pattern1)) {
int pos1=position(finder1);
String<char> mid=infixWithLength(seq1, (pos1-4),4);
String<char> uid=infixWithLength(seq1, (pos1-31),28);

//::std::cout << mid << ":" << uid << ::std::endl;
Pattern<CharString> pattern(uid);
if(find(finder, pattern)) {
        //std::cout << beginPosition(finder).i1 << std::endl;
 assignSeq(seq3, multiSeqFile3[beginPosition(finder).i1], format3);    // read sequence

Finder<String<char> > finder3(seq3);
Pattern<String<char>, MyersUkkonen> pattern3(seq2);
setScoreLimit(pattern3, -10);//Edit Distance = -1. -2 would be worse 
int mpos=-1;
int mscr=-100;
while (find(finder3, pattern3)) {
     int cscr=getScore(pattern3);
    if(cscr > mscr){mscr=cscr;
     mpos=position(finder3);
     }
}
     if(mpos >=0){
         mpos=mpos-length(seq2);
     ::std::cout << mid << " " << uid << " " << mpos << " " << mscr <<" "<< seq2 << ::std::endl;
}

}
clear(finder);

}
}
	std::cout << "Loading " << seqCount1 << " sequences took " << SEQAN_PROTIMEDIFF(loadTime);
	std::cout << " seconds." << std::endl << std::endl;
	return 0;
}

