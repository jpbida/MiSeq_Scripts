// FRAGMENT(includes)
#define SEQAN_PROFILE // enable time measurements
#include <seqan/file.h>
#include <iostream>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/basic.h>
using namespace seqan;

typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
typedef Value<TReadSeqStore>::Type TReadSeq;
typedef FragmentStore<>::TContigStore TContigStore;
typedef Value<TContigStore>::Type TContigStoreElement;
typedef TContigStoreElement::TContigSeq TContigSeq;
typedef Index<TReadSeqStore, IndexQGram<Shape<Dna, UngappedShape<11> >, OpenAddressing> > TIndex;
typedef Pattern<TIndex, Swift<SwiftSemiGlobal> > TPattern;
typedef Finder<TContigSeq, Swift<SwiftSemiGlobal> > TFinder;
typedef FragmentStore<>::TAlignedReadStore TAlignedReadStore;
typedef Value<TAlignedReadStore>::Type TAlignedRead;
typedef FragmentStore<>::TContigStore TContigStore;
const double EPSILON = 0.1;
typedef Iterator<String<char> >::Type TIterator;
// FRAGMENT(open_file)
int main (int argc, char const * argv[])
{
	SEQAN_PROTIMESTART(loadTime);

FragmentStore<> fragStore;
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
	StringSet<String<Dna5Q> > seqs1;
//	StringSet<CharString> seqIDs;

	reserve(seqs1, seqCount1, Exact());
//	reserve(seqIDs, seqCount, Exact());
	
    unsigned seqCount2 = length(multiSeqFile2);
	StringSet<String<Dna5Q> > seqs2;
//	StringSet<CharString> seqIDs;
	reserve(seqs2, seqCount2, Exact());
	
    unsigned seqCount3 = length(multiSeqFile3);

//	reserve(seqIDs, seqCount, Exact());

// FRAGMENT(read_sequences)
	String<char> seq1;
	String<char> seq2;
	String<char> seq3;
	CharString qual1;
	CharString qual2;
	CharString id1;
	CharString id2;

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
	for (unsigned j = 0; j < seqCount3; ++j)
	{
        assignSeq(seq3,multiSeqFile3[j], format3);
Finder<String<char> > finder(seq3);
Pattern<String<char>, WildShiftAnd> pattern(uid);
if(find(finder, pattern)) {
//::std::cout << mid << ":" << seq2 << ::std::endl;
Finder<String<char> > finder3(seq3);
Pattern<String<char>, MyersUkkonen> pattern3(seq2);
setScoreLimit(pattern3, -4);//Edit Distance = -1. -2 would be worse 
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
     ::std::cout << mid << " " << uid << " " << mpos << " " << mscr <<" "<< seq2 << "\n";
}

}
    }


//Get Library Sequence with UID 


//::std::FILE * fl = ::std::fopen("testfile.fa", "ab");
//write(fl, seq2, "frag", Fasta());
//close (fl);

/*
String<char> ndl("CCT");

Finder<String<char> > fnd(haystk);
Pattern<String<char>, MyersUkkonen> pat(ndl);

setScoreLimit(pat, -1);//Edit Distance = -1. -2 would be worse 
while (find(fnd, pat)) {
     ::std::cout << position(fnd) << ": " << getScore(pat) << "\n";
}

*/
}

//Get Library Sequence

//Align seq2 to library

//If alignment quality is > threshold 


		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
//		for (unsigned j = 0; j < length(qual1) && j < length(seq1); ++j)
//			assignQualityValue(seq1[j], (int)(ordValue(qual1[j]) - 33));

		// we use reserve and append, as assign is not supported
		// by StringSet<..., Owner<ConcatDirect<> > >
//		appendValue(seqs2, seq2, Generous());
//		appendValue(seqIDs, id, Generous());
	}

// FRAGMENT(output)
	std::cout << "Loading " << seqCount1 << " sequences took " << SEQAN_PROTIMEDIFF(loadTime);
	std::cout << " seconds." << std::endl << std::endl;

	return 0;
}

