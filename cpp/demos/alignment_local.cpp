///A tutorial about local alignments.
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
///Example 1: This program applies the Smith-Waterman algorithm to compute the best local alignment between two given sequences.
	Align< String<char> > ali;
	appendValue(rows(ali), "aphilologicaltheorem");
	appendValue(rows(ali), "bizarreamphibology");
    ::std::cout << "Score = " << localAlignment(ali, Score<int>(3,-3,-2, -2), SmithWaterman()) << ::std::endl;
	::std::cout << ali;
	::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0))-1) << "]";
	::std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1))-1) << "]" << ::std::endl << ::std::endl;


///Example 2: This program applies the Waterman-Eggert algorithm to compute all non-overlapping local alignments with score better or equal 2.
	Align< String<Dna> > ali2;
	appendValue(rows(ali2), "ataagcgtctcg");
	appendValue(rows(ali2), "tcatagagttgc");

	LocalAlignmentFinder<> finder(ali2);
	Score<int> scoring(2, -1, -2, 0);
	while (localAlignment(ali2, finder, scoring, 2, WatermanEggert())) {
		::std::cout << "Score = " << getScore(finder) << ::std::endl;
		::std::cout << ali2;
		::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali2, 0)) << ":" << (clippedEndPosition(row(ali2, 0))-1) << "]";
		::std::cout << " and Seq2[" << clippedBeginPosition(row(ali2, 1)) << ":" <<  (clippedEndPosition(row(ali2, 1))-1) << "]" << ::std::endl << ::std::endl;
	}

///Example 3
    Align< String<Dna> > ali3;
	appendValue(rows(ali3), "cccccc");
	appendValue(rows(ali3), "tttttggccccccgg");
	LocalAlignmentFinder<> finder3(ali3);
	Score<int> scoring3(1, -1, -1, -1);
    while (localAlignment(ali3, finder3, scoring3, 5, WatermanEggert())) {
        ::std::cout << "Score = " << getScore(finder3) << ::std::endl;
	    ::std::cout << ali3;
	    ::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali3, 0)) << ":" << (clippedEndPosition(row(ali3, 0))-1) << "]";
	    ::std::cout << " and Seq2[" << clippedBeginPosition(row(ali3, 1)) << ":" <<  (clippedEndPosition(row(ali3, 1))-1) << "]" << ::std::endl << ::std::endl;
    }

///Example 4: This program applies the banded Waterman-Eggert algorithm to compute all non-overlapping local alignments with score or equal 5
///           in the band from diagonal -1 to diagonal 8.
    Align< String<Dna5> > ali4;
    appendValue(rows(ali4), "AAAAAAANAAAGGGNGGGGGGGGNGGGGGANAA");
    appendValue(rows(ali4), "GGGGGGCGGGGGGGA");

    LocalAlignmentFinder<> finder4(ali4);
    Score<int> scoring4(1, -1, -1, -1);
    while (localAlignment(ali4, finder4, scoring4, 5, -1, 8, BandedWatermanEggert())) {
        ::std::cout << "Score = " << getScore(finder4) << ::std::endl;
        ::std::cout << ali4;
        ::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali4, 0)) << ":" << (clippedEndPosition(row(ali4, 0))-1) << "]";
        ::std::cout << " and Seq2[" << clippedBeginPosition(row(ali4, 1)) << ":" <<  (clippedEndPosition(row(ali4, 1))-1) << "]" << ::std::endl << ::std::endl;
    }

	return 0;
}
