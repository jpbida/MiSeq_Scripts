//FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
//FRAGMENT(init1)
	Align< String<char> > ali;
	appendValue(rows(ali), "aphilologicaltheorem");
	appendValue(rows(ali), "bizarreamphibology");

//FRAGMENT(ali1)
	::std::cout << "Score = " << localAlignment(ali, Score<int>(3,-3,-2, -2), SmithWaterman()) << ::std::endl;
	::std::cout << ali;
	::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0))-1) << "]";
	::std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1))-1) << "]" << ::std::endl << ::std::endl;

//FRAGMENT(init2)
	Align< String<Dna> > ali2;
	appendValue(rows(ali2), "ataagcgtctcg");
	appendValue(rows(ali2), "tcatagagttgc");

//FRAGMENT(ali2)
	LocalAlignmentFinder<> finder(ali2);
	Score<int> scoring(2, -1, -2, 0);
	while (localAlignment(ali2, finder, scoring, 4, WatermanEggert())) {
		::std::cout << "Score = " << getScore(finder) << ::std::endl;
		::std::cout << ali2;
		::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali2, 0)) << ":" << (clippedEndPosition(row(ali2, 0))-1) << "]";
		::std::cout << " and Seq2[" << clippedBeginPosition(row(ali2, 1)) << ":" <<  (clippedEndPosition(row(ali2, 1))-1) << "]" << ::std::endl << ::std::endl;
	}
	return 0;
}
