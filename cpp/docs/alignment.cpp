#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_align.h>

using namespace seqan;

int main()
{
	typedef String<Dna> TSequence;
	TSequence seq1 = "atcgaatgcgga";
	TSequence seq2 = "actcgttgca";
	Score<int> score(0, -1, -1, -2);
	Align<TSequence, ArrayGaps> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);

	::std::cout << "Score = " << globalAlignment(align, score) << ::std::endl;
	::std::cout << align << ::std::endl;
	::std::cout << "Score = " << globalAlignment(align, score, MyersHirschberg()) << ::std::endl;
	::std::cout << align << ::std::endl;
	typedef StringSet<TSequence, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

	TStringSet string_set;
	appendValue(string_set, seq1);
	appendValue(string_set, seq2);
	TAlignmentGraph alignment_graph(string_set);

	::std::cout << "Score = " << globalAlignment(alignment_graph, score, Gotoh()) << ::std::endl;
	::std::cout << alignment_graph << ::std::endl;
	return 0;
}
