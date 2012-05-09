// FRAGMENT(includes)
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;
// FRAGMENT(showAllLetterOfMyAlphabet)
template <typename TAlphabet>
void showAllLetterOfMyAlphabet(TAlphabet const &)
{
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	for (TSize i = 0; i < alphSize; ++i)
		std::cout << i << ',' << TAlphabet(i) << "  ";
	std::cout << std::endl;
}
// FRAGMENT(main)
int main()
{
	showAllLetterOfMyAlphabet(AminoAcid());
	showAllLetterOfMyAlphabet(Dna());
	showAllLetterOfMyAlphabet(Dna5());
	return 0;
}
