///A tutorial about the counts fibre of the q-gram index.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
///First, we create a @Class.StringSet@ of 4 @Class.String.Strings@.
	StringSet< String<char> > mySet;
	resize(mySet, 4);
	mySet[0] = "tobeornottobe";
	mySet[1] = "thebeeonthecomb";
	mySet[2] = "helloworld";
	mySet[3] = "beingjohnmalkovich";

///Then we create an @Class.Index@ of our @Class.StringSet@ and
///a @Class.Finder@ of the @Class.Index@.
	typedef Index< StringSet<String<char> >, IndexQGram<UngappedShape<2> > > TIndex;
	typedef Infix<Fibre<TIndex, QGramCounts>::Type const>::Type TCounts;

  TIndex myIndex(mySet);

///Finally we output how often $"be"$ occurs in each sequence.
  hash(indexShape(myIndex), "be");
	TCounts cnts = countOccurrencesMultiple(myIndex, indexShape(myIndex));
	for (unsigned i = 0; i < length(cnts); ++i)
    std::cout << cnts[i].i2 << " occurrences in sequence " << cnts[i].i1  << std::endl;

	return 0;
}
