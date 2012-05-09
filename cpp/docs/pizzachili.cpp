#include <iostream>
#include <seqan/index.h>

using namespace seqan;

template <typename TSpec>
void testPizzaChili() {
    typedef Index<String<char>, PizzaChili<TSpec> > index_t;
    index_t index_pc;
    indexText(index_pc) = "This is the best test with a bast jest.";

    ::std::cout << indexText(index_pc) << ::std::endl;

    Finder<index_t> finder(index_pc);
    while (find(finder, "est"))
        ::std::cout << "Hit at position " << position(finder) << ::std::endl;

    typename Fibre<index_t, PizzaChiliText>::Type text = indexText(index_pc);
    ::std::cout << "infix(text, 12, 21): " << infix(text, 12, 21) << ::std::endl;

    save(index_pc, "pizzachili");
    index_t index2;
    open(index2, "pizzachili");
    ::std::cout << indexText(index2) << ::std::endl;
}

int main() {
    ::std::cout << "Test the alphabet-friendly FM index:" << ::std::endl;
    testPizzaChili<PizzaChiliAF>();
    ::std::cout << ::std::endl << "Test the compressed compact suffix array index:" << ::std::endl;
    testPizzaChili<PizzaChiliCcsa>();
    return 0;
}
