/*=========================================================================
  Copyright (C) 2009 by Stephan Aiche

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================
  $Id: rep_sep.cpp 8349 2010-12-20 10:34:27Z reinert $
 ==========================================================================*/

#define SEQAN_PROFILE

#include <iostream>
#include <sstream>
#include <fstream>

// load command line parser
#include <seqan/misc/misc_cmdparser.h>

#include <seqan/store.h>

#include "utils.h"
#include "assembly_parser.h"
#include "column_scanner.h"
#include "rgraph.h"

using namespace seqan;
using namespace std;

void _makeParser(CommandLineParser & parser)
{
    // needed input file    
    //addOption(parser,CommandLineOption('a',"afg","filename of the afg assembly file",(OptionType::String | OptionType::Mandatory ) ));    
    
    addOption(parser, CommandLineOption("a", "assembly", "input assembly filename", (OptionType::String | OptionType::Mandatory ) ));
    addHelpLine(parser, "Currently supported input formats are afg (AMOS message format)"); 
    addOption(parser,CommandLineOption("c","contig","(zero based) number of the contig in the assembly that should be analyzed, default is 0",OptionType::Int ));
    addOption(parser,CommandLineOption("n","copy-number","number of compressed repeat copies in the given contig (default is 2)",OptionType::Int ));
    addHelpLine(parser, "");
    
    addOption(parser,CommandLineOption("","no-clean","do not clean the graph",OptionType::Boolean));
    addOption(parser,CommandLineOption("","dotfile","write constructed graph as dotfile to visualize in Graphviz",OptionType::Boolean));
    addHelpLine(parser, "");
    
    addOption(parser,CommandLineOption("p","output-prefix","filename prefix for the result files.",(OptionType::String | OptionType::Mandatory ) ));
    addHelpLine(parser, "Files for the ILP and the Result will be named PREFIX.(lp|rs|dot|heu)");
    // only needed for Gunnar
    //addHelpLine(parser, "");
    //addOption(parser,CommandLineOption(' ',"side-chain-format","reformulate the problem into the side-chain positioning problem",OptionType::Boolean));
    
    addSection(parser, "Column Detection Strategy");
    addOption(parser, CommandLineOption("d","dnp", "use DNP strategy for column detection [DEFAULT]",OptionType::Boolean));
    addOption(parser, CommandLineOption("s","simple-columns", "use the simple (normal distributed) column detection strategy",OptionType::Boolean));
    addOption(parser, CommandLineOption("e","error", "expected sequencing error (is only used in combination with simple strategy)", OptionType::Double));
    //addHelpLine(parser, "The default is DNP");

    addSection(parser, "Heuristic selection");
    addOption(parser, CommandLineOption("","hpcm", "try to solve the problem with the multi-component-expansion (mce) heuristic [DEFAULT]",OptionType::Boolean));
    addOption(parser, CommandLineOption("","hsce", "try to solve the problem with the single-component-expansion (sce) heuristic",OptionType::Boolean));
    //addHelpLine(parser, "The default is hpcm");

    addTitleLine(parser,"***********************************************");
    addTitleLine(parser,"***                 RepSep                  ***");
    addTitleLine(parser,"***          (c) Stephan Aiche 2009         ***");
    addTitleLine(parser,"***********************************************");

    addUsageLine(parser, "[OPTION]... --assembly <input file> --output-prefix <prefix>");
}

int main(int argc, const char *argv[])
{
    CommandLineParser parser;
    addVersionLine(parser, "RepSep version 0.1 20090723");
    _makeParser(parser);

    if (argc == 1)
    {
        shortHelp(parser, cerr);	// print short help and exit
        return 0;
    }

    bool stop = !parse(parser, argc, argv, cerr);

    if (isSetLong(parser, "help") || isSetLong(parser, "version") || stop) return 0;	// print help or version and exit

    // print some hello text
    version(parser); 

    // run starts
    SEQAN_PROTIMESTART(profileTime);

    typedef FragmentStore<> TFragmentStore;
    typedef Size<TFragmentStore>::Type TSize;
    TFragmentStore fragStore;

    CharString assembly_filename;
    getOptionValueLong(parser, "assembly", assembly_filename);
    
    cout << "loading data from " << assembly_filename << " .. " << flush;

    TSize numberOfContigs = 0;
    TSize numberOfReads = 0;

    FILE* strmReads = fopen(toCString(assembly_filename), "rb");
    read(strmReads, fragStore, Amos());
    fclose(strmReads);
    numberOfContigs = length(fragStore.contigStore);
    numberOfReads = length(fragStore.readStore);

    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;
    
    /*
    cout << endl << "------------------------------- " << endl;
    cout << "stats: " << endl;
    cout << "number of contigs: " << numberOfContigs << endl;
    cout << "number of reads:   " << numberOfReads << endl;
    cout << "------------------------------- " << endl;
    */
    int selected_contig = 0;
    getOptionValueLong(parser, "contig", selected_contig);

    if(! (static_cast<TSize>(selected_contig) < numberOfContigs) ) {
        cout << "You have selected an invalid contig! Only " << numberOfContigs << " different contig(s) were found" << endl;
        cout << "in the currently assembly" << endl;
        return 1;
    }

    cout << "#INFO: you have selected contig nr. " << selected_contig << endl;

    // construct matrix for parsing
    typedef char TAlpahbet;
    typedef Id<TFragmentStore::TAlignedReadStore>::Type TId;
    typedef Size<TFragmentStore::TReadSeq>::Type TReadPos;

    typedef Triple<TAlpahbet,TId,TReadPos> TMatrixValue;
    typedef String<TMatrixValue> TMatrix;

    TMatrix matrix;
    typedef Triple<char, TId, TReadPos> TColumnAlphabet;
    typedef String<TColumnAlphabet> TCandidateColumn;
    typedef Pair<TReadPos, TCandidateColumn> TAnnotatedCandidateColumn;
    String<TAnnotatedCandidateColumn> candidates;

    SEQAN_PROTIMEUPDATE(profileTime);
    cout << "parsing contig for candidate columns .. " << flush;
    
    // TODO: add code for selection of scanner
    parseContig(fragStore, selected_contig, candidates, SimpleColumn());
    
    for(TSize x = 0 ; x < length(candidates) ; ++x)
    {
      cout << candidates[x].i1 << endl;
    }
    
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;
    cout << "#INFO: <" << length(candidates) << "> possible candidate columns were identified" << endl; 

    
    ReadGraph<TColumnAlphabet, Value<TFragmentStore::TAlignedReadStore>::Type, TReadPos > rgraph;
    
    cout << "adding candidates to graph .. " << flush;
    construct(rgraph, candidates, fragStore, selected_contig);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    cout << "adding mate pairs to graph .. " << flush;
    add_mate_pairs(rgraph, fragStore, selected_contig);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;   

    GraphScoring scoring_scheme;
    
    cout << "compute scores for graph edges .. " << flush;
    scoreGraph_(rgraph, fragStore, static_cast<TSize>( selected_contig ), scoring_scheme);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;   

    if(!isSetLong(parser, "no-clean") )
    {
        cout << "cleaning graph .. " << flush;
        // TODO: cleaning needs to be implemented
        cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;   
    }
    else 
    {
        cout << "#INFO: skipped graph cleaning" << endl;
    }

    cout << "analyzing constructed graph .. " << flush;
    bool hasMultiComp = hasMultipleComponents(rgraph);
    cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;

    if(hasMultiComp) {
        cout << "This graph has multiple components!!" << endl;
    }

    // TODO: add some statistics for the graph
    cout << endl << "################ graph info ################" << endl;
    
    cout << " |edges| " << numEdges(rgraph) << endl;
    cout << " |vertices| " << numVertices(rgraph) << endl;

    cout << "############################################" << endl << endl << flush;

    // TODO: implement heurisitics 
    typedef Graph_< ReadGraph<TColumnAlphabet, Value<TFragmentStore::TAlignedReadStore>::Type, TReadPos > >::Type TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    typedef String<TVertexDescriptor> TComponent;
    typedef String<TComponent> TComponentList;

    TComponentList components;

    int copy_number = 0;
    getOptionValueLong(parser, "copy-number", copy_number);

    if(isSetLong(parser, "hpcm"))
    {
      cout << "separate the graph using (GuidedParallelComponentMerge) .. " << flush;
      solve(rgraph,components,copy_number,GuidedParallelComponentMerge());
      cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;       
    }
    else 
    {
      cout << "separate the graph using (SingleComponentExpansion) .. " << flush;
      solve(rgraph,components,copy_number,SingleComponentExpansion());
      cout << "done (" << SEQAN_PROTIMEUPDATE(profileTime) << " seconds)" << endl;        
    }
    
    for(TSize c = 0 ; c < length(components) ; ++c)
    {
      cout << "Component <" << c << "> contains " << length(components[c]) << " reads " << endl;
      for(TSize r = 0 ; r < length(components[c]) ; ++r)
      {
        cout << value(rgraph.vertexCargo,components[c][r]).alignedRead.readId << " (" << fragStore.readNameStore[value(rgraph.vertexCargo,value(components[c],r)).alignedRead.readId] << ")" << endl;
      }
    }
    // TODO: implement IO stuff
}
