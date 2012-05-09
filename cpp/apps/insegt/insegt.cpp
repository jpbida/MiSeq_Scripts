#include <fstream>
#include <iostream>
#include <sstream>

#define SEQAN_PROFILE
#ifndef RELEASE	
//#define SEQAN_DEBUG			
//#define SEQAN_TEST	
#endif

#include <string>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/map.h>
#include <seqan/refinement.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

#include "base.h"
#include "read_gff.h"
#include "read_gtf.h"
#include "create_gff.h"
#include "fusion.h"
#include "overlap_module.h"

using namespace seqan;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
////// main
///////////////////////////////////////////////////////////////////////////////

int main( int argc, const char *argv[] ) 
{
	CharString nameSAM;
	CharString nameGFF;
	CharString outputPath;
	unsigned nTuple;
	unsigned offsetInterval;
	unsigned thresholdGaps;
	unsigned thresholdCount;
	double thresholdRPKM = 0.0;
	bool maxTuple = 0;
	bool exact_nTuple = 0;
	bool unknownO = 0;
	bool fusion = 0;
	bool gtf = 0;
	
	CommandLineParser parser;
	
	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "**********************************************************************");
	addTitleLine(parser, "*** INSEGT                                                         ***");
	addTitleLine(parser, "*** INtersecting SEcond Generation sequencing daTa with annotation ***");
	addTitleLine(parser, "**********************************************************************");

	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("s", "sam", "Sam-file with aligned reads", (int)OptionType::String), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("g", "gff", "GFF_file with annotations", (int)OptionType::String), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("p", "outputPath", "path for output-files", (int)OptionType::String, ""), "<String>"));
	addOption(parser, addArgumentText(CommandLineOption("n", "nTuple", "nTuple", (int)OptionType::Int, 2), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("o", "offsetInterval", "offset to short alignment-intervals for search", (int)OptionType::Int, 5), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("t", "thresholdGaps", "threshold for allowed gaps in alignment (not introns)", (int)OptionType::Int, 5), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("c", "thresholdCount", "threshold for min. count of tuple for output", (int)OptionType::Int, 1), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("r", "thresholdRPKM", "threshold for min. RPKM of tuple for output", (int)OptionType::Double, 0.0), "<Double>"));
	addOption(parser, CommandLineOption("m", "maxTuple", "create only maxTuple", (int)OptionType::Boolean));
	addOption(parser, CommandLineOption("e", "exact_nTuple", "create only Tuple of exact length n", (int)OptionType::Boolean));
	addOption(parser, CommandLineOption("u", "unknown_orientation", "orientation of reads is unknown", (int)OptionType::Boolean));
	addOption(parser, CommandLineOption("f", "fusion_genes", "check for fusion genes and create separate output for matepair tuple", (int)OptionType::Boolean));
	addOption(parser, CommandLineOption("z", "gtf", "Gtf-format as input for annotation (instead of Gff-format)", (int)OptionType::Boolean));	

	if (argc == 1)
	{
		shortHelp(parser, cerr);	// print short help and exit
		return 0;
	}

	if (!parse(parser, argc, argv, ::std::cerr)) return 1;
	
	getOptionValueLong(parser, "sam", nameSAM);
	getOptionValueLong(parser, "gff", nameGFF);
	getOptionValueLong(parser, "outputPath", outputPath);
	getOptionValueLong(parser, "nTuple", nTuple);
	getOptionValueLong(parser, "offsetInterval", offsetInterval);
	getOptionValueLong(parser, "thresholdGaps", thresholdGaps);
	getOptionValueLong(parser, "thresholdCount", thresholdCount);
	getOptionValueLong(parser, "thresholdRPKM", thresholdRPKM);

	if (isSetLong(parser, "maxTuple")) maxTuple = 1;
	if (isSetLong(parser, "exact_nTuple")) exact_nTuple = 1;
	if (isSetLong(parser, "unknown_orientation")) unknownO = 1;
	if (isSetLong(parser, "fusion_genes")) fusion = 1;
	if (isSetLong(parser, "gtf")) gtf = 1;	

	if (maxTuple) 
	{
		nTuple = 0;		// sign for maxTuple
		exact_nTuple = 0;	// not both possible: maxTuple is prefered over exact_nTuple and n
	}

	
	ngsOverlapper(nameSAM, nameGFF, outputPath, nTuple, exact_nTuple, thresholdGaps, offsetInterval, thresholdCount, thresholdRPKM, unknownO, fusion, gtf);
	return 0;
}

