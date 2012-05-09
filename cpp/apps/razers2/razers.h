 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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
 ==========================================================================*/

#ifndef SEQAN_HEADER_RAZERS_H
#define SEQAN_HEADER_RAZERS_H

#include <iostream>
#include <fstream>

#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/store.h>

#ifdef RAZERS_PARALLEL
#include "tbb/spin_mutex.h"
#endif

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options

	template < bool DONT_VERIFY_ = false, bool DONT_DUMP_RESULTS_ = false >
	struct RazerSSpec 
	{
		enum { DONT_VERIFY = DONT_VERIFY_ };				// omit verifying potential matches
		enum { DONT_DUMP_RESULTS = DONT_DUMP_RESULTS_ };	// omit dumping results
	};

	template < typename TSpec = RazerSSpec<> >
	struct RazerSOptions
	{
	// main options
		TSpec		spec;
		bool		forward;			// compute forward oriented read matches
		bool		reverse;			// compute reverse oriented read matches
		double		errorRate;			// Criteria 1 threshold
		unsigned	maxHits;			// output at most maxHits many matches
		unsigned	distanceRange;		// output only the best, second best, ..., distanceRange best matches
		bool		purgeAmbiguous;		// true..remove reads with more than maxHits best matches, false..keep them
		CharString	output;				// name of result file
		int			_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
		bool		hammingOnly;		// no indels
		int			trimLength;			// if >0, cut reads to #trimLength characters
		
	// output format options
		unsigned	outputFormat;		// 0..Razer format
										// 1..enhanced Fasta
										// 2..ELAND format
		bool		dumpAlignment;		// compute and dump the match alignments in the result files
		unsigned	genomeNaming;		// 0..use Fasta id
										// 1..enumerate reads beginning with 1
		unsigned	readNaming;			// 0..use Fasta id
										// 1..enumerate reads beginning with 1
										// 2..use the read sequence (only for short reads!)
										// 3..use Fasta id, do not append /L and /R for mate pairs.
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		unsigned	positionFormat;		// 0..gap space
										// 1..position space
		const char	*runID;				// runID needed for gff output

	// filtration parameters
		::std::string shape;			// shape (e.g. 11111111111)
		int			threshold;			// threshold
		int			tabooLength;		// taboo length
		int			repeatLength;		// repeat length threshold
		double		abundanceCut;		// abundance threshold

	// mate-pair parameters
		int			libraryLength;		// offset between two mates
		int			libraryError;		// offset tolerance
		unsigned	nextPairMatchId;	// use this id for the next mate-pair

	// verification parameters
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];

	// statistics
		__int64		countFiltration;	// matches returned by the filter
		__int64		countVerification;	// matches returned by the verifier
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results
		
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		bool		maqMapping;
		int		maxMismatchQualSum;
		int		mutationRateQual;
		int		absMaxQualSumErrors;
		unsigned	artSeedLength;
		bool		noBelowIdentity;
#endif

#ifdef RAZERS_MICRO_RNA
		bool		microRNA;
		unsigned	rnaSeedLength;
		bool 		exactSeed;
#endif			

		bool		lowMemory;		// set maximum shape weight to 13 to limit size of q-gram index
		bool		fastaIdQual;		// hidden option for special fasta+quality format we use


	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh

	// multi-threading

#ifdef RAZERS_PARALLEL
		typedef ::tbb::spin_mutex	TMutex;

		TMutex		*patternMutex;
		TMutex		optionsMutex;
		TMutex		matchMutex;
#endif

		RazerSOptions() 
		{
			forward = true;
			reverse = true;
			errorRate = 0.08;
			maxHits = 100;
			distanceRange = 0;	// disabled
			purgeAmbiguous = false;
			output = "";
			_debugLevel = 0;
			printVersion = false;
			hammingOnly = false;
			trimLength = 0;
			
			outputFormat = 0;
			dumpAlignment = false;
			genomeNaming = 0;
			readNaming = 0;
			sortOrder = 0;
			positionFormat = 0;
			runID = "s"; 	//

			matchN = false;

			shape = "11111111111";
			threshold = 1;
			tabooLength = 1;
			repeatLength = 1000;
			abundanceCut = 1;

			libraryLength = 220;
			libraryError = 50;
			nextPairMatchId = 0;
			
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;

			compactThresh = 1024;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			maqMapping = false;
			maxMismatchQualSum = 70;
			mutationRateQual = 30;
			artSeedLength = 28;	// the "artificial" seed length that is used for mapping quality assignment 
						// (28bp is maq default)
			absMaxQualSumErrors = 100; // maximum for sum of mism qualities in total readlength
			noBelowIdentity = false;
#endif

#ifdef RAZERS_MICRO_RNA
			microRNA = false;
			rnaSeedLength = 16;
			exactSeed = true;
#endif			

			lowMemory = false;		// set maximum shape weight to 13 to limit size of q-gram index
			fastaIdQual = false;

		}
	};

struct MicroRNA{};	

#ifdef RAZERS_MICRO_RNA
#define RAZERS_EXTENDED_MATCH
#endif

#ifdef RAZERS_DIRECT_MAQ_MAPPING 
#define RAZERS_EXTENDED_MATCH
#endif
	
	
//////////////////////////////////////////////////////////////////////////////
// Typedefs
/*
	// definition of a Read match
	template <typename TGPos_>
	struct ReadMatch 
	{
		typedef typename MakeSigned_<TGPos_>::Type TGPos;

		unsigned		gseqNo;			// genome seqNo
		unsigned		rseqNo;			// read seqNo
		TGPos			beginPos;			// begin position of the match in the genome
		TGPos			gEnd;			// end position of the match in the genome
#ifdef RAZERS_MATEPAIRS
		unsigned		pairId;			// unique id for the two mate-pair matches (0 if unpaired)
		int				mateDelta:24;	// outer coordinate delta to the other mate 
		int				pairScore:8;	// combined score of both mates
#endif
		unsigned short	editDist;		// Levenshtein distance
#ifdef RAZERS_EXTENDED_MATCH
		short	 		mScore;
		short			seedEditDist;
#endif
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};
*/	
	enum RAZERS_ERROR 
	{
		RAZERS_INVALID_OPTIONS = -1,
		RAZERS_READS_FAILED    = -2,
		RAZERS_GENOME_FAILED   = -3,
		RAZERS_INVALID_SHAPE   = -4
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef Dna5String									TGenome;
	typedef StringSet<TGenome>							TGenomeSet;
//	typedef Dna5String									TRead;
	typedef String<Dna5Q>								TRead;
/*#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead>							TReadSet;
#endif
*/
/*	typedef ReadMatch<Difference<TGenome>::Type>		TMatch;		// a single match
	typedef String<TMatch>								TMatches;	// array of matches
*/

	template <typename TReadSet, typename TShape, typename TSpec>
	struct Cargo< Index<TReadSet, IndexQGram<TShape, TSpec> > > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};

//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, IndexQGram<TShape, TSpec> > > 
	{
		typedef Pair<
			unsigned,				
			unsigned,
			BitCompressed<22, 10>	// max. 4M reads of length < 1024
		> Type;
	};
	
#else

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, IndexQGram<TShape, TSpec> > > 
	{
		typedef Pair<
			unsigned,			// many reads
			unsigned,			// of arbitrary length
			Compressed
		> Type;
	};

#endif

	template <>
	struct Size<Dna5String>
	{
		typedef unsigned Type;
	};

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, IndexQGram<TShape> > >
	{
		typedef unsigned Type;
	};

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, IndexQGram<TShape, OpenAddressing> > >
	{
		typedef unsigned Type;
	};
	

#ifdef RAZERS_PRUNE_QGRAM_INDEX

	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TReadSet, typename TShape, typename TSpec>
	inline bool _qgramDisableBuckets(Index<TReadSet, IndexQGram<TShape, TSpec> > &index) 
	{
		typedef Index<TReadSet, IndexQGram<TShape>	>		TReadIndex;
		typedef typename Fibre<TReadIndex, QGramDir>::Type	TDir;
		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename Value<TDir>::Type					TSize;

		TDir &dir    = indexDir(index);
		bool result  = false;
		unsigned counter = 0;
		TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
		if (thresh < 100) thresh = 100;

		TDirIterator it = begin(dir, Standard());
		TDirIterator itEnd = end(dir, Standard());
		for (; it != itEnd; ++it)
			if (*it > thresh) 
			{
				*it = (TSize)-1;
				result = true;
				++counter;
			}

		if (counter > 0 && cargo(index)._debugLevel >= 1)
			::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

		return result;
	}

#endif


	template <
		typename TFragmentStore_, 
		typename TRazerSOptions_,
		typename TPreprocessing_,
		typename TSwiftPattern_,
		typename TCounts_
	>
	struct MatchVerifier
	{
		typedef TFragmentStore_									TFragmentStore;
		typedef TRazerSOptions_									TOptions;
		typedef TPreprocessing_									TPreprocessing;
		typedef TSwiftPattern_									TSwiftPattern;
		typedef TCounts_										TCounts;
		
		typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
		typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
		typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
		typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
		typedef typename Size<TGenome>::Type					TSize;

		TFragmentStore	&store;
		TOptions		&options;			// RazerS options
		TPreprocessing	&preprocessing;
		TSwiftPattern	&swiftPattern;
		TCounts			&cnts;
		
		TAlignedRead	m;
		TAlignQuality	q;
		bool			onReverseComplement;
		TSize			genomeLength;
		bool			oneMatchPerBucket;
		
		MatchVerifier(TFragmentStore_ &_store, TOptions &_options, TPreprocessing &_preprocessing, TSwiftPattern &_swiftPattern, TCounts &_cnts):
			store(_store),
			options(_options),
			preprocessing(_preprocessing),
			swiftPattern(_swiftPattern),
			cnts(_cnts)
		{
			onReverseComplement = false;
			genomeLength = 0;
			oneMatchPerBucket = false;
		}
		
		inline void push()
		{
			if (onReverseComplement) 
			{
				// transform coordinates to the forward strand
				m.beginPos = genomeLength - m.beginPos;
				m.endPos = genomeLength - m.endPos;
			}
			if (!options.spec.DONT_DUMP_RESULTS)
			{
				m.id = length(store.alignedReadStore);
				appendValue(store.alignedReadStore, m, Generous());
				appendValue(store.alignQualityStore, q, Generous());
				if (length(store.alignedReadStore) > options.compactThresh)
				{
					typename Size<TAlignedReadStore>::Type oldSize = length(store.alignedReadStore);
#ifndef RAZERS_DONTMASKDUPLICATES
					maskDuplicates(store);	// overlapping parallelograms cause duplicates
#endif
#ifdef RAZERS_DIRECT_MAQ_MAPPING
					if (options.maqMapping)
						compactMatches(store, cnts, options, false, swiftPattern, true);
					else	
#endif
						compactMatches(store, cnts, options, false, swiftPattern);
					if (length(store.alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
						options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
					if (options._debugLevel >= 2)
						::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
				}
			}
			++options.countVerification;
		}
	};
 


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences with or w/o quality values
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
	FragmentStore<TFSSpec, TFSConfig> &store,
	const char *fileName, 
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA) countN = false;
#endif

	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);

	unsigned seqCount = length(multiFasta);

	String<Dna5Q>	seq;
	CharString		qual;
	CharString		id;
	
	unsigned kickoutcount = 0;
	unsigned maxReadLength = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0 || options.readNaming == 3
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			|| options.fastaIdQual
#endif
			)
			assignSeqId(id, multiFasta[i], format);	// read Fasta id
		assignSeq(seq, multiFasta[i], format);					// read Read sequence
		assignQual(qual, multiFasta[i], format);				// read ascii quality values  
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		if(options.fastaIdQual)
		{
			qual = suffix(id, length(id) - length(seq));
			if (options.readNaming == 0 || options.readNaming == 3)
				id = prefix(id,length(id) - length(seq));
			else 
				clear(id);
		}
#endif
		if (countN)
		{
			int count = 0;
			int cutoffCount = (int)(options.errorRate * length(seq));
			for (unsigned j = 0; j < length(seq); ++j)
				if (getValue(seq, j) == 'N')
					if (++count > cutoffCount)
					{
						clear(seq);
						clear(id);
						++kickoutcount;
						break;
					}
// low qual. reads are empty to output them and their id later as LQ reads
//			if (count > cutoffCount) continue;
		}

		// store dna and quality together
		assignQualities(seq, qual); 
		if (options.trimLength > 0 && length(seq) > (unsigned)options.trimLength)
			resize(seq, options.trimLength);

		appendRead(store, seq, id);
		if (maxReadLength < length(seq))
			maxReadLength = length(seq);
	}
	// memory optimization
	reserve(store.readSeqStore.concat, length(store.readSeqStore.concat), Exact());
//	reserve(store.readNameStore.concat, length(store.readNameStore.concat), Exact());

	typedef Shape<Dna, SimpleShape> TShape;
	typedef typename SAValue< Index<StringSet<TRead>, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
	TSAValue sa(0, 0);
	sa.i1 = ~sa.i1;
	sa.i2 = ~sa.i2;
	
	if ((unsigned)sa.i1 < length(store.readSeqStore) - 1)
	{
		::std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}
	if ((unsigned)sa.i2 < maxReadLength - 1)
	{
		::std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality reads.\n";
	return (seqCount > 0);
}

//////////////////////////////////////////////////////////////////////////////
// Read the first sequence of a multi-sequence file
// and return its length
inline int estimateReadLength(char const *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY))	// open the whole file
		return RAZERS_READS_FAILED;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);					// guess file format
	split(multiFasta, format);								// divide into single sequences

	if (length(multiFasta) == 0)
		return 0;

	Dna5String firstRead;
	assignSeq(firstRead, multiFasta[0], format);			// read the first sequence
	return length(firstRead);
}

#ifdef RAZERS_MICRO_RNA
	template <typename TReadMatch>
	struct LessRNoGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// genome position and orientation
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;
			if (a.beginPos < b.beginPos) return true;
			if (a.beginPos > b.beginPos) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;


			if(a.editDist < b.editDist) return true;
			if(a.editDist > b.editDist) return false;
			return (a.mScore > b.mScore);
		}
	};

	template <typename TReadMatch>
	struct LessRNoEdistHLen : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			if(a.editDist < b.editDist) return true;
			if(a.editDist > b.editDist) return false;
			return (a.mScore > b.mScore);

		}
	};
#else
	
	
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessRNoGPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TAlignedReadQualityStore &qualStore;
		
		LessRNoGPos(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (TAlignedRead const &a, TAlignedRead const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// beginning position
			typename TAlignedRead::TPos ba = _min(a.beginPos, a.endPos);
			typename TAlignedRead::TPos bb = _min(b.beginPos, b.endPos);
			if (ba < bb) return true;
			if (ba > bb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			// qualities
			if (a.id == TAlignedRead::INVALID_ID) return false;
			if (b.id == TAlignedRead::INVALID_ID) return true;
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			if (qa.score > qb.score) return true;
			if (qb.score > qa.score) return false;
			
			// prefer reads that support more of the reference
			return _max(a.beginPos, a.endPos) > _max(b.beginPos, b.endPos);
		}
	};

	// ... to sort matches and remove duplicates with equal gEnd
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessRNoGEndPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TAlignedReadQualityStore &qualStore;
		
		LessRNoGEndPos(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// end position
			typename TAlignedRead::TPos ea = _max(a.beginPos, a.endPos);
			typename TAlignedRead::TPos eb = _max(b.beginPos, b.endPos);
			if (ea < eb) return true;
			if (ea > eb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			// qualities
			if (a.id == TAlignedRead::INVALID_ID) return false;
			if (b.id == TAlignedRead::INVALID_ID) return true;
			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			if (qa.score > qb.score) return true;
			if (qb.score > qa.score) return false;
			
			// prefer reads that support more of the reference
			return _min(a.beginPos, a.endPos) < _min(b.beginPos, b.endPos);
		}
	};

#endif
		
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore>
	struct LessErrors : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		TAlignedReadQualityStore &qualStore;
		
		LessErrors(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// quality
			if (a.id == TAlignedRead::INVALID_ID) return false;
			if (b.id == TAlignedRead::INVALID_ID) return true;
			return qualStore[a.id].errors < qualStore[b.id].errors;
		}
	};
	
#ifdef RAZERS_DIRECT_MAQ_MAPPING

	struct QualityBasedScoring{};

	template <typename TReadMatch>
	struct LessRNoMQ : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;
			
			// quality
			if (a.mScore < b.mScore) return true; // sum of quality values of mismatches (the smaller the better)
			if (a.mScore > b.mScore) return false;
			
			return (a.editDist < b.editDist); // seedEditDist?
			// genome position and orientation
	/*		if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;
			if (a.beginPos < b.beginPos) return true;
			if (a.beginPos > b.beginPos) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;
	*/		
		}
	};
#endif

//////////////////////////////////////////////////////////////////////////////
// Mark duplicate matches for deletion
template <typename TFragmentStore>
void maskDuplicates(TFragmentStore &store)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

	sortAlignedReads(store.alignedReadStore, LessRNoGEndPos<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));

	TContigPos	beginPos = -1;
	TContigPos	endPos = -1;
	unsigned	contigId = TAlignedRead::INVALID_ID;
	unsigned	readId = TAlignedRead::INVALID_ID;
	bool		orientation = false;

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).pairMatchId != TAlignedRead::INVALID_ID) continue;
		TContigPos itEndPos = _max((*it).beginPos, (*it).endPos);
		if (endPos == itEndPos && orientation == ((*it).beginPos < (*it).endPos) &&
			contigId == (*it).contigId && readId == (*it).readId) 
		{
			(*it).id = TAlignedRead::INVALID_ID;	// mark this alignment for deletion
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		endPos = itEndPos;
		orientation = (*it).beginPos < (*it).endPos;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

	sortAlignedReads(store.alignedReadStore, LessRNoGPos<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));

	contigId = TAlignedRead::INVALID_ID;
	it = begin(store.alignedReadStore, Standard());
	itEnd = end(store.alignedReadStore, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID || (*it).pairMatchId != TAlignedRead::INVALID_ID) continue;

		TContigPos itBeginPos = _min((*it).beginPos, (*it).endPos);
		if (beginPos == itBeginPos && readId == (*it).readId &&
			contigId == (*it).contigId && orientation == ((*it).beginPos < (*it).endPos))
		{
			(*it).id = TAlignedRead::INVALID_ID;	// mark this alignment for deletion
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		beginPos = itBeginPos;
		orientation = (*it).beginPos < (*it).endPos;
	}

	sortAlignedReads(store.alignedReadStore, LessErrors<TAlignedReadStore, TAlignQualityStore>(store.alignQualityStore));
}

//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template <typename TFragmentStore, typename TCounts>
void countMatches(TFragmentStore &store, TCounts &cnt)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TCounts>::Type							TRow;
	typedef typename Value<TRow>::Type								TValue;
	
	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	
	unsigned readId = TAlignedRead::INVALID_ID;
	short errors = -1;
	__int64 count = 0;
	__int64 maxVal = MaxValue<TValue>::VALUE;

	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID) continue;
		if (readId == (*it).readId && errors == store.alignQualityStore[(*it).id].errors)
			++count;
		else
		{
			if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < length(cnt))
				cnt[errors][readId] = (maxVal < count)? (TValue)maxVal : (TValue)count;
			readId = (*it).readId;
			errors = store.alignQualityStore[(*it).id].errors;
			count = 1;
		}
	}
	if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < length(cnt))
		cnt[errors][readId] = (TValue)count;
}

//////////////////////////////////////////////////////////////////////////////

template < typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(Nothing &, TReadNo, TMaxErrors)
{
}

template < typename TSwift, typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(TSwift &swift, TReadNo readNo, TMaxErrors maxErrors)
{
	int minT = _qgramLemma(swift, readNo, maxErrors);
	if (minT > 1)
	{
		if (maxErrors < 0) minT = MaxValue<int>::VALUE;
//		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
		setMinThreshold(swift, readNo, (unsigned)minT);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TFragmentStore, typename TCounts, typename TSpec, typename TSwift >
void compactMatches(TFragmentStore &store, TCounts & 
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		cnts
#endif
	, RazerSOptions<TSpec> &options
	, bool compactFinal,
	TSwift &
#if defined RAZERS_DIRECT_MAQ_MAPPING || defined RAZERS_MASK_READS
		swift
#endif
	)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;

#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) compactMatches(store, cnts,options,compactFinal,swift,true);
#endif
	
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA)
		::std::sort(
			begin(matches, Standard()),
			end(matches, Standard()), 
			LessRNoEdistHLen<TMatch>());
	int bestMScore = 0;
#endif
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
#ifdef RAZERS_MICRO_RNA
	if(options.microRNA && options.purgeAmbiguous)
		++hitCountCutOff;	// we keep one more match than we actually want, so we can later decide
							// whether the read mapped more than maxhits times 
#endif
	unsigned errorsCutOff = MaxValue<unsigned>::VALUE;

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID) continue;
		unsigned errors = store.alignQualityStore[(*it).id].errors;
		if (readNo == (*it).readId && (*it).pairMatchId == TAlignedRead::INVALID_ID)
		{ 
			if (errors >= errorsCutOff) continue;
#ifdef RAZERS_MICRO_RNA
			if ((*it).mScore < bestMScore) continue;
#endif
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = errors - 1;
					if (options.purgeAmbiguous && (options.distanceRange == 0 || errors < options.distanceRange))
						maxErrors = -1;

					setMaxErrors(swift, readNo, maxErrors);

					if (maxErrors == -1 && options._debugLevel >= 2)
						::std::cerr << "(read #" << readNo << " disabled)";

					if (options.purgeAmbiguous)
					{
						if (options.distanceRange == 0 || errors < options.distanceRange || compactFinal)
							dit = ditBeg;
						else {
							*dit = *it;
							++dit;
						}

					}
				}
#endif
				continue;
			}
		}
		else
		{
			readNo = (*it).readId;
			hitCount = 0;
			if (options.distanceRange > 0)
				errorsCutOff = errors + options.distanceRange;
#ifdef RAZERS_MICRO_RNA
			bestMScore = (*it).mScore;
#endif
			ditBeg = dit;
		}
		*dit = *it;
		++dit;
	}
	resize(store.alignedReadStore, dit - begin(store.alignedReadStore, Standard()));
	compactAlignedReads(store);
}


#ifdef RAZERS_DIRECT_MAQ_MAPPING
//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec, typename TSwift >
void compactMatches(TMatches &matches, TCounts &cnts, RazerSOptions<TSpec> &, bool, TSwift &, bool dontCountFirstTwo)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoMQ<TMatch>());
	
	unsigned readNo = -1;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;

	//number of errors may not exceed 31!
	bool second = true;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).readId)
		{
			//second best match
			if (second)
			{
				second = false;
				if((cnts[1][(*it).readId] & 31)  > (*it).editDist)
				{
					//this second best match is better than any second best match before
					cnts[1][(*it).readId] = (*it).editDist; // second best dist is this editDist
										// count is 0 (will be updated if countFirstTwo)
				}
				if(!dontCountFirstTwo) 
					if((cnts[1][(*it).readId]>>5) != 2047) cnts[1][(*it).readId] += 32;
			}
			else
			{
				if ((*it).editDist <= (cnts[0][(*it).readId] & 31))
					if(cnts[0][(*it).readId]>>5 != 2047)
						cnts[0][(*it).readId] +=32;
				if ((*it).editDist <= (cnts[1][(*it).readId] & 31))
					if((cnts[1][(*it).readId]>>5) != 2047)
						cnts[1][(*it).readId] +=32;
				continue;
			}
		} else
		{	//best match
			second = true;
			readNo = (*it).readId;
			//cnts has 16bits, 11:5 for count:dist
			if((cnts[0][(*it).readId] & 31)  > (*it).editDist)
			{
				//this match is better than any match before
				cnts[1][(*it).readId] = cnts[0][(*it).readId]; // best before is now second best 
									       // (count will be updated when match is removed)
				cnts[0][(*it).readId] = (*it).editDist; // best dist is this editDist
									// count is 0 (will be updated if countFirstTwo)
			}
			if(!dontCountFirstTwo) 
				if((cnts[0][(*it).readId]>>5) != 2047) cnts[0][(*it).readId] += 32;	// shift 5 to the right, add 1, shift 5 to the left, and keep count
		}
		*dit = *it;
		++dit;
	}

	resize(matches, dit - begin(matches, Standard()));
}
#endif


#ifdef RAZERS_MICRO_RNA

template < typename TMatches, typename TSpec >
void purgeAmbiguousRnaMatches(TMatches &matches, RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type                          TMatch;
	typedef typename Iterator<TMatches, Standard>::Type             TIterator;

	::std::sort(begin(matches, Standard()),end(matches, Standard()),LessRNoEdistHLen<TMatch>());
	int bestMScore = 0;

	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int errorsCutOff = MaxValue<int>::VALUE;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	for (; it != itEnd; ++it)
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).readId)
		{
			if ((*it).editDist >= errorsCutOff) continue;
			if ( (*it).mScore < bestMScore) continue;
			if (++hitCount >= hitCountCutOff)
			{
				if (hitCount == hitCountCutOff)
					dit = ditBeg;
				continue;
			}
		}
		else
		{
			readNo = (*it).readId;
			hitCount = 0;
			if (options.distanceRange > 0)
				errorsCutOff = (*it).editDist + options.distanceRange;
			bestMScore = (*it).mScore;
			ditBeg = dit;
		}
		*dit = *it;
		++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
}

/*

//////////////////////////////////////////////////////////////////////////////
// Hamming verification recording sum of mismatching qualities in m.mScore
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,					// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,				// read number
	TReadSet &readSet,				// reads
	RazerSOptions<TSpec> const &options,		// RazerS options
	MicroRNA)					// MaqMapping
{

	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read		= readSet[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git	= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);

	// this is max number of errors the seed should have
	unsigned maxErrorsSeed = (unsigned)(options.rnaSeedLength * options.errorRate);	
	unsigned minSeedErrors = maxErrorsSeed + 1;
	unsigned bestHitLength = 0;

	for (; git < gitEnd; ++git)
	{
		bool hit = true;
		unsigned hitLength = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		TGenomeIterator g = git;	//maq would count errors in the first 28bp only (at least initially. not for output)
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
		{
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
				if (count < options.rnaSeedLength)		// the maq (28bp-)seed
				{
					if(++seedErrors > maxErrorsSeed)
					{
						hit = false;
						break;
					}
				}
				else 
					break;
			}
			++count;
		}
		if (hit) hitLength = count;
		if (hitLength > bestHitLength) //simply take the longest hit
		{
			minSeedErrors = seedErrors;
			bestHitLength = hitLength;
			m.beginPos = git - begin(host(inf), Standard());
		}
	}

//	std::cout  << "options.absMaxQualSumErrors" << options.absMaxQualSumErrors << std::endl;
//	std::cout  << "maxSeedErrors" << maxErrorsSeed << std::endl;
//	std::cout  << "minErrors" << minSeedErrors << std::endl;
//	if(derBesgte) ::std::cout << minErrors <<"minErrors\n";
	if (minSeedErrors > maxErrorsSeed) return false;
	
	m.endPos = m.beginPos + bestHitLength;
	m.editDist = minSeedErrors;			// errors in seed or total number of errors?
	m.mScore = bestHitLength;
	m.seedEditDist = minSeedErrors;
	return true;
}	
*/
#endif

struct SemiGlobalHamming_;
struct SemiGlobalEdit_;


//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned readId,						// read number
	TReadSet &readSet,						// reads
	SwiftSemiGlobalHamming)					// Hamming only
{
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) 
		return matchVerify(m,inf,rseqNo,readSet,options,QualityBasedScoring());
#endif

#ifdef RAZERS_MICRO_RNA
	if(options.microRNA) 
		return matchVerify(m,inf,rseqNo,readSet,options,MicroRNA());
#endif

	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId] << ::std::endl;
#endif

	// verify
	TRead &read				= readSet[readId];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	unsigned ndlLength		= ritEnd - ritBeg;

	if (length(inf) < ndlLength) return false;
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	
	unsigned maxErrors = (unsigned)(ndlLength * verifier.options.errorRate);
	unsigned minErrors = maxErrors + 1;
	unsigned errorThresh = (verifier.oneMatchPerBucket)? MaxValue<unsigned>::VALUE: maxErrors;

	for (; git < gitEnd; ++git)
	{
		unsigned errors = 0;
		TGenomeIterator g = git;
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
			if ((verifier.options.compMask[ordValue(*g)] & verifier.options.compMask[ordValue(*r)]) == 0)
				if (++errors > maxErrors)
					break;
		
		if (errors < minErrors)
		{
			minErrors = errors;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (errorThresh < errors)
		{
			if (minErrors <= maxErrors)
			{
				verifier.m.endPos = verifier.m.beginPos + ndlLength;
				verifier.q.pairScore = verifier.q.score = -(int)minErrors;
				verifier.q.errors = minErrors;
				verifier.push();
				minErrors = maxErrors + 1;
			}
		}
	}

	if (minErrors <= maxErrors)
	{
		verifier.m.endPos = verifier.m.beginPos + ndlLength;
		verifier.q.pairScore = verifier.q.score = -(int)minErrors;
		verifier.q.errors = minErrors;
		if (!verifier.oneMatchPerBucket)
			verifier.push();
		return true;
	}
	return false;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned readId,						// read number
	TReadSet &readSet,						// reads
	SwiftSemiGlobal)						// Mismatches and Indels
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Prefix<TRead>::Type					TReadPrefix;
	typedef typename Position<TGenomeInfix>::Type			TPosition;

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
	typedef typename TMatchVerifier::TPreprocessing			TPreprocessing;
	typedef typename Value<TPreprocessing>::Type			TMyersPattern;

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef Finder<TGenomeInfixRev>							TMyersFinderRev;

#ifdef RAZERS_NOOUTERREADGAPS
	typedef ModifiedString<TReadPrefix, ModReverse>			TReadPrefixRev;
	typedef Pattern<TReadPrefixRev, MyersUkkonenGlobal>		TMyersPatternRev;
#else
	typedef ModifiedString<TRead, ModReverse>				TReadRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;
#endif

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = verifier.preprocessing[readId];

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId]<<::std::endl;
#endif

    unsigned ndlLength = sequenceLength(readId, readSet);
	int maxScore = MinValue<int>::VALUE;
	int minScore = -(int)(ndlLength * verifier.options.errorRate);
	TPosition maxPos = 0;
	TPosition lastPos = length(inf);
	unsigned minDistance = (verifier.oneMatchPerBucket)? lastPos: 1;

#ifdef RAZERS_NOOUTERREADGAPS
	TGenomeInfix origInf(inf);
	setEndPosition(inf, endPosition(inf) - 1);
	--ndlLength;
#endif
	
	// find end of best semi-global alignment
	while (find(myersFinder, myersPattern, minScore))
	{
		TPosition pos = position(hostIterator(myersFinder));
		int score = getScore(myersPattern);
		
#ifdef RAZERS_NOOUTERREADGAPS
		// Manually align the last base of the read
		//
		// In this case myersPattern contains the whole read without the
		// last base. We compare the bases and adjust the score.
		// We also have to adjust inf and remove the last base of the
		// genomic region that has to be verified.
		SEQAN_ASSERT_LT(pos + 1, length(origInf));
		if ((verifier.options.compMask[ordValue(origInf[pos + 1])] & verifier.options.compMask[ordValue(back(readSet[readId]))]) == 0)
			if (--score < minScore) continue;
#endif		
		if (lastPos + minDistance < pos)
		{
			if (minScore <= maxScore)
			{
				TPosition infBeginPos = beginPosition(inf);
				TPosition infEndPos = endPosition(inf);
				verifier.q.pairScore = verifier.q.score = maxScore;
				verifier.q.errors = -maxScore;

				// find beginning of best semi-global alignment
				setEndPosition(inf, verifier.m.endPos = (beginPosition(inf) + maxPos + 1));

#ifdef RAZERS_NOOUTERREADGAPS
				// The best score must be corrected to hold the score of the prefix w/o the last read base
				if ((verifier.options.compMask[ordValue(origInf[maxPos + 1])] & verifier.options.compMask[ordValue(back(readSet[readId]))]) == 0)
					++maxScore;

				TReadPrefixRev		readRev(prefix(readSet[readId], ndlLength));
				TMyersPatternRev	myersPatternRev(readRev);
#else
				TReadRev			readRev(readSet[readId]);
				TMyersPatternRev	myersPatternRev(readRev);
#endif

//				// limit the beginning to needle length plus errors (== -maxScore)
//				if (length(inf) > ndlLength - maxScore)
//					setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);

				// we eventually have to search before the beginning of our parallelogram
				// otherwise alignments of an island in the previous parallelogram
				// could be cut and prevent that an island in this parallelgram is found
				if (endPosition(inf) > (unsigned)(ndlLength - maxScore))
					setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
				else
					setBeginPosition(inf, 0);
				
				TGenomeInfixRev		infRev(inf);
				TMyersFinderRev		myersFinderRev(infRev);

				_patternMatchNOfPattern(myersPatternRev, verifier.options.matchN);
				_patternMatchNOfFinder(myersPatternRev, verifier.options.matchN);
				while (find(myersFinderRev, myersPatternRev, /*score*/maxScore))
					verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);

#ifdef RAZERS_NOOUTERREADGAPS
				// The match end position must be increased by the omitted base.
				++verifier.m.endPos;
#endif
				verifier.push();
				setBeginPosition(inf, infBeginPos);
				setEndPosition(inf, infEndPos);
				maxScore = minScore - 1;
			}
		}
		if (score >= maxScore) 
		{
			maxScore = score;
			maxPos = pos;
		}
		lastPos = pos;
	}

	if (minScore <= maxScore)
	{
		verifier.q.pairScore = verifier.q.score = maxScore;
		verifier.q.errors = -maxScore;
		setEndPosition(inf, verifier.m.endPos = (beginPosition(inf) + maxPos + 1));

#ifdef RAZERS_NOOUTERREADGAPS
		// The best score must be corrected to hold the score of the prefix w/o the last read base
		if ((verifier.options.compMask[ordValue(origInf[maxPos + 1])] & verifier.options.compMask[ordValue(back(readSet[readId]))]) == 0)
			++maxScore;

		TReadPrefixRev		readRev(prefix(readSet[readId], ndlLength));
		TMyersPatternRev	myersPatternRev(readRev);
#else
		TReadRev			readRev(readSet[readId]);
		TMyersPatternRev	myersPatternRev(readRev);
#endif

//		// limit the beginning to needle length plus errors (== -maxScore)
//		if (length(inf) > ndlLength - maxScore)
//			setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);

		// we eventually have to search before the beginning of our parallelogram
		// otherwise alignments of an island in the previous parallelogram
		// could be cut and prevent that an island in this parallelgram is found
		if (endPosition(inf) > (unsigned)(ndlLength - maxScore))
			setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
		else
			setBeginPosition(inf, 0);
		
		// find beginning of best semi-global alignment
		TGenomeInfixRev		infRev(inf);
		TMyersFinderRev		myersFinderRev(infRev);

		_patternMatchNOfPattern(myersPatternRev, verifier.options.matchN);
		_patternMatchNOfFinder(myersPatternRev, verifier.options.matchN);
		while (find(myersFinderRev, myersPatternRev, maxScore))
			verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);

#ifdef RAZERS_NOOUTERREADGAPS
		// The match end position must be increased by the omitted base.
		++verifier.m.endPos;
#endif
		if (!verifier.oneMatchPerBucket)
			verifier.push();
		return true;
	}
	return false;
}


#ifdef RAZERS_DIRECT_MAQ_MAPPING
/*
//////////////////////////////////////////////////////////////////////////////
// Hamming verification recording sum of mismatching qualities in m.mScore
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,					// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,				// read number
	TReadSet &readSet,				// reads
	RazerSOptions<TSpec> const &options,		// RazerS options
	QualityBasedScoring)					// MaqMapping
{
	
	typedef Segment<TGenome, InfixSegment>				TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

// #ifdef RAZERS_DEBUG
// 	cout<<"Verify: "<<::std::endl;
// 	cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
// 	cout<<"Read:   "<<host(myersPattern)<<::std::endl;
// #endif

//	bool derBesgte = false;
	//if(rseqNo == 2) derBesgte = true;
//	if(derBesgte) ::std::cout << "der besagte\n";
	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read		= readSet[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git	= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	TGenomeIterator bestIt	= begin(inf, Standard());

	// this is max number of errors the 28bp 'seed' should have
	//assuming that maxErrors-1 matches can be found with 100% SN 
	unsigned maxErrorsSeed = (unsigned)(options.artSeedLength * options.errorRate) + 1;	
	unsigned maxErrorsTotal = (unsigned)(ndlLength * 0.25); //options.maxErrorRate);
	unsigned minErrors = maxErrorsTotal + 1;
	int minQualSumErrors = options.absMaxQualSumErrors + 10;
	unsigned minSeedErrors = maxErrorsSeed + 1;

	for (; git < gitEnd; ++git)
	{
		bool hit = true;
		unsigned errors = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		int qualSumErrors = 0;
		TGenomeIterator g = git;	//maq would count errors in the first 28bp only (at least initially. not for output)
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
		{
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
			//	::std::cout << count << "<-";
				if (++errors > maxErrorsTotal)
				{
					hit = false;
					break;
				}
				//int qualityValue = (int)((unsigned char)*r >> 3);
				//qualSumErrors += (qualityValue < options.mutationRateQual) ? qualityValue : options.mutationRateQual;
				qualSumErrors += (getQualityValue(*r) < options.mutationRateQual) ? getQualityValue(*r) : options.mutationRateQual;
				if(qualSumErrors > options.absMaxQualSumErrors || qualSumErrors > minQualSumErrors)
				{
					hit = false;
					break;
				}
				if (count < options.artSeedLength)		// the maq (28bp-)seed
				{
					if(++seedErrors > maxErrorsSeed)
					{
						hit = false;
						break;
					}
					if(qualSumErrors > options.maxMismatchQualSum)
					{
						hit = false;
						break;
					}// discard match, if 'seed' is bad (later calculations are done with the quality sum over the whole read)
				}
			}
			++count;
		}
//		if (hit && (qualSumErrors < minQualSumErrors && seedErrors <=maxErrorsSeed)) //oder (seedErrors < minSeedErrors)
		if (hit && (qualSumErrors < minQualSumErrors))
		{
			minErrors = errors;
			minSeedErrors = seedErrors;
			minQualSumErrors = qualSumErrors;
			m.beginPos = git - begin(host(inf), Standard());
			bestIt = git;
		}
	}
	if (minQualSumErrors > options.absMaxQualSumErrors || minSeedErrors > maxErrorsSeed || minErrors > maxErrorsTotal) return false;
	
	m.endPos = m.beginPos + ndlLength;
	m.editDist = minErrors;			// errors in seed or total number of errors?
	m.mScore = minQualSumErrors;
	m.seedEditDist = minSeedErrors;
	return true;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet,
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,	    				// reads
	TMyersPatterns const & pat,				// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobal const &swiftsemi,				// Hamming only
	QualityBasedScoring)						// Swift specialization
{
	//if(!options.maqMapping) 
		return matchVerify(m,inf,rseqNo,readSet,pat,options,swiftsemi);
	//else
	//	return matchVerify(m,inf,rseqNo,readSet,pat,options,swiftsemi); // todo!
}
*/
#endif




#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TFragmentStore, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions >
void mapSingleReads(
	TFragmentStore							& store,
	TGenome									& genome,				// genome ...
	unsigned								  contigId,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPattern,
	TPreprocessing							& preprocessing,
	TCounts									& cnts,
	char									  orientation,				// q-gram index of reads
	TRazerSOptions							& options)
{
	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >								TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >							TSwiftPattern;

	// VERIFICATION
	typedef MatchVerifier <
		TFragmentStore, 
		TRazerSOptions, 
		TPreprocessing, 
		TSwiftPattern,
		TCounts >															TVerifier;
	typedef typename Fibre<TReadIndex, FibreText>::Type					TReadSet;
	
	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]";
		else                    ::std::cerr << "[rev]";
	}

	TReadSet		&readSet = host(host(swiftPattern));
	TSwiftFinder	swiftFinder(genome, options.repeatLength, 1);
	TVerifier		verifier(store, options, preprocessing, swiftPattern, cnts);

	verifier.onReverseComplement = (orientation == 'R');
	verifier.genomeLength = length(genome);
	verifier.m.contigId = contigId;

	// iterate all verification regions returned by SWIFT
#ifdef RAZERS_MICRO_RNA
	while (find(swiftFinder, swiftPattern, 0.2)) 
#else
	while (find(swiftFinder, swiftPattern, options.errorRate)) 
#endif
	{
		verifier.m.readId = (*swiftFinder.curHit).ndlSeqNo;
		//-i 98 -vv -id -rr 100    -of 4 -o razers100_02.sam data/saccharomyces/genome.fasta.cut    data/saccharomyces/reads_454/SRR001317.1k.fasta
//		if (store.readNameStore[verifier.m.readId] == "SRR001317.770.1")
//			std::cout<<"gotit"<<std::endl;
		if (!options.spec.DONT_VERIFY)
			matchVerify(verifier, infix(swiftFinder), verifier.m.readId, readSet, TSwiftSpec());
		++options.countFiltration;
	}
}
#endif



#ifdef RAZERS_MICRO_RNA

	// multiple sequences
	template <
		typename TSA, 
		typename TString, 
		typename TSpec, 
		typename TShape, 
		typename TDir, 
		typename TValue, 
		typename TWithConstraints
	>
	inline void
	_qgramFillSuffixArray(
		TSA &sa, 
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape, 
		TDir &dir,
		TWithConstraints const,
		TValue prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type				TSize;
		typedef typename Value<TShape>::Type				THash;

		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = prefixLen - length(shape) + 1;

			typename Value<TSA>::Type localPos;
			assignValueI1(localPos, seqNo);
			assignValueI2(localPos, 0);

			TIterator itText = begin(sequence, Standard());
			if (TWithConstraints::VALUE) {
				THash h = hash(shape, itText) + 1;						// first hash
				if (dir[h] != (TSize)-1) sa[dir[h]++] = localPos;		// if bucket is enabled
			} else
				sa[dir[hash(shape, itText) + 1]++] = localPos;			// first hash

			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				assignValueI2(localPos, i);
				if (TWithConstraints::VALUE) {
					THash h = hashNext(shape, itText) + 1;				// next hash
					if (dir[h] != (TSize)-1) sa[dir[h]++] = localPos;	// if bucket is enabled
				} else
					sa[dir[hashNext(shape, itText) + 1]++] = localPos;	// next hash
			}
		}
	}


	template < typename TDir, typename TString, typename TSpec, typename TShape, typename TValue >
	inline void
	_qgramCountQGrams(TDir &dir, StringSet<TString, TSpec> const &stringSet, TShape &shape, TValue prefixLen)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Value<TDir>::Type							TSize;
	
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = prefixLen - length(shape) + 1;

			TIterator itText = begin(sequence, Standard());
			++dir[hash(shape, itText)];
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				++dir[hashNext(shape, itText)];
			}
		}
	}
	
	
	template < typename TIndex, typename TValue>
	void createQGramIndex(TIndex &index, TValue prefixLen)
	{
	SEQAN_CHECKPOINT
		typename Fibre<TIndex, QGramText>::Type const &text  = indexText(index);
		typename Fibre<TIndex, QGramSA>::Type         &sa    = indexSA(index);
		typename Fibre<TIndex, QGramDir>::Type        &dir   = indexDir(index);
		typename Fibre<TIndex, QGramShape>::Type      &shape = indexShape(index);
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		_qgramCountQGrams(dir, text, shape, prefixLen);

		if (_qgramDisableBuckets(index))
		{
			// 3. cumulative sum
			_qgramCummulativeSum(dir, True());

			// 4. fill suffix array
			_qgramFillSuffixArray(sa, text, shape, dir, True(), prefixLen);

			// 5. correct disabled buckets
			_qgramPostprocessBuckets(dir);
		}
		else
		{
			// 3. cumulative sum
			_qgramCummulativeSum(dir, False());
			
			// 4. fill suffix array
			_qgramFillSuffixArray(sa, text, shape, dir, False(),prefixLen);
		} 
	}


#endif



//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TGNoToFile,
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString>				& genomeFileNameList,
	String<TGNoToFile>					& gnoToFileMap,
	TCounts								& cnts,
	RazerSOptions<TSpec>				& options,
	TShape const						& shape,
	Swift<TSwiftSpec> const)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	typedef Index<TReadSeqStore, IndexQGram<TShape,OpenAddressing> >	TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

/*	// try opening each genome file once before running the whole mapping procedure
	int filecount = 0;
	int numFiles = length(genomeFileNameList);
	while(filecount < numFiles)
	{
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;
		file.close();
		++filecount;
	}
	*/

	// configure q-gram index
	TIndex swiftIndex(store.readSeqStore, shape);
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPattern(swiftIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;
	swiftPattern.params.printDots = options._debugLevel > 0;

	// init edit distance verifiers
	unsigned readCount = countSequences(swiftIndex);
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
#ifdef RAZERS_NOOUTERREADGAPS
			if (!empty(indexText(swiftIndex)[i]))
				setHost(forwardPatterns[i], prefix(indexText(swiftIndex)[i], length(indexText(swiftIndex)[i]) - 1));
#else
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
#endif
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}
#ifdef RAZERS_MICRO_RNA
	typename Size<TIndex>::Type qgram_count = 0;
	if(options.microRNA)
	{
		for(unsigned i = 0; i < countSequences(swiftIndex); ++i)
			if (sequenceLength(i, swiftIndex) >= options.rnaSeedLength)
				qgram_count += options.rnaSeedLength - (length(shape) - 1);
		resize(indexSA(swiftIndex), qgram_count, Exact());
		resize(indexDir(swiftIndex), _fullDirLength(swiftIndex), Exact());
		createQGramIndex(swiftIndex,options.rnaSeedLength);
	}
#endif

	
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping)
	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			resize(cnts[i], readCount, 31); //initialize with maxeditDist, 11:5 for count:dist
	}
#endif

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	unsigned filecount = 0;
	unsigned numFiles = length(genomeFileNameList);
	unsigned gseqNo = 0;

	// open genome files, one by one	
	while (filecount < numFiles)
	{
		// open genome file	
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;

		// remove the directory prefix of current genome file
		::std::string genomeFile(toCString(genomeFileNameList[filecount]));
		size_t lastPos = genomeFile.find_last_of('/') + 1;
		if (lastPos == genomeFile.npos) lastPos = genomeFile.find_last_of('\\') + 1;
		if (lastPos == genomeFile.npos) lastPos = 0;
		::std::string genomeName = genomeFile.substr(lastPos);
		

		CharString	id;
		Dna5String	genome;
		unsigned gseqNoWithinFile = 0;
		// iterate over genome sequences
		SEQAN_PROTIMESTART(find_time);
		for(; !_streamEOF(file); ++gseqNo)
		{
			if (options.genomeNaming == 0)
			{
				//readID(file, id, Fasta());			// read Fasta id
				readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
				appendValue(store.contigNameStore, id, Generous());
			}
			read(file, genome, Fasta());			// read Fasta sequence
			
			appendValue(gnoToFileMap, TGNoToFile(genomeName, gseqNoWithinFile));
			
			if (options.forward)
				mapSingleReads(store, genome, gseqNo, swiftPattern, forwardPatterns, cnts, 'F', options);

			if (options.reverse)
			{
				reverseComplement(genome);
				mapSingleReads(store, genome, gseqNo, swiftPattern, forwardPatterns, cnts, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		++filecount;
	}

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (given as StringSet)
template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet,
	typename TCounts, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type				TRead;
	typedef Index<TReadSet, IndexQGram<TShape,OpenAddressing> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

	// configure q-gram index
	TIndex swiftIndex(readSet, shape);
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPattern(swiftIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;

	// init edit distance verifiers
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		unsigned readCount = countSequences(swiftIndex);
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
#ifdef RAZERS_NOOUTERREADGAPS
			if (!empty(indexText(swiftIndex)[i]))
				setHost(forwardPatterns[i], prefix(indexText(swiftIndex)[i], length(indexText(swiftIndex)[i]) - 1));
#else
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
#endif
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	CharString	id;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
 	if(options.maqMapping)
 	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			resize(cnts[i], length(readSet), 31);
	}
#endif
	
	
	
	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; gseqNo < length(genomeSet); ++gseqNo)
	{
		if (options.forward)
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, cnts, 'F', options);

		if (options.reverse)
		{
			reverseComplement(genomeSet[gseqNo]);
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, cnts, 'R', options);
			reverseComplement(genomeSet[gseqNo]);
		}

	}
	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for single/mate-pair mapping
template <
	typename TFSSpec, 
	typename TFSConfig,
	typename TGNoToFile,
	typename TCounts,
	typename TSpec,
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString> &	genomeFileNameList,
	String<TGNoToFile> & gnoToFileMap,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(store, genomeFileNameList, gnoToFileMap, cnts, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(store, genomeFileNameList, gnoToFileMap, cnts, options, shape, Swift<TSwiftSpec>());
}

template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TFSSpec, typename TFSConfig, typename TGNoToFile, typename TCounts, typename TSpec>
int mapReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	StringSet<CharString> &	genomeFileNameList,
	String<TGNoToFile> & gnoToFileMap,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, cnts, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, cnts, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, cnts, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, cnts, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, cnts, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(store, genomeFileNameList, gnoToFileMap, cnts, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TCounts, typename TSpec>
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

}

#endif
