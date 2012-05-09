 /*==========================================================================
  SNP Calling routine of RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by Anne-Katrin Emde

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

#ifndef SEQAN_HEADER_CALLSNPS_H
#define SEQAN_HEADER_CALLSNPS_H


#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options

struct SnpStoreSpec_;
template <>
struct FragmentStoreConfig<SnpStoreSpec_ > 
{
	typedef String<Dna5Q>		TReadSeq;
	typedef String<Dna5Q>		TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void			TReadStoreElementSpec;
	typedef Owner<>			TReadSeqStoreSpec;
	typedef void			TMatePairStoreElementSpec;
	typedef void			TLibraryStoreElementSpec;
	typedef void			TContigStoreElementSpec;
	typedef void			TContigFileSpec;
	typedef void			TAlignedReadStoreElementSpec;
	typedef Owner<>			TAlignedReadTagStoreSpec;
	typedef void			TAnnotationStoreElementSpec;
};


struct SnpStoreGroupSpec_;
template <>
struct FragmentStoreConfig<SnpStoreGroupSpec_ > 
{
	typedef String<Dna5Q>		TReadSeq;
	typedef String<Dna5Q>		TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void			TReadStoreElementSpec;
	typedef Dependent<>		TReadSeqStoreSpec;
	typedef void			TMatePairStoreElementSpec;
	typedef void			TLibraryStoreElementSpec;
	typedef void			TContigStoreElementSpec;
	typedef void			TContigFileSpec;
	typedef void			TAlignedReadStoreElementSpec;
	typedef Owner<>			TAlignedReadTagStoreSpec;
	typedef void			TAnnotationStoreElementSpec;
};



	template <typename TGPos_>
	struct SimplePosition
	{
		typedef typename MakeSigned_<TGPos_>::Type TGPos;
		
		TGPos			gBegin;			// begin position in the genome 
		unsigned		gseqNo;


	};

	template < bool _HAMMING_ONLY = true >
	struct SNPCallingSpec 
	{
		enum { HAMMING_ONLY = _HAMMING_ONLY };				// omit verifying potential matches
	};
	
	
	
	template < typename TSpec = SNPCallingSpec<> >
	struct SNPCallingOptions
	{
		// main options
		TSpec		spec;
		const char	*output;		// name of snp output file
		int		_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
		bool		hammingOnly;		// no indels
		std::stringstream	programCall;
	
		// output format options
		unsigned	outputFormat;		// 0 (detailed output of all candidate positios) 
							// or 1 (only successful candidates)
		unsigned	inputFormat;		// 0 = razers, 1=eland, 2 = maq
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		unsigned	positionFormat;		// 0..gap space
										// 1..position space

		// for compring nucleotides
		unsigned char compMask[5];
		String<char> toIupac;

		// statistics
		double		timeLoadFiles;		// time for loading input files
		double		timeDumpResults;	// time for dumping the results

		unsigned	forceReadLength;	// specify read length, trim reads even if they are longer in mapped read file and 
							// 
		unsigned	readLength;		// specify read length, for direct trimming
		std::string	runID;			// runID needed for gff output

		const char 	*positionFile;
		const char 	*outputPositionAnalysis;
		const char	*outputCNV;			// name of result file
	
		unsigned 	method;
		unsigned	maxPile;
		bool		laneSpecificMaxPile;			// 
		const char	*outputSNP;			// name of result file
		const char	*errorPrbFileName;			// name of result file
		const char 	*outputCompletePositionAnalysis;
		bool 		ocp_SuperCompact;
		const char	*coverageFile;			// name of coverage file (wiggle format)
		unsigned 	stepInterval;			//step size for coverage file

		bool		qualityFile;			// name of .qual input file
		const char	*tabFile;			// name of .tab.txt result file
		const char	*regionFile;			//
		const char	*regionCoverageFile;			//
		
		bool 		useBaseQuality;		// use base quality instead of min{base,mapping}
		float		avgQualT;
		float		percentageT;
		unsigned	minMutT;
#ifdef RAZERS_MAQ_MAPPING
		bool 		maqMapping;
#endif
		bool		keepSuboptimalReads;
		bool		keepMultiReads;
		bool		dontClip;
		double		hetRate;
		int		numHaplotypes;
		int		minMapQual;
		double 		theta;
		double		eta;
		String<long double> cnks;
		String<long double> fks;
		String<long double> hetTable;
		double 		priorHetQ;
		bool 		realign;
		int 		forceCallCount;
		int		realignAddBorder;

	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh
		unsigned 	maxHitLength;
		unsigned	minCoord;
		unsigned	maxCoord;

		unsigned	minCoverage;
		unsigned	extractMinCov;
		unsigned	extractMaxCov;
		
		bool 		orientationAware;
		bool		showQualityStrings;

		const char	*outputIndel;			// name of result file
		bool 		doIndelCalling;
		int		maxIndelLen;
		int		slopeSmooth;
		double 		slopeTolerance;
		
		float		indelPercentageT;
		unsigned	indelCountThreshold;
		unsigned 	indelWindow;
		unsigned 	rLen; //temporary
		
		unsigned 	minDifferentReadPos;
		unsigned 	excludeBorderPos;
		
		unsigned	windowSize;
		unsigned	windowBuff;
		unsigned	expectedReadsPerBin;
		unsigned	expectedReadsSD;
		unsigned	cnvWindowSize;
		int			minClippedLength;
		bool		clipTagsInFile;

		const char 	*outputLog;
						
		SNPCallingOptions() 
		{
			output = "";
			qualityFile = false;
			tabFile = "";
			regionFile = "";
			regionCoverageFile = "";
			programCall << "";
	
			_debugLevel = 0;
			printVersion = false;
			hammingOnly = true;
			
			outputFormat = 0;
			sortOrder = 0;
			positionFormat = 1;
			dontClip = false;
			
			keepSuboptimalReads = false;
			keepMultiReads = false;
			
			for (unsigned i = 0; i < 5; ++i)
				compMask[i] = 1 << i;
			//compMask[4] = 0;
			toIupac = "AMRWMCSYRSGKWYKT";
			
			windowSize = 1000000;
			windowBuff = 70;

			forceReadLength = 0;
			readLength = 0; //-> use entire read length
			runID = ""; //
			
			rLen = 32;
			method = 1;
			outputSNP = "";
			positionFile = "";
			outputPositionAnalysis = "";
			errorPrbFileName = "";
			inputFormat = 0;
			outputCompletePositionAnalysis = "";
			ocp_SuperCompact = false;
			
			forceCallCount = 0;
			maxPile = 4;
			laneSpecificMaxPile = true;
			avgQualT = 10;
			percentageT = (float)0.3;
			minMutT = 5;
			useBaseQuality = true;
			
			minCoverage = 5;
			extractMinCov = 0;
			extractMaxCov = 1024*4;	//sollte reichen
			orientationAware = false;
			minMapQual = 1;
			priorHetQ = 0;
			realign = false;
			realignAddBorder = 20;

			stepInterval = 1;
			coverageFile = "";
			maxHitLength = 1;			

			hetRate = 0.001;
			numHaplotypes = 2;
			theta = 0.85;
			showQualityStrings = true;
			eta = 0.03;
#ifdef RAZERS_MAQ_MAPPING

			maqMapping = false;
#endif

			indelPercentageT = (float) 0.3;
			indelCountThreshold = 5;
			indelWindow = 0;
			outputIndel = "";
			outputCNV = "";
			doIndelCalling = true;
			maxIndelLen = 10; // names are misleading, needs to be rethought
			slopeSmooth = 32; // both are < rLen-2
			slopeTolerance = 10.0;
			
			minDifferentReadPos = 0;
			excludeBorderPos = 0;
			expectedReadsPerBin = 125;
			expectedReadsSD = 25;
			cnvWindowSize = 1000;
			minClippedLength = 10; //20;
			clipTagsInFile = false;
			outputLog = "";
		}
	};
	
	
	
	
	
//////////////////////////////////////////////////////////////////////////////
// Typedefs

	// definition of a Read match
	template <typename TGPos_>
	struct MappedReadMatch 
	{
		typedef typename MakeSigned_<TGPos_>::Type TGPos;
		
// 		TGPos		Batch	gBegin;			// begin position of the match in the genome			--> beginPos
		TGPos			gEnd;			// end position of the match in the genome				--> endPos
		unsigned		rseqNo;			// read seqNo											--> readId
		
		unsigned 		gseqNo:15;		// genome seqNo		<32K sequences						--> contigId
		unsigned 		hasIndel:1;		// is 1 if read match contains indels, 0 else			--> gaps
		
		unsigned 		editDist:3;		// Levenshtein distance <8								--> errors
		unsigned 		mScore:7;		// mapping quality	<128								--> currently not in use anyway
		unsigned 		avgQuality:6;	// avg read quality	<64									--> score

		char			orientation;		// 'F'..forward strand, 'R'..reverse comp. strand	--> endPos > beginPos ?	

	};
	


	enum CALLSNPS_ERROR {
		CALLSNPS_GFF_FAILED = 1,
		CALLSNPS_GENOME_FAILED = 2,
		CALLSNPS_QUALITY_FAILED = 3,
		CALLSNPS_OUT_FAILED = 4
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions


// sort operators	
	template <typename TMatches, typename TMatchQualities>
	struct LessGStackMQ : 
		public ::std::binary_function < typename Value<TMatches>::Type, typename Value<TMatchQualities>::Type, bool >
	{
		TMatchQualities &qualStore;
		
		LessGStackMQ(TMatchQualities &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (
			typename Value<TMatches>::Type const &a, 
			typename Value<TMatches>::Type const &b) const 
		{
			typedef typename Value<TMatches>::Type TMatch;
			
		
			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;
	
			// begin position
			typename TMatch::TPos ba = _min(a.beginPos, a.endPos);
			typename TMatch::TPos bb = _min(b.beginPos, b.endPos);
	
			if (ba < bb) return true;
			if (ba > bb) return false;
		
			// end position
			typename TMatch::TPos ea = _max(a.beginPos, a.endPos);
			typename TMatch::TPos eb = _max(b.beginPos, b.endPos);

			if (ea < eb) return true;
			if (ea > eb) return false;
			
			// quality
			if (a.id == TMatch::INVALID_ID) return false;
			if (b.id == TMatch::INVALID_ID) return true;

			if (qualStore[a.id].score > qualStore[b.id].score) return true;
			if (!(qualStore[a.id].score >= qualStore[b.id].score)) return false;
			return (qualStore[a.id].errors < qualStore[b.id].errors);


		}
	};



    template <typename TPosLen>
    struct LessPosLen : public ::std::binary_function < TPosLen, TPosLen, bool >
    {
        inline bool operator() (TPosLen const &a, TPosLen const &b) const 
        {
            // read number
            if (a.i1 < b.i1) return true;
            if (a.i1 > b.i1) return false;

            return (a.i2 < b.i2);

        }
    };

	template <typename TMatches, typename TMatchQualities>
	struct LessGStackOaMQ : 
		public ::std::binary_function < typename Value<TMatches>::Type, typename Value<TMatchQualities>::Type, bool >
	{
		TMatchQualities &qualStore;
		
		LessGStackOaMQ(TMatchQualities &_qualStore):
			qualStore(_qualStore) {}
		
		inline bool operator() (
			typename Value<TMatches>::Type const &a, 
			typename Value<TMatches>::Type const &b) const 
		{
			typedef typename Value<TMatches>::Type TMatch;

			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// begin position
			typename TMatch::TPos ba = _min(a.beginPos, a.endPos);
			typename TMatch::TPos bb = _min(b.beginPos, b.endPos);
			if (ba < bb) return true;
			if (ba > bb) return false;

			// end position
			typename TMatch::TPos ea = _max(a.beginPos, a.endPos);
			typename TMatch::TPos eb = _max(b.beginPos, b.endPos);
			if (ea < eb) return true;
			if (ea > eb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			// quality
			if (a.id == TMatch::INVALID_ID) return false;
			if (b.id == TMatch::INVALID_ID) return true;
			if (qualStore[a.id].score > qualStore[b.id].score) return true;
			if (!(qualStore[a.id].score >= qualStore[b.id].score)) return false;
			return (qualStore[a.id].errors < qualStore[b.id].errors);

		}
	};


	
	template <typename TReadMatch>
	struct LessId : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome sequence
			return (a.readId < b.readId);

		}
	};


	// ... to sort matches according to gBegin 
	template <typename TReadMatch>
	struct LessGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome sequence
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// begin position
			return (_min(a.endPos,a.beginPos) < _min(b.endPos,b.beginPos));
		}
	};



	// ... to sort matches according to gEnd 
	template <typename TReadMatch>
	struct LessGPosEnd : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome sequence
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// end position
			return (_max(a.endPos,a.beginPos) < _max(b.endPos,b.beginPos));
		}
	};

	
	
	template <typename TReadMatch>
	struct LessGPosEndOa : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome sequence
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// end position
			if (_max(a.endPos,a.beginPos) < _max(b.endPos,b.beginPos)) return true;
			if (_max(a.endPos,a.beginPos) > _max(b.endPos,b.beginPos)) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			return oa;

		}
	};
	
		// ... to sort quality values // 
	template <typename TQual>
	struct HigherQ : public ::std::binary_function < TQual, TQual, bool >
	{
		inline bool operator() (TQual const &a, TQual const &b) const 
		{
			// quality
			return ordValue(a) > ordValue(b); // 
		}
	};

//_____________________________________________________________________________________//
/////////////////////////////////////////////////////////////////////////////////////////	





#ifdef PLATFORM_WINDOWS

template<typename TVal>
double
lgamma(TVal x)
{
	// TODO: replace
	x -= 1;
	if(x < 2) return log((double)1);
	else
	{
		double f = log((double)2);
		for(int s = 3; s <= x;++s) f += log((double)s);
		return f;
	}
  
}

#endif

template<typename TValue>
inline bool my_isnan(TValue value)
{
	return value != value;
}


// this is basically maq's source code translated into seqan
// see Li, H., Ruan, J. & Durbin, R. Mapping short DNA sequencing reads and calling variants
// using mapping quality scores. Genome Res. 2008.
// and http://maq.sourceforge.net
template <typename THomoTable, typename TDependencies, typename TOptions>
void computeCnks(THomoTable & cnks, TDependencies & fks, TOptions & options)
{
	typedef typename Value<THomoTable>::Type TValue;

	String<TValue> sum_a, beta, q_c, temp, lFks, lC;
	resize(sum_a,257);	
	resize(beta,256); 
	resize(q_c,256);
	resize(temp,256);
	resize(lFks,256); 
	resize(lC, 256*256); 
	resize(fks,256);
	resize(cnks, 256*256*64,0.0); // n<256, k<256, q<64 
	
	
	
	fks[0] = lFks[0] = 1.0;
	for (int n = 1; n < 256; ++n) 
	{
		fks[n] = pow(options.theta, n) * (1.0 - options.eta) + options.eta;
		lFks[n] = fks[n>>1]; //
	}
	
	for (int n = 1; n < 256; ++n)
		for (int k = 0; k <= n; ++k)  //
			lC[n<<8|k] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1); // for every combination of k errors in n reads,
										// (n and k share 16bit in lC)
	
	for (int q = 1; q < 64; ++q)
	{
		// for all quality values up to 64
		// these are the 'average' values computed from the mapped reads
		double e = pow(10.0, -q/10.0);
		double le = log(e);				
		double le1 = log(1.0-e);	

		for (int n = 1; n < 256; ++n) 
		{
			sum_a[n+1] = 0.0;
			for (int k = n; k >= 0; --k) { // a_k = \sum_{i=k}^n C^n_k \epsilon^k (1-\epsilon)^{n-k}
				sum_a[k] = sum_a[k+1] + expl(lC[n<<8|k] + k*le + (n-k)*le1); 
				beta[k] = sum_a[k+1] / sum_a[k];		
				if (beta[k] > 0.99) beta[k] = 0.99; 		
			}
			
			for (int k = 0; k < n; ++k) 					// c_k
				q_c[k] = -4.343 * lFks[k] * logl(beta[k] / e);	
			for (int k = 1; k < n; ++k)
				q_c[k] += q_c[k-1]; // \prod_{i=0}^k c_i
			
			for (int k = 0; k <= n; ++k) 
			{ 
				temp[k] = -4.343 * logl(1.0 - expl(lFks[k] * logl(beta[k])));
				cnks[q<<16|n<<8|k] = (k > 0 ? q_c[k-1] : 0);// + temp[k]; 
				cnks[q<<16|n<<8|k] += (my_isnan(temp[k]) ? 0 : temp[k]);
			}
			
		}
	}

// 	std::fstream testf;
// 	testf.open("cnks", std::ios_base::out);
// 	if(testf.is_open())
// 	{
// 		for(unsigned f = 0; f < length(cnks); ++f)
//  			testf << cnks[f] << " ";
//  		testf <<std::endl;
// 	}
// 	testf.close();

}


//maqs heterozygote probabilites
//returns het prior in phred scale
template<typename THeteroTable, typename TOptions>
double computeHetTable(THeteroTable & hetTable, TOptions & options)
{
	typedef typename Value<THeteroTable>::Type TValue;
	double poly_rate;
	int numHaplotypes = 2;//options.numHaplotypes;

	resize(hetTable,256*256); // for n,k < 256
	TValue sum_harmo = 0.0;
	for (int k = 1; k <= numHaplotypes - 1; ++k)
		sum_harmo += 1.0 / k;
	for (int n1 = 0; n1 < 256; ++n1) 
	{
		for (int n2 = 0; n2 < 256; ++n2) 
		{
			long double sum = 0.0;
			double lC = lgamma(n1+n2+1) - lgamma(n1+1) - lgamma(n2+1); // \binom{n1+n2}{n1}
			for (int k = 1; k <= numHaplotypes - 1; ++k) 
			{
				double pk = 1.0 / k / sum_harmo;
				double log1 = log((double)k/numHaplotypes);
				double log2 = log(1.0 - (double)k/numHaplotypes);
				sum += pk * 0.5 * (expl(log1*n2) * expl(log2*n1) + expl(log1*n1) * expl(log2*n2));
			}
			hetTable[n1<<8|n2] = lC + logl(sum);
		}
	}
	poly_rate = options.hetRate * sum_harmo;
	double hetPriorQ = -4.343 * log(2.0 * poly_rate / (1.0 - poly_rate));

	return hetPriorQ;
}



template<typename THomoTable, typename TDependencies, typename TQStrings,typename TVal>
void
getHomoProbs(THomoTable & cnks,
			TDependencies & fks,
			TQStrings & qualitiesForward,
			TQStrings & qualitiesReverse,
			int & best,
			int & secondBest,
			long double & probQ1,
			long double & probQ2, TVal
#ifdef SNPSTORE_DEBUG_CANDPOS
			candidatePos
#endif
			)
{

	typedef typename Value<THomoTable>::Type TValue;
	typedef typename Value< typename Value<TQStrings>::Type >::Type TQuality;
	
#ifdef SNPSTORE_DEBUG_CANDPOS
	bool extraV = false; //|| candidatePos < 9335310
	if(candidatePos==118487871) extraV = true; 
#endif	
	
	String<double> sumE, sumF;
	resize(sumE,4,0.0);
	resize(sumF,4,0.0);

	for(unsigned i = 0; i < 4; ++i)
	{
		sort(begin(qualitiesForward[i],Standard()),end(qualitiesForward[i],Standard()),HigherQ<TQuality>());
		sort(begin(qualitiesReverse[i],Standard()),end(qualitiesReverse[i],Standard()),HigherQ<TQuality>());
		//compute average q
		double fk = 0.0;
		double qual = 0.0;
#ifdef SNPSTORE_DEBUG_CANDPOS
		if(extraV) std::cout << "F base"<<i<<": " << std::flush;
#endif		
		for(unsigned j = 0; j < length(qualitiesForward[i]); ++j)
		{
			qual = ordValue(qualitiesForward[i][j])-33;
			//qual = rescale into regular log
			if(j>=256) fk = fks[255];
			else fk = fks[j];
			sumE[i] += fk * qual;  
			sumF[i] += fk;
#ifdef SNPSTORE_DEBUG_CANDPOS
			if(extraV)
			{
				std::cout << sumE[i] << " " << std::flush;
				std::cout << sumF[i] << " " << std::flush;
			}
#endif
		}
#ifdef SNPSTORE_DEBUG_CANDPOS
		if(extraV) std::cout << std::endl;
		if(extraV) std::cout << "R base"<<i<<": " << std::flush;
#endif
		for(unsigned j = 0; j < length(qualitiesReverse[i]); ++j)
		{
			qual = ordValue(qualitiesReverse[i][j])-33;
			if(j>=256) fk = fks[255];
			else fk = fks[j];
			sumE[i] += fk * qual;
			sumF[i] += fk;
#ifdef SNPSTORE_DEBUG_CANDPOS
			if(extraV)
			{
				std::cout << sumE[i] << " " << std::flush;
				std::cout << sumF[i] << " " << std::flush;
			}
#endif
		}
	
	}
#ifdef SNPSTORE_DEBUG_CANDPOS
	if(extraV)
	{
		for(unsigned j = 0; j < 256; ++j)
			std::cout << fks[j] << " " << std::flush;
		std::cout << std::endl;
	
	}
#endif
		
	best = -1;
	secondBest = -1;
	double bestSum = 0.0;
	double secondBestSum = 0.0;
	for (int j = 0; j < 4; ++j) 
	{
		if (sumE[j] > bestSum) 
		{
			secondBestSum = bestSum; 
			secondBest = best;
			bestSum = sumE[j]; 
			best = j;
		}
		else 
		{
			if (sumE[j] > secondBestSum) 
			{
				secondBestSum = sumE[j];
				secondBest = j;
			}
		}
	}
#ifdef SNPSTORE_DEBUG_CANDPOS
	if(extraV) std::cout <<"best="<<best <<" secondbest="<<secondBest << std::flush << std::endl;
#endif
	
	int qAvgBest = 0, qAvgSecondBest = 0;
	int countBest = 0, countSecondBest = 0;
	if(best != -1)
	{
		// normalized quality of the base with the best weighted sum of qualities
		qAvgBest = (int)(sumE[best]/sumF[best] + 0.5); 
		countBest = length(qualitiesForward[best]) + length(qualitiesReverse[best]);
	}
	else countBest = 0;
	if(secondBest != -1)
	{
		qAvgSecondBest = (int)(sumE[secondBest]/sumF[secondBest] + 0.5);
		countSecondBest = length(qualitiesForward[secondBest]) + length(qualitiesReverse[secondBest]);	
	}
	else countSecondBest = 0;

	//runter skalieren
	int countTotal = countBest + countSecondBest;
	if (countTotal > 255) 
	{
		countBest = int(255.0 * countBest / countTotal + 0.5);
		countSecondBest = int(255.0 * countSecondBest / countTotal + 0.5);
	//	if(extraV) ::std::cout <<  "countBest = " << countBest << "\tcountSecondBest = " << countSecondBest << ::std::endl;
		countTotal = 255;
	}
#ifdef SNPSTORE_DEBUG_CANDPOS
	if(extraV)std::cout << "qAvgBest" <<  qAvgBest << " qAvgSecond"<< qAvgSecondBest << "\n";
	if(extraV)std::cout << "totalCount" <<  countTotal << " countBest"<< countBest << " countSecondBest"<< countSecondBest << "\n";
#endif
	probQ1 = ((countSecondBest > 0) ? sumE[secondBest] : 0);
	probQ1 += (my_isnan(cnks[qAvgSecondBest<<16|countTotal<<8|countSecondBest])) ? 0 : cnks[qAvgSecondBest<<16|countTotal<<8|countSecondBest];
	probQ2 = ((countBest > 0) ? sumE[best] : 0);
	probQ2 += (my_isnan(cnks[qAvgBest<<16|countTotal<<8|countBest])) ? 0 : cnks[qAvgBest<<16|countTotal<<8|countBest];
	
#ifdef SNPSTORE_DEBUG_CANDPOS
	if(extraV)std::cout << "cnkBest" <<  cnks[qAvgBest<<16|countTotal<<8|countBest] << "  bei cnkindex " <<(qAvgBest<<16|countTotal<<8|countBest)<<"\n";
	if(extraV)std::cout << "cnkSecondBest" <<  cnks[qAvgSecondBest<<16|countTotal<<8|countSecondBest] <<  "  bei cnkindex " <<(qAvgSecondBest<<16|countTotal<<8|countSecondBest)<< "\n";
	if(extraV)std::cout << "probQ1" <<  probQ1 << "\n";
	if(extraV)std::cout << "probQ2" <<  probQ2 << "\n";
#endif
	
// 	if(extraV)
 //	{
 //		for(unsigned f = 0; f < length(cnks); ++f)
 //			std::cout << cnks[f] << " ";
 //		std::cout <<std::endl;
 //	}


}


///////////////////////////////////////////////////////////////////////////////////////////////////




//looks for mismatches in alignemnt and returns positions with respect to 2nd row (read sequence)
template<typename TAlign, typename TString>
void
getMismatchMutations(TAlign & align, TString & mutations)
{
	
	typedef typename Source<TAlign>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;


	TAlignIterator ali_it0_stop = iter(row(align,0),endPosition(cols(align)));
	TAlignIterator ali_it1_stop = iter(row(align,1),endPosition(cols(align)));
	TAlignIterator ali_it0 = iter(row(align,0),beginPosition(cols(align)));
	TAlignIterator ali_it1 = iter(row(align,1),beginPosition(cols(align)));   


	//std::cout << "getting cigar line\n";//ali0 len = " <<ali_it0_stop-ali_it0 << " \t ali1 len = "<<ali_it1_stop-ali_it1<<"\n";
	int refPos = 0;
	//int readPos = 0;
	
	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
	{
		while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
		{
			if(*ali_it1 != *ali_it0)
				appendValue(mutations,refPos);
			++refPos;
			++ali_it0;
			++ali_it1;
		}
		while(ali_it0!=ali_it0_stop && isGap(ali_it0))
		{
			++refPos;
			++ali_it0;
			++ali_it1;
		}
		while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
		{
			++ali_it0;
			++ali_it1;
		}
	}
	
}



//looks for position in source sequence of row1 and returns aligned position in row0
template<typename TAlign, typename TPosition>
int
getReadPos(TAlign & align, TPosition pos_row1, bool extraV = false)
{
	typedef typename Iterator<typename Row<TAlign>::Type, Rooted>::Type TAlignIterator;
	
	TAlignIterator ali_it0_stop = iter(row(align,0),endPosition(cols(align)));
	TAlignIterator ali_it1_stop = iter(row(align,1),endPosition(cols(align)));
	TAlignIterator ali_it0 = iter(row(align,0),beginPosition(cols(align))); // walks over read
	TAlignIterator ali_it1 = iter(row(align,1),beginPosition(cols(align))); // walks over ref
	
	int refPos = 0;
	int readPos = 0;
	if(extraV) std::cout << align ;
	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop && refPos < pos_row1)
	{
		while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1) &&  refPos < pos_row1)
		{
			++refPos;
			++readPos;
			++ali_it0;
			++ali_it1;
		}
		while(ali_it0!=ali_it0_stop && isGap(ali_it0) &&  refPos < pos_row1)
		{
			++refPos;
			++ali_it0;
			++ali_it1;
		}
		while(isGap(ali_it1)&& ali_it1!=ali_it1_stop &&  refPos <= pos_row1)
		{
			++readPos;
			++ali_it0;
			++ali_it1;
		}
	}
	if(isGap(ali_it0)) return -1;
	else return readPos;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


struct SingleBaseVariant{
	bool called;  // did this variant pass calling criteria?
	int genotype; // called diploid genotype (allele1 << 2 | allele2) 
	int count;	  // number of non-ref observations (only counting most frequent mutational base)
	int quality;  // a quality value associated with the call
	int coverage; // totalCoverage at position
};


/*
struct IndelVariant{
	bool called;   // did this variant pass calling criteria?
	int indelSize; // called diploid genotype (allele1 << 2 | allele2) 
	int count;	   // number of supporting reads
	int quality;   // a quality value associated with the call
	int coverage;  // totalCoverage at position
	DnaString sequence;
};
*/





struct TagMaqMethod_;
typedef Tag<TagMaqMethod_> const MaqMethod;

struct TagThresholdMethod_;
typedef Tag<TagThresholdMethod_> const ThresholdMethod;


template<typename TFragmentStore, typename TGroupStore, typename TMatchIterator>
void
copyFragmentStore(TGroupStore &fragStoreGroup,
				  TFragmentStore 			&fragmentStore,
				  TMatchIterator 			matchItBatchBegin,
				  TMatchIterator 			matchItBatchEnd,
				  typename TFragmentStore::TContigPos	groupStartPos,
				  typename TFragmentStore::TContigPos	groupEndPos)
{
	//TFragmentStore fragStoreGroup = fragmentStore; //clear(fragStoreGroup.alignedReadStore); resize; arrayCopy(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
	
	// pointers are enough
	fragStoreGroup.readSeqStore = fragmentStore.readSeqStore;
	fragStoreGroup.readStore = fragmentStore.readStore;
	fragStoreGroup.readNameStore = fragmentStore.readNameStore;
	fragStoreGroup.alignQualityStore = fragmentStore.alignQualityStore;
	

	// need to be copied / moved
	resize(fragStoreGroup.alignedReadStore,matchItBatchEnd-matchItBatchBegin,Exact());
	arrayCopy(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
	//arrayMoveForward(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore

	// shorten reference sequence to the current region (groupStartPos to groupEndPos)
	// has to be copied because it will be overwritten
	// fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartPos,groupEndPos);
	typedef typename TGroupStore::TContigStore		TContigStore;			// TGenomeSet
	typedef typename Value<TContigStore>::Type		TContig;

	TContig conti;
	conti.seq = infix(fragmentStore.contigStore[0].seq,groupStartPos,groupEndPos);
	appendValue(fragStoreGroup.contigStore, conti, Generous() );
	appendValue(fragStoreGroup.contigNameStore, fragmentStore.contigNameStore[0], Generous() );

}


///////////////////////////////////////////////////////////////////////////////////////
// Output SNPs, do realignment if a certain number of indels is observed in the reads
template <
	typename TFragmentStore,
	typename TReadCigars,
	typename TReadCounts,
	typename TGenomeName,
	typename TFile,
	typename TOptions
>
void dumpVariantsRealignBatchWrap(
	TFragmentStore				&fragmentStore,			// forward/reverse matches
	TReadCigars				&readCigars,
	TReadCounts const			&readCounts,
	TGenomeName const 			genomeID,			// genome name
	typename TFragmentStore::TContigPos	startCoord,			// startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
	typename TFragmentStore::TContigPos	currWindowBegin,
	typename TFragmentStore::TContigPos	currWindowEnd,
	TFile					&fileSNPs,
	TFile					&fileIndels,
	TOptions 				&options)
{

	typedef typename TFragmentStore::TAlignedReadStore 	TMatches;
	typedef typename Value<TMatches>::Type 			TMatch;
	typedef typename TFragmentStore::TAlignQualityStore 	TMatchQualities;
	typedef typename Value<TMatchQualities>::Type 		TMatchQuality;
	typedef typename TFragmentStore::TReadSeqStore	 	TReads;
	typedef typename Value<TReads>::Type 			TRead;
	typedef typename TFragmentStore::TContigStore 		TContigStore;
	typedef typename TFragmentStore::TContigPos 		TContigPos;
	typedef typename Value<TContigStore>::Type	 	TContig;
	typedef typename TFragmentStore::TContigSeq 		TContigSeq;
	typedef typename Iterator<TMatches,Standard>::Type	TMatchIterator;
	
	
	TMatches &matches		= fragmentStore.alignedReadStore;
	TMatchQualities &matchQualities = fragmentStore.alignQualityStore;

	::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPos<TMatch>());		

	TMatchIterator matchIt = begin(matches,Standard());
	TMatchIterator matchItEnd = end(matches,Standard());
	
	unsigned minNumIndels = options.indelCountThreshold;
	
//	std::fstream tmpfile;
//	tmpfile.open("Z:/seqan071010/projects/library/apps/chr4.beforeTotal.sam", ::std::ios_base::out);
//	write(tmpfile, fragmentStore, Sam());
//	tmpfile.close();

#ifdef SNPSTORE_DEBUG
	CharString strstr = "test";
	_dumpMatches(fragmentStore,strstr);
#endif

	// now find connected subsets, i.e. groups of reads that overlap
	// dont realign regions unworthy of realignment (no indel reads)
	while(matchIt != matchItEnd)
	{
		TMatchIterator matchItBatchBegin = matchIt;
		TContigPos groupEndPos = _max((*matchIt).endPos,(*matchIt).beginPos);
		TContigPos groupStartPos = _min((*matchIt).endPos,(*matchIt).beginPos);

		TContigPos groupStartCoordLocal = _max(0,(int)groupStartPos-options.realignAddBorder);
	
		int indelReadCount = 0; // how many reads have indels in the current group
		while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < groupEndPos)
		{
			groupEndPos = (_max((*matchIt).beginPos,(*matchIt).endPos) > groupEndPos) ? _max((*matchIt).beginPos,(*matchIt).endPos) : groupEndPos;
			// reads wont be needed anymore! (make sure this is the case!!!)
			(*matchIt).beginPos -= groupStartCoordLocal;
			(*matchIt).endPos -= groupStartCoordLocal;
			if(matchQualities[(*matchIt).id].pairScore == 1 ) ++indelReadCount;
			++matchIt;
		}
		TMatchIterator matchItBatchEnd = matchIt;
		unsigned numMatches = matchItBatchEnd -matchItBatchBegin;

		TContigPos groupEndCoordLocal = _min(groupEndPos+(TContigPos)options.realignAddBorder,(TContigPos)length(fragmentStore.contigStore[0].seq));
		
		if(numMatches >= options.minCoverage)
		{
			//make temporary fragstore for group
			// shorten reference sequence to the current region (groupStartPos to groupEndPos)

			//FragmentStore<SnpStoreGroupSpec_> fragStoreGroup;
			//copyFragmentStore(fragStoreGroup,fragmentStore,matchItBatchBegin,matchItBatchEnd,groupStartPos,groupEndPos);

			TFragmentStore fragStoreGroup = fragmentStore; //clear(fragStoreGroup.alignedReadStore); resize; arrayCopy(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
			resize(fragStoreGroup.alignedReadStore,numMatches,Exact());
			arrayMoveForward(matchItBatchBegin,matchItBatchEnd,begin(fragStoreGroup.alignedReadStore,Standard())); // reads wont be needed anymore
//			fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartPos,groupEndPos);
			fragStoreGroup.contigStore[0].seq = infix(fragmentStore.contigStore[0].seq,groupStartCoordLocal,groupEndCoordLocal);
	
#ifdef SNPSTORE_DEBUG
			std::cout << "in realign wrap: groupEndPos = " <<  groupEndPos << " groupStartPos=" <<  groupStartPos << std::endl;
			std::cout << "genomeLength= " <<  length(fragmentStore.contigStore[0].seq) << std::endl;

			CharString strstre = "testgroup";
			_dumpMatches(fragStoreGroup,strstre);

#endif
			groupStartPos += startCoord;
			groupEndPos += startCoord;
			//groupStartCoord = groupStartPos;
			TContigPos groupStartCoord = startCoord + groupStartCoordLocal;
			groupStartPos = _max(groupStartPos,currWindowBegin);
			groupEndPos = _min(groupEndPos,currWindowEnd);
		
			//the current group is formed by all reads from matchItBatchBegin until matchItBatchEnd
			if(indelReadCount >= (int)minNumIndels && options.realign)
			{
				//do realignment
				dumpVariantsRealignBatch(fragStoreGroup,readCigars,
					readCounts,genomeID,
					groupStartCoord,groupStartPos,groupEndPos,
					fileSNPs,fileIndels,options);
			}
			else
			{
				// todo: switch between with or without realignment in dumpSNPsBatch.. make global in any case
/*				dumpVariantsRealignBatch(fragStoreGroup,readCigars,
					readCounts,genomeID,
					groupStartCoord,groupStartPos,groupEndPos,
					fileSNPs,fileIndels,options);*/
					dumpSNPsBatch(fragStoreGroup,readCigars,
						readCounts,genomeID,
						groupStartCoord,groupStartPos,groupEndPos,
						fileSNPs,options); 
			}
		}
	
	}
	
	
}

///////////////////////////////////////////////////////////////////////
// SNP calling Maq style
template<typename TCounts, typename TQualities, typename TOptions>
inline bool
_doSnpCall(TCounts & countF,
		  TCounts & countR,
		  TQualities & qualF,
		  TQualities & qualR,
		  int &refAllele,
		  TOptions & options,
		  SingleBaseVariant &snp,
		  MaqMethod&
#ifdef SNPSTORE_DEBUG_CANDPOS
			, int candPos
#endif
		  )
{


		// the diploid reference genotype
		int genotypeRef = (refAllele<<2) | refAllele; 
		int genotypeCalled = genotypeRef, qCall1 = 0; 
		int genotypeCalled2 = genotypeRef, qCall2 = 0; 
		

#ifdef SNPSTORE_DEBUG_CANDPOS
		bool extraV = false;	
		if(candPos == 118487871) extraV = true;
		if(extraV)
		{
			::std::cout << "Forward qualities:\n" << std::flush;
			for(unsigned x = 0; x < length(qualF); ++x)
				::std::cout << qualF[x] << "\t";
			::std::cout << "\nReverse qualities:\n" << std::flush;
			for(unsigned x = 0; x < length(qualR); ++x)
				::std::cout << qualR[x] << "\t";
			::std::cout << "\n" << std::flush;
		}
#endif

		// do the Maq statistics
		// 
		// argmax P(g|D)=P(D|g)*P(g)/P(D)
		//    g
		// 

		// get pHomo for best and second best nucleotide
		int best, secondBest;
		long double pHet = 0, pHomo1 = 0, pHomo2 = 0;
		getHomoProbs(options.cnks,options.fks,qualF,qualR,best,secondBest,pHomo1,pHomo2,
#ifdef SNPSTORE_DEBUG_CANDPOS			
			candPos
#else
			0
#endif		
			);
		if(secondBest == -1)
		{
			if(best==refAllele) // shouldnt happen
				return false;
			secondBest = refAllele;
		}

		//get pHet
		int n = countF[best] + countR[best] + countF[secondBest] + countR[secondBest];
#ifdef SNPSTORE_DEBUG_CANDPOS
		if(extraV)
		{
			std::cout << " n = " <<n << std::endl;
			std::cout << "(countF[secondBest] + countR[secondBest]) = " << (countF[secondBest] + countR[secondBest]) << std::endl;
		}
#endif				
		if (n > 255) 
		{
			int temp2 = (int)((countF[secondBest] + countR[secondBest])*255.0/n + 0.5);
			int temp1 = (int)((countF[best] + countR[best])*255.0/n + 0.5);
#ifdef SNPSTORE_DEBUG_CANDPOS
		if(extraV)
			std::cout << "temp1 = " << temp1 << std::endl;
#endif				

			pHet = options.priorHetQ - 4.343 * options.hetTable[temp2<<8|temp1]; 
//			pHet = options.priorHetQ - 4.343 * options.hetTable[255<<8|temp]; 
		}
		else 
			pHet = options.priorHetQ - 4.343 * options.hetTable[(countF[secondBest] + countR[secondBest])<<8|(countF[best] + countR[best])]; 
//			pHet = options.priorHetQ - 4.343 * options.hetTable[n<<8|(countF[secondBest] + countR[secondBest])]; 
				
#ifdef SNPSTORE_DEBUG_CANDPOS
		if(extraV)
		{
			std::cout << "best = " << best << " secondBest = " << secondBest << std::endl;
			std::cout << "pHet = " << pHet << std::endl;
			std::cout << "pHomo1 = " << pHomo1 << " pHomo2 = " << pHomo2 << std::endl;
		}
#endif				

		pHet = (pHet > 0.0) ? pHet : 0.0;
		pHomo1 = (pHomo1 > 0.0) ? pHomo1 : 0.0;
		pHomo2 = (pHomo2 > 0.0) ? pHomo2 : 0.0;
			
		int het,homo1,homo2; //0,1,2
			
		//rank them and create the genotype
		if(pHet < pHomo1)
		{
			if(pHet < pHomo2)
			{
				het = 0; //het is best
				if(best==refAllele)
					genotypeCalled = (best<<2) | secondBest;
				else
					genotypeCalled = (secondBest<<2) | best;
					
				if(pHomo1<=pHomo2)    //(1)
				{
					qCall1 = (int)(pHomo1 - pHet  + 0.5);
					homo1 = 1; //second best
					genotypeCalled2 = (best<<2)| best;
					qCall2 = (int)(pHomo2 - pHomo1 + 0.5);
					homo2 = 2; // last 
				}
				else				//(2)
				{
					qCall1 = (int)(pHomo2 - pHet + 0.5);
					homo2 = 1;
					genotypeCalled2 = (secondBest<<2)| secondBest;
					qCall2 = (int)(pHomo1 - pHomo2 + 0.5);
					homo1 = 2;
				}
			}
			else
			{						//(3)
				// shouldnt happen
				homo2 = 0;
				qCall1 = (int)(pHet - pHomo2 + 0.5);
				genotypeCalled = (secondBest<<2)| secondBest;
				het = 1;
				qCall2 = (int)(pHomo1 - pHet + 0.5);
				if(best==refAllele)
					genotypeCalled2 = (best<<2)| secondBest;
				else
					genotypeCalled2 = (secondBest<<2) | best;
				homo1 = 2;
			}
		}
		else
		{
			if(pHomo2 < pHomo1)	//(4)
			{
				//this case shouldnt happen 
				homo2 = 0;
				qCall1 = (int)(pHomo1 - pHomo2 + 0.5);
				genotypeCalled = (secondBest<<2)| secondBest;
				homo1 = 1;
				qCall2 = (int)(pHet - pHomo1 + 0.5);
				genotypeCalled2 = (best<<2)| best;
				het = 2;
			}
			else 
			{
				homo1 = 0;
				genotypeCalled = (best<<2)| best;
				if(pHet<pHomo2)		//(5)
				{
					qCall1 = (int)(pHet - pHomo1 + 0.5);
					het = 1;
					if(best==refAllele)
						genotypeCalled2 = (best<<2)| secondBest;
					else
						genotypeCalled2 = (secondBest<<2) | best;
					qCall2 = (int)(pHomo2 - pHet + 0.5);
					homo2 = 2;
				}
				else		// (6)
				{
					qCall1 = (int)(pHomo2 - pHomo1+ 0.5);
					homo2 = 1;
					qCall2 = (int)(pHet - pHomo2 + 0.5);
					genotypeCalled2 = (secondBest<<2)| secondBest;
					het = 2;
				}
			}
		}

		if (het != 0 && homo2 == 0) // 
		{
#ifdef SNPSTORE_DEBUG_CANDPOS
			std::cout << "Second best is best homozygote?!" << std::endl;
#endif
			//return false;
			genotypeCalled2 = genotypeCalled;
			genotypeCalled = genotypeRef;	// disable call
			homo2 = 1;
			het = 2; 
			qCall1 = 0;
			qCall2 = 0;
		}
		
	//}

	unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
						   + countR[0] + countR[1] +countR[2] +countR[3] +countR[4];

	snp.genotype = genotypeCalled;
	snp.count    = countF[best] + countR[best];
	snp.quality  = qCall1;
	snp.coverage = totalCoverage;
	if (genotypeCalled == genotypeRef) snp.called = false;
	else snp.called = true;

	return true;
}


template<typename TCounts, typename TQualities, typename TOptions>
inline bool
_doSnpCall(TCounts & countF,
		  TCounts & countR,
		  TQualities & qualF,		// columnQualityF
		  TQualities & qualR,		// columnQualityR
		  int &refAllele,
		  TOptions & options,
		  SingleBaseVariant &snp,
		  ThresholdMethod&)
{
		
	// find potential mutation allele
	int allele1 = -1;	// most frequent allele
	int allele2 = -1;	// second most frequent allele
				
	unsigned maxCount=0;
	for(int k=0; k < 5; ++k)
	{
		if(countF[k]+countR[k] > maxCount)
		{
			maxCount = countF[k]+countR[k];
			allele1 = k;
		}
	}
	maxCount = 0;
	for(int k=0; k < 5; ++k)
	{
		if(k != allele1 && countF[k]+countR[k] >= maxCount)
		{
			maxCount = countF[k]+countR[k];
			allele2 = k;
		}
	}

	// No evidence of non-ref bases left... (used to happen with onthefly-pileupcorrection
	// cannot happen anymore as these positions would never be inspected)
	if(allele1==refAllele && allele2==refAllele)
	{
		::std::cout << "No non-ref base observed. Correct??\n";
		return false;
	}

	// get the mutational allele
	int mutAllele=allele1;
	if(allele1==refAllele) mutAllele=allele2;

	unsigned mutCoverage   = countF[mutAllele] + countR[mutAllele];
	unsigned totalCoverage = countF[0] + countF[1] +countF[2] +countF[3] +countF[4]
						   + countR[0] + countR[1] +countR[2] +countR[3] +countR[4];

	// the diploid reference genotype
	int genotypeRef = (refAllele<<2) | refAllele; 
	int genotypeCalled = genotypeRef; 

	// threshold method
	if( mutCoverage >= options.minMutT 
		&& (float)mutCoverage/totalCoverage >= options.percentageT 
		&& (float)(qualF[mutAllele]+qualR[mutAllele])/mutCoverage >= options.avgQualT)
	{
		genotypeCalled = (mutAllele<<2)|mutAllele; // we dont attempt real genotype calling here
	}

	
	snp.genotype = genotypeCalled;
	snp.count    = mutCoverage;
	snp.quality  = qualF[mutAllele]+qualR[mutAllele]/ mutCoverage;
	snp.coverage = totalCoverage;
	if (genotypeCalled == genotypeRef) snp.called = false;
	else snp.called = true;

	return true;
}


// write to file
template<typename TFile, typename TString, typename TQualities, typename TPos, typename TOptions>
inline bool
_write(TFile &file, 
	   SingleBaseVariant &snp, 
	   TQualities &qualityStringF, 
	   TQualities &qualityStringR,
	   int refAllele, 
	   TString &genomeID, 
	   TPos candPos, 
	   unsigned realCoverage,
	   TOptions &options)
{
	if (!file.is_open()) 
	{
		::std::cerr << "SNP output file is not open" << ::std::endl;
		return false;
	}

	//chromosome
	file << genomeID << '\t';
	file << candPos + options.positionFormat<< '\t';
	file << (Dna5)refAllele <<'\t';
	if(options.orientationAware)
	{
		if(options.showQualityStrings)
		{
			file << "["<<qualityStringF[0] <<"]\t";
			file << "["<<qualityStringF[1] <<"]\t";
			file << "["<<qualityStringF[2] <<"]\t";
			file << "["<<qualityStringF[3] <<"]\t";
			file << "["<<qualityStringR[0] <<"]\t";
			file << "["<<qualityStringR[1] <<"]\t";
			file << "["<<qualityStringR[2] <<"]\t";
			file << "["<<qualityStringR[3] <<"]\t";
		}
		else
		{
			file << length(qualityStringF[0]) <<"\t";
			file << length(qualityStringF[1]) <<"\t";
			file << length(qualityStringF[2]) <<"\t";
			file << length(qualityStringF[3]) <<"\t";
			file << length(qualityStringR[0]) <<"\t";
			file << length(qualityStringR[1]) <<"\t";
			file << length(qualityStringR[2]) <<"\t";
			file << length(qualityStringR[3]) <<"\t";
		}
	}
	else
	{
		if(options.showQualityStrings)
		{
			file << "["<<qualityStringF[0]<<qualityStringR[0] <<"]\t";
			file << "["<<qualityStringF[1]<<qualityStringR[1] <<"]\t";
			file << "["<<qualityStringF[2]<<qualityStringR[2] <<"]\t";
			file << "["<<qualityStringF[3]<<qualityStringR[3] <<"]\t";
		}
		else
		{
			file << length(qualityStringF[0])+length(qualityStringR[0]) <<"\t";
			file << length(qualityStringF[1])+length(qualityStringR[1]) <<"\t";
			file << length(qualityStringF[2])+length(qualityStringR[2]) <<"\t";
			file << length(qualityStringF[3])+length(qualityStringR[3]) <<"\t";
		}
	}
	file << realCoverage;
	
	if(options.method == 1)
	{
		//genotypeCalled to string
		if(snp.called)//genotypeCalled != genotypeRef) 
			file << '\t' << (char)options.toIupac[(unsigned)(snp.genotype&15)]<< '\t' << snp.quality;
		else 
			file << "\t\t"; 
	}
	else
	{
		if(snp.called) 
			file  << '\t' << (Dna)(snp.genotype & 3) << '\t' << snp.quality;// mutAllele;// << "/" << (Dna)mutAllele;
		else file << "\t\t";
	}
	file << std::endl;
	return true;
}



//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
	typename TFragmentStore,
	typename TReadCounts,
	typename TReadCigars,
	typename TGenomeName,
	typename TFile,
	typename TOptions
>
void dumpVariantsRealignBatch(
	TFragmentStore				&fragmentStore,				// forward/reverse matches
	TReadCigars				&,
	TReadCounts const			&,
	TGenomeName const	 		genomeID,					// genome name
	typename TFragmentStore::TContigPos	startCoord,			// startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
	typename TFragmentStore::TContigPos	currStart,
	typename TFragmentStore::TContigPos	currEnd,
	TFile					&file,
	TFile					&indelfile,
	TOptions 				&options)
{

	typedef typename TFragmentStore::TAlignedReadStore 	TMatches;
	typedef typename Value<TMatches>::Type 				TMatch;
	typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
	typedef typename Value<TMatchQualities>::Type 		TMatchQuality;
	typedef typename Iterator<TMatches,Standard>::Type	TMatchIterator;

	typedef typename TFragmentStore::TReadSeqStore	 	TReads;
	typedef typename Value<TReads>::Type	 			TRead;

	typedef typename TFragmentStore::TContigStore		TContigStore;
	typedef typename Value<TContigStore>::Type			TContig;
	typedef typename TFragmentStore::TContigPos 		TContigPos;
	typedef typename TFragmentStore::TContigSeq 		TContigSeq;

	typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
	typedef Gaps<TRead, AnchorGaps<typename TMatch::TGapAnchors> >			TReadGaps;
	typedef typename Iterator<TContigGaps>::Type							TContigGapIter;
	typedef typename Iterator<TReadGaps>::Type								TReadGapIter;
	
	
	SEQAN_PROTIMESTART(dump_time);
	
	// matches need to be ordered according to genome position
	TReads &reads					= fragmentStore.readSeqStore;
	TMatches &matches				= fragmentStore.alignedReadStore;
	TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
	TContigPos genomeLen				= (TContigPos)length(fragmentStore.contigStore[0].seq);

	::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPos<TMatch>());		

	// make sure both output files are open
	if (!file.is_open() && !indelfile.is_open()) 
	{
		::std::cerr << "Neither SNP nor indel output file is open" << ::std::endl;
		return;
	}

	// log file business
	::std::ofstream logfile;
	if(*options.outputLog != 0)
	{
		logfile.open(options.outputLog, ::std::ios_base::out | ::std::ios_base::app);
		if (!logfile.is_open()) 
			::std::cerr << "Failed to write to log file" << ::std::endl;
		logfile << "#stats for window " << currStart << " " << currEnd << " of " << genomeID << std::endl;
	}

	if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << genomeID << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;

#ifdef SNPSTORE_DEBUG
	bool extraV = true;
	
	std::cout << genomeLen << " <-length genome \n";
	std::cout << length(fragmentStore.alignedReadStore) << " <-nummatches \n";
	std::cout << length(fragmentStore.readSeqStore) << " <-numreads \n";
	CharString str = "realignBatch";
//	_dumpMatches(fragmentStore, str);
	std::cout << "startcoord=" << startCoord << std::endl;
#endif

//	std::fstream tmpfile;
//  	tmpfile.open("Z:/seqan071010/projects/library/apps/chr4.before.sam", ::std::ios_base::out);
//  	write(tmpfile, fragmentStore, Sam());
//  	tmpfile.close()
;
	// convert matches in fragmentstore into a globally consistent alignment
	//Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	Score<int> scoreType = Score<int>(0, -1, -2, -10);	// (match, mismatch,gapExtend,gapOpen) 
	convertMatchesToGlobalAlignment(fragmentStore, scoreType, Nothing());

//	std::fstream tmpfile2;
//	tmpfile2.open("Z:/seqan071010/projects/library/apps/chr4.after.sam", ::std::ios_base::out);
//	write(tmpfile2, fragmentStore, Sam());
//	tmpfile2.close();

#ifdef SNPSTORE_DEBUG
	if(extraV)
	{
		TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
		TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
		maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
		std::cout << "maxPos visual = " << maxPos;
		std::cout << " genomeLen = " << genomeLen << std::endl;

		AlignedReadLayout layout;
		layoutAlignment(layout, fragmentStore);
		printAlignment(std::cout, Raw(), layout, fragmentStore, 0, (TContigPos)0, (TContigPos)maxPos, 0, 150);
	}
	::std::cout << "done.\n" << std::flush;
	if(extraV)
	{
		CharString strstr = "befReal";
		_dumpMatches(fragmentStore,strstr);
	}
#endif

	//todo: do realignment here if options.realign, else work on global alignment
//	if(options.realign)
//	{
	// TODO: check out if a different scoring scheme makes more sense
	Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > consScore;
	int bandWidth = 10; // ad hoc

	if(options._debugLevel > 1)
	{
		::std::cout << "Realigning "<< length(matches)<<" reads to genome of length " <<genomeLen << std::flush;
//		::std::cout << " StartCoord="<< startCoord << std::endl;
	}


#ifndef TRACE_PIPELINE
	bandWidth = 25; // allow for some more bandWidth in general
#endif

#ifdef READS_454
	reAlign(fragmentStore,consScore,0,1,bandWidth,true);
#else
	reAlign(fragmentStore,consScore,0,1,bandWidth,true);
#endif



#ifdef SNPSTORE_DEBUG
	if(extraV)
	{
		CharString strstr = "befRefRal";
		_dumpMatches(fragmentStore,strstr);
		TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
		TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
		maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
		std::cout << "maxPos visual = " << maxPos;
		std::cout << " genomeLen = " << genomeLen << std::endl;
			AlignedReadLayout layout;
			layoutAlignment(layout, fragmentStore);
			printAlignment(std::cout, Raw(), layout, fragmentStore, 0, (TContigPos)0, (TContigPos)maxPos, 0, 150);
		}
#endif		

		
	if(options._debugLevel > 1)::std::cout << "Realigning reads including reference..." << std::flush;

#ifdef READS_454
	reAlign(fragmentStore,consScore,0,1,/*bandWidth*/5,false);
#else
	reAlign(fragmentStore,consScore,0,1,/*bandWidth*/5,false);
#endif

	if(options._debugLevel > 1) ::std::cout << "Finished realigning." << std::endl;	



#ifdef SNPSTORE_DEBUG
	::std::cout << "Realignment done.\n";
	if(extraV)
	{
		TContigGaps contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
		TContigPos maxPos = positionSeqToGap(contigGaps,length(fragmentStore.contigStore[0].seq)-1)+1;
		maxPos = _max(maxPos,(TContigPos)length(fragmentStore.contigStore[0].seq));
		std::cout << "maxPos visual = " << maxPos;
		std::cout << " genomeLen = " << genomeLen << std::endl;
		AlignedReadLayout layout;
		layoutAlignment(layout, fragmentStore);
		printAlignment(std::cout, Raw(), layout, fragmentStore, 0, (TContigPos)0, (TContigPos)maxPos, 0, 150);
		std::fstream tmpfile3;
		tmpfile3.open("test.realigned.sam", ::std::ios_base::out);
		write(tmpfile3, fragmentStore, Sam());
		tmpfile3.close();
	}
	::std::cout << "done." << std::flush;

	//	std::fstream tmpfile2;
	//	tmpfile2.open("tmpfile_realigned.sam", ::std::ios_base::out);
	//	write(tmpfile2, fragmentStore, Sam());
	//	tmpfile2.close();

#endif
	//}

	
	// forward match qualities
	String<int> columnQualityF;			resize(columnQualityF,5);
	String<unsigned> countF;			resize(countF,5);
	String<CharString> qualityStringF;	resize(qualityStringF,5);

	// reverse match qualities
	String<int> columnQualityR;			resize(columnQualityR,5);
	String<unsigned> countR;			resize(countR,5);
	String<CharString> qualityStringR;	resize(qualityStringR,5);

	FunctorComplement<Dna5> f;
	
	// sort reads according to begin position
	sortAlignedReads(fragmentStore.alignedReadStore, SortBeginPos());
	TMatchIterator matchIt		= begin(matches, Standard());
	TMatchIterator matchItEnd	= end(matches, Standard());	
	unsigned numReads = length(matches)-1; // exclude reference sequence
	unsigned refId = length(matchQualities); // reference id (there may be more matchQs than matches due to pile up correction)

	// look for reference sequence and move it to the end of alignedreads
	// todo: only do this when realignment was done
	bool refFound = false;
	TMatchIterator matchItKeep = matchIt;
	TMatch tempRef;
	while(matchIt != matchItEnd)
	{
		if((*matchIt).readId == refId) // this is the reference
		{
			refFound = true;
			tempRef = *matchIt;
			matchItKeep = matchIt;
			++matchIt;
			continue;
		}
		if(refFound)
		{
			*matchItKeep = *matchIt; // matchItKeep lags behind by one match 
			++matchIt;++matchItKeep;
		}
		else ++matchIt;
	}
	*matchItKeep = tempRef;
	SEQAN_ASSERT_TRUE(refFound);
		
#ifdef SNPSTORE_DEBUG
	if(!refFound) ::std::cout << "ref not Found!\n";
	::std::cout << "done looking for ref." << std::flush << std::endl;
#endif

	matchIt		= begin(matches, Standard());
	matchItEnd	= end(matches, Standard());	
	matchItEnd--; // exclude reference sequence
	
	TRead		&reference = fragmentStore.readSeqStore[fragmentStore.alignedReadStore[numReads].readId]; // last read is reference sequence
	TReadGaps	referenceGaps(reference, fragmentStore.alignedReadStore[numReads].gaps);
	TContigPos      refStart = (TContigPos)fragmentStore.alignedReadStore[numReads].beginPos;
	TContigGaps	contigGaps(fragmentStore.contigStore[0].seq, fragmentStore.contigStore[0].gaps);
	SingleBaseVariant snp = {0,0,0,0,0};

#ifdef SNPSTORE_DEBUG
	if(extraV)
	{
		::std::cout << "lengthrefgaps=" << length(referenceGaps)<< std::endl;
		::std::cout << "length(genome) = " << genomeLen << " length(ref)=" << length(reference) << std::endl;
	}
#endif


	// for indels:
	// i1 keeps track of consensus character
	// i2 keeps track of coverage (last 8 bits) and indelcount (first 8 bits)
	String<Pair<short unsigned,short unsigned> > indelConsens;
	resize(indelConsens,refStart + length(referenceGaps));
	for(unsigned i = 0; i < refStart + length(referenceGaps) ; ++i)
	{
		indelConsens[i].i1 = 6;
		indelConsens[i].i2 = 0;
	}

	if(options._debugLevel>1) std::cout << "Start inspecting alignment..." << std::endl;
	// now walk through the reference sequence in gaps view space,
	// i.e. position may be a gap
	// example:
	// Ref      ACCGTGCACTAGCATCATT--ACTAGCATCATA
	// Reads    ACCGTACA--AGCATCAT
	//              TACA--AGCATCATT--ACT
	//                          ATTTTACTAGCATCATA
	for(TContigPos candidateViewPos = refStart; candidateViewPos < refStart + (TContigPos)length(referenceGaps); ++candidateViewPos)
	{
		// first check if reference has a gap (potential insertion in reads) at this position
		TContigGapIter refIt = iter(referenceGaps,candidateViewPos-refStart);
		bool refGap = false; 
		if(isGap(refIt)) refGap = true;
			
		//get position in sequence space
		TContigPos candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);
	
		// not in the current window yet
		if(candidatePos + startCoord < currStart) continue;
		// not in the current window anymore
		if(candidatePos + startCoord >= currEnd) break;

		Dna5 refBase = reference[candidatePos];	// what happens if refGap==true, esp. for leading gaps?
		if(refBase=='N') continue;
		
#ifdef SNPSTORE_DEBUG
		std::cout << "candidateViewPos = " << candidateViewPos << std::endl;
		std::cout << "candidatePos = " << candidatePos << std::endl;
		std::cout << "candidatePosMitStart = " << candidatePos + startCoord << " refBase = " << refBase << std::endl;
		if(refGap) std::cout << "refGap!" << std::endl;
		bool extraVVVV = false;
		if(candidatePos + startCoord == 19388258) extraVVVV=true;
#endif		

		//find range of relevant read matches 
		// CHECK: remove unnecessarily walking through same matches multiple times
		while(matchIt != matchItEnd &&  _max((*matchIt).endPos,(*matchIt).beginPos) <= candidateViewPos)
			++matchIt;
		TMatchIterator matchRangeBegin = matchIt;
		while(matchIt != matchItEnd &&  _min((*matchIt).endPos,(*matchIt).beginPos)  <= candidateViewPos)
			++matchIt;	
		TMatchIterator matchRangeEnd = matchIt; // could remember this for next round
		matchIt = matchRangeBegin;
		
		int coverage = matchRangeEnd-matchRangeBegin;
#ifdef SNPSTORE_DEBUG
		if(extraVVVV) std::cout <<"cov=" << coverage << std::endl;
#endif
		if(coverage<(int)options.minCoverage) 
			continue; // coverage too low
		
		// start checking reads for this position, prepare some helpers
		Dna5 candidateBase;
		int quality;
		std::set<int> readPosMap;

		for(unsigned t=0;t<5;++t) 
		{
			countF[t] = 0;
			columnQualityF[t] = 0;
			clear(qualityStringF[t]);

			countR[t] = 0;
			columnQualityR[t] = 0;
			clear(qualityStringR[t]);
		}

		bool observedAtLeastOneMut = false;
		int numIndelsObserved = 0;  // if refGap then this counts the number of insertions
						// else it counts the number of deletions
		int positionCoverage = 0;   // how many reads actually span the position?

		// now check reads
		while(matchIt != matchRangeEnd)
		{
			TContigPos currViewBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
			TContigPos currViewEnd = _max((*matchIt).beginPos,(*matchIt).endPos);

			// make sure this match is really spanning the position 
			if(!(currViewBegin <= candidateViewPos && candidateViewPos < currViewEnd))
			{
				++matchIt;
				continue;
			}
			++positionCoverage;

			char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';

			TReadGaps readGaps(reads[(*matchIt).readId],(*matchIt).gaps);
			TReadGapIter rgIt = iter(readGaps,candidateViewPos - currViewBegin);

			// check out which position is hit in this read
			int readPos;
			if(isGap(rgIt)) readPos = -1; //potential deletion in reads (insertion in reference)
			else
			{
				readPos = positionGapToSeq(readGaps,candidateViewPos - currViewBegin);
				if(orientation == 'R')
					readPos = length(reads[(*matchIt).readId]) - readPos - 1;
			}

#ifdef SNPSTORE_DEBUG
			std::cout << "ReadPos = " << readPos << std::endl;
#endif

			if(readPos != -1) //-1 indicates gap in read
			{
				if(options.minDifferentReadPos > 0)
					if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
							(unsigned) readPos >= options.excludeBorderPos )
						readPosMap.insert(readPos);

				if(orientation == 'R') candidateBase = f((Dna5)reads[(*matchIt).readId][readPos]);
				else candidateBase = (Dna5)reads[(*matchIt).readId][readPos];
				
				if(refGap) ++numIndelsObserved; // count insertions
				else if(candidateBase != refBase) observedAtLeastOneMut = true;
				quality = getQualityValue(reads[(*matchIt).readId][readPos]);

				if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)
				{	// dont trust the quality of this position more
					// than the average quality of this read
					quality = (int) matchQualities[(*matchIt).id].score;
				}
				
				if(orientation == 'F')
				{
					columnQualityF[ordValue(candidateBase)] += quality;
					++countF[ordValue(candidateBase)];
					appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33),Generous());
				}
				else
				{
					columnQualityR[ordValue(candidateBase)] += quality;
					++countR[ordValue(candidateBase)];
					appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
				}
			}
			else
			{	//potential deletions
#ifdef SNPSTORE_DEBUG
				if(extraVVVV) std::cout <<"readPos = -1 " << std::endl;
#endif

				if(!refGap)
					++numIndelsObserved;
				
			}

			++matchIt;
		}
		matchIt = matchRangeBegin; //set iterator back to where we started from, same matches might be involved in next cand pos
		
#ifdef SNPSTORE_DEBUG
		if(extraVVVV)
		{
			std::cout << "posCov=" << positionCoverage << "numIndels = " << numIndelsObserved << std::endl;
			if(observedAtLeastOneMut) std::cout << "observed at least one mut " << std::endl;
		}
#endif

		// too few reads actually cover the position
		if(positionCoverage < (int)options.minCoverage) 
			continue;

		//all observed bases match the reference allele or there were too few indels
		//if(!observedAtLeastOneMut && numIndelsObserved< options.indelCountThreshold) 
		//	continue; 

	
		// do SNP calling
		if(file.is_open() && observedAtLeastOneMut) 
		{
			bool isSnp = true;

			// coverage depth	
			int refAllele = ordValue(reference[candidatePos]);
			unsigned realCoverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
			unsigned realCoverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
			unsigned realCoverage  = realCoverageF + realCoverageR;

			// Coverage too low after discarding Ns and gaps from alignment column
			if(realCoverage<options.minCoverage) isSnp = false;

			// is the min. number of different read positions supporting the mutation met?
			if(isSnp && options.minDifferentReadPos > 0 && readPosMap.size() < options.minDifferentReadPos)
				isSnp = false;

			// do genotype calling
			if(isSnp && options.method == 1)
				isSnp = _doSnpCall(countF,countR,qualityStringF,qualityStringR,refAllele,options,snp,MaqMethod()
#ifdef SNPSTORE_DEBUG_CANDPOS
				,(int) candidatePos + startCoord
#endif
				);
			else if(isSnp && options.method == 0)
				isSnp = _doSnpCall(countF,countR,columnQualityF,columnQualityR,refAllele,options,snp,ThresholdMethod() );
			
			// write SNP to file
			if(isSnp && (snp.called || options.outputFormat == 0))
				_write(file,snp,qualityStringF,qualityStringR,refAllele,genomeID,candidatePos+startCoord,realCoverage,options);
			
		}

		// do indel calling
		if (indelfile.is_open() && numIndelsObserved >= (int)options.indelCountThreshold 
			&& ((float)numIndelsObserved/(float)positionCoverage) >= options.indelPercentageT)
		{
			char mostCommonBase = 5; // 5 represents gap char "-", potential deletion
			if(refGap) // potential insertion
			{
#ifdef SNPSTORE_DEBUG
				if(extraVVVV) std::cout << "potential insertion" << std::endl;
#endif
				SEQAN_ASSERT_TRUE(!observedAtLeastOneMut);
				mostCommonBase = 0;
				unsigned maxCount = countF[0] + countR[0];
				for(unsigned j = 0; j < length(countF); ++j)
					if(countF[j] + countR[j] > maxCount)
					{
						maxCount = countF[j] + countR[j];
						mostCommonBase = j;
					}
			}
			indelConsens[candidateViewPos].i1 = mostCommonBase;
#ifdef SNPSTORE_DEBUG
				if(extraVVVV) std::cout << "mosCommonBase = " << (int)mostCommonBase << std::endl;
#endif
			if(positionCoverage > 255) //downscaling if numbers get too large
			{
				numIndelsObserved *= (int)((float)255.0/(float)positionCoverage);
				positionCoverage = 255;
#ifdef SNPSTORE_DEBUG
				if(extraVVVV) std::cout << "downscaled to " << numIndelsObserved << std::endl;
#endif

			}
			indelConsens[candidateViewPos].i2 = numIndelsObserved << 8 | positionCoverage;
		}

	
	}

	CharString chrPrefix = "";

#ifdef TRACE_PIPELINE
	chrPrefix = "chr";
#endif

	::std::string runID = options.runID;

#ifdef SNPSTORE_DEBUG
	// write out indels
	for(unsigned i = refStart; i < refStart + length(referenceGaps); ++i)
		std::cout << indelConsens[i].i1;
	std::cout << std::endl;
#endif

	if(indelfile.is_open())	//indelcalling
	{
		if(options._debugLevel > 1) std::cout << "Calling indels..." << std::endl;
		TContigPos candidateViewPos = refStart;
		Dna5String insertionSeq;
		while(candidateViewPos < refStart + (TContigPos)length(referenceGaps))
		{

			if(candidateViewPos < refStart + (TContigPos)length(referenceGaps) &&
				indelConsens[candidateViewPos].i1==6) // not a relevant position
			{
#ifdef SNPSTORE_DEBUG
				::std::cout << candidateViewPos << "not relevant for indels" <<  std::endl;
#endif
				++candidateViewPos;
				continue;
			}
		
			//get position in sequence space
			TContigPos candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);
			int indelSize = 0;
			unsigned depth = 0;
			float percentage = 0.0;
			// gap position
#ifdef SNPSTORE_DEBUG
			::std::cout << candidateViewPos << " indel?" <<  std::endl;
#endif
			while(candidateViewPos < refStart + (TContigPos)length(referenceGaps) && // shouldnt happen actually
				(indelConsens[candidateViewPos].i1==5 || 		// deletion in consens
				(indelConsens[candidateViewPos].i1==6 && isGap(referenceGaps, candidateViewPos-refStart)))	// position in consensus is the same as in reference
				) 													// and reference is a gap (same candidatePosition as before)
			{
#ifdef SNPSTORE_DEBUG
				::std::cout << startCoord + candidateViewPos << " del!!" <<  std::endl;
#endif
				if(indelConsens[candidateViewPos].i1==5) 
				{
					++indelSize;
					depth += (indelConsens[candidateViewPos].i2 & 255);
					percentage += (float)((indelConsens[candidateViewPos].i2 >> 8) & 255);
				}
				++candidateViewPos;
				
			}
			if(indelSize>0)
			{
				percentage = percentage/(float)depth; // low coverage positions get a lower weight here
				depth = (unsigned)round(depth/indelSize);   // coverage is spread over all positions
				//print deletion
				indelfile << chrPrefix << genomeID << '\t' << runID << "\tdeletion\t";
				indelfile << candidatePos + startCoord + options.positionFormat  << '\t';
				indelfile << candidatePos + startCoord + options.positionFormat + indelSize - 1;
				indelfile << "\t" << percentage;
				indelfile << "\t+\t.\tID=" << candidatePos + startCoord + options.positionFormat ;
				indelfile << ";size=" << indelSize;
				indelfile << ";depth=" << depth;
				//if(splitSupport>0) indelfile << ";splitSupport=" << splitSupport;
				indelfile << std::endl;
		
		
				//reset		
				candidatePos = positionGapToSeq(referenceGaps, candidateViewPos-refStart);
				indelSize = 0;
				depth = 0;
				percentage = 0.0;
			}
			clear(insertionSeq);
			while(candidateViewPos < refStart + (TContigPos)length(referenceGaps) && // shouldnt happen actually
				(indelConsens[candidateViewPos].i1<5 ||		     // insertion in consensus
				(indelConsens[candidateViewPos].i1==6 && candidatePos == positionGapToSeq(referenceGaps, candidateViewPos-refStart)))	// position in consensus is the same as in reference
				) 													// and reference is a gap (same candidatePosition as before)
			{
#ifdef SNPSTORE_DEBUG
				::std::cout << candidateViewPos << " ins!!!" <<  std::endl;
#endif
				if(indelConsens[candidateViewPos].i1<5)
				{
					--indelSize;
					depth += (indelConsens[candidateViewPos].i2 & 255);
					percentage += (float)((indelConsens[candidateViewPos].i2 >> 8) & 255);
					appendValue(insertionSeq,(Dna5)indelConsens[candidateViewPos].i1);
				}
				++candidateViewPos;
			}
			if(indelSize<0)
			{
				percentage = percentage/(float)depth; // low coverage positions get a lower weight here
				depth = (unsigned) round(depth/-indelSize);   // coverage is spread over all positions
		
				//print insertion
				indelfile << chrPrefix <<genomeID << '\t' << runID << "\tinsertion\t";
				indelfile << candidatePos + startCoord + options.positionFormat - 1 << '\t';
				indelfile << candidatePos + startCoord;// + options.positionFormat; //VORSICHT!!!
				indelfile << "\t" << percentage;
				indelfile << "\t+\t.\tID=" << candidatePos + startCoord + options.positionFormat;
				indelfile << ";size=" << indelSize;
				indelfile << ";seq="<< insertionSeq;
				indelfile << ";depth=" << depth;
				//if(splitSupport>0) indelfile << ";splitSupport=" << splitSupport;
				indelfile << std::endl;
		
				//resetting will be done in next round
			}
		
		}
		if(options._debugLevel > 1) std::cout << "Finished calling indels..." << std::endl;

	}
	
	if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;

	if((*options.outputLog != 0) && logfile.is_open())
		logfile.close();
		
	
}




//////////////////////////////////////////////////////////////////////////////
// Output SNPs
template <
	typename TFragmentStore,
	typename TReadCigars,
	typename TReadCounts,
	typename TGenomeName,
	typename TFile,
	typename TOptions
>
void dumpSNPsBatch(
	TFragmentStore				&fragmentStore,				// forward/reverse matches
	TReadCigars				&,
	TReadCounts const			&readCounts,
	TGenomeName const 			genomeID,					// genome name
	typename TFragmentStore::TContigPos	startCoord,			// startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
	typename TFragmentStore::TContigPos	currStart,
	typename TFragmentStore::TContigPos	currEnd,
	TFile				&file,
	TOptions 			&options)
{

	typedef typename TFragmentStore::TAlignedReadStore 	TMatches;
	typedef typename Value<TMatches>::Type 				TMatch;
	typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
	typedef typename Value<TMatchQualities>::Type 		TMatchQuality;
	typedef typename TFragmentStore::TReadSeqStore	 	TReads;
	typedef typename Value<TReads>::Type 				TRead;
	typedef typename TFragmentStore::TContigPos 		TContigPos;
	typedef typename TFragmentStore::TContigSeq 		TContigSeq;
	typedef typename Iterator<TMatches,Standard>::Type	TMatchIterator;
	
	SEQAN_PROTIMESTART(dump_time);
	//options._debugLevel = 2;
	String<char> toIupac = "AMRWMCSYRSGKWYKT";
	//std::cout << "Hier\n";
	// matches need to be ordered accordign to genome position
	TReads &reads					= fragmentStore.readSeqStore;
	TMatches &matches				= fragmentStore.alignedReadStore;
	TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
	TContigSeq &genome				= fragmentStore.contigStore[0].seq;
	
	::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPos<TMatch>());		
	
	Align<String<Dna5>, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	resize(rows(align), 2);

	if (!file.is_open()) 
	{
		::std::cerr << "Output file is not open" << ::std::endl;
		return;
	}

	::std::ofstream logfile;
	if(*options.outputLog != 0)
	{
		logfile.open(options.outputLog, ::std::ios_base::out | ::std::ios_base::app);
		if (!logfile.is_open()) 
			::std::cerr << "Failed to write to log file" << ::std::endl;
		logfile << "#stats for window " << currStart << " " << currEnd << " of " << genomeID << std::endl;
	}
	
	TMatchIterator matchIt	= begin(matches, Standard());
	TMatchIterator matchItEnd	= end(matches, Standard());	
	//matchItEnd--;
	unsigned countLowerMQ = 0, countHigherMQ = 0;
	
	if(options._debugLevel > 1) ::std::cout << "Scanning chromosome " << genomeID << " window (" << currStart<<","<< currEnd << ") for SNPs..." << ::std::endl;
	
	// forward match qualities
	String<int> columnQualityF;			resize(columnQualityF,5);
	String<unsigned> countF;			resize(countF,5);
	String<CharString> qualityStringF;	resize(qualityStringF,5);
	
	// reverse match qualities
	String<int> columnQualityR;			resize(columnQualityR,5);
	String<unsigned> countR;			resize(countR,5);
	String<CharString> qualityStringR;	resize(qualityStringR,5);
	
	// both
	String<unsigned> count;				resize(count,5);
	String<unsigned> columnQuality;		resize(columnQuality,5);
	
#ifdef SNPSTORE_DEBUG
	bool extraV = false;
#endif
	SingleBaseVariant snp = {0,0,0,0,0};

	for(TContigPos candidatePos = 0; candidatePos < (TContigPos)length(genome); ++candidatePos)
	{
		if(options._debugLevel > 1) ::std::cout << "Next pos\n";
		
		if(candidatePos + startCoord < currStart) continue;

		// not in the current window anymore
		if(candidatePos + startCoord >= currEnd)
			break;

		Dna5 refBase = genome[candidatePos];
		if(refBase=='N') continue;

#ifdef SNPSTORE_DEBUG
		::std::cout << "candPos=" << candidatePos + startCoord << ::std::endl;
		if(candidatePos + startCoord == 798634) 
			::std::cout << "ab jetzt.." << ::std::flush;
#endif
		
		Dna5 candidateBase;
		int quality;
		
		if(options._debugLevel > 1)std::cout << candidatePos+startCoord << "<-candidatePos\n";
		for(unsigned t=0;t<5;++t) 
		{
			countF[t] = 0;
			columnQualityF[t] = 0;
			clear(qualityStringF[t]);

			countR[t] = 0;
			columnQualityR[t] = 0;
			clear(qualityStringR[t]);
		}

		//find range of relevant read matches 
		while(matchIt != matchItEnd &&  _max((*matchIt).endPos,(*matchIt).beginPos) <= candidatePos)
			++matchIt;
		TMatchIterator matchRangeBegin = matchIt;
		while(matchIt != matchItEnd &&  _min((*matchIt).endPos,(*matchIt).beginPos)  <= candidatePos)
			++matchIt;
		TMatchIterator matchRangeEnd = matchIt;
		matchIt = matchRangeBegin;

		int coverage = matchRangeEnd-matchRangeBegin;
		if(coverage<(int)options.minCoverage) continue; // coverage too low

		if(options._debugLevel > 1)::std::cout << "Match range:" << matchRangeEnd - matchRangeBegin << ::std::endl;
#ifdef SNPSTORE_DEBUG
		if(extraV)
		{
			for (TMatchIterator tempIt = matchRangeBegin; tempIt != matchRangeEnd; ++tempIt)
				::std::cout << reads[(*tempIt).readId]<<"\n";
		}
#endif
		std::set<unsigned> readPosMap;
		bool observedAtLeastOneMut = false;
	
		while(matchIt != matchRangeEnd)
		{
			TContigPos currentBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
			TContigPos currentEnd	= _max((*matchIt).beginPos,(*matchIt).endPos);
			char orientation = ((*matchIt).beginPos > (*matchIt).endPos) ? 'R' : 'F';
			
#ifdef SNPSTORE_DEBUG
			if(extraV)
			{
				::std::cout <<"currentBegin = "<<currentBegin << "\n";
				::std::cout <<"currentEnd = "<<currentEnd << "\n";
			}
#endif			
			if(!(currentBegin <= candidatePos && candidatePos < currentEnd))// this match is not really spanning the position 
			{																// (can happen because of indels or variable-length reads)
				++matchIt;													
				continue;
			}
			/*if(!empty(readCigars[(*matchIt).readId]))//splitRead, dont use for snp calling for now
			{ 
				++matchIt; 
				continue;
			}*/
			// or do use and only :
			// if(candidatePos lies in split gap) continue;
			
			// do edit alignment
			if((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/) // splitReads: hamming: pairScore=0
			{
				Dna5String gInf = infix(genome, currentBegin, currentEnd);
				if (orientation == 'R')
					reverseComplement(gInf);

				assignSource(row(align, 0), reads[(*matchIt).readId]);
				assignSource(row(align, 1), gInf);
				globalAlignment(align, scoreType);	//splitReads: get alignment from cigar string
			}

			if (orientation == 'R')
			{
				FunctorComplement<Dna5> f;
				
				int readPos = currentEnd - candidatePos - 1;
				if ((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/)
					readPos = getReadPos(align,readPos,false); // 
				
#ifdef SNPSTORE_DEBUG
				if(extraV) std::cout << "readPosNacher = " << readPos << std::endl;
#endif				
				if(readPos != -1) //-1 indicates gap
				{
					if(options.minDifferentReadPos > 0)
						if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
							(unsigned) readPos >= options.excludeBorderPos )
							readPosMap.insert(readPos);
					candidateBase = f((Dna5)reads[(*matchIt).readId][readPos]);
#ifdef SNPSTORE_DEBUG
					if(extraV) std::cout << candidateBase << "candBase\n";
#endif
					quality = getQualityValue(reads[(*matchIt).readId][readPos]);
					if(candidateBase != refBase) observedAtLeastOneMut = true;

					if(!options.useBaseQuality && quality > (int)matchQualities[(*matchIt).id].score)	// dont trust the quality of this position more
					{																				// than the average quality of this read
						quality = (int) matchQualities[(*matchIt).id].score;
						++countLowerMQ;
					}
					else ++countHigherMQ;
					//if(quality < 0 || quality > 40)::std::cout << "falschQ candPos = " << candidatePos + startCoord << std::endl;

					unsigned tmpCount = 1;
					if(!empty(readCounts)) tmpCount = readCounts[(*matchIt).readId];
					for (unsigned k = 0; k < tmpCount; ++k)
					{
						columnQualityR[ordValue(candidateBase)] += quality;
						++countR[ordValue(candidateBase)];
						appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33),Generous());
					}
				}
			}
			else
			{
				int readPos = candidatePos - currentBegin;
				
				if((int)matchQualities[(*matchIt).id].pairScore == 1 /*!empty((*matchIt).gaps)*/)
					readPos = getReadPos(align,readPos);
			
				if(readPos != -1) //-1 indicates gap
				{
					if(options.minDifferentReadPos > 0)
						if((unsigned)(length(reads[(*matchIt).readId]) - readPos) > options.excludeBorderPos  &&
							(unsigned) readPos >= options.excludeBorderPos )
							readPosMap.insert(readPos);
							
					candidateBase = (Dna5)reads[(*matchIt).readId][readPos];
					quality = getQualityValue(reads[(*matchIt).readId][readPos]);
					if(candidateBase != refBase) observedAtLeastOneMut = true;

					if(!options.useBaseQuality && quality > (int) matchQualities[(*matchIt).id].score) 
					{
						quality = (int) matchQualities[(*matchIt).id].score;
						++countLowerMQ;
					}
					else ++countHigherMQ;
					//if(quality < 0 || quality > 40)::std::cout << "falschQ candPos = " << candidatePos + startCoord << std::endl;

					unsigned tmpCount = 1;
					if(!empty(readCounts)) tmpCount = readCounts[(*matchIt).readId];
					for (unsigned k = 0; k < tmpCount; ++k)
					{
						++countF[ordValue(candidateBase)];
						columnQualityF[ordValue(candidateBase)] += quality;
						appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33),Generous());
					}
				}
			}
			++matchIt;
		}
		matchIt = matchRangeBegin; //set iterator back to where we started from, same matches might be involved in next cand pos

		if(!observedAtLeastOneMut) continue; //all observed bases match the reference allele
	
		// do SNP calling
		bool isSnp = true;

		// coverage depth	
		int refAllele = ordValue(genome[candidatePos]);
		unsigned realCoverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
		unsigned realCoverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
		unsigned realCoverage  = realCoverageF + realCoverageR;

		// Coverage too low after discarding Ns and gaps from alignment column
		if(realCoverage<options.minCoverage) isSnp = false;

		// is the min. number of different read positions supporting the mutation met?
		if(isSnp && options.minDifferentReadPos > 0 && readPosMap.size() < options.minDifferentReadPos)
			isSnp = false;

		// do genotype calling
		if(isSnp && options.method == 1)
			isSnp = _doSnpCall(countF,countR,qualityStringF,qualityStringR,refAllele,options,snp,MaqMethod()
#ifdef SNPSTORE_DEBUG_CANDPOS
				,(int) candidatePos + startCoord
#endif
			);
		else if(isSnp && options.method == 0)
			isSnp = _doSnpCall(countF,countR,columnQualityF,columnQualityR,refAllele,options,snp,ThresholdMethod() );
		
		// write SNP to file
		if(isSnp && (snp.called || options.outputFormat == 0))
			_write(file,snp,qualityStringF,qualityStringR,refAllele,genomeID,candidatePos+startCoord,realCoverage,options);
			
	}

	if(options._debugLevel>1) std::cout <<"Finished scanning window.\n"<<std::flush;

	if((*options.outputLog != 0) && logfile.is_open())
		logfile.close();
		

	return;

}




template<typename TAlign, typename TString, typename TPosition>
void
getIndels(TAlign & align,TString &insertions, TString &deletions, TPosition begin_, TPosition end_)
{
	
	typedef typename Source<TAlign>::Type TSource;
	typedef typename Iterator<TSource, Rooted>::Type TStringIterator;
	
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;
	
	TAlignIterator ali_it0 = iter(row(align,0),begin_);
	TAlignIterator ali_it1 = iter(row(align,1),begin_);					
	TAlignIterator ali_it0_stop = iter(row(align,0),end_);
	TAlignIterator ali_it1_stop = iter(row(align,1),end_);
	TStringIterator readBase = begin(source(row(align,0))); 
	//std::cout << "getting cigar line\n";//ali0 len = " <<ali_it0_stop-ali_it0 << " \t ali1 len = "<<ali_it1_stop-ali_it1<<"\n";
	int readPos = 0;
	int refPos = 0;
	while(ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
	{
		int inserted = 0;
		int deleted = 0;
		while(ali_it0!=ali_it0_stop && ali_it1!=ali_it1_stop && !isGap(ali_it0)&& !isGap(ali_it1))
		{
			++readPos;
			++refPos;
			++ali_it0;
			++ali_it1;
		}
		while(ali_it0!=ali_it0_stop && isGap(ali_it0))
		{
			++refPos;
			++ali_it0;
			++ali_it1;
			++deleted;
		}
		if(deleted>0)
		{
			appendValue(deletions,Pair<int,Pair<int,int> >(refPos-deleted,Pair<int,int>(readPos,deleted)));
//			appendValue(deletions,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,deleted)));
		}
		while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
		{
			++ali_it0;
			++ali_it1;
			++readPos;
			++inserted;
		}
		if(inserted>0)
		{
			appendValue(insertions,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos-inserted,inserted)));
		}
	}
	
}



template <
	typename TFragmentStore,
	typename TReadCigars,
	typename TGenome,
	typename TGenomeName,
	typename TFile,
	typename TOptions
>
void dumpShortIndelPolymorphismsBatch(
	TFragmentStore 				&fragmentStore,				// forward/reverse matches
	TReadCigars				&readCigars,
	TGenome 				&genome,				// genome sequence
	TGenomeName const 			genomeID,				// genome name
	typename TFragmentStore::TContigPos	startCoord,			// startCoordinate + posOnGenomeInfix = real coordinate on whole chromosome
	typename TFragmentStore::TContigPos	currStart,
	typename TFragmentStore::TContigPos	currEnd,
	TFile					&indelfile,
	TOptions 				&options)
{
	
	typedef typename TFragmentStore::TAlignedReadStore 		TMatches;
	typedef typename TFragmentStore::TAlignQualityStore 		TMatchQualities;
	typedef typename Value<TMatches>::Type 				TMatch;
	typedef typename TFragmentStore::TReadSeqStore	 		TReads;
	typedef typename Value<TReads>::Type 				TRead;
	typedef typename Infix<TRead>::Type				TReadInf;
	typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;
	typedef typename TFragmentStore::TContigPos TContigPos;
	// matches need to be ordered accordign to genome position
	TReads &reads = fragmentStore.readSeqStore;
	TMatches &matches = fragmentStore.alignedReadStore;
	TMatchQualities &matchQualities = fragmentStore.alignQualityStore;
	::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPos<TMatch>());		

	if (!indelfile.is_open()) 
	{
		::std::cerr << "Failed to open indel output file" << ::std::endl;
		return;
	}

	TMatchIterator matchIt = begin(matches, Standard());
	TMatchIterator matchItEnd = end(matches, Standard());	
	
	Align<Dna5String, ArrayGaps> align;
//	Score<int> scoreType = Score<int>(1, -3, -11, -1);	//
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	resize(rows(align), 2);

	::std::string runID = options.runID;

	typedef Pair<unsigned,int> TPosLen;
	typedef typename std::map<TPosLen,Pair<unsigned,TReadInf>, LessPosLen<TPosLen> > TIndelMap;
	typedef typename TIndelMap::iterator TIndelIt;
	typedef typename std::map<TPosLen,unsigned,LessPosLen<TPosLen> >  TSplitMap;
	typedef typename TSplitMap::iterator TSplitIt;
	
	// position,length and count,sequence
	TIndelMap indels; //readinf empty for deletions

	//remember how many split reads supported the indel
	TSplitMap splitCounts;
	
	TIndelIt indelIt;
	TSplitIt splitCountIt;
	CharString chrPrefix = "";

#ifdef TRACE_PIPELINE
	chrPrefix = "chr";
#endif

	TReadInf dummyInf;

#ifdef SNPSTORE_DEBUG
	bool extraV = true;
#endif
	
	// collect potential indels	
	for(;matchIt != matchItEnd; ++matchIt) 
	{
		if(matchQualities[(*matchIt).id].pairScore == 0 && empty(readCigars[(*matchIt).readId])) //if(length(cigar)>0)  dont skip!
			continue;
			
		String<Pair<int,Pair<int,int> > > readInserts;
		String<Pair<int,Pair<int,int> > > readDeletes;
		
		TRead& read = reads[(*matchIt).readId];
		int readLen = length(read);
		if(empty(readCigars[(*matchIt).readId]))// if this is not a split read --> do edit alignment
		{
#ifdef SNPSTORE_DEBUG
			if(extraV) ::std::cout << "read is edit indel mapped" << std::endl;
			if(extraV) ::std::cout << "read=" << read << " beg,end="<<(*matchIt).beginPos << ","<<(*matchIt).endPos <<::std::endl;
#endif
			assignSource(row(align, 0), reads[(*matchIt).readId]);
			assignSource(row(align, 1), infix(genome, _min((*matchIt).beginPos,(*matchIt).endPos), _max((*matchIt).beginPos,(*matchIt).endPos)));
			if ((*matchIt).beginPos > (*matchIt).endPos)
				reverseComplement(source(row(align, 0))); // check if reversing read is better for gap placement
			
			globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());
//			globalAlignment(align, scoreType, AlignConfig<false,false,false,false>(), Gotoh());
#ifdef SNPSTORE_DEBUG
			if(extraV) ::std::cout << align << std::endl;
#endif
			// transform first and last read character to genomic positions
			unsigned viewPosReadFirst  = toViewPosition(row(align, 0), 0);
			unsigned viewPosReadLast   = toViewPosition(row(align, 0), length(reads[(*matchIt).readId]) - 1);
			
			getIndels(align,readInserts,readDeletes, viewPosReadFirst,viewPosReadLast+1);

#ifdef SNPSTORE_DEBUG
			for (unsigned i = 0; i < length(readInserts); ++i)
				if(extraV) ::std::cout <<"ins: "<< readInserts[i].i1 << ","<<readInserts[i].i2.i1 <<","<< readInserts[i].i2.i2 << ::std::endl;
			for (unsigned i = 0; i < length(readDeletes); ++i)
				if(extraV) ::std::cout <<"del: "<<  readDeletes[i].i1 << ","<<readDeletes[i].i2.i1 <<","<< readDeletes[i].i2.i2 << ::std::endl;
#endif
		}
		else
		{
			
	//		if(extraV) ::std::cout << "read is split mapped" << std::endl;
			//this is where i have to get rid of adjacent insertions/deletions in edit-split-mapped reads
			typename Value<TReadCigars>::Type &cigar = readCigars[(*matchIt).readId];
			int readPos = 0;    
			int refPos = 0;
			if((*matchIt).endPos > (*matchIt).beginPos)
			{
				for(unsigned i = 0; i < length(cigar); ++i)
				{
					if(cigar[i].i1 == 'D') //deletion
					{
						appendValue(readDeletes,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
						//::std::cout << " "<<cigar[i].i2 << " d at refPos " << refPos ;
						refPos += cigar[i].i2;
					}
					if(cigar[i].i1 == 'I') //deletion
					{
						appendValue(readInserts,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
						readPos += cigar[i].i2;
						//::std::cout << " "<<cigar[i].i2<< " i at refPos " << refPos ;
					}
					if(cigar[i].i1 == 'M') //matches
					{
						refPos += cigar[i].i2;
						readPos += cigar[i].i2;
						//::std::cout <<"  "<< cigar[i].i2<< " m " ;
					}
				}
//				::std::cout << std::endl;::std::cout << std::endl;
			}
			else{
				for(int i = length(cigar)-1; i >= 0; --i)
				{
					if(cigar[i].i1 == 'D') //deletion
					{
						appendValue(readDeletes,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
//						::std::cout << " "<<cigar[i].i2 << " d at refPos " << refPos ;
						refPos += cigar[i].i2;
					}
					if(cigar[i].i1 == 'I') //deletion
					{
						appendValue(readInserts,Pair<int,Pair<int,int> >(refPos,Pair<int,int>(readPos,cigar[i].i2)));
//						::std::cout << " "<<cigar[i].i2<< " i at refPos " << refPos << " readPos" << readPos << std::endl;
						readPos += cigar[i].i2;
					}
					if(cigar[i].i1 == 'M') //matches
					{
						refPos += cigar[i].i2;
						readPos += cigar[i].i2;
//						::std::cout <<"  "<< cigar[i].i2<< " m " ;
					}
				}
//				::std::cout << std::endl;::std::cout << std::endl;
			}
			//cigar[read][1].i1 => I/D? cigar[read][1].i2 = indelLen, cigar[read][0].i2 = readPos
		}
		
		//::std::cout << align;
		//unsigned refStart = min((*matchIt).beginPos,(*matchIt).endPos);
		for (unsigned i = 0; i < length(readInserts); ++i)
		{
			// go to genomic indel position
			unsigned indelCandPos;//
			if((*matchIt).beginPos > (*matchIt).endPos) 
				indelCandPos = (*matchIt).endPos + readInserts[i].i1;
			
			else indelCandPos = (*matchIt).beginPos + readInserts[i].i1;
#ifdef SNPSTORE_DEBUG
			if(extraV)  //62
				std::cout << "Pos=" << indelCandPos  + startCoord << " len=" <<  (readInserts[i].i2).i2 << std::endl;
#endif

			//TODO: make use of i2
			indelIt = indels.find(TPosLen((unsigned)indelCandPos,-(int)(readInserts[i].i2).i2));
			if(indelIt == indels.end())
			{
				// this is the first insertion with this length found at this genomic position 
				// --> remember with count 1 and also store indelsize and readInf
	//			if(extraV) ::std::cout << "new inds pos" << std::endl;

				TReadInf rInf;
				if((*matchIt).beginPos < (*matchIt).endPos)
					rInf = infix(read,
					(readInserts[i].i2).i1,
					(readInserts[i].i2).i1+(readInserts[i].i2).i2);
				else
					rInf = infix(read,
					readLen - (readInserts[i].i2).i1-(readInserts[i].i2).i2,
					readLen - (readInserts[i].i2).i1);
				if((*matchIt).beginPos > (*matchIt).endPos) reverseComplement(rInf);
				
				indels.insert(std::make_pair<Pair<unsigned,int>,Pair<unsigned,TReadInf> >
					(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2),
					 Pair<unsigned,TReadInf>(1,rInf)));
				
	//			if(extraV)std::cout << rInf << " <-" << (*matchIt).id<<std::endl;
			}
			else
			{
	//			if(extraV) ::std::cout << "increase counter ins pos" << std::endl;
				++(indelIt->second.i1);
/*				if (storeall insertion sequences)
				{
					
					TReadInf rInf;
					if((*matchIt).beginPos < (*matchIt).endPos)
						rInf = infix(read,
						(readInserts[i].i2).i1,
						(readInserts[i].i2).i1+(readInserts[i].i2).i2);
					else
						rInf = infix(read,
						readLen - (readInserts[i].i2).i1-(readInserts[i].i2).i2,
						readLen - (readInserts[i].i2).i1);
					if((*matchIt).beginPos > (*matchIt).endPos) 
						reverseComplement(rInf);
					if(extraV)std::cout << rInf << " <-" << (*matchIt).id<<std::endl;
					indelIt->second.i2 = rInf;
				}*/
				
			}
			if(!empty(readCigars[(*matchIt).readId]))// if this is a split read --> increase counter
			{
				splitCountIt = splitCounts.find(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2));
				if(splitCountIt == splitCounts.end())
				{
					splitCounts.insert(std::make_pair<Pair<unsigned,int>,unsigned>
					(Pair<unsigned,int>(indelCandPos,-(int)(readInserts[i].i2).i2),1));
				
				}
				else
					++(splitCountIt->second);
			}
			
		}
		for (unsigned i = 0; i < length(readDeletes); ++i)
		{
			unsigned indelCandPos;// = refStart + readDeletes[i].i1;
			if((*matchIt).beginPos > (*matchIt).endPos) 
			//	indelCandPos = (*matchIt).beginPos - readDeletes[i].i1 - (readDeletes[i].i2).i2;
				indelCandPos = (*matchIt).endPos + readDeletes[i].i1;
			else indelCandPos = (*matchIt).beginPos + readDeletes[i].i1;
#ifdef SNPSTORE_DEBUG
			if(extraV) 
				std::cout << "Pos=" << indelCandPos  + startCoord << " len=" <<  (readDeletes[i].i2).i2;
#endif

			//TODO: make use of i2
			indelIt = indels.find(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2));
			if(indelIt == indels.end())
			{
				indels.insert(std::make_pair<Pair<unsigned,int>,Pair<unsigned,TReadInf> >
					(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2),
					 Pair<unsigned,TReadInf>(1,dummyInf)));
			}
			else
				++(indelIt->second.i1);
			if(!empty(readCigars[(*matchIt).readId]))// if this is a split read --> increase counter
			{
				splitCountIt = splitCounts.find(Pair<unsigned,int>(indelCandPos,(readDeletes[i].i2).i2));
				if(splitCountIt == splitCounts.end())
				{
					splitCounts.insert(std::make_pair<Pair<unsigned,int>,unsigned>
					(Pair<unsigned,int>(indelCandPos,(int)(readDeletes[i].i2).i2),1));
				
				}
				else
					++(splitCountIt->second);
			}
			
		}
	}
	
	// now output all indels that meet the filter criteria 
	matchIt = begin(matches, Standard());
	while(matchIt != matchItEnd) 
	{
		unsigned currSeqNo = (*matchIt).contigId;
		TMatchIterator currSeqMatchItBegin = matchIt;
		while(matchIt != matchItEnd)
		{
			if ((*matchIt).contigId != currSeqNo) break;
			++matchIt;
		}
		TMatchIterator currSeqMatchItEnd = matchIt;
		
		matchIt = currSeqMatchItBegin;
		indelIt = indels.begin();
		splitCountIt = splitCounts.begin();
		TIndelIt endIt = indels.end();
		TSplitIt splitEndIt = splitCounts.end();
		
		//indel-merging, possibly suboptimal
		if(options.indelWindow > 0)
		{
			while(indelIt != endIt)
			{
				unsigned currPos = indelIt->first.i1;
				unsigned oriCurrPos = currPos;
				if(indelIt->second.i1 == 0) {++indelIt;continue;}
				TIndelIt nextIt = indelIt;
				++nextIt;
				//for all positions that are 
				while(nextIt != endIt && currPos + options.indelWindow > nextIt->first.i1 && oriCurrPos + 2*options.indelWindow > nextIt->first.i1 )
				{
					//add the number of found indel-reads to it
					if(indelIt->first.i2 == nextIt->first.i2)
					{
						if(indelIt->second.i1 < nextIt->second.i1 ) //nextIT has a higher count for that position --> add counts of indelIT to nextIT
						{
							nextIt->second.i1 += indelIt->second.i1;
							indelIt->second.i1 = 0;
							currPos = indelIt->first.i1;
						}
						else
						{
							indelIt->second.i1 += nextIt->second.i1;
							nextIt->second.i1 = 0;
						}
					}
					++nextIt;
				}
				++indelIt;
			}
			indelIt = indels.begin();
		}
		
		
/*		for(splitCountIt = splitCounts.begin(); splitCountIt != splitEndIt; ++splitCountIt)
		{
			std::cout << splitCountIt->first.i1 << ","  << splitCountIt->first.i2 << "," << splitCountIt->second << std::endl;
		}
		splitCountIt = splitCounts.begin();*/
		for(; indelIt != endIt; ++indelIt)
		{
			bool debug = false;
#ifdef SNPSTORE_DEBUG
			debug=true;
#endif
			int splitSupport = 0;
			if(splitCountIt != splitEndIt && 
				(splitCountIt->first.i1 == indelIt->first.i1) && (splitCountIt->first.i2 == indelIt->first.i2))
			{
				splitSupport = splitCountIt->second;
				++splitCountIt;
			}
			
			if(indelIt->second.i1 < options.indelCountThreshold)
			{
				if(debug)::std::cout << "indel: count too low "<<indelIt->second.i1<<"\n";
				continue;
			}
			if((TContigPos)indelIt->first.i1 + startCoord < currStart || (TContigPos)indelIt->first.i1 + startCoord >= currEnd)
			{
				if(debug)::std::cout << "indel: pos outside range "<<indelIt->first.i1<<"\n";
				continue;
			}
			unsigned candidatePos = indelIt->first.i1;
			while(matchIt != currSeqMatchItEnd && _max((*matchIt).endPos,(*matchIt).beginPos) <= (TContigPos) candidatePos)
				++matchIt;
			
			TMatchIterator matchRangeBegin = matchIt;
			while(matchIt != currSeqMatchItEnd && _min((*matchIt).endPos,(*matchIt).beginPos) <= (TContigPos) candidatePos)
				++matchIt;
			TMatchIterator matchRangeEnd = matchIt;

			int coverage = matchRangeEnd-matchRangeBegin;
			if(coverage<(int)options.minCoverage)
			{
				matchIt = matchRangeBegin;
				continue;
			
			}
			matchIt = matchRangeBegin;
			
			Dna5 refBase = genome[candidatePos];
			if(refBase=='N') continue;
			
			unsigned covF = 0; 
			unsigned covR = 0;
			
			while(matchIt != matchRangeEnd)
			{
				
				if(!(_min((*matchIt).beginPos,(*matchIt).endPos) <= (TContigPos)candidatePos 
					&& (TContigPos)candidatePos < _max((*matchIt).beginPos,(*matchIt).endPos)))
				{
					++matchIt;
					continue;
				}
			
				if((*matchIt).beginPos < (*matchIt).endPos) ++covF;
				else ++covR;
				++matchIt;
			}
			
			if((covF+covR) < options.minCoverage)
			{
				if(options._debugLevel > 0)
					::std::cout << "Coverage " << covF+covR << " after applying max pile filter and discarding Ns" << ::std::endl;
				matchIt = matchRangeBegin;
				continue;
			}
			
			if((float)indelIt->second.i1/(covF+covR) < options.indelPercentageT)
			{
				matchIt = matchRangeBegin;
				continue;
			}
			
			int indelSize=indelIt->first.i2;
		
				
			if(options.outputFormat < 2) //
			{
				if(indelSize > 0 )indelfile << chrPrefix << genomeID << '\t' << runID << "\tdeletion\t";
				else indelfile << chrPrefix <<genomeID << '\t' << runID << "\tinsertion\t";
				if(indelSize > 0 )indelfile << candidatePos + startCoord + options.positionFormat << '\t';
				else indelfile << candidatePos + startCoord + options.positionFormat - 1 << '\t';
				if(indelSize > 0 )indelfile << candidatePos + startCoord + options.positionFormat  + indelSize - 1;
				else indelfile << candidatePos + startCoord;// + options.positionFormat; //VORSICHT!!!
				indelfile << "\t" << (float)indelIt->second.i1/(covF+covR);
				indelfile << "\t+\t.\tID=" << candidatePos + startCoord + options.positionFormat;
				indelfile << ";size=" << indelSize;
				if(indelSize < 0)indelfile << ";seq="<<indelIt->second.i2;
				indelfile << ";depth=" << covF+covR;
				if(splitSupport>0) indelfile << ";splitSupport=" << splitSupport;
				indelfile << std::endl;
			}
			else
			{
				
				//chromosome
				indelfile << genomeID << '\t';
				indelfile <<  candidatePos + startCoord + options.positionFormat << '\t';
				if(options.orientationAware)
				{
					indelfile << covF  <<'\t';
					indelfile << covR  <<'\t';
				}
				else
				{
					indelfile << covF+covR  <<'\t';
				}
				indelfile << indelIt->second.i1 << std::endl;
			}

			matchIt = matchRangeBegin;
		}

		matchIt = currSeqMatchItEnd;
		if(options._debugLevel>1) std::cout <<"Finished scanning window for deletions.\n"<<std::flush;
	}


	return;

}





//do indel calling based on pronounced drops in coverage
template <
	typename TFragmentStore,
	typename TGenomeName,
	typename TFile,
	typename TOptions
>
void dumpCopyNumberPolymorphismsBatch(
	TFragmentStore 				&fragmentStore,				// forward/reverse matches
	TGenomeName const 			genomeID,				// genome name
	typename TFragmentStore::TContigPos	startCoord,
	typename TFragmentStore::TContigPos	currStart,
	typename TFragmentStore::TContigPos	currEnd,
	TFile			&file,
	TOptions 		&options)
{
	typedef typename TFragmentStore::TAlignedReadStore 		TMatches;
	typedef typename Value<TMatches>::Type 				TMatch;
	typedef typename TFragmentStore::TReadSeqStore	 		TReads;
	typedef typename Value<TReads>::Type 				TRead;
	typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;

	// matches need to be ordered accordign to genome position
	TMatches &matches = fragmentStore.alignedReadStore;
	::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPos<TMatch>());		

	if (!file.is_open()) 
	{
		::std::cerr << "Failed to open cnv output file" << ::std::endl;
		return;
	}
	
	typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;

	TMatchIterator matchIt = begin(matches, Standard());
	TMatchIterator matchItEnd = end(matches, Standard());	

	::std::string runID = options.runID;

	
	//VORSICHT!! windowsize must be a multiple of cnvwindowsize !!!!!
	// bin matches according to their start! positions
	for(unsigned currBinStart = currStart;  currBinStart < currEnd; currBinStart += options.cnvWindowSize)
	{
		
		while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < currBinStart ) // havent reached bin begin yet
			++matchIt;
		
		unsigned count = 0;
		while(matchIt != matchItEnd && _min((*matchIt).beginPos,(*matchIt).endPos) < currBinStart + options.cnvWindowSize ) //count reads
		{
			++count;
			++matchIt;
		}
		CharString guess = "normal";
		if (count > options.expectedReadsPerBin + 3 *options.expectedReadsSD)
			guess = "insertion";
		else if (count < options.expectedReadsPerBin - 3 *options.expectedReadsSD)
			guess = "deletion";
		//if(guess != "normal" )
		//{
			file << genomeID << "\tcoverage\t"<< guess << "\t";
			file << currBinStart+ startCoord+1  << "\t" << currBinStart + startCoord + options.cnvWindowSize << "\t";
			file << count << "\t+\t.\t.\n";
		//}
	}


	return;

}




}

#endif
