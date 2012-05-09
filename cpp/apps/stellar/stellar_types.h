 /*==========================================================================
                     STELLAR - Fast Local Alignment

 ============================================================================
  Copyright (C) 2010 by Birte Kehr

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

#ifndef SEQAN_HEADER_STELLAR_TYPES_H
#define SEQAN_HEADER_STELLAR_TYPES_H

#include <seqan/align.h>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// Options for Stellar
struct StellarOptions {
	// i/o options
	CharString databaseFile;		// name of database file
	CharString queryFile;			// name of query file
	CharString outputFile;			// name of result file
	CharString disabledQueriesFile;	// name of result file containing disabled queries
	CharString outputFormat;		// Possible formats: gff, text

	// main options
	unsigned qGram;				// length of the q-grams
	double epsilon;				// maximal error rate
	int minLength;				// minimal length of an epsilon-match
	double xDrop;				// maximal x-drop

	// more options
	bool forward;				// compute matches to forward strand of database
	bool reverse;				// compute matches to reverse complemented database
	CharString fastOption;		// verification strategy: exact, bestLocal, bandedGlobal
	unsigned disableThresh;		// maximal number of matches allowed per query before disabling verification of hits for that query
	unsigned compactThresh;		// number of matches after which removal of overlaps and duplicates is started
	unsigned numMatches;		// maximal number of matches per query and database
	unsigned maxRepeatPeriod;	// maximal period of low complexity repeats to be filtered
	unsigned minRepeatLength;	// minimal length of low complexity repeats to be filtered
	double qgramAbundanceCut;
	char verbose;				// verbosity mode: 0 - low, 1 - medium (, 2 - high)


	StellarOptions() {
		outputFile = "stellar.gff";
		disabledQueriesFile = "stellar.disabled.fasta";
		outputFormat = "gff";

		qGram = (unsigned)-1;
		epsilon = 0.05;
		minLength = 100;
		xDrop = 5;

		forward = true;
		reverse = true;
		fastOption = "exact";		// exact verification
		disableThresh = (unsigned)-1;
		compactThresh = 500;
		numMatches = 50;
		maxRepeatPeriod = 1;
		minRepeatLength = 1000;
		qgramAbundanceCut = 1;
		verbose = 0;
	}
}; 


///////////////////////////////////////////////////////////////////////////////
// Container for storing local alignment matches of one query sequence
template<typename TMatch_>
struct QueryMatches {
	String<TMatch_> matches;
	bool disabled;

	QueryMatches() {
		disabled = false;
	}
};

///////////////////////////////////////////////////////////////////////////////
// Container for storing a local alignment match
template<typename TSequence_, typename TId_>
struct StellarMatch {
	typedef TSequence_							TSequence;
	typedef TId_								TId;
	typedef typename Position<TSequence>::Type	TPos;

	typedef Align<TSequence>					TAlign;
	typedef typename Row<TAlign>::Type			TRow;

	static const TId INVALID_ID;

	TPos begin1;
	TPos end1;
	TRow row1;

	TId id;
	TPos begin2;
	TPos end2;
	TRow row2;

	StellarMatch() {}

	template<typename TAlign, typename TId>
	StellarMatch(TAlign & _align, TId _id) {
		id = _id;
		row1 = row(_align, 0);
		row2 = row(_align, 1);

		begin1 = clippedBeginPosition(row1);
		end1 = clippedEndPosition(row1);

		begin2 = clippedBeginPosition(row2);
		end2 = clippedEndPosition(row2);
	}
};

template <typename TSequence, typename TId> 
const TId
StellarMatch<TSequence, TId>::INVALID_ID = "###########";


///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// to sort matches by position and remove overlapping matches
template <typename TMatch>
struct LessPos : public ::std::binary_function <TMatch, TMatch, bool> {		
	LessPos() {}
	
	inline int compare(TMatch const & a, TMatch const & b) const {
		// query number
		if ((a.id) < (b.id)) return -1;
		if ((a.id) > (b.id)) return 1;

		// database begin position
		typename TMatch::TPos aBegin1 = _min(a.begin1, a.end1);
		typename TMatch::TPos bBegin1 = _min(b.begin1, b.end1);
		if (aBegin1 < bBegin1) return -1;
		if (aBegin1 > bBegin1) return 1;

		// database end position
		typename TMatch::TPos aEnd1 = _max(a.begin1, a.end1);
		typename TMatch::TPos bEnd1 = _max(b.begin1, b.end1);
		if (aEnd1 < bEnd1) return -1;
		if (aEnd1 > bEnd1) return 1;

		// query begin position
		typename TMatch::TPos aBegin2 = _min(a.begin2, a.end2);
		typename TMatch::TPos bBegin2 = _min(b.begin2, b.end2);
		if (aBegin2 < bBegin2) return -1;
		if (aBegin2 > bBegin2) return 1;

		// query end position
		typename TMatch::TPos aEnd2 = _max(a.begin2, a.end2);
		typename TMatch::TPos bEnd2 = _max(b.begin2, b.end2);
		if (aEnd2 < bEnd2) return -1;
		if (aEnd2 > bEnd2) return 1;

		//// orientation
		//bool oa = a.begin1 < a.end1;
		//bool ob = b.begin1 < b.end1;
		//if (oa != ob) return oa;

		return 0;
	}
		
	inline bool operator() (TMatch const & a, TMatch const & b) const {
		return compare(a, b) == -1;
	}
};

///////////////////////////////////////////////////////////////////////////////
// to sort matches by length
template <typename TMatch>
struct LessLength : public ::std::binary_function <TMatch, TMatch, bool> {		
	LessLength() {}

	inline int compare(TMatch const & a, TMatch const & b) const {
		typename TMatch::TPos aLength = abs((int)a.end1 - (int)a.begin1);
		typename TMatch::TPos bLength = abs((int)b.end1 - (int)b.begin1);

		if (a.id == TMatch::INVALID_ID) return 1;
		if (b.id == TMatch::INVALID_ID) return -1;

		if (aLength < bLength) return 1;
		if (aLength > bLength) return -1;

		return 0;
	}

	inline bool operator() (TMatch const & a, TMatch const & b) const {
		return compare(a, b) == -1;
	}
};

///////////////////////////////////////////////////////////////////////////////
// Determines whether match1 is upstream of match2 in specified row.
//  If matches overlap, the non-overlapping parts have to be longer than minLenght.
template<typename TSequence, typename TId, typename TRowNo, typename TSize>
inline bool
_isUpstream(StellarMatch<TSequence, TId> & match1, StellarMatch<TSequence, TId> & match2, TRowNo row, TSize minLength) {
SEQAN_CHECKPOINT 
	typedef typename StellarMatch<TSequence, TId>::TPos TPos;

	TPos e1, b2;
	if (row == 0) {
		e1 = match1.end1;
		b2 = match2.begin1;
	} else {
		e1 = match1.end2;
		b2 = match2.begin2;
	}

    if (e1 <= b2) return true;
    
	TPos b1, e2;
	if (row == 0) {
		e2 = match2.end1;
		b1 = match1.begin1;
	} else {
		e2 = match2.end2;
		b1 = match1.begin2;
	}

    if (b1 < b2 && (b2 - b1 >= minLength)) {
        if ((e1 < e2) && (e2 - e1 >= minLength)) return true;
    }
    
    return false;
}

///////////////////////////////////////////////////////////////////////////////
// sorts StellarMatchees by specified functor
template <typename TMatches, typename TFunctorLess>
inline void
sortMatches(TMatches & stellarMatches, TFunctorLess const & less) {
	std::stable_sort(
		begin(stellarMatches, Standard()), 
		end(stellarMatches, Standard()), 
		less);
}

///////////////////////////////////////////////////////////////////////////////
// returns the length of the longer row from StellarMatch
template <typename TSequence, typename TId>
inline typename Size<TSequence>::Type
length(StellarMatch<TSequence, TId> & match) {
	return _max(length(match.row1), length(match.row2));
}

#endif
