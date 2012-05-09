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

#ifndef SEQAN_HEADER_STELLAR_OUTPUT_H
#define SEQAN_HEADER_STELLAR_OUTPUT_H

#include <iostream>
#include <seqan/align.h>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// Computes a CIGAR string and mutations from rows of StellarMatch.
template<typename TRow, typename TString>
void
_getCigarLine(TRow const & row0, TRow const & row1, TString & cigar, TString & mutations) { 
SEQAN_CHECKPOINT
    typedef typename Size<TRow>::Type TSize;

    TSize pos = 0;

    TSize dbEndPos = endPosition(row0);
    TSize queryEndPos = endPosition(row1);

    bool first = true;
    TSize readBasePos = pos;
    TSize readPos = 0;
	while (pos != dbEndPos && pos != queryEndPos) {
		int matched = 0;
		int inserted = 0;
		int deleted = 0;
		while (pos != dbEndPos && pos != queryEndPos &&
               !isGap(row0, pos) && !isGap(row1, pos)) {
            ++readPos;
			if (value(row0, pos) != value(row1, pos)) {
				if (first) first = false;
				else mutations << ",";
				mutations << readPos << value(source(row1), readBasePos);
			}
			++readBasePos;
			++pos;
			++matched;
		}
		if (matched > 0) cigar << matched << "M" ;
		while (pos != queryEndPos && isGap(row1, pos)) {
			++pos;
			++deleted;
		}
		if (deleted > 0) cigar << deleted << "D";
		while (pos != dbEndPos && isGap(row0, pos)) {
			++pos;
			++readPos;
			if (first) first = false;
			else mutations << ",";
			mutations << readPos << value(source(row1), readBasePos);
			++readBasePos;
			++inserted;
		}
		if (inserted > 0) cigar << inserted << "I";
	}
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the identity of two alignment rows (percentage of matching positions).
template<typename TRow>
double
_calculateIdentity(TRow const & row0, TRow const & row1) {
SEQAN_CHECKPOINT
    typedef typename Size<TRow>::Type TSize;
    TSize matches = 0;
    TSize len = _max(length(row0), length(row1));

    TSize pos = 0;

    TSize end0 = endPosition(row0);
    TSize end1 = endPosition(row1);

    while ((pos < end0) && (pos < end1)) {
        if (!isGap(row0, pos) && !isGap(row1, pos)) {
            if (value(row0, pos) == value(row1, pos)) {
                ++matches;
            }
        }
        ++pos;
    }

    return ((double)matches/(double)len)*100.0;
}

///////////////////////////////////////////////////////////////////////////////
// Writes rows of a StellarMatch in gff format to a file.
template<typename TId, typename TRow, typename TFile>
void
_writeMatchGff(TId const & databaseID,
              TId const & patternID,
              bool const databaseStrand,
              TRow const & row0,
              TRow const & row1,
              TFile & file) {
SEQAN_CHECKPOINT    
    for (typename Position<TId>::Type i = 0; i < length(databaseID) && value(databaseID, i) > 32; ++i) {
        file << value(databaseID, i);
    }

    file << "\tStellar";
    file << "\teps-matches";

    if (databaseStrand) {
        file << "\t" << 
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0)) + 1;
        file << "\t" << 
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
    } else {
        file << "\t" << length(source(row0)) - 
            (toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0))) + 1;
        file << "\t" << length(source(row0)) - 
            (toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0)));
    }

    file << "\t" << _calculateIdentity(row0, row1);

    file << "\t" << (databaseStrand ? '+' : '-');

    file << "\t.\t";
    for (typename Position<TId>::Type i = 0; i < length(patternID) && value(patternID, i) > 32; ++i) {
        file << value(patternID, i);
    }

	//file << ";seq2Length=" << length(source(row1));

    file << ";seq2Range=" << 
		toSourcePosition(row1, beginPosition(row1)) + beginPosition(source(row1)) + 1;
    file << "," << 
		toSourcePosition(row1, endPosition(row1)) + beginPosition(source(row1));

    std::stringstream cigar, mutations;
    _getCigarLine(row0, row1, cigar, mutations);
    file << ";cigar=" << cigar.str();
    file << ";mutations=" << mutations.str();
    file << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// Writes rows of a StellarMatch in human readable format to file.
template<typename TId, typename TRow, typename TFile>
void
_writeMatch(TId const & databaseID,
            TId const & patternID,
            bool const databaseStrand,
            TRow const & row0,
            TRow const & row1,
            TFile & file) {
SEQAN_CHECKPOINT
	// write database ID
	file << "Database sequence: " << databaseID;
	if (!databaseStrand) file << " (complement)" << std::endl;
	else file << std::endl;

	// write database positions
	file << "Database positions: ";
	if (databaseStrand) {
		file << toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0));
		file << ".." << toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
	} else {
		file << length(source(row0)) - toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0));
		file << ".." << length(source(row0)) - toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
	}
	file << std::endl;

	// write query ID
	file << "Query sequence: " << patternID << std::endl;

	// write query positions
	file << "Query positions: ";
	file << toSourcePosition(row1, beginPosition(row1)) + beginPosition(source(row1));
	file << ".." << toSourcePosition(row1, endPosition(row1)) + beginPosition(source(row1));
	file << std::endl;

	file << std::endl;

	// write match
	Align<typename Source<TRow>::Type> align;
	appendValue(align.data_rows, row0);
	appendValue(align.data_rows, row1);
	file << align;
	file << "----------------------------------------------------------------------\n" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Calls _writeMatchGff for each match in StringSet of String of matches.
//   = Writes matches in gff format to a file.
template<typename TInfix, typename TQueryId, typename TId, typename TIds, typename TMode, typename TFile,
         typename TString>
int
_outputMatches(StringSet<QueryMatches<StellarMatch<TInfix, TQueryId> > > const & matches,
			   TId const & databaseID,
			   bool const databaseStrand,
			   TIds const & ids,
			   TMode verbose,
			   TFile & fileName,
			   TString & format) {
SEQAN_CHECKPOINT
	typedef StellarMatch<TInfix, TQueryId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;

	TSize numMatches = 0;
	TSize totalLength = 0;
    TSize maxLength = 0;

	std::ofstream file;
	file.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::app);
	if (!file.is_open()) {
		std::cerr << "Could not open output file." << std::endl;
		return 1;
	}

	for (TSize i = 0; i < length(matches); i++) {
		QueryMatches<TMatch> queryMatches = value(matches, i);
		if (length(queryMatches.matches) == 0) continue;

		for (TSize j = 0; j < length(queryMatches.matches); j++) {
			TMatch m = value(queryMatches.matches, j);
			if (format == "gff")
				_writeMatchGff(databaseID, ids[i], databaseStrand, m.row1, m.row2, file);
			else 
				_writeMatch(databaseID, ids[i], databaseStrand, m.row1, m.row2, file);

			TSize len = _max(length(m.row1), length(m.row2));
			totalLength += len;
			if(len > maxLength) maxLength = len;
		}
		numMatches += length(queryMatches.matches);
	}

	file.close();

	if (verbose > 0) {
		if (numMatches > 0) {
			std::cout << "    Longest eps-match : " << maxLength << std::endl;
			std::cout << "    Avg match length  : " << totalLength / numMatches << std::endl;
		}
		std::cout << "    # Eps-matches     : " << numMatches << std::endl;
	}

	return numMatches;
}

///////////////////////////////////////////////////////////////////////////////
// Calls _writeMatchGff for each match in StringSet of String of matches.
//   = Writes matches in gff format to a file.
// Writes disabled query sequences to disabledFile.
template<typename TInfix, typename TQueryId, typename TId, typename TQueries, typename TIds,
         typename TMode, typename TFile, typename TString>
int
_outputMatches(StringSet<QueryMatches<StellarMatch<TInfix, TQueryId> > > const & matches, 
			   TId const & databaseID,
			   bool const databaseStrand,
			   TQueries & queries,
			   TIds const & ids,
			   TMode verbose,
			   TFile & fileName,
			   TString & format,
			   TString & disabledFile) {
SEQAN_CHECKPOINT
	typedef StellarMatch<TInfix, TQueryId> TMatch;
	typedef typename Size<typename TMatch::TAlign>::Type TSize;

    TSize maxLength = 0;
    TSize totalLength = 0;
    TSize numMatches = 0;
    TSize numDisabled = 0;

    std::ofstream daFile, file;
	file.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::app);
	if (!file.is_open()) {
		std::cerr << "Could not open output file." << std::endl;
		return 1;
	}

	daFile.open(toCString(disabledFile), ::std::ios_base::out | ::std::ios_base::app);
	if (!daFile.is_open()) {
		std::cerr << "Could not open file for disabled queries." << std::endl;
		return 1;
	}

    for (TSize i = 0; i < length(matches); i++) {
		QueryMatches<TMatch> queryMatches = value(matches, i);
		if (queryMatches.disabled) {
			if (numDisabled == 0) {
				daFile << "Database sequence " << databaseID << ":\n\n";
			}
			daFile << ">" << ids[i] << "\n";
			daFile << queries[i] << "\n\n";
			++numDisabled;
		}
        if (length(queryMatches.matches) == 0) continue;
        //std::cout << "Pattern sequence: " << ids[i] << "\n";
        for (TSize j = 0; j < length(queryMatches.matches); j++) {
            TMatch m = value(queryMatches.matches, j);

            TSize len = _max(length(m.row1), length(m.row2));
            totalLength += len;
            if(len > maxLength) maxLength = len;

			if (format == "gff")
	            _writeMatchGff(databaseID, ids[i], databaseStrand, m.row1, m.row2, file);
			else 
				_writeMatch(databaseID, ids[i], databaseStrand, m.row1, m.row2, file);
        }
        numMatches += length(queryMatches.matches);
        //std::cout << "  # Eps-matches: " << length(queryMatches.matches) << std::endl;
    }
	daFile.close();
	file.close();

	if (verbose > 0 ) {
		if (numMatches > 0) {
			std::cout << "    Longest eps-match : " << maxLength << std::endl;
			std::cout << "    Avg match length  : " << totalLength / numMatches << std::endl;
		}
		std::cout << "    # Eps-matches     : " << numMatches << std::endl;
		std::cout << "    # Disabled queries: " << numDisabled << std::endl;
	}

	return numMatches;
}

#endif
