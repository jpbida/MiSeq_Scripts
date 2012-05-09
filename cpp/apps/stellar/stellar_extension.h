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

#ifndef SEQAN_HEADER_STELLAR_EXTENSION_H
#define SEQAN_HEADER_STELLAR_EXTENSION_H

#include <seqan/seeds2.h>


///////////////////////////////////////////////////////////////////////////////
// Container for storing possible end positions in extension of eps-core
template<typename TPos_>
struct ExtensionEndPosition {
	typedef TPos_			TPosition;
	typedef Pair<TPosition>	TCoordinate;

	TPosition length;
	TCoordinate coord;

	ExtensionEndPosition():
		length(0), coord(TCoordinate(0,0)) {}

	ExtensionEndPosition(TPosition len, TPosition row, TPosition col):
		length(len), coord(TCoordinate(row, col)) {}
};

///////////////////////////////////////////////////////////////////////////////
// returns true if align has a match at pos, otherwise false
template<typename TSource, typename TSize>
inline bool
isMatch(Align<TSource> const & align, TSize pos) {
SEQAN_CHECKPOINT
    if(isGap(row(align, 0), pos)) {
        return false;
    } else if(isGap(row(align, 1), pos)) {
        return false;
    } else if(row(align, 0)[pos] != row(align, 1)[pos]) {
        return false;
    } else {
        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Computes possible end positions of an eps-match in a given alignment.
template<typename TSource, typename TPos>
void
_fillGapsString(Align<TSource> const & align,
                String<Triple<TPos, TPos, TPos> > & gaps) {
SEQAN_CHECKPOINT
    typedef Triple<TPos, TPos, TPos> TGapInfo;
    TPos totalErrors = 0;
	typename Row<Align<TSource> >::Type row0 = row(align, 0);
    TPos i = 0;
	TPos endPos = endPosition(row0);
    TPos gapBegin = i;

    // append gap starting at beginPosition (also if its length is 0!)
    while(i < endPos && !isMatch(align, i)) {
        ++i;
        ++totalErrors;
    }
    appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));

    // iterate over alignment and append gaps
    while (i < endPos) {
        // skip matches
        while(i < endPos && isMatch(align, i)) {
            ++i;
        }
        gapBegin = i;
        // skip and count mismatches/indels
        while(i < endPos && !isMatch(align, i)) {
            ++i;
            ++totalErrors;
        }
        appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));
    }
    /*for(unsigned l = 0; l < length(gaps); ++l) {
        std::cout << gaps[l].i1 << "  " << gaps[l].i2 << "  " << gaps[l].i3 << std::endl;
    }*/
}

///////////////////////////////////////////////////////////////////////////////
// Checks the error rate of the fragment between end of left and start of right.
template<typename TPos, typename TFloat>
inline bool
_isEpsMatch(Triple<TPos, TPos, TPos> const & left,
           Triple<TPos, TPos, TPos> const & right,
           TFloat eps) {
SEQAN_CHECKPOINT
    // compute mismatches/indels and length
    TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
    TPos len = right.i1 - left.i2;

    // check error rate
    return errors/(TFloat)(len) <= eps;
}

///////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TSource, typename TSize, typename TFloat>
void 
longestEpsMatch(Align<TSource> & align,
                TSize matchMinLength,
                TFloat epsilon) {
SEQAN_CHECKPOINT
    // Preprocessing: compute and store gaps and lengths
    // A gap is a triple of gap begin position, gap end position, and total number of errors in sequence from begin
    //   to end position of this gap.
    typedef typename Position<Align<TSource> >::Type TPosition;
    typedef String<Triple<TPosition, TPosition, TPosition> > TGapsString;
    TGapsString gaps;
    _fillGapsString(align, gaps);

    // Identify longest eps match by iterating over combinations of left and right positions
    typename Iterator<TGapsString >::Type rightIt = end(gaps) - 1;
    typename Iterator<TGapsString >::Type leftIt = begin(gaps);

    TPosition beginPos = 0;
    TPosition endPos = 0;
    TSize minLength = matchMinLength - 1;
    
    while ((*leftIt).i2 + minLength < (*rightIt).i1) {
        while ((*leftIt).i2 + minLength < (*rightIt).i1) {
            if(_isEpsMatch(*leftIt, *rightIt, epsilon)) {
                beginPos = (*leftIt).i2;
                endPos = (*rightIt).i1;
                minLength = endPos - beginPos;
                break;
            }
            --rightIt;
        }
        rightIt = end(gaps) - 1;
        ++leftIt;
    }

    // Set view positions to the eps-match
	setClippedBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	setClippedBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	setBeginPosition(row(align, 0), beginPos);
	setBeginPosition(row(align, 1), beginPos);
	setClippedEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	setClippedEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
}

///////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align from possEndsLeft and possEndsRight and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TLength, typename TSize, typename TEps>
Pair<typename Iterator<String<ExtensionEndPosition<TLength> > >::Type>
longestEpsMatch(String<ExtensionEndPosition<TLength> > const & possEndsLeft,
				String<ExtensionEndPosition<TLength> > const & possEndsRight,
				TLength const alignLen,
				TLength const alignErr,
				TSize const matchMinLength,
				TEps const epsilon) {
SEQAN_CHECKPOINT
	typedef ExtensionEndPosition<TLength> TEnd;
	typedef typename Iterator<String<TEnd> >::Type	TIterator;

    // Identify longest eps match by iterating over combinations of left and right positions
    TIterator rightIt = end(possEndsRight) - 1;
    TIterator leftIt = end(possEndsLeft) - 1;
	TIterator right, left;

	/*for (int i = 0; i < length(possEndsRight); ++i) {
		std::cout << possEndsRight[i].length << "  " << possEndsRight[i].coord.i1 << "," << possEndsRight[i].coord.i2 << std::endl;
	}
	for (int i = 0; i < length(possEndsLeft); ++i) {
		std::cout << possEndsLeft[i].length << "  " << possEndsLeft[i].coord.i1 << "," << possEndsLeft[i].coord.i2 << std::endl;
	}*/

	TSize leftErr = length(possEndsLeft) - 1;

    TSize minLength = matchMinLength;
	bool found = false;

	while (leftIt >= begin(possEndsLeft)) {
		TSize totalLen = (*leftIt).length + alignLen + (*rightIt).length;
		if (totalLen < minLength) break;
		TSize totalErr = leftErr + alignErr + length(possEndsRight) - 1;
		while (rightIt >= begin(possEndsRight)) {
			totalLen = (*leftIt).length + alignLen + (*rightIt).length;
			if (totalLen < minLength) break;
			if ((TEps)totalErr/(TEps)totalLen <= epsilon) {
				right = rightIt;
				left = leftIt;
				//std::cout << totalLen << std::endl;
				minLength = totalLen;
				found = true;
				break;
			}
			--rightIt;
			--totalErr;
		}
		rightIt = end(possEndsRight) - 1;
		--leftIt;
		--leftErr;
	}

	if (found) return Pair<TIterator>(left, right);
	else return Pair<TIterator>(0,0);
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix and additionally a string with the best
//   alignment end point for each alignment length.
template <typename TTrace, typename TEnd, typename TStringSet, typename TScore, typename TDiagonal>
inline void
_align_banded_nw_best_ends(TTrace& trace,
						   String<TEnd> & bestEnds,
						   TStringSet const & str,
						   TScore const & sc,
						   TDiagonal const diagL,
						   TDiagonal const diagU)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TTrace>::Type TSize;

	SEQAN_ASSERT_GEQ(diagU, diagL);

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	TSize hi_diag = diagonalWidth;
	TSize lo_diag = 0;
	if (diagL > 0) lo_diag = 0;
	else lo_diag = (diagU < 0) ? hi_diag : (TSize) (1 - diagL); 
	TSize lo_row = (diagU <= 0) ? -diagU : 0;
	TSize hi_row = len2;
	if (len1 - diagL < hi_row) hi_row = len1 - diagL;
	TSize height = hi_row - lo_row;

	typedef String<TScoreValue> TRow;
	TRow mat, len;
	resize(mat, diagonalWidth);
	resize(len, diagonalWidth);
	resize(trace, height * diagonalWidth);

	// Classical DP with affine gap costs
	typedef typename Iterator<TRow, Standard>::Type TRowIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TSize actualCol = 0;
	TSize actualRow = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	TScoreValue hori_len = len1+len2+1;
	TSize errors;

	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if ((TDiagonal)actualRow >= (TDiagonal)len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TRowIter matIt = begin(mat, Standard()) + lo_diag;
		TRowIter lenIt = begin(len, Standard()) + lo_diag;
		hori_val = MinValue<TScoreValue>::VALUE;
		hori_len = len1+len2+1;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++matIt, ++traceIt, ++lenIt) {
			actualCol = col + diagL + actualRow;
			if (actualCol >= len1) break;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
				*traceIt = Diagonal;
				++(*lenIt);
				if ((verti_val = (col < diagonalWidth - 1) ? *(matIt+1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : MinValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = verti_val;
					*traceIt = Vertical;
					*lenIt = *(lenIt+1) + 1;
				}						
				if ((hori_val = (col > 0) ? hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : MinValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = hori_val;
					*traceIt = Horizontal;
					*lenIt = hori_len + 1;
				}
				hori_val = *matIt;
				hori_len = *lenIt;
			} else {			
				// Usual initialization for first row and column
				if (actualRow == 0) {
					*matIt = actualCol * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2);
					*lenIt = actualCol;
				}
				else {
					*matIt = actualRow * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2);
					*lenIt = actualRow;
					hori_val = *matIt;
					hori_len = actualRow;
				}
			}
			errors = (*matIt - (*lenIt * scoreMatch(const_cast<TScore&>(sc)))) /
						(scoreGap(const_cast<TScore&>(sc)) - scoreMatch(const_cast<TScore&>(sc)));
			SEQAN_ASSERT_LEQ(errors, length(bestEnds));
			if (errors == length(bestEnds)) {
					appendValue(bestEnds, TEnd(*lenIt, row, col));
			} else if (*lenIt > static_cast<TScoreValue>(value(bestEnds, errors).length))
				value(bestEnds, errors) = TEnd(*lenIt, row, col);
			//std::cerr << row << ',' << col << ':' << *matIt << std::endl;
		}
	}
	TSize newLength = length(bestEnds) - 1;
	while (newLength > 0 && bestEnds[newLength].length <= bestEnds[newLength-1].length) {
		--newLength;
	}
	resize(bestEnds, newLength + 1);
}

///////////////////////////////////////////////////////////////////////////////
// Reverses the infixes of the left extension in place in hosts of a and b.
template<typename TInfixA, typename TInfixB, typename TSeed>
void
_reverseLeftExtension(TInfixA const & a,
					  TInfixB const & b,
					  TSeed & seed,
					  TSeed & seedOld) {
SEQAN_CHECKPOINT
	TInfixB infixA = infix(host(a), getBeginDim0(seed), getBeginDim0(seedOld));
	TInfixB infixB = infix(host(b), getBeginDim1(seed), getBeginDim1(seedOld));
	reverse(infixA);
	reverse(infixB);
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix for the left extension and 
//   returns a string with possible start positions of an eps-match.
// Caution: The infixes of the left extension is reversed in place in hosts of a and b!
template<typename TMatrix, typename TPossEnd, typename TInfixA, typename TInfixB, typename TSeed, typename TScore>
void
_fillMatrixBestEndsLeft(TMatrix & matrixLeft,
							String<TPossEnd> & possibleEndsLeft,
							TInfixA const & a,
							TInfixB const & b,
							TSeed & seed,
							TSeed & seedOld,
							TScore const & scoreMatrix) {
SEQAN_CHECKPOINT
	TInfixB infixA = infix(host(a), getBeginDim0(seed), getBeginDim0(seedOld));
	TInfixB infixB = infix(host(b), getBeginDim1(seed), getBeginDim1(seedOld));

	reverse(infixA);
	reverse(infixB);

	StringSet<TInfixB> str;
	appendValue(str, infixA);
	appendValue(str, infixB);

	_align_banded_nw_best_ends(matrixLeft, possibleEndsLeft, str, scoreMatrix, 
							   getUpperDiagonal(seedOld) - getUpperDiagonal(seed),
							   getUpperDiagonal(seedOld) - getLowerDiagonal(seed));
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix for the right extension and 
//   returns a string with possible end positions of an eps-match.
template<typename TMatrix, typename TPossEnd, typename TInfixA, typename TInfixB, typename TSeed, typename TScore>
void
_fillMatrixBestEndsRight(TMatrix & matrixRight,
							String<TPossEnd> & possibleEndsRight,
							TInfixA const & a,
							TInfixB const & b,
							TSeed & seed,
							TSeed & seedOld,
							TScore const & scoreMatrix) {
SEQAN_CHECKPOINT
	StringSet<TInfixB> str;
	appendValue(str, infix(host(a), getEndDim0(seedOld), getEndDim0(seed)));
	appendValue(str, infix(host(b), getEndDim1(seedOld), getEndDim1(seed)));

	_align_banded_nw_best_ends(matrixRight, possibleEndsRight, str, scoreMatrix, 
							   getLowerDiagonal(seedOld) - getUpperDiagonal(seed),
							   getLowerDiagonal(seedOld) - getLowerDiagonal(seed));
}

///////////////////////////////////////////////////////////////////////////////
// Traceback from an arbitrary point (coordinate) in the banded alignment trace matrix (trace).
template <typename TAlign, typename TStringSet, typename TTrace, typename TCoord, typename TDiagonal>
inline void
_alignBandedNeedlemanWunschTrace(TAlign& align,
					   TStringSet const& str,
					   TTrace const& trace,
					   TCoord const& coordinate,
					   TDiagonal const diagL,
					   TDiagonal const diagU)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize lo_row = (diagU <= 0) ? -1 * diagU : 0;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);

	// Start the trace from the cell with the max value
	TSize row = coordinate.i1;
	TSize col = coordinate.i2;

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;

	if ((actualRow != 0) && (actualCol != 0)) {
		// Find initial direction
		TTraceValue tv = trace[row * diagonalWidth + col];
		if (tv == Horizontal) --col;
		else if (tv == Vertical) {--row; ++col;} 
		else --row;
	
		// Walk until we hit a border
		TSize seqLen = 1;
		TTraceValue newTv = tv;
		while(true) {
			actualRow = row + lo_row;
			actualCol = col + diagL + actualRow;
			newTv = trace[row * diagonalWidth + col];

			// Check if we hit a border
			if ((actualRow == 0) || (actualCol == 0)) break;
			else {
				//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
				if (tv == Diagonal) {
					if (newTv == Horizontal) {
						_alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
						--col; seqLen = 1;
					} else if (newTv == Vertical) {
						_alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
						--row; ++col; seqLen = 1;
					} else {
						--row; ++seqLen;
					}
				} else {
					if (tv == Horizontal) { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Vertical) {
							_alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
							--row; ++col; seqLen = 1;
						} else {
							--col; ++seqLen;
						}
					} else { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Horizontal) {
							_alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
							--col; seqLen = 1;
						} else {
							--row; ++col; ++seqLen;
						}
					}
				}
				tv = newTv;
			}
		}
	
		// Align left overs
		if (seqLen) _alignTracePrint(align, str, id1, actualCol, id2, actualRow, seqLen, tv);
	}

	// Handle the remaining sequence
	if (actualCol != 0) _alignTracePrint(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _alignTracePrint(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}

///////////////////////////////////////////////////////////////////////////////
// Conducts the traceback on the extension to the left from best start position
//   and writes the result into align.
template<typename TMatrix, typename TCoord, typename TInfixA, typename TInfixB, typename TSeed, typename TPos, typename TAlign>
void
_tracebackLeft(TMatrix const & matrixLeft,
			   TCoord const & coordinate,
			   TInfixA const & a,
			   TInfixB const & b,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const endLeftA,
			   TPos const endLeftB,
			   TAlign & align) {
SEQAN_CHECKPOINT
	typedef StringSet<TInfixB>							TInfixSet;
	typedef typename Size<TInfixSet>::Type				TSize;
	typedef typename Iterator<String<TraceBack> >::Type	TIterator;

	TInfixSet str;
	TInfixB infixA = infix(host(a), getBeginDim0(seed), getBeginDim0(seedOld));
	TInfixB infixB = infix(host(b), getBeginDim1(seed), getBeginDim1(seedOld));
	appendValue(str, infixA);
	appendValue(str, infixB);

	AlignTraceback<TPos> traceBack;
	_alignBandedNeedlemanWunschTrace(traceBack, str, matrixLeft, coordinate,
				   getUpperDiagonal(seedOld) - getUpperDiagonal(seed), getUpperDiagonal(seedOld) - getLowerDiagonal(seed));
	
	reverse(traceBack.sizes);
	reverse(traceBack.tvs);

	Align<TInfixB> infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(str[0], length(str[0]) - endLeftA, length(str[0])));
	assignSource(row(infixAlign, 1), infix(str[1], length(str[1]) - endLeftB, length(str[1])));

	_pumpTraceToAlign(infixAlign, traceBack);
	integrateAlign(align, infixAlign);
}


///////////////////////////////////////////////////////////////////////////////
// Conducts the traceback on the extension to the right from best end position
//   and writes the result into align.
template<typename TMatrix, typename TCoord, typename TInfixA, typename TInfixB, typename TSeed, typename TPos, typename TAlign>
void
_tracebackRight(TMatrix const & matrixRight,
			   TCoord const & coordinate,
			   TInfixA const & a,
			   TInfixB const & b,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const endRightA,
			   TPos const endRightB,
			   TAlign & align) {
SEQAN_CHECKPOINT
	typedef StringSet<TInfixB>							TInfixSet;
	typedef typename Size<TInfixSet>::Type				TSize;
	typedef typename Iterator<String<TraceBack> >::Type	TIterator;
		
	TInfixSet str;
	appendValue(str, infix(host(a), getEndDim0(seedOld), getEndDim0(seed)));
	appendValue(str, infix(host(b), getEndDim1(seedOld), getEndDim1(seed)));

	AlignTraceback<TPos> traceBack;
	_alignBandedNeedlemanWunschTrace(traceBack, str, matrixRight, coordinate,
				   getLowerDiagonal(seedOld) - getUpperDiagonal(seed), getLowerDiagonal(seedOld) - getLowerDiagonal(seed));

	Align<TInfixB> infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(str[0], 0, endRightA));
	assignSource(row(infixAlign, 1), infix(str[1], 0, endRightB));

	_pumpTraceToAlign(infixAlign, traceBack);
	integrateAlign(align, infixAlign);
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix and fills a string with possible start 
//   and end positions of an eps-match. Determines the optimal start and end
//   position for the longest eps-match and writes the trace into align.
template<typename TInfixA, typename TInfixB, typename TSeed, typename TPos, typename TDir, typename TScore,
		 typename TSize, typename TEps, typename TAlign>
bool
_bestExtension(TInfixA const & a,
			   TInfixB const & b,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const alignLen,
			   TPos const alignErr,
			   TScore const & scoreMatrix,
			   TDir const direction,
			   TSize const minLength,
			   TEps const eps,
			   TAlign & align) {
SEQAN_CHECKPOINT
	typedef String<TraceBack>							TAlignmentMatrix;
	typedef ExtensionEndPosition<TPos>					TEndInfo;
	typedef typename Iterator<String<TEndInfo> >::Type	TEndIterator;
	typedef typename Value<TScore>::Type				TScoreValue;
	typedef StringSet<TInfixB>							TInfixSet;

	// variables for banded alignment and possible ends of match
	TAlignmentMatrix matrixRight, matrixLeft;
	String<TEndInfo> possibleEndsLeft, possibleEndsRight;

	// fill banded matrix and gaps string for ...
	if (direction == EXTEND_BOTH || direction == EXTEND_LEFT) { // ... extension to the left
		_fillMatrixBestEndsLeft(matrixLeft, possibleEndsLeft, a, b, seed, seedOld, scoreMatrix);
		// Caution: left extension infix is now reversed in host(a and b) !!!
	} else appendValue(possibleEndsLeft, TEndInfo());
	if (direction == EXTEND_BOTH || direction == EXTEND_RIGHT) { // ... extension to the right
		_fillMatrixBestEndsRight(matrixRight, possibleEndsRight, a, b, seed, seedOld, scoreMatrix);
	} else appendValue(possibleEndsRight, TEndInfo());

	// longest eps match on poss ends string
	Pair<TEndIterator> endPair = longestEpsMatch(possibleEndsLeft, possibleEndsRight, alignLen, alignErr, minLength, eps);

	if (endPair == Pair<TEndIterator>(0, 0)) { // no eps-match found
		if (direction != 1) 
			_reverseLeftExtension(a, b, seed, seedOld); // back to original orientation
		return false;
	}

	// determine end positions of maximal eps-match in ...
	TPos endLeftA = 0, endLeftB = 0;
	TPos endRightA = 0, endRightB = 0;
	if((*endPair.i1).length != 0) { // ... extension to the left
		endLeftB = (*endPair.i1).coord.i1;
		// correction for banded coordinates to unbanded:
		if (getUpperDiagonal(seedOld) - getLowerDiagonal(seed) <= 0)
			endLeftB -= (TPos)(getUpperDiagonal(seedOld) - getLowerDiagonal(seed));
		endLeftA = (TPos)((*endPair.i1).coord.i2 + endLeftB + getUpperDiagonal(seedOld) - getUpperDiagonal(seed));
	}
	if((*endPair.i2).length != 0) { // ... extension to the right
		endRightB = (*endPair.i2).coord.i1;
		// correction for banded coordinates to unbanded:
		if (getLowerDiagonal(seedOld) - getLowerDiagonal(seed) <= 0)
			endRightB -= (TPos)(getLowerDiagonal(seedOld) - getLowerDiagonal(seed));
		endRightA = (TPos)((*endPair.i2).coord.i2 + endRightB + getLowerDiagonal(seedOld) - getUpperDiagonal(seed));
	}

	// set begin and end positions of align
	setClippedBeginPosition(row(align, 0), getBeginDim0(seedOld) - endLeftA);
	setClippedBeginPosition(row(align, 1), getBeginDim1(seedOld) - endLeftB);
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setClippedEndPosition(row(align, 0), getEndDim0(seedOld) + endRightA);
	setClippedEndPosition(row(align, 1), getEndDim1(seedOld) + endRightB);

	// traceback through matrix from begin/end pos on ...
	if((*endPair.i1).length != 0) { // ... extension to the left
		_tracebackLeft(matrixLeft, (*endPair.i1).coord, a, b, seed, seedOld, endLeftA, endLeftB, align);
	}
	if((*endPair.i2).length != 0) { // ... extension to the right
		_tracebackRight(matrixRight, (*endPair.i2).coord, a, b, seed, seedOld, endRightA, endRightB, align);
	}

	if (direction == EXTEND_BOTH || direction == EXTEND_LEFT) 
		_reverseLeftExtension(a, b, seed, seedOld); // back to original orientation

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Conducts best X-drop extension and calls _bestExtension.
//  After the call align contains the longest eps-Match that spans the eps-core (localAlign).
template<typename TScoreValue, typename TScore, typename TInfixA, typename TInfixB, typename TSize, typename TEps, typename TAlign>
bool
_extendAndExtract(Align<TInfixB> const & localAlign,
				  TScoreValue scoreDropOff,
				  TScore const & scoreMatrix,
				  TInfixA const & a,
				  TInfixB const & b,
				  ExtensionDirection direction,
				  TSize minLength,
				  TEps eps,
				  TAlign & align) {
SEQAN_CHECKPOINT
    typedef typename Position<TInfixB>::Type TPos;
    typedef Seed<Simple> TSeed;

	integrateAlign(align, localAlign);

	// begin and end position of local alignment (seed)
	TPos seedBeginA = clippedBeginPosition(row(localAlign, 0)) + beginPosition(a);
	TPos seedBeginB = clippedBeginPosition(row(localAlign, 1)) + beginPosition(b);
	TPos seedEndA = clippedEndPosition(row(localAlign, 0)) + beginPosition(a);
	TPos seedEndB = clippedEndPosition(row(localAlign, 1)) + beginPosition(b);

	if (direction == EXTEND_NONE) {
		// set begin and end positions of align
		setClippedBeginPosition(row(align, 0), seedBeginA);
		setClippedBeginPosition(row(align, 1), seedBeginB);
		setBeginPosition(row(align, 0), 0);
		setBeginPosition(row(align, 1), 0);
		setClippedEndPosition(row(align, 0), seedEndA);
		setClippedEndPosition(row(align, 1), seedEndB);

		if ((TSize)length(row(align, 0)) < minLength)
			return false;

		longestEpsMatch(align, minLength, eps);
	} else {
		// gapped X-drop extension of local alignment (seed)
		TSeed seed(seedBeginA, seedBeginB, seedEndA, seedEndB);
		TSeed seedOld(seed);
		extendSeed(seed, host(a), host(b), direction, scoreMatrix, scoreDropOff, GappedXDrop());

		if (static_cast<__int64>(getSeedSize(seed)) < minLength - (int)floor(minLength*eps))
			return false;

		// determine length and number of error columns of local alignment (seed)
		TPos alignLen = _max(length(row(localAlign, 0)), length(row(localAlign, 1)));
		TPos alignErr = 0;
		for (TPos i = 0; i < alignLen; ++i) {
			if (!isMatch(localAlign, i)) ++alignErr;
		}

		// determine best extension lengths and write the trace into align
		if (!_bestExtension(a, b, seed, seedOld, alignLen, alignErr, scoreMatrix, direction, minLength, eps, align))
			return false;
	}
	return true;
}

#endif
