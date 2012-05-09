// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef APPS_RABEMA_EVALUATION_H_
#define APPS_RABEMA_EVALUATION_H_

#include <seqan/basic.h>
#include <seqan/store.h>
#include <seqan/find.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/misc/misc_interval_tree.h>           // For interval trees.

#include "rabema.h"

#include "wit_store.h"
#include "find_hamming_simple_ext.h"
#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"
#include "find_approx_dp_quality.h"
#include "find_hamming_simple_quality.h"
#include "evaluation_options.h"

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// Counters for the comparison result.
struct ComparisonResult {
    // Total number of intervals in golden standard.
    size_t totalIntervalCount;

    // Number of intervals in golden standard that were found by the
    // read mapper.
    size_t foundIntervalCount;

    // Number of intervals with a too high distance that were found by
    // the read mappers, i.e. "junk output".
    size_t superflousIntervalCount;

    // Number of intervals with a good score that the read mapper
    // found which were not in our golden standard.
    size_t additionalIntervalCount;

    size_t totalReadCount;
    double normalizedIntervals;
    
    ComparisonResult()
            : totalIntervalCount(0), foundIntervalCount(0),
              superflousIntervalCount(0), additionalIntervalCount(0),
              totalReadCount(0), normalizedIntervals(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Resetting all counters to 0 for ComparisonResult.
void clear(ComparisonResult & result) {
    result.totalIntervalCount = 0;
    result.foundIntervalCount = 0;
    result.superflousIntervalCount = 0;
    result.additionalIntervalCount = 0;
    result.totalReadCount = 0;
    result.normalizedIntervals = 0;
}

// Output-to-stream operator for ComparisonResult.
template <typename TStream>
TStream & operator<<(TStream & stream, ComparisonResult const & result) {
    stream << "{\"total_intervals\": " << result.totalIntervalCount
           << ", \"found_intervals\": " << result.foundIntervalCount
           << ", \"superflous_intervals\": " << result.superflousIntervalCount
           << ", \"additional_intervals\": " << result.additionalIntervalCount
           << ", \"total_reads\": " << result.totalReadCount
           << ", \"normalized_intervals\": " << result.normalizedIntervals
           << "}";
    return stream;
}

// Copy-and-paste from reweight_wit.h
//
// Compute quality-based alignment score.  The read has to be given
// since we do not have qualities in the alignment object.
template <typename TAlign>
int computeQualityAlignmentScore(TAlign const & align,
                                 Score<int, ScoreMatrix<Dna5> > const & scoreMatrix,
                                 String<Dna5Q> const & read) {
    // TODO(holtgrew): Maybe convert to iterators for performance?
    typedef typename Row<TAlign const>::Type TRow;
    typedef typename Value<TRow>::Type TAlignChar;
    typedef typename Size<TRow>::Type TSize;

    TRow & rowContig = row(align, 0);
    TRow & rowRead = row(align, 1);

    int result = 0;
    for (TSize i = 0; i < length(rowContig); ++i) {
        TAlignChar contigChar = rowContig[i];
        TAlignChar readChar = rowRead[i];
        if (isGap(rowContig, i)) {
            result -= getQualityValue(read[toSourcePosition(rowRead, i)]);
        } else if (isGap(rowRead, i)) {
            if (toSourcePosition(rowRead, i) == 0) {
                result -= getQualityValue(read[0]);
            } else if (toSourcePosition(rowRead, i) == length(read)) {
                result -= getQualityValue(read[length(read) - 1]);
            } else {
                int x = 0;
                x += getQualityValue(read[toSourcePosition(rowRead, i) - 1]);
                x += getQualityValue(read[toSourcePosition(rowRead, i)]);
                result -= static_cast<int>(ceil(static_cast<double>(x) / 2.0));
            }
        } else {
            result += score(scoreMatrix, readChar, contigChar) * getQualityValue(read[toSourcePosition(rowRead, i)]);
        }
    }

    return result;
}


// Returns the best score for the alignment of the aligned read from
// the given fragment store.  The maximum error is given, to be able
// to limit the interval in the contig we are looking for.
template <typename TFragmentStore, typename TContigSeq2, typename TAlignedRead, typename TScore, typename TPatternSpec>
int bestScoreForAligned(TFragmentStore & fragments,
                        TContigSeq2 & contig2,
                        bool const & isForward,
                        TAlignedRead const & alignedRead,
                        int /*maxError*/,
                        Options<EvaluateResults> const & options,
                        TScore const & scoringScheme,
                        TPatternSpec const &) {
    typedef size_t TContigId;  // TODO(holtgrew): Better type.
    typedef size_t TAlignedReadPos;  // TODO(holtgrew): Better type.
    typedef typename TFragmentStore::TContigStore      TContigStore;
    typedef typename TFragmentStore::TContigSeq        TContigSeq;
    typedef typename Value<TContigStore>::Type         TContig;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
    typedef typename TFragmentStore::TReadSeqStore     TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq        TReadSeq;
    typedef typename Size<TReadSeq>::Type TSize;
    TContigStore & contigs = fragments.contigStore;
    TContigId contigId = alignedRead.contigId;
    TContigSeq & contig = contigs[contigId].seq;
    TContigGaps contigGaps(contigs[contigId].seq, contigs[contigId].gaps);
    TAlignedReadPos endPos = positionGapToSeq(contigGaps, alignedRead.endPos);
    TAlignedReadPos beginPos = positionGapToSeq(contigGaps, alignedRead.beginPos);
    TReadSeqStore & readSeqs = fragments.readSeqStore;
    TReadSeq read = readSeqs[alignedRead.readId];

    // If we are aligning on the reverse strand then we have to compute the
    // begin and end position on this strand.
    if (!isForward) {
        SEQAN_ASSERT_GT(beginPos, endPos);
        beginPos = length(contig) - beginPos;
        endPos = length(contig) - endPos;
    }

    // Initialize finder and pattern, configure to match N with none or all,
    // depending on configuration.
    Finder<TContigSeq2> finder(contig2);
    Pattern<TReadSeq, TPatternSpec> pattern(read, -static_cast<int>(length(read)) * 1000);
    _patternMatchNOfPattern(pattern, options.matchN);
    _patternMatchNOfFinder(pattern, options.matchN);
    bool ret = setEndPosition(finder, pattern, endPos);
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(endPos, endPosition(finder));
    
    // No explicit alignment is required if distances are not to be weighted.
    if (!options.weightedDistances)
        return getScore(pattern);

    // Otherwise, we need to build an alignment and compute the score from it.
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);

    // Prepare alignment datastructures.
    Align<String<Dna5>, ArrayGaps> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(finder));
    assignSource(row(align, 1), read);

    // Perform banded Needleman-Wunsch alignment.
    // TODO(holtgrew): Use wrapper for trace pumping, once it is in place.
    StringSet<String<Dna5> > stringSet;
    appendValue(stringSet, infix(finder));
    appendValue(stringSet, read);
    int alignmentScore = globalAlignment(align, stringSet, scoringScheme, getScore(pattern), -getScore(pattern), BandedNeedlemanWunsch());
    (void)alignmentScore; // Supress warning in non-debug mode.
    SEQAN_ASSERT_EQ(alignmentScore, getScore(pattern));

    // Compute quality-based score of alignment.  We pass the
    // score matrix to allow for N-is-wildcard mode.
    int qualityValue = computeQualityAlignmentScore(align, scoringScheme, read);
    return qualityValue;
}


// Return maximum error count for maximum error rate
inline int maxErrorRateToMaxErrors(int maxErrorRate, size_t len) {
    return (int)floor(maxErrorRate / 100.0 * len);
}


template <typename TFragmentStore, typename TAlignedReadIter, typename TWitRecordIter, typename TContigSeq, typename TContigGaps, typename TPatternSpec>
void
compareAlignedReadsToReferenceOnContigForOneRead(Options<EvaluateResults> const & options,
                                                 TFragmentStore & fragments,
                                                 size_t const & contigId,
                                                 TContigSeq /*const*/ & contig,
                                                 TContigGaps /*const*/ & contigGaps,
                                                 bool const & isForward,
                                                 size_t const & readId,
                                                 TAlignedReadIter const & alignedReadsBegin,
                                                 TAlignedReadIter const & alignedReadsEnd,
                                                 TWitRecordIter const & witRecordsBegin,
                                                 TWitRecordIter const & witRecordsEnd,
                                                 String<size_t> & result,
                                                 TPatternSpec const &) {
    typedef size_t TPos;
    typedef typename TFragmentStore::TContigStore TContigStore;

//     std::cout << "wit records are" << std::endl;
//     for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
//         std::cout << value(it) << std::endl;
//     }
//     std::cout << "XXX" << std::endl;

    // Build scoring matrix that allows N to match with all.
    int gapExtensionScore = -1;
    int gapOpenScore = -1;
    if (IsSameType<TPatternSpec, HammingSimple>::VALUE) {
        // No gaps for hamming distance.
        gapOpenScore = -static_cast<int>(length(fragments.readSeqStore[readId]));
        gapExtensionScore = -static_cast<int>(length(fragments.readSeqStore[readId]));
    }
    // Build scoring matrix.
    Score<int, ScoreMatrix<Dna5> > matrixScore(gapExtensionScore, gapOpenScore);
    for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
        for (int y = 0; y < ValueSize<Dna5>::VALUE; ++y)
            setScore(matrixScore, Dna5(x), Dna5(y), -1);
        setScore(matrixScore, Dna5(x), Dna5(x), 0);
    }
    // Update score matrix if N works as a wildcard character.
    if (options.matchN) {
        for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
            setScore(matrixScore, Dna5(x), Dna5('N'), 0);
            setScore(matrixScore, Dna5('N'), Dna5(x), 0);
        }
    }
        
    //if (readId == 71)
    //    std::cout << "stop here" << std::endl;

    // Build interval tree.
    //std::cerr << ">>> readId == " << readId << ", contigId == " << contigId << std::endl;
    //std::cerr << "-----------" << std::endl;
    typedef IntervalAndCargo<size_t, size_t> TInterval;
    String<TInterval> intervals;
    for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
        // Skip aligned reads on other strand.
        //if (it->isForward != isForward) continue;
        // Skip intervals with too high distance, ignore distance in oracle wit mode.
        //std::cerr << "it->distance == " << it->distance << std::endl;
        if (!options.oracleWitMode && static_cast<int>(it->distance) > options.maxError) continue;
        //std::cerr << "insert(intervals, TInterval(" << value(it).firstPos << ", " << value(it).lastPos + 1 << ", " << value(it).id << "))" << std::endl;

        SEQAN_ASSERT_LEQ(value(it).firstPos, value(it).lastPos);
        appendValue(intervals, TInterval(value(it).firstPos, value(it).lastPos + 1, value(it).id));
    }
    IntervalTree<size_t, size_t> intervalTree(intervals, ComputeCenter());
    //std::cerr << "-----------" << std::endl;

    // Now, try to hit all entries in interval tree with the aligned reads on this strand.
    std::set<size_t> intervalsInResult;
    for (TAlignedReadIter it = alignedReadsBegin; it != alignedReadsEnd; ++it) {
        SEQAN_ASSERT_EQ(it->readId, readId);
        // Skip aligned reads on other strand.
        if (isForward && it->beginPos > it->endPos)
            continue;
        if (!isForward && it->beginPos <= it->endPos)
            continue;

        // Convert from gap space to sequence space and maybe into reverse strand position.
        TPos endPos = positionGapToSeq(contigGaps, it->endPos);
        TPos beginPos = positionGapToSeq(contigGaps, it->beginPos);
        if (!isForward) {
            endPos = length(contig) - endPos;
            beginPos = length(contig) - beginPos;
        }

        // Skip if aligning too far to the left.
        if (endPos < length(fragments.readSeqStore[it->readId])) 
            continue;

        // Skip reads that aligned with a too bad score.  Ignore score if in wit-oracle mode, i.e. comparing against simulated data.
        int bestScore = 1;  // Marker for "not computed, oracle wit mode."
        if (!options.oracleWitMode) {
            int bestScore = bestScoreForAligned(fragments, contig, isForward, *it, maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])), options, matrixScore, TPatternSpec());
            if (bestScore < -maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId]))) {
                if (options.showSuperflousIntervals) {
                    std::cerr << "log> {\"type\": \"log.superflous_hit"
                              << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                              << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                              << "\", \"distance\": " << -bestScore
                              << ", \"strand\": \"" << (isForward ? "forward" : "reverse")
                              << "\", \"alignment_begin\": " << beginPos
                              << ", \"alignment_end\": " << endPos
                              << ", \"read_seq\": \"" << fragments.readSeqStore[it->readId]
                              << "\", \"contig_infix_seq\": \"" << infix(contig, beginPos, endPos);
                    std::cerr << "\", \"qualities\": [";
                    for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                        if (i > 0)
                            std::cerr << ", ";
                        std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]);
                    }
                    std::cerr << "]}" << std::endl;
                }
                appendValue(result, IntervalOfReadOnContig::superflousId());
                continue;
            }
        }

        // Compute last position of alignment and search for intervals this
        // position falls into.
        TPos lastPos = endPos - 1;
        SEQAN_ASSERT_LEQ(beginPos, lastPos);
        if (options.showTryHitIntervals) {
            std::cerr << "log> {\"type\": \"log.try_hit"
                      << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                      << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                      << "\", \"strand\": \"" << (isForward ? "forward" : "reverse")
                      << "\", \"alignment_first\": " << beginPos
                      << ", \"alignment_end\": " << endPos
                      << ", \"read_seq\": \"" << fragments.readSeqStore[it->readId]
                      << "\", \"contig_infix_seq\": \"" << infix(contig, beginPos, endPos);
            std::cerr << "\", \"qualities\": [";
            for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                if (i > 0)
                    std::cerr << ", ";
                std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]);
            }
            std::cerr << "]}" << std::endl;
        }

        // Query interval tree with lastPos.
        String<size_t> foundIntervalIds;
        //if (lastPos == 470875) {
        //  std::cout << "stop here" << std::endl;
        //  std::cout << "readId == " << readId << ", getMateNo() == " << getMateNo(fragments, readId) << std::endl;
        //}
        if (length(intervals) > 0u)
          findIntervals(intervalTree, lastPos, foundIntervalIds);

        // Handle alignment out of target intervals.
        if (length(foundIntervalIds) == 0) {
            if (options.weightedDistances) {
                if (options.showAdditionalIntervals) {
                    std::cerr << "log> {\"type\": \"log.additional_hit"
                              << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                              << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                              << "\", \"distance\": " << -bestScore
                              << ", \"strand\": \"" << (isForward ? "forward" : "reverse")
                              << "\", \"begin_pos\": " << beginPos
                              << ", \"end_pos\": " << endPos
                              << ", \"mate_no\": " << getMateNo(fragments, readId)
                              << ", \"read_seq\": \"" << fragments.readSeqStore[it->readId]
                              << "\", \"contig_infix_seq\": \"" << infix(contig, beginPos, endPos);
                    std::cerr << "\", \"qualities\": [";
                    for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                        if (i > 0)
                            std::cerr << ", ";
                        std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]);
                    }
                    std::cerr << "]}" << std::endl;
                }
                appendValue(result, IntervalOfReadOnContig::additionalId());
                continue;
            } else { // if (!options.weightedDistances) {
                if (options.dontPanic)
                    std::cerr << "WARNING: ";
                else
                    std::cerr << "PANIC: ";
                std::cerr << "A read in the Sam file aligns out of all target intervals for this read in the WIT file." << std::endl;
                std::cerr << "bestScore = " << bestScore << std::endl;
                std::cerr << "read name = " << fragments.readNameStore[it->readId] << std::endl;
                std::cerr << "read is = " << fragments.readSeqStore[it->readId] << std::endl;
                Dna5String rcRead(fragments.readSeqStore[it->readId]);
                reverseComplement(rcRead);
                std::cerr << "          " << rcRead << std::endl;
                std::cerr << "mate no is = " << static_cast<int>(getMateNo(fragments, it->readId)) << std::endl;
                std::cerr << "on forward strand? " << isForward << std::endl;
                std::cerr << "original begin pos = " << it->beginPos << std::endl;
                std::cerr << "original end pos = " << it->endPos << std::endl;
                std::cerr << "begin pos = " << beginPos << std::endl;
                std::cerr << "end pos = " << endPos << std::endl;
                std::cerr << "last pos = " << lastPos << std::endl;
                std::cerr << "contigId = " << it->contigId << std::endl;
                std::cerr << "mateNo " << getMateNo(fragments, it->readId) << std::endl;
                std::cerr << "max error rate is " << options.maxError << std::endl;
                std::cerr << "contig_infix_seq = " << infix(contig, beginPos, endPos) << std::endl;
                std::cerr << "read length is " << length(fragments.readSeqStore[it->readId]) << std::endl;
                std::cerr << "max errors is " <<  maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])) << std::endl;
                if (options.oracleWitMode)
                    std::cerr << "NOTE: Oracle WIT mode is enabled!" << std::endl;
                else if (!options.dontPanic)
                    exit(1);
                appendValue(result, IntervalOfReadOnContig::additionalId());
            }
        }

        // Record found intervals.
        typedef typename Iterator<String<size_t>, Standard>::Type TFoundIdsIterator;
        for (TFoundIdsIterator iter = begin(foundIntervalIds, Standard()); iter != end(foundIntervalIds, Standard()); ++iter) {
            // Count interval if hit for the first time.
            if (intervalsInResult.find(value(iter)) == intervalsInResult.end()) {
                appendValue(result, value(iter));
                intervalsInResult.insert(value(iter));
            }
        }
    }

    // If configured so, show hit and missed intervals.
    if (options.showHitIntervals) {
        for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
            bool found = intervalsInResult.find(value(it).id) != intervalsInResult.end();
			if (found && options.showHitIntervals) {
                std::cerr << "log> {\"type\": \"log.hit_interval"
                          << "\", \"interval_id\": " << value(it).id
                          << ", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                          << "\", \"strand\": \"" << (isForward ? "forward" : "reverse")
                          << "\", \"read_id\": \"" << fragments.readNameStore[readId]
                          << "\", \"interval_first\": " << value(it).firstPos
                          << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
            }
        }
    }
}


// Compare the aligned reads in [alignedReadsBegin, alignedReadsEnd)
// to the intervals in [witRecordsBegin, witRecordsEnd) on the given
// contig on the forward strand iff isForward.
//
// foundIntervalCount is incremented by the number of hit intervals,
// the number of relevant (on the selected strand) intervals is added
// to relevantIntervalCount.
template <typename TFragmentStore, typename TAlignedReadsIter, typename TWitRecordsIter, typename TContigSeq, typename TPatternSpec>
void
compareAlignedReadsToReferenceOnContig(Options<EvaluateResults> const & options,
                                       TFragmentStore & fragments,
                                       size_t const & contigId,
                                       TContigSeq /*const*/ & contig,
                                       bool const & isForward,
                                       TAlignedReadsIter const & alignedReadsBegin,
                                       TAlignedReadsIter const & alignedReadsEnd,
                                       TWitRecordsIter const & witRecordsBegin,
                                       TWitRecordsIter const & witRecordsEnd,
                                       String<size_t> & result,
                                       TPatternSpec const &) {
    typedef size_t TPos;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContigStoreElement::TGapAnchors> > TContigGaps;

    // Build contig gaps datastructure for gap space to sequence space conversion.
    TContigGaps contigGaps(fragments.contigStore[contigId].seq, fragments.contigStore[contigId].gaps);

    // The aligned reads and wit records are filtered to this contig and
    // sorted by read id.
    TAlignedReadsIter alignedReadsBeginForRead = alignedReadsBegin;
    TAlignedReadsIter alignedReadsEndForRead = alignedReadsBeginForRead;
    TWitRecordsIter witRecordsBeginForRead = witRecordsBegin;
    TWitRecordsIter witRecordsEndForRead = witRecordsBegin;

    while (alignedReadsBeginForRead != alignedReadsEnd && witRecordsBeginForRead != witRecordsEnd) {
        // Get current contigId.
        size_t readId = _min(alignedReadsBeginForRead->readId, witRecordsBeginForRead->readId);
        // Get aligned reads iterator for the next contigId.
        while (alignedReadsEndForRead != alignedReadsEnd && alignedReadsEndForRead->readId <= readId)
            ++alignedReadsEndForRead;
        // Get wit records iterator for the next contigId.
        while (witRecordsEndForRead != witRecordsEnd && witRecordsEndForRead->readId <= readId)
            ++witRecordsEndForRead;

        // Actually compare the aligned reads for this contig on forward and backwards strand.
//         std::cout << "read id " << fragments.readNameStore[readId] << std::endl;
        compareAlignedReadsToReferenceOnContigForOneRead(options, fragments, contigId, contig, contigGaps, isForward, readId, alignedReadsBeginForRead, alignedReadsEndForRead, witRecordsBeginForRead, witRecordsEndForRead, result, TPatternSpec());
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBeginForRead = alignedReadsEndForRead;
        witRecordsBeginForRead = witRecordsEndForRead;
    }
}


// Compare aligned reads in fragment store to the intervals specified
// in witStore, results is used for the counting statistics.
template <typename TFragmentStore, typename TPatternSpec>
void
compareAlignedReadsToReference(String<size_t> & result,
                               TFragmentStore & fragments,  // non-const so reads can be sorted
                               WitStore & witStore,  // non-const so it is sortable
                               Options<EvaluateResults> const & options,
                               TPatternSpec const &) {
    // Type aliases.
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename WitStore::TIntervalStore TIntervalStore;
    typedef typename Iterator<TIntervalStore>::Type TIntervalIter;
    
    // Initialization.
    clear(result);

    // Sort aligned reads and wit records by contig id.
    sortAlignedReads(fragments.alignedReadStore, SortEndPos());
    sortAlignedReads(fragments.alignedReadStore, SortReadId());
    sortAlignedReads(fragments.alignedReadStore, SortContigId());
    sortWitRecords(witStore, SortReadId());
    sortWitRecords(witStore, SortContigId());

    TAlignedReadsIter alignedReadsBegin = begin(fragments.alignedReadStore, Standard());
    TAlignedReadsIter alignedReadsEnd = alignedReadsBegin;
    TIntervalIter witRecordsBegin = begin(witStore.intervals, Standard());
    TIntervalIter witRecordsEnd = witRecordsBegin;

    while (alignedReadsBegin != end(fragments.alignedReadStore, Standard()) &&
           witRecordsBegin != end(witStore.intervals, Standard())) {
        // Get current contigId.
        size_t contigId = _min(alignedReadsBegin->contigId, witRecordsBegin->contigId);
        // Get aligned reads iterator for the next contigId.
        while (alignedReadsEnd != end(fragments.alignedReadStore, Standard()) && alignedReadsEnd->contigId <= contigId)
            ++alignedReadsEnd;
        // Get wit records iterator for the next contigId.
        while (witRecordsEnd != end(witStore.intervals, Standard()) && witRecordsEnd->contigId <= contigId)
            ++witRecordsEnd;

        // Actually compare the aligned reads for this contig on forward and backwards strand.
        TContigSeq & contig = fragments.contigStore[contigId].seq;
        compareAlignedReadsToReferenceOnContig(options, fragments, contigId, contig, true, alignedReadsBegin, alignedReadsEnd, witRecordsBegin, witRecordsEnd, result, TPatternSpec());
        TContigSeq rcContig(contig);
        reverseComplement(rcContig);
        compareAlignedReadsToReferenceOnContig(options, fragments, contigId, rcContig, false, alignedReadsBegin, alignedReadsEnd, witRecordsBegin, witRecordsEnd, result, TPatternSpec());
        
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBegin = alignedReadsEnd;
        witRecordsBegin = witRecordsEnd;
    }
}


// Tags for selecting an evaluateFoundIntervals_compareToIntervals specialization.
struct CategoryAll_;
typedef Tag<CategoryAll_> CategoryAll;
struct CategoryAnyBest_;
typedef Tag<CategoryAnyBest_> CategoryAnyBest;
struct CategoryAllBest_;
typedef Tag<CategoryAllBest_> CategoryAllBest;


// result must contain all ids of the intervals in
// [beginIntervalsForRead, endIntervalsForRead) that have the error
// rate options.maxError.  The intervals are sorted by (read id,
// error rate).
template <typename TFragmentStore>
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               TFragmentStore const & /*fragments*/,
                                               Options<EvaluateResults> const & options,
                                               CategoryAll const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // Now, get begin and end iterators to the intervals for the error rate
    // configured in options, 0 in oracle wit mode.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = options.oracleWitMode ? 0 : options.maxError;
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    // Intervals to be found for this read and counter for found intervals.
    size_t intervalsForRead = boundsForDistance.second - boundsForDistance.first;
    size_t intervalsFound = 0;

    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        SEQAN_ASSERT_EQ(static_cast<int>(value(it).distance), options.maxError);
        bool found = std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id);
        intervalsFound += found;
        if (!found && options.showMissedIntervals) {
            std::cerr << "log> {\"type\": \"log.missed_interval"
                      << "\", \"interval_id\": " << value(it).id
                      << ", \"contig_id\": \"" << value(witStore.contigNames)[value(it).contigId]
                      << "\", \"strand\": \"" << (value(it).isForward ? "forward" : "reverse")
                      << "\", \"read_id\": " << value(it).readId
                      << ", \"read_name\": \"" << value(witStore.readNames)[value(it).readId]
                      << "\", \"interval_first\": " << value(it).firstPos
                      << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
        }
    }
    if (intervalsFound > intervalsForRead)
      exit(-1);

    // Update comparison results.
    comparisonResult.totalIntervalCount += intervalsForRead;
    comparisonResult.foundIntervalCount += intervalsFound;
    comparisonResult.totalReadCount += 1;
    comparisonResult.normalizedIntervals += 1.0 * intervalsFound / intervalsForRead;
}


template <typename TFragmentStore>
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               TFragmentStore const & /*fragments*/,
                                               Options<EvaluateResults> const & options,
                                               CategoryAnyBest const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // The read mapper has to find one of the intervals with the smallest
    // error rate in [beginIntervalsForRead, endIntervalsForRead).
    //
    // The intervals in this range are sorted by error rate.  We get the
    // smallest error rate first, then get the subrange with the smallest
    // range.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = value(beginIntervalsForRead).distance;
    if (static_cast<int>(refInterval.distance) > options.maxError)
        return;  // Guard: Skip if best has too large distance.
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    // One interval is to be found.
    comparisonResult.totalIntervalCount += 1;
    comparisonResult.totalReadCount += 1;

    // Now, try to find one of the intervals.
    bool found = false;
    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        if (std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id)) {
            found = true;
            comparisonResult.foundIntervalCount += 1;
            comparisonResult.normalizedIntervals += 1;
            break;
        }
    }
    if (!found && options.showMissedIntervals) {
      for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
          std::cerr << "log> {\"type\": \"log.missed_interval"
                    << "\", \"interval_id\": " << value(it).id
                    << ", \"contig_id\": \"" << value(witStore.contigNames)[value(it).contigId]
                    << "\", \"strand\": \"" << (value(it).isForward ? "forward" : "reverse")
                    << "\", \"read_id\": \"" << value(witStore.readNames)[value(it).readId]
                    << "\", \"distance\": " << value(beginIntervalsForRead).distance
                    << ", \"interval_first\": " << value(it).firstPos
                    << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
      }
    }
}


template <typename TFragmentStore>
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               TFragmentStore const & /*fragments*/,
                                               Options<EvaluateResults> const & options,
                                               CategoryAllBest const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // The read mapper has to find one of the intervals with the smallest
    // error rate in [beginIntervalsForRead, endIntervalsForRead).
    //
    // The intervals in this range are sorted by error rate.  We get the
    // smallest error rate first, then get the subrange with the smallest
    // range.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = value(beginIntervalsForRead).distance;  // Difference to all: Distance of smallest here not options.maxError.
    if (static_cast<int>(refInterval.distance) > options.maxError)
        return;  // Guard: Skip if best has too large distance.
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    comparisonResult.totalReadCount += 1;
    size_t foundReads = 0;
    size_t bestReadCount = boundsForDistance.second - boundsForDistance.first;

    // Now, try to find one of the intervals.
    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        comparisonResult.totalIntervalCount += 1;
        bool found = std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id);
        comparisonResult.foundIntervalCount += found;
        foundReads += found;
        if (!found && options.showMissedIntervals) {
            std::cerr << "log> {\"type\": \"log.missed_interval"
                      << "\", \"interval_id\": " << value(it).id
                      << ", \"contig_id\": \"" << value(witStore.contigNames)[value(it).contigId]
                      << "\", \"strand\": \"" << (value(it).isForward ? "forward" : "reverse")
                      << "\", \"read_id\": \"" << value(witStore.readNames)[value(it).readId]
                      << "\", \"interval_first\": " << value(it).firstPos
                      << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
        }
    }

    comparisonResult.normalizedIntervals += 1.0 * foundReads / bestReadCount;
}


// Evaluate the ids in result pointing to intervals in witStore.
// Result is written to comparisonResult.
template <typename TFragmentStore>
void evaluateFoundIntervals(ComparisonResult & comparisonResult,
                            WitStore & witStore,
                            String<size_t> & result,
                            TFragmentStore const & fragments,
                            Options<EvaluateResults> const & options)
{
    typedef Iterator<String<size_t>, Standard>::Type TIdIterator;
    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // Initialization.
    clear(comparisonResult);
    sortWitRecords(witStore, SortDistance());
    sortWitRecords(witStore, SortReadId());
    std::sort(begin(result, Standard()), end(result, Standard()));

    // Count ids of superflous and additional intervals (at end of result).
    std::pair<TIdIterator, TIdIterator> superflousBounds = std::equal_range(begin(result, Standard()), end(result, Standard()), IntervalOfReadOnContig::superflousId());
    comparisonResult.superflousIntervalCount = superflousBounds.second - superflousBounds.first;
    std::pair<TIdIterator, TIdIterator> additionalBounds = std::equal_range(begin(result, Standard()), end(result, Standard()), IntervalOfReadOnContig::additionalId());
    comparisonResult.additionalIntervalCount = additionalBounds.second - additionalBounds.first;

    // For each read: Get list of intervals that are to be found,
    // depending on options.  Then, look whether they are in result.
    for (size_t readId = 0; readId < length(value(witStore.readNames)); ++readId) {
        // We need a read id to search for in a sequence of
        // IntervalOfReadOnContig objects.  To do this, we create a
        // reference interval with the current read id.
        IntervalOfReadOnContig refInterval;
        refInterval.readId = readId;
        std::pair<TIntervalIterator, TIntervalIterator> boundsForRead;
        boundsForRead = std::equal_range(begin(witStore.intervals, Standard()), end(witStore.intervals, Standard()), refInterval, WitStoreLess<SortReadId>(witStore));
        if (boundsForRead.first == boundsForRead.second)
            continue;  // Skip if there are no intervals to be found for read.

        // Now, depending on the configured benchmark category, perform the
        // evaluation.
        if (options.benchmarkCategory == "all") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, fragments, options, CategoryAll());
        } else if (options.benchmarkCategory == "any-best") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, fragments, options, CategoryAnyBest());
        } else if (options.benchmarkCategory == "all-best") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, fragments, options, CategoryAllBest());
        } else {
            SEQAN_ASSERT_FAIL("Invalid benchmark category '%s'.", toCString(options.benchmarkCategory));
        }
    }
}

// Entry point for the read mapper evaluation subprogram.
int evaluateReadMapperResult(Options<EvaluateResults> const & options)
{
    // =================================================================
    // Load FASTA Sequence And Sam File Into FragmentStore.
    // =================================================================
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;

    // Load Contigs.
    double startTime = sysTime();
    std::cerr << "Reading FASTA contigs sequence file " << options.seqFileName << " ..." << std::endl;
    if (!loadContigs(fragments, options.seqFileName)) {
        std::cerr << "Could not read contigs." << std::endl;
        return 1;
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // Load Sam File.
    std::cerr << "Reading Sam file file " << options.samFileName << " ..." << std::endl;
    startTime = sysTime();
    {
        std::fstream fstrm(toCString(options.samFileName),
                           std::ios_base::in | std::ios_base::binary);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open Sam file." << std::endl;
            return 1;
        }
        read(fstrm, fragments, Sam());
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    //for (unsigned i = 0; i < length(fragments.readNameStore); ++i) {
      //std::cerr << ">>>" << fragments.readNameStore[i] << " " << i << std::endl;
    //}

    // =================================================================
    // Load WIT file.
    // =================================================================
    std::cerr << "Loading intervals from " << options.witFileName << std::endl;
    startTime = sysTime();
    WitStore witStore;
    loadWitFile(witStore, fragments, options.witFileName);
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Compare The Sam Hits Against WIT Intervals.
    // =================================================================
    std::cerr << "Compare reader hits from Sam file against WIT file." << std::endl;
    startTime = sysTime();
    typedef Position<WitStore::TIntervalStore>::Type TPos;
    // The result will be a list of ids to entries in witStore.
    String<size_t> result;
    if (options.distanceFunction == "edit")
        compareAlignedReadsToReference(result, fragments, witStore, options, MyersUkkonenReads());
    else  // options.distanceFunction == "hamming"
        compareAlignedReadsToReference(result, fragments, witStore, options, HammingSimple());
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Perform Counting On Result Indices, Yield ComparisonResult.
    // =================================================================
    ComparisonResult comparisonResult;
    evaluateFoundIntervals(comparisonResult, witStore, result, fragments, options);

    // =================================================================
    // Write Output.
    // =================================================================
    startTime = sysTime();
    // The output consists of one line that describes the total and
    // found intervals as a JSON record with the entries
    // "total_intervals", "found_itervals", "superflous_intervals",
    // "additional_intervals".
    if (options.outFileName == "-") {
        // Print to stdout.
        std::cout << comparisonResult << std::endl;
    } else {
        // Write output to file.
        std::fstream fstrm(toCString(options.outFileName), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open output JSON file." << std::endl;
            return 1;
        }
        fstrm << comparisonResult << std::endl;
    }

    return 0;
}

#endif  // #ifndef APPS_RABEMA_EVALUATION_H_
