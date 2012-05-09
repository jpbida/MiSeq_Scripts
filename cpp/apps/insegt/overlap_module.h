 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_OVERLAP_MODULE_H
#define SEQAN_HEADER_OVERLAP_MODULE_H
//#define DEBUG_OVERLAP_MODULE

namespace SEQAN_NAMESPACE_MAIN
{ 

//////////////////////////////////////////////////////////////////////////////
// getIdsFroRead
//////////////////////////////////////////////////////////////////////////////

template<typename TAnnoIds, typename TSpec, typename TConfig, typename TIntervalTree, typename TIntervals>
inline void
getIdsForRead(TAnnoIds & ids, FragmentStore<TSpec, TConfig> & me, TIntervalTree & intervalTree, TIntervals & alignIntervals, unsigned offsetInterval)
{
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	typedef typename Value<TIntervals>::Type 				TInterval;
	typedef 	 String<TId> 						TResult;
	typedef typename Iterator<TIntervals >::Type 				TIntervalIter;
	typedef typename Iterator<StringSet<TResult > >::Type			TResultIter;
	typedef typename Iterator<TResult >::Type				TIdIter;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(ids, length(alignIntervals));
	
	TIntervalIter itI = begin(alignIntervals);
	TIntervalIter itIEnd = end(alignIntervals);
	TResultIter itR = begin(ids);
	TResultIter itREnd = end(ids);

	// search mapped annotations for each interval of the aligned read and store them in the corresponding list 'ids'
	for ( ; itI != itIEnd; goNext(itI), goNext(itR))
	{
		findIntervalsForInterval(value(itR), intervalTree, getValue(itI), offsetInterval);	
	}
	
	// check for each aligment-interval, if the inner interval-borders fit to the borders of the annotation id:
	itI = begin(alignIntervals);
	itR = begin(ids);
	TId currentId;
	TContigPos beginPos;
	TContigPos endPos;
	TContigPos help;
	
	for ( ; itI != itIEnd; goNext(itI), goNext(itR))
	{
		for (unsigned i = 0; i < length(getValue(itR)); ++i)
		{
			currentId = getValue(getValue(itR), i);
			beginPos = getValue(me.annotationStore, currentId).beginPos;
			endPos = getValue(me.annotationStore, currentId).endPos;
			
			if (beginPos > endPos)
			{
				help = beginPos;
				beginPos = endPos;
				endPos = beginPos;
			}
			// begin of read
			if (itR == begin(ids) && length(ids) > 1)
			{
				if (static_cast<TContigPos>(getValue(itI).i2 + offsetInterval) < endPos) // if the borders don't fit: delete annotation-id
				{
					erase(value(itR), i);
					--i;
				}
			}
			// end of read
			else if (position(itR, ids) == endPosition(ids) - 1 && length(ids) > 1u)
			{
				if (static_cast<TContigPos>(getValue(itI).i1 - offsetInterval) > beginPos)
				{
					erase(value(itR), i);
					--i;
				}
			}
			// in the middle of the read
			else if (length(ids) > 2)
			{
				if (static_cast<TContigPos>(getValue(itI).i2 + offsetInterval) < endPos)
				{
					erase(value(itR), i);
					--i;
				}
				else if (static_cast<TContigPos>(getValue(itI).i1 - offsetInterval) > beginPos)
				{
					erase(value(itR), i);
					--i;
				}
			}
		}
		if (empty(getValue(itR)) )  // if aligment-interval doesn't fit to any annotation, append INVALID_ID to mark this
			appendValue(value(itR), INVALID_ID, Generous());
	}
}


//////////////////////////////////////////////////////////////////////////////
////// ReadAnnoStoreELement
//////////////////////////////////////////////////////////////////////////////
template <typename TId>
struct ReadAnnoStoreElement
{
	typedef StringSet<String<TId> > TAnnoIds;
	
	TAnnoIds	annoIds;
	String<TId> 	parentIds;     // not only for exon-annotations -> more than one parentId possible, only if whole read mapped in parent
	TId		contigId;
};


//////////////////////////////////////////////////////////////////////////////
////// assign Ids to ReadAnnoStore
//////////////////////////////////////////////////////////////////////////////
template<typename TReadAnnoStore, typename TSpec, typename TConfig, typename TId, typename TAnnoIds>
inline void
assignToReadAnnoStore(TReadAnnoStore &readAnnoStore, FragmentStore<TSpec, TConfig> & me, TId readId, TAnnoIds &annoIds)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename Iterator<TAnnoIds>::Type				TAnnoIdsIter;	
	typedef typename Value<TAnnoIds>::Type					TIds;
	typedef typename Iterator<TIds>::Type					TIdsIter;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	// assign annotationIds:
	value(readAnnoStore, readId).annoIds = annoIds;	
		
	// assign coresponding parentIds:
	clear (value(readAnnoStore, readId).parentIds);
	if(!empty(annoIds))
	{
		TIdsIter itId = begin(front(annoIds));
		TIdsIter itIdEnd = end(front(annoIds));
		for ( ; itId != itIdEnd; goNext(itId))		// read maps in gene, if all intervals map in gene: at least one exon of the gene has to occur in the id-list of the first interval
			if (getValue(itId) != INVALID_ID && !isElement_unsorted(getValue(me.annotationStore, getValue(itId)).parentId, getValue(readAnnoStore, readId).parentIds))	
				appendValue(value(readAnnoStore, readId).parentIds, getValue(me.annotationStore, getValue(itId)).parentId, Generous() );
	
		if (!empty(getValue(readAnnoStore, readId).parentIds))
		{
			TAnnoIdsIter itA = begin(annoIds);
			TAnnoIdsIter itAEnd = end(annoIds);
			goNext(itA);
			for ( ; itA != itAEnd; goNext(itA))	// not only for exon-annotations -> more than one parentId possible
			{
				itId = begin(getValue(itA));		// for each interval of read:
				itIdEnd = end(getValue(itA));
				for (unsigned i = 0; i < length(getValue(readAnnoStore, readId).parentIds); ++i) // check if at least one child of the parentId occurs 
				{						
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) != INVALID_ID && getValue(me.annotationStore, getValue(itId)).parentId == getValue(getValue(readAnnoStore, readId).parentIds, i) )
							break;
					}
					if (itId == itIdEnd)			 // if not, delete parentId
					{
						erase(value(readAnnoStore, readId).parentIds, i);
						--i;
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// buildTupleCountStore
//////////////////////////////////////////////////////////////////////////////
template <typename TId>
struct TupleCountStoreElement
{
	typedef String<TId>			TTuple;				
	typedef String<TTuple > 		TTupleList;
	typedef String<unsigned>		TTupleCounts;
	typedef String<double>			TTupleNorm;

	TTupleList		readConnections;
	TTupleCounts		readConnectionCounts;
	TTupleNorm		readConnectionNorm;
	TTupleList		matePairConnections; 
	TTupleCounts		matePairConnectionCounts;
	TTupleNorm		matePairConnectionNorm;
};


//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore, typename TSpec, typename TConfig, typename TReadAnnoStore>
inline void
buildTupleCountStore(TTupleCountStore & tupleCountStore, 
		     FragmentStore<TSpec, TConfig> &  me, 
		     TReadAnnoStore & readAnnoStore, 
		     unsigned n, 
		     bool exact_nTuple)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore 		TReadStore;
	typedef typename Value<TReadStore>::Type 				TReadStoreElement;
	typedef typename TReadStoreElement::TId					TReadId;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId; 
	typedef typename Iterator<TReadAnnoStore>::Type 			TReadIter;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIds			TAnnoIds;
	typedef typename Iterator<TAnnoIds>::Type 				TAnnoIdsIter;
	typedef typename Value<TAnnoIds>::Type					TIds;
	typedef typename Iterator<TIds>::Type					TIdsIter;
	
	static const TReadId INVALID_READ_ID = TReadStoreElement::INVALID_ID;
	static const TId INVALID_ANNO_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(tupleCountStore, length(me.annotationStore)); 
	
	bool validMate;
	TReadIter itRead = begin(readAnnoStore);
	TReadIter itReadEnd = end(readAnnoStore);
	TIdsIter itP;
	TIdsIter itPEnd;
	TAnnoIds annoIds;
	TAnnoIds tupleSet;
	TReadId readId;
	TReadId matePairId;
	TReadId secReadId;
	TAnnoIds secTupleSet;
	TId firstAnnoId1;
	TId firstAnnoId2;
	TAnnoIdsIter itTuple;
	TAnnoIdsIter itTupleEnd;
	TAnnoIdsIter itSecTuple;
	TAnnoIdsIter itSecTupleEnd;
	TAnnoIdsIter itAnnoIds;
	TAnnoIdsIter itAnnoIdsEnd;
	TIds matePairTuple;
	TPos beginPos1;	
	TPos endPos1;
	TPos beginPos2;
	TPos endPos2;
	unsigned pos;
	
	for ( ; itRead != itReadEnd; goNext(itRead))
	{
		if (!empty(getValue(itRead).parentIds) )
		{ 
			itP = begin(getValue(itRead).parentIds);
			itPEnd = end(getValue(itRead).parentIds);
			for ( ; itP != itPEnd; goNext(itP) )
			{
				validMate = false;
				// create list of all possible tuples for current read:
				annoIds = getValue(itRead).annoIds;
				clear(tupleSet);
				// create all Tuple of length n:
				if (exact_nTuple && n <= length(annoIds)) create_nTuple(tupleSet, me, annoIds, getValue(itP), n);
				// create all max-Tuple (whole read) for current parentId:
				else if (!exact_nTuple && n == 0) create_nTuple(tupleSet, me, annoIds, getValue(itP), length(annoIds));	
				// create all tuple >= n for current parentId:
				else if (!exact_nTuple) create_Tuple(tupleSet, me, annoIds, getValue(itP), n);	
				if (!empty(tupleSet))
				{
					// create if necessary list of all possible tuples for second matepair-read:
					readId = position(itRead, readAnnoStore);
					matePairId = getValue(me.readStore, readId).matePairId;
					clear(secTupleSet);
					if (matePairId != INVALID_READ_ID)
					{
						if(getValue(getValue(me.matePairStore, matePairId).readId, 0) == readId)
							secReadId = getValue(getValue(me.matePairStore, matePairId).readId, 1);
						else
							secReadId = getValue(getValue(me.matePairStore, matePairId).readId, 0);
				
						if ( secReadId != INVALID_READ_ID )	
						{
							//if (!empty(getValue(readAnnoStore, secReadId).annoIds)) 
							if ( isElement_unsorted(getValue(itP), getValue(readAnnoStore, secReadId).parentIds) )	// p in parents of matepair? -> annoIds is not empty
							{
								validMate = true;
								annoIds = getValue(readAnnoStore, secReadId).annoIds;
								firstAnnoId1 = front(front(tupleSet));	// ids necessary to check positions in aligment  
								firstAnnoId2 = front(front(annoIds));	// can't be INVALID_ID, because parents was checked
							
								// check if current read-position is smaller than the position of the second read -> tuple are ordered by position
								if ( (getValue(me.annotationStore, firstAnnoId1).beginPos <= 
									getValue(me.annotationStore,firstAnnoId1).endPos && 
								      getValue(me.annotationStore, firstAnnoId1).beginPos < 
								      	getValue(me.annotationStore, firstAnnoId2).endPos) ||
								     (getValue(me.annotationStore, firstAnnoId1).beginPos > 
								     	getValue(me.annotationStore, firstAnnoId1).endPos && 
								      getValue(me.annotationStore, firstAnnoId1).endPos < 
								      	getValue(me.annotationStore, firstAnnoId2).beginPos)  ) 
								{	
									if (exact_nTuple && n <= length(annoIds)) create_nTuple(secTupleSet, me, annoIds, getValue(itP), n);
									else if (!exact_nTuple && n == 0) create_nTuple(secTupleSet, me, annoIds, getValue(itP), length(annoIds));		
									else if (!exact_nTuple) create_Tuple(secTupleSet, me, annoIds, getValue(itP), n);
								}
							}
						}
					}
					else validMate = true;
			
					// access to tupleCountStore for all tuple of current read:
					if (validMate)
					{
						itTuple = begin(tupleSet);
						itTupleEnd = end(tupleSet);
						for ( ; itTuple != itTupleEnd; goNext(itTuple))	
						{
							firstAnnoId1 = front(getValue(itTuple));
							erase(value(itTuple), 0);			// first id is not stored; is know by position in tupleCountStore

							// readConnections:
							if (!empty(getValue(itTuple)))
							{
								if (searchValue(pos, getValue(itTuple), getValue(tupleCountStore, firstAnnoId1).readConnections)) 
									++value(value(tupleCountStore, firstAnnoId1).readConnectionCounts, pos);
								else 
								{
									if (pos != endPosition(getValue(tupleCountStore, firstAnnoId1).readConnections) )
									{
										resizeSpace(value(tupleCountStore, firstAnnoId1).readConnections, 1, pos, pos, Generous());
										assignValue(value(tupleCountStore, firstAnnoId1).readConnections, pos, getValue(itTuple));
										insertValue(value(tupleCountStore, firstAnnoId1).readConnectionCounts, pos, 1, Generous());
									}
									else
									{
										appendValue(value(tupleCountStore, firstAnnoId1).readConnections, getValue(itTuple), Generous());
										appendValue(value(tupleCountStore, firstAnnoId1).readConnectionCounts, 1, Generous());
									}
								}
							}
							// matePairConnections: 
							if (!empty(secTupleSet))
							{
								itSecTuple = begin(secTupleSet);
								itSecTupleEnd = end(secTupleSet);
								for ( ; itSecTuple != itSecTupleEnd; goNext(itSecTuple) )
								{
									matePairTuple = getValue(itTuple);
									// INVALID_ID: sign for connection by matepair (apart from that, there are no INVALID_IDs in the list)
									appendValue(matePairTuple, INVALID_ANNO_ID, Generous());				
									if (!empty(getValue(itTuple)) && back(getValue(itTuple)) == front(getValue(itSecTuple)) )		// no id 2x allowed
									{	
										if (exact_nTuple == 0 && n == 0) erase(value(itSecTuple), 0);
										else continue;							// tupel would be created double or tupel wouldn't have the length n anymore
									}
									append(matePairTuple, getValue(itSecTuple), Generous());
					
									if (empty(getValue(itTuple))) 							
									{
										beginPos1 = getValue(me.annotationStore, firstAnnoId1).beginPos;
										endPos1 = getValue(me.annotationStore, firstAnnoId1).endPos;
									}
									else
									{
										beginPos1 = getValue(me.annotationStore, back(getValue(itTuple))).beginPos;
										endPos1 = getValue(me.annotationStore, back(getValue(itTuple))).endPos;
									}
									// begin position of first annotation in tuple of second read
									beginPos2 = getValue(me.annotationStore, front(getValue(itSecTuple))).beginPos; 
									endPos2 = getValue(me.annotationStore, front(getValue(itSecTuple))).endPos;
									if ( (beginPos1 <= endPos1 && endPos1 < beginPos2) ||			// no overlapping annotations allowed
									     (endPos1 < beginPos1 && beginPos1 < endPos2) )
									{
										if (searchValue(pos, matePairTuple, getValue(tupleCountStore, firstAnnoId1).matePairConnections))
											++value(value(tupleCountStore, firstAnnoId1).matePairConnectionCounts, pos);
										else 
										{
											if (pos != endPosition(getValue(tupleCountStore, firstAnnoId1).matePairConnections) )
											{
												resizeSpace(value(tupleCountStore, firstAnnoId1).matePairConnections, 1, pos, pos, Generous());
												assignValue(value(tupleCountStore, firstAnnoId1).matePairConnections, pos, matePairTuple);
												insertValue(value(tupleCountStore, firstAnnoId1).matePairConnectionCounts, pos, 1, Generous());
											}
											else
											{
												appendValue(value(tupleCountStore, firstAnnoId1).matePairConnections, matePairTuple, Generous());
												appendValue(value(tupleCountStore, firstAnnoId1).matePairConnectionCounts, 1, Generous());
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// buildAnnoCountStoreg++  -I../seqan/projects/library/ -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -O3 -pedantic -lrt  main.cpp   -o main
//////////////////////////////////////////////////////////////////////////////
template<typename TAnnoCountStore, typename TSpec, typename TConfig, typename TReadAnnoStore>
inline void
buildAnnoCountStore(TAnnoCountStore & annoCountStore, FragmentStore<TSpec, TConfig> & me, TReadAnnoStore & readAnnoStore)
{
	typedef typename Iterator<TReadAnnoStore>::Type 			TReadIter;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore 		TReadStore;
	typedef typename Value<TReadStore>::Type 				TReadStoreElement;
	typedef typename Value<TReadAnnoStore>::Type				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIds			TAnnoIds;
	typedef typename Iterator<TAnnoIds>::Type				TAnnoIdsIter;
	typedef typename Value<TAnnoIds>::Type					TIds;
	typedef typename Iterator<TIds>::Type 					TIdsIter;
	typedef typename Value<TIds>::Type					TId;
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	
	static const TId INVALID_READ_ID = TReadStoreElement::INVALID_ID;
	static const TId INVALID_ANNO_ID = TAnnotationStoreElement::INVALID_ID;
	
	resize(annoCountStore, length(me.annotationStore), 0);
	
	TReadIter itRead = begin(readAnnoStore);
	TReadIter itReadEnd = end(readAnnoStore);
	TId readId;
	TId matePairId;
	TId secReadId;
	TIds interSecIds;
	TIdsIter itP;
	TIdsIter itPEnd;
	TAnnoIdsIter itAnnoIds;
	TAnnoIdsIter itAnnoIdsEnd;
	TIdsIter itId;
	TIdsIter itIdEnd;
	
	// increment for each read respective to its mapped ids the count in the annoCountStore 
	for ( ; itRead != itReadEnd; goNext(itRead))
	{
		if (!empty(getValue(itRead).annoIds) )
		{
			readId = position(itRead, readAnnoStore);
			matePairId = getValue(me.readStore, readId).matePairId;
			if (matePairId != INVALID_READ_ID)
			{
				if (getValue(getValue(me.matePairStore, matePairId).readId, 0) == readId)
					secReadId = getValue(getValue(me.matePairStore, matePairId).readId, 1);
				else
					secReadId = getValue(getValue(me.matePairStore, matePairId).readId, 0);
			}
			// for each parentId: we just want to count annotations, in which the read mapped
			if (!empty(getValue(itRead).parentIds))
			{
				itP = begin(getValue(itRead).parentIds);
				itPEnd = end(getValue(itRead).parentIds);
				for (; itP != itPEnd; goNext(itP) )
				{
					// check mate-read to prevent double counts
					if (matePairId != INVALID_READ_ID)
					{
						if (!isElement_unsorted(getValue(itP), getValue(readAnnoStore, secReadId).parentIds) )
							continue; // if matepair read doesn't map in same parentId: no count (go to next parentId)
					
						// count annotations, which occur in both reads, shouldn't be increment for the read with the bigger readId 
						if (secReadId < readId )	 	
						{
							clear(interSecIds);								
							// just check the periphery annotations
							// if the end of the current read mapped in a same annotation as the start of the second read:  
							if ( interSec(interSecIds, back(getValue(itRead).annoIds), front(getValue(readAnnoStore, secReadId).annoIds)) ) 
							{
								for (unsigned i = 0; i < length(interSecIds); ++i)
								{
									if (getValue(me.annotationStore, getValue(interSecIds, i) ).parentId ==  getValue(itP))
										--value(annoCountStore, getValue(interSecIds, i));
									// decrement the corresponding count
								}
							}
							// or if the start of the current read mapped in a same annotation as the end of the second read: 
							else if ( interSec(interSecIds, front(getValue(itRead).annoIds), 
								  back(getValue(readAnnoStore, secReadId).annoIds)) )	
							{
								for (unsigned i = 0; i < length(interSecIds); ++i)
								{
									if (getValue(me.annotationStore, getValue(interSecIds, i) ).parentId ==  getValue(itP))
										--value(annoCountStore, getValue(interSecIds, i));
								}
							}	
							if (getValue(itP) != INVALID_ANNO_ID) --value(annoCountStore, getValue(itP));
						}
					}
			
					// count for all annoIds
					itAnnoIds = begin(getValue(itRead).annoIds);
					itAnnoIdsEnd = end(getValue(itRead).annoIds);
					for ( ; itAnnoIds != itAnnoIdsEnd; goNext(itAnnoIds))
					{
						itId = begin(getValue(itAnnoIds));
						itIdEnd = end(getValue(itAnnoIds));
						for ( ; itId != itIdEnd; goNext(itId))
							if (getValue(itId) != INVALID_ANNO_ID && getValue(me.annotationStore, getValue(itId)).parentId == getValue(itP) )
								++value(annoCountStore, getValue(itId));
					}
		
					// count for parentIds (already selected)
					if (getValue(itP) != INVALID_ANNO_ID) ++value(annoCountStore, getValue(itP));
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
////// Overlap Module
//////////////////////////////////////////////////////////////////////////////
template<typename TReadAnnoStore, typename TAnnoCountStore, typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
getResults(TReadAnnoStore & readAnnoStore,
	   TAnnoCountStore & annoCountStore,
	   TTupleCountStore & tupleCountStore, 
	   FragmentStore<TSpec, TConfig> & me, 
	   unsigned tupelSize,
	   bool exact_nTuple,
	   unsigned offsetInterval,
	   unsigned thresholdGaps,
	   bool unknownO)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	typedef typename FragmentStore<TSpec, TConfig>::TIntervalTreeStore 	TIntervalTreeStore;
	typedef typename Iterator<TIntervalTreeStore>::Type			TIntervalTree;
	typedef typename Value<TReadAnnoStore>::Type 				TReadAnnoStoreElement;
	typedef typename TReadAnnoStoreElement::TAnnoIds 			TAnnoIds;
	
	typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore	TAlignedReadStore;
	typedef typename Position<TAlignedReadStore>::Type 			TAlignPos;
	typedef 	 String<AlignIntervalsStoreElement<> > 			TAlignIntervalsStore;
	typedef typename Iterator<TAlignIntervalsStore>::Type 			TAlignIntervalsStoreIter;
	
	resize(readAnnoStore, length(me.readStore));
	
	TIntervalTree intervalTree;
	
	// extract intervals from alignedReadStore and store them in AlignIntervalsStore:
	TAlignIntervalsStore alignIntervalsStore;
	buildAlignIntervalsStore(alignIntervalsStore, me, thresholdGaps);

	if (!empty(alignIntervalsStore))
	{
		TAlignPos alignPos;
		TId contigId;
		TId readId;
		TAnnoIds  ids;
	
		TAlignIntervalsStoreIter it = begin(alignIntervalsStore);
		TAlignIntervalsStoreIter itEnd = end(alignIntervalsStore);
		// for each item in alignIntervalsStore:
		for ( ; it != itEnd; goNext(it))
		{
			// get ids from alignedReadStore (same position as in alignIntervalsStore):
			alignPos = position(it, alignIntervalsStore);
			contigId = getValue(me.alignedReadStore, alignPos).contigId;
			readId = getValue(me.alignedReadStore, alignPos).readId;
			// get respective intervalTree
			if (unknownO ||  getValue(me.alignedReadStore, alignPos).beginPos <= getValue(me.alignedReadStore, alignPos).endPos)
				intervalTree = begin(me.intervalTreeStore_F, Standard()) + contigId; 	//getValue(me.intervalTreeStore_F, contigId);
			else 
				intervalTree = begin(me.intervalTreeStore_R, Standard()) + contigId;        //getValue(me.intervalTreeStore_R, contigId);
			
			// get annotationStore-Ids for these intervals:
			clear(ids);
			if ((*intervalTree).interval_counter != 0)
				getIdsForRead(ids, me, *intervalTree, getValue(it).intervals, offsetInterval);
			// assign Ids from mapped annotations to readAnnoStore:
			value(readAnnoStore, readId).contigId = contigId;
			assignToReadAnnoStore(readAnnoStore, me, readId, ids);
		}
	}
	buildAnnoCountStore(annoCountStore, me, readAnnoStore);
	buildTupleCountStore(tupleCountStore, me, readAnnoStore, tupelSize, exact_nTuple);
}


//////////////////////////////////////////////////////////////////////////////
/// get normalized values for annotations And get Map for Gene orientations (necessary for annotation Output)
//////////////////////////////////////////////////////////////////////////////
template<typename TAnnoNormStore, typename TMapO, typename TAnnoCountStore, typename TSpec, typename TConfig>
inline void
normalizeAnnoCounts(TAnnoNormStore &annoNormStore, TMapO &mapO, TAnnoCountStore &annoCountStore, FragmentStore<TSpec, TConfig> &me)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId				TId;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore		TReadStore;
	typedef typename Size<TReadStore>::Type					TReadStoreSize; 
	typedef typename Iterator<TAnnotationStore>::Type			TAnnoIter;
	typedef typename Iterator<TAnnoCountStore>::Type 			TCountIter;
	typedef typename Iterator<TAnnoNormStore>::Type				TNormIter;
	typedef typename Size<TPos>::Type					TSize;
	typedef 	 String<TSize>						TChildrenLengths;
	typedef typename Iterator<String<TSize> >::Type				TLengthIter;	
	typedef 	 Pair<TId, TChildrenLengths>				TPair;
	typedef typename Value<TMapO>::Type					TPairO;
	typedef 	 Map<TPair>						TMap;
	typedef typename Iterator<TMap>::Type					TMapIter;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	static const TPos INVALID_POS = TAnnotationStoreElement::INVALID_POS;
	
	resize(annoNormStore, length(annoCountStore), 0);
	
	TReadStoreSize readNo = length(me.readStore) - length(me.matePairStore);
	
	TMap map;
	clear(map);
	clear(mapO);
	
	if(!empty(me.annotationStore))
	{
		TAnnoIter itA = begin(me.annotationStore);
		TAnnoIter itAEnd = end(me.annotationStore);
		TCountIter itC = begin(annoCountStore);
		TNormIter itN = begin(annoNormStore);
		TChildrenLengths childrenLengths;
		TPair pair;
		TPairO pairO;
		TSize length;
		
		for ( ; itA != itAEnd; goNext(itA), goNext(itC), goNext(itN) )
		{
			if (getValue(itA).beginPos == INVALID_POS && getValue(itA).parentId == INVALID_ID) 	// make entry for each gene/parent in map:
			{
				clear(childrenLengths);
				pair.i1 = position(itA, me.annotationStore);
				pair.i2 = childrenLengths;
				insert(map, pair);								// for lengths
		
				pairO.i1 = position(itA, me.annotationStore);
				pairO.i2 = 0;
				insert(mapO, pairO);								// for orientation
			}
			else if (getValue(itA).beginPos != INVALID_POS)						// for each exon/child: 
			{
				if (getValue(itA).beginPos <= getValue(itA).endPos)
					length = getValue(itA).endPos - getValue(itA).beginPos + 1;
				else
					length = getValue(itA).beginPos - getValue(itA).endPos + 1;
		
				value(itN) = ((double)1000000000 * (double)getValue(itC))/((double)readNo * (double)length);		// calculate normalized expression-value
		
				if (getValue(itA).parentId != INVALID_ID)					// append length to gene/parent lengths
				{
					appendValue(mapValue(map, getValue(itA).parentId), length, Generous());
					if (getValue(itA).beginPos > getValue(itA).endPos)
						mapValue(mapO, getValue(itA).parentId) = 1;
				}
			}
		}
	}
	
	if (!empty(map))
	{
		TMapIter itM = begin(map);
		TMapIter itMEnd = end(map);
		TSize length;
		TLengthIter itL;
		TLengthIter itLEnd;
		for ( ; itM != itMEnd; goNext(itM))	// calculate normalized gene/parent expression-values
		{
			length = 0;
			itL = begin(value(itM).i2);
			itLEnd = end(value(itM).i2);
			for ( ; itL != itLEnd; goNext(itL))
				length += getValue(itL);
			value(annoNormStore, value(itM).i1) = ((double)1000000000 * (double)getValue(annoCountStore, value(itM).i1) )/((double)readNo * (double)length);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
/// get normalized values for tuple
//////////////////////////////////////////////////////////////////////////////
template<typename TTupleCountStore, typename TSpec, typename TConfig>
inline void
normalizeTupleCounts(TTupleCountStore &tupleCountStore, FragmentStore<TSpec, TConfig> &me)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId				TId;
	typedef typename Value<TTupleCountStore>::Type				TTupleCountStoreElement;
	typedef typename TTupleCountStoreElement::TTupleList			TTupleList;
	typedef typename TTupleCountStoreElement::TTupleCounts			TTupleCounts;
	typedef typename TTupleCountStoreElement::TTupleNorm			TTupleNorm;
	typedef typename TTupleCountStoreElement::TTuple			TTuple;
	typedef typename Iterator<TTupleCountStore>::Type 			TStoreIter;
	typedef typename Iterator<TTupleList>::Type				TTupleListIter;
	typedef typename Iterator<TTupleCounts>::Type				TCountIter;
	typedef typename Iterator<TTupleNorm>::Type				TNormIter;
	typedef typename Iterator<TTuple>::Type					TTupleIter;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos		TPos;
	typedef typename Size<TPos>::Type					TSize;
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore		TReadStore;
	typedef typename Size<TReadStore>::Type					TReadStoreSize; 
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	TReadStoreSize readNo = length(me.readStore) - length(me.matePairStore);
	
	if (!empty(tupleCountStore))
	{
		TStoreIter itS = begin(tupleCountStore);
		TStoreIter itSEnd = end(tupleCountStore);
		TTupleListIter itT;
		TTupleListIter itTEnd;
		TCountIter itC;
		TNormIter itN;
		TSize tupleLength;
		TTupleIter itId;
		TTupleIter itIdEnd;
		for ( ; itS != itSEnd; goNext(itS))
		{
			// readConnections:
			resize(value(itS).readConnectionNorm, length(getValue(itS).readConnections));
			if (!empty(getValue(itS).readConnections))
			{
				itT = begin(getValue(itS).readConnections);
				itTEnd = end(getValue(itS).readConnections);
				itC = begin(getValue(itS).readConnectionCounts);
				itN = begin(getValue(itS).readConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{	
						if (getValue(me.annotationStore, getValue(itId)).beginPos <= getValue(me.annotationStore, getValue(itId)).endPos)
							tupleLength += getValue(me.annotationStore, getValue(itId)).endPos - getValue(me.annotationStore, getValue(itId)).beginPos + 1;
						else
							tupleLength += getValue(me.annotationStore, getValue(itId)).beginPos - getValue(me.annotationStore, getValue(itId)).endPos + 1;
					}
					value(itN) = ((double)1000000000 * (double)getValue(itC)) / ((double)readNo * (double)tupleLength);
				}
			}
			// matePairConnections:
			resize(value(itS).matePairConnectionNorm, length(getValue(itS).matePairConnections));
			if (!empty(getValue(itS).matePairConnections))
			{
				itT = begin(getValue(itS).matePairConnections);
				itTEnd = end(getValue(itS).matePairConnections);
				itC = begin(getValue(itS).matePairConnectionCounts);
				itN = begin(getValue(itS).matePairConnectionNorm);
				for ( ; itT != itTEnd; goNext(itT), goNext(itC), goNext(itN))
				{
					tupleLength = 0;
					itId = begin(getValue(itT));
					itIdEnd = end(getValue(itT));
					for ( ; itId != itIdEnd; goNext(itId))
					{
						if (getValue(itId) != INVALID_ID)
						{
							if (getValue(me.annotationStore, getValue(itId)).beginPos <= getValue(me.annotationStore, getValue(itId)).endPos)
								tupleLength += getValue(me.annotationStore, getValue(itId)).endPos - getValue(me.annotationStore, getValue(itId)).beginPos + 1;
							else
								tupleLength += getValue(me.annotationStore, getValue(itId)).beginPos - getValue(me.annotationStore, getValue(itId)).endPos + 1;
						}
					}
					value(itN) = ((double)1000000000 * (double)getValue(itC)) / ((double)readNo * (double)tupleLength);
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
/// NGS Overlapper main function
//////////////////////////////////////////////////////////////////////////////

inline bool
ngsOverlapper(CharString const &nameSAM, CharString const &nameGFF, CharString const &outputPath,
	      unsigned tupleValue, bool exact_nTuple,
	      unsigned thresholdGaps, unsigned offsetInterval,
	      unsigned thresholdCount, double thresholdRPKM,
	      bool unknownO,
	      bool fusion,
	      bool gtf)
{	
	FragmentStore<> me;
#ifdef DEBUG_OVERLAP_MODULE
	SEQAN_PROTIMESTART(find1_time);	
#endif 
	// build contigStore from FASTA file
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "load Sam..." << std::endl;
#endif 
	// read aligned reads in FragmentStore from Sam files
	std::fstream fileSAM;
	fileSAM.open(toCString(nameSAM), std::ios_base::in | std::ios_base::binary);
	if(!fileSAM.is_open()) return false;
	read(fileSAM, me, Sam());
	fileSAM.close();

	// readAnnotations from Gff or Gtf:
	if (gtf == 0)
	{
#ifdef DEBUG_OVERLAP_MODULE
		SEQAN_PROTIMESTART(find2_time);	
		std::cout << "load Gff..." << std::endl;
#endif 
		readAnnotationsFromGFF(me, toCString(nameGFF));
	}
	else
	{
#ifdef DEBUG_OVERLAP_MODULE
		SEQAN_PROTIMESTART(find2_time);	
		std::cout << "load Gff..." << std::endl;
#endif 
		readAnnotationsFromGTF(me, toCString(nameGFF));
	}

	// create IntervalTreeStore:
#ifdef DEBUG_OVERLAP_MODULE
	SEQAN_PROTIMESTART(find3_time);
#endif 
	createIntervalTreeStore(me, unknownO);
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "create intervalTreeStores from annotationStore took: \t" << SEQAN_PROTIMEDIFF(find3_time) << " seconds" << std::endl;
#endif 
	
	// build stores for results:
	String<ReadAnnoStoreElement<unsigned> >			readAnnoStore; 	
	String<unsigned> 					annoCountStore;
	String<TupleCountStoreElement<unsigned> >		tupleCountStore;
	String<TupleCountStoreElement_Fusion<unsigned> >	tupleCountStore_Fusion;			// additional Store, if fusion genes should be checked



	// get results with additional check for transfusion genes (will be changed later additionally) 
	if (fusion == 1)
		getResults_Fusion(readAnnoStore, annoCountStore, tupleCountStore, tupleCountStore_Fusion, me, tupleValue, exact_nTuple, offsetInterval, thresholdGaps, unknownO);
	else // get normal results:
		getResults(readAnnoStore, annoCountStore, tupleCountStore, me, tupleValue, exact_nTuple, offsetInterval, thresholdGaps, unknownO);
	
	
	// normalize:
	String<double>			annoNormStore;
	Map<Pair<unsigned, bool> > 	mapO;
	normalizeAnnoCounts(annoNormStore, mapO, annoCountStore, me);
	normalizeTupleCounts(tupleCountStore, me);
	if (fusion == 1)
		normalizeTupleCounts_Fusion(tupleCountStore_Fusion, me);


	// output:
	CharString pathReadOutput = outputPath;
	append(pathReadOutput, "readOutput.gff");
	std::fstream readOutput;
	readOutput.open(toCString(pathReadOutput), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
	createReadCountGFF(readOutput, readAnnoStore, me);
	readOutput.close();
	
	CharString pathAnnoOutput = outputPath;
	append(pathAnnoOutput, "annoOutput.gff");
	std::fstream annoOutput;
	annoOutput.open(toCString(pathAnnoOutput), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
	createAnnoCountGFF(annoOutput, annoCountStore, annoNormStore, me, mapO);
	annoOutput.close();
	
	CharString pathTupleOutput = outputPath;
	append(pathTupleOutput, "tupleOutput.gff");
	std::fstream tupleOutput;
	tupleOutput.open(toCString(pathTupleOutput), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
	createTupleCountGFF(tupleOutput, tupleCountStore, me, thresholdCount, thresholdRPKM);
	tupleOutput.close();

	// additional output, if fusion genes were checked
	if (fusion == 1)
	{
		CharString pathTupleOutput_Fusion = outputPath;
		append(pathTupleOutput_Fusion, "tupleOutput_Fusion.gff");
		std::fstream tupleOutput_Fusion;
		tupleOutput_Fusion.open(toCString(pathTupleOutput_Fusion), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
		createTupleCountGFF_Fusion(tupleOutput_Fusion, tupleCountStore_Fusion, me, thresholdCount, thresholdRPKM);
		tupleOutput_Fusion.close();
	}

	
#ifdef DEBUG_OVERLAP_MODULE
	std::cout << "ngsOverlapper-function took: \t" << SEQAN_PROTIMEDIFF(find1_time) << " seconds" << std::endl;
	std::cout << "ngsOverlapper-function without reading Sam took: \t" << SEQAN_PROTIMEDIFF(find2_time) << " seconds" << std::endl;
	std::cout << "ngsOverlapper-function and create IntervalTreeStore without reading Sam took:\t" << SEQAN_PROTIMEDIFF(find2_time) - SEQAN_PROTIMEDIFF(find3_time) << " seconds" << std::endl;
#endif 
	return true;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
