/*==========================================================================
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

//#define RAZERS_DUMP_SNPS
//#define COVBASED_INDELCALLING

#define TRACE_PIPELINE
//#define READS_454

//#define SNPSTORE_DEBUG
//#define SNPSTORE_DEBUG_CANDPOS

#ifdef SNPSTORE_DEBUG
#define SNPSTORE_DEBUG_CANDPOS
#endif

#include "seqan/platform.h"
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/consensus.h>

#ifdef PLATFORM_WINDOWS
#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
#define SEQAN_DEFAULT_TMPDIR "./"
#endif

//#include "/home/takifugu2/emde/seqan/seqan-trunk/projects/library/apps/razers/outputFormat.h"
#include "snp_store.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using namespace std;
using namespace seqan;

// get reference file names
template<typename TOptions>
int getGenomeFileNameList(char const * filename, StringSet<CharString> & genomeFileNames, TOptions &options)
{
	::std::ifstream file;
	file.open(filename,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GENOME_FAILED;
	
	CharString nameStr;
	char c = _streamGet(file);
	if (c != '>' && c != '@')	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		if(options._debugLevel >=1)
			::std::cout << ::std::endl << "Reading multiple genome files:" <<::std::endl;
		/*		//locations of genome files are relative to list file's location
		 ::std::string tempGenomeFile(filename);
		 size_t lastPos = tempGenomeFile.find_last_of('/') + 1;
		 if (lastPos == tempGenomeFile.npos) lastPos = tempGenomeFile.find_last_of('\\') + 1;
		 if (lastPos == tempGenomeFile.npos) lastPos = 0;
		 ::std::string filePrefix = tempGenomeFile.substr(0,lastPos);*/
		unsigned i = 1;
		while(!_streamEOF(file))
		{
			clear(nameStr);
			_parseSkipWhitespace(file, c);
			while ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r'))
			{
				appendValue(nameStr,c,Generous());
				c = _streamGet(file);
			}
			appendValue(genomeFileNames,nameStr,Generous());
			if(options._debugLevel >=2)
				::std::cout <<"Genome file #"<< i <<": " << genomeFileNames[length(genomeFileNames)-1] << ::std::endl;
			++i;
			_parseSkipWhitespace(file, c);
		}
		if(options._debugLevel >=1)
			::std::cout << i-1 << " genome files total." <<::std::endl;
	}
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames,filename,Generous());
	file.close();
	return 0;
	
}

template<typename TFragmentStore, typename TStr>
void
_dumpMatches(TFragmentStore &fragmentStore, TStr str)
{
	typedef typename TFragmentStore::TAlignedReadStore 			TMatches;
	typedef typename Value<TMatches>::Type 						TMatch;
	typedef typename TFragmentStore::TAlignQualityStore 		TMatchQualities;
	typedef typename Value<TMatchQualities>::Type 				TMatchQuality;
	typedef typename TFragmentStore::TReadSeqStore	 			TReads;
	typedef typename Value<TReads>::Type 						TRead;
	typedef typename Iterator<TReads,Standard>::Type			TReadIt;
	typedef typename Iterator<TMatchQualities,Standard>::Type	TMatchQIt;
	typedef typename Iterator<TMatches,Standard>::Type			TMatchIt;
	
	std::cout << "Length of matches = " << length(fragmentStore.alignedReadStore)  << "\n";
	std::cout << "Length of reads   = " << length(fragmentStore.readSeqStore)  << "\n";
	std::cout << "Length of matchqs = " << length(fragmentStore.alignQualityStore)  << "\n";
	
	for(unsigned i = 0 ; i < length(fragmentStore.alignedReadStore); ++i)
	{
		char ori = (fragmentStore.alignedReadStore[i].beginPos < fragmentStore.alignedReadStore[i].endPos) ? 'F' : 'R';
		std::cout << "--"<<str<<"Match number " << i << ":\n";
		std::cout << "--"<<str<<"MatchId  = " << fragmentStore.alignedReadStore[i].id << "\n";
		std::cout << "--"<<str<<"ReadId   = " << fragmentStore.alignedReadStore[i].readId << "\n";
		std::cout << "--"<<str<<"ContigId = " << fragmentStore.alignedReadStore[i].contigId << std::flush << "\n";
		std::cout << "--"<<str<<"gBegin   = " << _min(fragmentStore.alignedReadStore[i].beginPos, fragmentStore.alignedReadStore[i].endPos) << "\n";
		std::cout << "--"<<str<<"gEnd     = " << _max(fragmentStore.alignedReadStore[i].beginPos, fragmentStore.alignedReadStore[i].endPos) << "\n";
		std::cout << "--"<<str<<"orient   = " << ori << std::flush << std::endl;
		if(length(fragmentStore.alignQualityStore) > fragmentStore.alignedReadStore[i].id)
		{
			std::cout << "--"<<str<<"EditDist = " << (int) fragmentStore.alignQualityStore[fragmentStore.alignedReadStore[i].id].errors << "\n";
			std::cout << "--"<<str<<"AvgQ     = " << (int)fragmentStore.alignQualityStore[fragmentStore.alignedReadStore[i].id].score << "\n";
		}
		std::cout << "--"<<str<<"Readseq  = " << fragmentStore.readSeqStore[fragmentStore.alignedReadStore[i].readId] << std::flush << "\n";
		
	}
}



// load entire genome into memory
template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList, ::std::map<CharString,unsigned> &gIdStringToIdNumMap)
{
	unsigned gSeqNo = 0;
	unsigned filecount = 0;
	CharString temp;
	while(filecount < length(fileNameList))
	{
		clear(temp);
		MultiFasta multiFasta;
		if (!open(multiFasta.concat, toCString(fileNameList[filecount]), OPEN_RDONLY)) return false;
		split(multiFasta, Fasta());
		
		unsigned seqCount = length(multiFasta);
		if(length(genomes) < gSeqNo+seqCount) 
			resize(genomes,gSeqNo+seqCount);
		for(unsigned i = 0; i < seqCount; ++i)
		{
			assignSeq(genomes[gSeqNo+i], multiFasta[i], Fasta());		// read Genome sequence
			assignSeqId(temp, multiFasta[i], Fasta());
			for (unsigned pos = 0; pos < length(temp); ++pos)
			{
				if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
				{
					resize(temp,pos);
					break;
				}
			}
			gIdStringToIdNumMap.insert(::std::make_pair<CharString,unsigned>(temp,gSeqNo+i)); //TODO shortID
		}
		gSeqNo += seqCount;
		++filecount;
	}
	resize(genomes,gSeqNo);
	return (gSeqNo > 0);
}





//read filename (read line and trim trailing whitespaces)
template<typename TFile, typename TChar, typename TString>
void
_parseReadWordUntilWhitespace(TFile& file, TString& str, TChar& c)
{
	append(str,c);
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
		append(str, c);
	}
	return;
}



//read filename (read line and trim trailing whitespaces)
template<typename TFile, typename TChar>
void
_parse_skipUntilWhitespace(TFile& file, TChar& c)
{
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
	}
	return;
}




// transform global cooridnates to coordinates relative to chromosomal segment
template<typename TFragmentStore, typename TContigPos, typename TOptions>
void 
transformCoordinates(TFragmentStore &fragmentStore, TContigPos startCoord, TOptions&)
{
	typedef typename TFragmentStore::TAlignedReadStore 			TMatches;
	typedef typename Value<TMatches>::Type 						TMatch;
	typedef typename Iterator<TMatches,Standard>::Type			TMatchIt;
	
	TMatchIt mIt		= begin(fragmentStore.alignedReadStore,Standard());
	TMatchIt mItEnd 	= end(fragmentStore.alignedReadStore,Standard());
	
	while(mIt != mItEnd)
	{
		(*mIt).endPos -= startCoord;
		(*mIt).beginPos -= startCoord;
		++mIt;
	}
	
}

// copy matches relevant for next window
template<typename TFragmentStore, typename TReadCigars, typename TReadCounts, typename TReadClips, typename TSize, typename TContigPos, typename TOptions>
void 
copyNextWindowMatchesAndReads(TFragmentStore &fragmentStore,
							  TReadCounts &readCounts,
							  TReadCigars &readCigars,
							  TReadCounts &tmpReadCounts,
							  typename TFragmentStore::TReadSeqStore &tmpReads,
							  typename TFragmentStore::TReadStore &tmpRs,
							  typename TFragmentStore::TAlignedReadStore &tmpMatches,
							  typename TFragmentStore::TAlignQualityStore &tmpQualities,
							  TReadClips &tmpReadClips,
							  TReadCigars &tmpReadCigars,
							  TSize ,
							  TContigPos currentWindowEnd,
							  TOptions &options)
{
	
	typedef typename TFragmentStore::TAlignedReadStore 			TMatches;
	typedef typename Value<TMatches>::Type 						TMatch;
	typedef typename Iterator<TMatches,Standard>::Type			TMatchIt;
	typedef typename Id<TFragmentStore>::Type					TId;
	typedef typename Value<TReadClips>::Type 					TPair;
	
	SEQAN_ASSERT(length(fragmentStore.readSeqStore) == length(fragmentStore.alignQualityStore))
	
	::std::sort(begin(fragmentStore.alignedReadStore, Standard()), end(fragmentStore.alignedReadStore, Standard()), LessGPos<TMatch>());	
	
	if(options._debugLevel > 1 )::std::cout << "Copying matches overlapping more than one window ... \n";
	
	TMatchIt mIt		= end(fragmentStore.alignedReadStore,Standard());
	TMatchIt mItBegin 	= begin(fragmentStore.alignedReadStore,Standard());
	--mIt;
	
	options.minCoord = MaxValue<unsigned>::VALUE;
	options.maxCoord = 0;
	//CharString str = "discBef";
	//_dumpMatches(fragmentStore, str);
	
	// look for matches that are inside our window of interest, copy corresponding matches,reads,qualities
	while(mIt >= mItBegin && _min((*mIt).beginPos,(*mIt).endPos) + (TContigPos)options.maxHitLength + (TContigPos)options.windowBuff >= currentWindowEnd )
	{
		if( _max((*mIt).beginPos,(*mIt).endPos) + (TContigPos)options.windowBuff > currentWindowEnd )
		{
			TId id = length(tmpMatches);
			appendValue(tmpMatches,*mIt);
			tmpMatches[id].id = id;
			tmpMatches[id].readId = id;
			appendValue(tmpReads,fragmentStore.readSeqStore[(*mIt).id]);
			appendValue(tmpRs,fragmentStore.readStore[(*mIt).id]);
			if(!empty(readCounts))appendValue(tmpReadCounts,readCounts[(*mIt).id]);
			appendValue(tmpQualities,fragmentStore.alignQualityStore[(*mIt).id]);
			appendValue(tmpReadCigars,readCigars[(*mIt).id]);
			appendValue(tmpReadClips,TPair(0,0));
			if(_max((*mIt).beginPos,(*mIt).endPos) > (TContigPos)options.maxCoord) options.maxCoord = (unsigned) _max((*mIt).beginPos,(*mIt).endPos);
			if(_min((*mIt).beginPos,(*mIt).endPos) < (TContigPos)options.minCoord) options.minCoord = (unsigned) _min((*mIt).beginPos,(*mIt).endPos);
			
		}
		--mIt;
	}
	
	if(options._debugLevel > 1)
		std::cout << length(tmpMatches)<<" matches left over from previous window.\n";
	
	
}


// little helper
template<typename TMatch>
char
orientation(TMatch & match)
{
	if(match.endPos > match.beginPos) 
		return 'F';
	else
		return 'R';
}

// discard reads stacking up, give preference to high quality reads
template<typename TFragmentStore, typename TSize, typename TOptions>
void
applyPileupCorrection(TFragmentStore	&fragmentStore, 
					  TSize							arrayBeginPos, 
					  TSize							arrayEndPos, 
					  TOptions						&options)
{
	typedef StringSet<String<Dna5Q>, Owner<ConcatDirect<> > >	 TFalseType;
	typedef typename TFragmentStore::TAlignedReadStore 	TMatches;
	typedef typename Value<TMatches>::Type 				TMatch;
	typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
	typedef typename Value<TMatchQualities>::Type 		TMatchQuality;
	typedef typename TFragmentStore::TReadSeqStore	 	TReads;
	typedef typename Value<TReads>::Type 				TRead;
	typedef typename TFragmentStore::TContigStore	 	TGenomeSet;
	typedef typename Value<TGenomeSet>::Type 			TGenome;
	typedef typename TFragmentStore::TContigPos 		TContigPos;
	typedef typename Iterator<TMatches,Standard>::Type	TMatchIterator;
	
	if(IsSameType<TReads, TFalseType >::VALUE)
		std::cout << "Hier stimmt was nciht. strinsetspec concat direct\n";
	
	if(options._debugLevel > 0) std::cout << arrayEndPos - arrayBeginPos  << " matches subject to pile up correction." << std::endl;
	
	//CharString str = "pileBef";
	//_dumpMatches(fragmentStore, str);
	
	if(options.orientationAware) 
		::std::sort(iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard()),
					iter(fragmentStore.alignedReadStore, arrayEndPos, Standard()), 
					LessGStackOaMQ<TMatches,TMatchQualities>(fragmentStore.alignQualityStore));	
	else
		::std::sort(iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard()),
					iter(fragmentStore.alignedReadStore, arrayEndPos, Standard()), 
					LessGStackMQ<TMatches,TMatchQualities>(fragmentStore.alignQualityStore));
	
	
	TMatchIterator matchIt			= iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
	TMatchIterator matchRangeEnd	= iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
	TMatchIterator matchItKeep		= matchIt;
	
	while(matchIt != matchRangeEnd)
	{
		TContigPos currentBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
		TContigPos currentEnd = _max((*matchIt).beginPos,(*matchIt).endPos);
		unsigned currentSeqno = (*matchIt).contigId;
		char currentOrientation = orientation(*matchIt);
		unsigned currPile = 0;
		while(matchIt != matchRangeEnd 
			  && (*matchIt).contigId == currentSeqno 
			  && _min((*matchIt).beginPos,(*matchIt).endPos) == currentBegin 
			  && _max((*matchIt).beginPos,(*matchIt).endPos) == currentEnd 
			  && (!options.orientationAware || orientation(*matchIt) == currentOrientation) 
			  && currPile < options.maxPile)
		{
			*matchItKeep = *matchIt;
			++currPile;
			++matchIt;
			++matchItKeep;
		}
		//if(matchRangeEnd > matchItEnd) ::std::cerr <<"neeeeeeee\n";
		while(matchIt != matchRangeEnd 
			  && (*matchIt).contigId == currentSeqno 
			  && _min((*matchIt).beginPos,(*matchIt).endPos) == currentBegin 
			  && _max((*matchIt).beginPos,(*matchIt).endPos) == currentEnd 
			  && (!options.orientationAware || orientation(*matchIt) == currentOrientation))
			++matchIt;
		
	}
	
	if(options._debugLevel > 0) std::cout << matchIt - matchItKeep << " matches discarded." << std::endl;
	resize(fragmentStore.alignedReadStore,matchItKeep - begin(fragmentStore.alignedReadStore,Standard()));
	
	//	::std::cout << "numMatches = " << length(fragmentStore.alignedReadStore) << ::std::endl;
	SEQAN_ASSERT(length(fragmentStore.alignedReadStore) <= length(fragmentStore.alignQualityStore))
	SEQAN_ASSERT(length(fragmentStore.readSeqStore) == length(fragmentStore.alignQualityStore))
	
	//	str="pileAft";
	//	_dumpMatches(fragmentStore,str);
	
}

// average quality of read is kept as extra info for each match
// used for prioritization in pile up correction
template<typename TFragmentStore, typename TSize, typename TOptions>
void
addReadQualityToMatches(TFragmentStore	&fragmentStore, 
						TSize							arrayBeginPos,
						TSize							arrayEndPos,
						TOptions &)
{
	typedef typename TFragmentStore::TAlignedReadStore		TMatches;
	typedef typename Value<TMatches>::Type 					TMatch;
	typedef typename TFragmentStore::TReadSeqStore	 		TReads;
	typedef typename Value<TReads>::Type 					TRead;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	TIterator it = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
	TIterator itEnd = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
	
	int avgRQ;
	for (; it != itEnd; ++it) 
	{
		TRead const &read = fragmentStore.readSeqStore[(*it).readId];
		avgRQ = 0;
		for(unsigned i = 0; i < length(read); ++i)
			avgRQ += (int) getQualityValue(read[i]);
		(fragmentStore.alignQualityStore[(*it).id]).score = (char)(avgRQ/length(read));
	}
	
}

// checks whether an alignment has indels
template<typename TValue, typename TSpec>
bool alignmentHasIndel(Align<TValue,TSpec> &align)
{
	typedef Align<TValue,TSpec>  TAlign;
	typedef typename Row<TAlign>::Type 		TAlignRow;
	typedef typename Iterator<TAlignRow>::Type 	TAlignIterator;
	
	bool hasIndel = false;
	for(unsigned i = 0; !hasIndel && i < length(rows(align)); ++i)
	{
		TAlignIterator it = begin(row(align,i));
		TAlignIterator itEnd = end(row(align,i));
		while (it != itEnd)
		{
			if(isGap(it))
			{
				hasIndel = true;
				break;
			}
			++it;
		}
	}
	return hasIndel;
}

// perform read clipping
template<typename TFragmentStore, typename TReadClips, typename TSize, typename TOptions>
void
clipReads(TFragmentStore 	&fragmentStore,
		  TReadClips 	&readClips,
		  TSize		arrayBeginPos,
		  TSize		arrayEndPos,
		  TOptions 	&options)
{
	typedef typename TFragmentStore::TAlignedReadStore		TMatches;
	typedef typename Value<TMatches>::Type 				TMatch;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;	 	// TMatchQualities
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
	typedef typename TFragmentStore::TReadSeqStore	 		TReads;
	typedef typename Value<TReads>::Type 				TRead;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	typedef typename TFragmentStore::TContigSeq			TContigSeq;				// TGenome
	typedef typename Position<TContigSeq>::Type			TContigPos;
	
	TIterator it = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
	TIterator itEnd = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
	Align<TRead, ArrayGaps> align;
	resize(rows(align), 2);
#ifdef SNPSTORE_DEBUG
	bool extraV = false;
#endif
	
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);
	if(length(readClips) < (arrayEndPos-arrayBeginPos)) ::std::cout << length(readClips) << " readclips but " << (arrayEndPos-arrayBeginPos) << " many reads.\n";
	TContigSeq &currGenome = fragmentStore.contigStore[0].seq;
	int kickout = 0;
	for (; it != itEnd; ++it) 
	{
		TMatch &m = *it;
		TRead &read = fragmentStore.readSeqStore[m.readId];
		int clipLeft = readClips[m.readId].i1;
		int clipRight = readClips[m.readId].i2;
		TContigPos beginPos = (m.beginPos < m.endPos ) ? m.beginPos : m.endPos;
		TContigPos endPos = (m.beginPos > m.endPos ) ? m.beginPos : m.endPos;
		TAlignQuality &aliQ = fragmentStore.alignQualityStore[m.id];
		
#ifdef SNPSTORE_DEBUG
		TContigPos beginBefore = beginPos;
#endif
		if(m.id != m.readId) ::std::cout << "match id != readId \n";
		if(clipLeft+clipRight > (int)length(read) || clipLeft > (int)length(read) || clipRight > (int)length(read))
		{
			if(options._debugLevel > 0)::std::cout << "WARNING: clipLeft+clipRight > readLen!!\n";
#ifdef SNPSTORE_DEBUG
			::std::cout << "readlength = "<<length(read)<< " \n";
			::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
			::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
			::std::cout << "read=" << read << std::endl;
			::std::cout << "beginPos=" << beginPos << std::endl;
			::std::cout << "endPos=" << endPos << std::endl; 
#endif
			clipLeft = length(read);
			clipRight = 0;
		}
		if(clipLeft > 0 || clipRight > 0)
		{
			//	if(extraV) ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << std::endl;
			if((int)length(read)-clipLeft-clipRight < options.minClippedLength)
			{
				if(options._debugLevel > 1 )
					::std::cout << "Discarded: "<<read<<" at position "<< beginPos <<"\n";
				m.endPos = m.beginPos = 0; // "mask" read
				read = "";
				++kickout;
				continue;
			}
			// adjust read sequence
			read = infix(read,clipLeft,length(read)-clipRight);
			
			// upate avg read quality
			int avgRQ = 0;
			for(unsigned i = 0; i < length(read); ++i)
				avgRQ += (int) getQualityValue(read[i]);
			aliQ.score = (char)(avgRQ/length(read));
			//		if(extraV) ::std::cout << "aliQ.score = " << (int)aliQ.score << ::std::endl;
			
			
			//do semi-global alignment
			assignSource(row(align, 0), read);
			assignSource(row(align, 1), infix(currGenome, beginPos, endPos));
			if ((*it).endPos < (*it).beginPos)
				reverseComplement(source(row(align, 1)));
			
			int score = globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());		
			aliQ.errors = (unsigned char) round(-score/1000);
			
#ifdef SNPSTORE_DEBUG
			if(extraV) ::std::cout << align << std::endl;
			if(extraV) ::std::cout << "aliQ.errors = " << aliQ.errors << std::endl;	
#endif
			
			// transform first and last read character to genomic positions
			if(aliQ.pairScore == 1)
			{
				unsigned viewPosReadFirst = toViewPosition(row(align, 0), 0);
				unsigned viewPosReadLast  = toViewPosition(row(align, 0), length(read) - 1);
				unsigned genomePosReadFirst = toSourcePosition(row(align, 1), viewPosReadFirst);
				unsigned genomePosReadLast  = toSourcePosition(row(align, 1), viewPosReadLast);
				//				if(isGap(row(align,1),viewPosReadFirst))
				//				{
				//					std::cout << "bgein gap --> do nothing " << std::endl;
				//					
				//				}
				
				if(isGap(row(align,1),viewPosReadLast))
				{
					genomePosReadLast--;					
				}
#ifdef SNPSTORE_DEBUG
				if(extraV)
				{	::std::cout << "viewPosReadFirst = " << viewPosReadFirst << ::std::endl;
					::std::cout << "viewPosReadLast = " << viewPosReadLast << ::std::endl;
					::std::cout << "genomePosReadFirst = " << genomePosReadFirst << ::std::endl;
				}
#endif
				if(m.beginPos < m.endPos) // forward
				{
					endPos = beginPos + (genomePosReadLast + 1);
					beginPos += genomePosReadFirst;
				}
				else
				{
					beginPos = endPos - genomePosReadLast - 1;
					endPos = endPos - genomePosReadFirst;
				}
				
				if(!alignmentHasIndel(align)) aliQ.pairScore = 0;
			}
			else
			{
				if(m.beginPos < m.endPos) // forward
				{
					endPos -= clipRight;
					beginPos += clipLeft;
				}
				else
				{
					endPos -= clipLeft;
					beginPos += clipRight;
				}
			}
			
			// set genomic positions
			if(m.beginPos < m.endPos) // forward
			{
				m.endPos = endPos;
				m.beginPos = beginPos;
			}
			else
			{
				m.endPos = beginPos;
				m.beginPos = endPos;
			}
			if(beginPos > 300000000 || endPos > 300000000)
			{
				::std::cout << "bgein groesser 300mio neu beginPos = "<< beginPos << " endpos=" << endPos << ::std::endl;
#ifdef SNPSTORE_DEBUG
				::std::cout << "WARNING: clipLeft+clipRight > readLen!!??\n";
				::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
				::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
				::std::cout << "read=" << read << std::endl;
				::std::cout << "beginPos before=" << beginBefore << std::endl;
				::std::cout << "beginPos=" << beginPos << std::endl;
				::std::cout << "endPos=" << endPos << std::endl;
#endif
			}
			
			
		}
	}
	if(options._debugLevel > 0) 
		::std::cout << kickout <<" reads too short after clipping, discarding!\n";
	
}






/////////////////////////////////////////////////////////////
// read sorted(!) Gff input file containing mapped reads
template <
typename TFile,
typename TFragmentStore,
typename TReadCounts,
typename TCigarStr,
typename TGenome,
typename TGenomeIdMap,
typename TContigPos,
typename TSize,
typename TValue,
typename TOptions
>
int readMatchesFromGFF_Batch(
							 TFile			 		&file,
							 TFragmentStore 				&fragmentStore,				// forward/reverse fragmentStore.alignedReadStore
							 TReadCounts				&readCounts,
							 String<Pair<int,int> >			&readClips,
							 StringSet<TCigarStr>			&readCigars,
							 TGenome					&genome,
							 TGenomeIdMap				&gIdStringToIdNumMap,
							 TSize					currSeqNo,
							 TContigPos				currentBegin,
							 TContigPos				currentEnd,
							 TValue					&highestChrId,
							 TOptions				&options)
{
	
	
	typedef typename TFragmentStore::TAlignedReadStore 	TMatches;
	typedef typename Value<TMatches>::Type 			TMatch;
	typedef typename TFragmentStore::TAlignQualityStore 	TMatchQualities;
	typedef typename Value<TMatchQualities>::Type 		TMatchQuality;
	typedef typename TFragmentStore::TReadSeqStore		TReads;
	typedef typename Value<TReads>::Type 			TRead;
	typedef typename TFragmentStore::TReadStore		TReadStore;
	typedef typename Value<TReadStore>::Type		TReadStoreElement;
	typedef typename Value<TCigarStr>::Type			TCigar;
	typedef typename Value<TReads>::Type 			TRead;
	typedef typename TFragmentStore::TContigStore	 	TGenomeSet;
	typedef typename Id<TFragmentStore>::Type 		TId;
	typedef typename Iterator<TMatches,Standard>::Type 	TMatchIterator;
	
	
	if(length(fragmentStore.readSeqStore)!=length(fragmentStore.alignQualityStore))
	{
		::std::cerr << "Lengths need to be equal!!\n";
		return 10;
	}
	int readCount = length(fragmentStore.readSeqStore);
	TContigPos genomeLen = length(genome);
	
	// general stuff that is needed
	typename TGenomeIdMap::const_iterator it;
	unsigned rSeq = readCount;
	Dna5String gInf;
	String<Dna5Q> curr_read;
	TCigarStr tmpCigarStr;
	CharString readTemplate, temp_read;
	CharString readName, temp_str;
	String<int> gAliPos;
	
	//	bool test = true;
	
	char c = _streamGet(*file);
	while (!_streamEOF(*file))
	{
		// our razers gff output looks like this:
		//X       razers      read            100919085       100919120       2       +       .       ID=s_3_1_3;unique;mutations=34A;quality=I)IEIIII-7IA>IIIIII07,-%I>)&#029.2-.
		
		typename std::ifstream::pos_type lineStart = (*file).tellg();
		lineStart = lineStart - (std::ifstream::pos_type) 1;
		
		TId contigId;
		
		// clear temporary variables
		clear(temp_str);
		clear(temp_read);
		
		unsigned pos= 0;
		unsigned pos2 = 0;
		int clipLeft = 0;
		int clipRight = 0;
		unsigned readCount = 0;
		bool qualityFound = false;
		bool readFound = false;
		char orientation = 'F';
		bool hasIndel = false;
		int editDist = 0;
		int mScore = 100;
		_parseSkipWhitespace(*file, c);
		
		// skip whitespaces just in case (actually there shouldnt be a whitespace at the beginning of a line)
		// and read entry in column 1  --> genomeID
		_parseReadWordUntilWhitespace(*file,temp_str,c); 
		
		//check if the genomeID is in our map of relevant genomeIDs, otherwise skip match
		it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) contigId = it->second;
		else
		{
			_parseSkipLine(*file,c);
			continue;
		}
		if((int)contigId < (int)highestChrId)
		{
			std::cerr << "Read files need to be sorted according to chromosomes in genome file.\n";
			return CALLSNPS_GFF_FAILED;
		}
		
		highestChrId = contigId;
		if(contigId < currSeqNo)	// havent reached the sequence of interest yet
		{
			_parseSkipLine(*file,c);
			continue;
		}
		
		if(contigId > currSeqNo)	// have passed the seq of interest
		{
			(*file).seekp(lineStart);
			break;
		}
		contigId = 0; // we only store one chromosome at a time
		
		
		// skip whitespaces and read entry in column 2
		_parseSkipWhitespace(*file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(*file,temp_str,c); 
		
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		
		// skip whitespaces and read entry in column 3
		_parseSkipWhitespace(*file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(*file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		
		// skip whitespaces and read entry in column 4  --> genomic begin position
		_parseSkipWhitespace(*file, c);
		TContigPos beginPos = _parseReadNumber(*file,c) - options.positionFormat;
		if(beginPos > currentEnd + (TContigPos)options.windowBuff)	// we have passed the relevant match positions
		{
			if(options._debugLevel > 1)
				std::cout  << "gBegin "<< beginPos<<"  of match is too large, seeking "<<lineStart<<"\n";
			(*file).seekp(lineStart);
			break;
		}
		if(options._debugLevel > 1) 
			::std::cout << beginPos << "\t";
		
		// skip whitespaces and read entry in column 5  --> genomic end position
		_parseSkipWhitespace(*file, c);
		TContigPos endPos = _parseReadNumber(*file,c);
		
		if(options._debugLevel > 1) 
			::std::cout << endPos << "\t";
		if(endPos + (TContigPos)options.windowBuff < currentBegin)	//we havent reached a relevant read yet
		{
			_parseSkipLine(*file,c);
			continue;
		}
		
		int gMatchLen = endPos - beginPos;	
		int rLen = gMatchLen;
		if(endPos > genomeLen)
		{
			(*file).seekp(lineStart);
			break;	
		}
		
		// skip whitespaces and read entry in column 6  --> score (percent identity or mapping quality) or a '.'
		_parseSkipWhitespace(*file, c);
		if(c=='.')
		{
			mScore = 100;  //not used, but needs to be >= options.minMapQual (default 1)
			c = _streamGet(*file);
		}
		else mScore = (int)_parseReadDouble(*file,c); 
		if(options._debugLevel > 1) 
			::std::cout << mScore << "\t";
		
		// skip whitespaces and read entry in column 7  --> strand information: '+' or '-' 
		_parseSkipWhitespace(*file, c);
		if (c=='+')
			orientation = 'F';
		else
			orientation = 'R';
		c = _streamGet(*file);
		
		// skip whitespaces and read entry in column 8  --> in razers output this is always a '.' 
		_parseSkipWhitespace(*file, c);
		c = _streamGet(*file);
		
		// skip whitespaces and read entry in column 9  --> tags, extra information. in razers output first tag is always "ID"
		_parseSkipWhitespace(*file, c);
		clear(temp_str);
		_parseReadIdentifier(*file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if(temp_str!="ID") ::std::cout << "first feature field should be 'ID'"<<::std::endl;
		
		// skip the "="
		c = _streamGet(*file);
		
		// read the readID
		clear(temp_str);
		CharString readName;
		_parseReadIdentifier(*file,readName,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<<readName << "\n";
#ifdef SNPSTORE_DEBUG
		bool extraV = false;
#endif
		// cut out the read template from the genomic coordinates
		gInf = infix(genome, beginPos, endPos);
		if (orientation == 'R')
			reverseComplement(gInf);
		readTemplate = gInf;
		if(options._debugLevel > 1) cout << readTemplate << "\n";
		
		// process tags in a loop
		CharString current_tag;
		_parseSkipWhitespace(*file,c); 
		bool multi = false, suboptimal = false, unique = true, splitRead = false;
		clipLeft = 0; clipRight = 0; 
		clear(gAliPos);
		clear(tmpCigarStr);
		while(!_streamEOF(*file) && !(c == '\n' || (c == '\r' && _streamPeek(*file) != '\n'))) // while in same line
		{
			// different tags are separated by ';'  ATTENTION: ascii qualities can contain ';', therefore the tag "quality"
			while(c != ';') 						// MUST be the last tag in a line!!!!!!!!
			{
				if(c == '\n' || (c == '\r' && _streamPeek(*file) != '\n')) // end of line
					break;
				c = _streamGet(*file);
			}
			if(c == '\n' || (c == '\r' && _streamPeek(*file) != '\n')) // end of line
				break;
			// get the current tag
			clear(current_tag);
			c = _streamGet(*file);
			_parseReadIdentifier(*file,current_tag,c);
#ifdef SNPSTORE_DEBUG
			if(options._debugLevel > 1) 
				::std::cout << current_tag << " in features\n";
#endif
			if(!options.qualityFile && current_tag=="quality")
			{
				// add the quality to the read
				qualityFound = true;
				if(!readFound)	//read fragmentStore.alignedReadStore with 0 errors --> read == genomeInfix
				{
					temp_read = infix(readTemplate,0,gMatchLen);  // without quality values
					curr_read = temp_read;				// initialized with q40
					readFound = true;
				}
				for(int i = 0; i < rLen ; ++i) //vorsicht! rLen muss hier schon bekannt sein! (i.e. quality tag is last tag!)
				{
					c = _streamGet(*file);
					if (!(c == '\n' || (c == '\r' && _streamPeek(*file) != '\n')))
					{
						assignQualityValue(curr_read[i],(int)ordValue(c)-33);
#ifdef SNPSTORE_DEBUG
						if(extraV) ::std::cout << (char)(getQualityValue(curr_read[i])+33);
#endif
					}
					else
					{
						// shouldnt happen
						if(i != rLen-1 && options._debugLevel > 1) std::cout << curr_read << "macht probleme\n";
						break;
					}
				}
			}
			else
			{
				// parse other tags
				pos = 0;
				pos2 = 0;
				if (current_tag == "unique") {unique = true;}
				else if (current_tag == "multi") {multi = true;}
				else if (current_tag == "suboptimal") {suboptimal = true;}
				else if (current_tag == "split") {splitRead = true;}
				else if (current_tag == "clip")
				{
					options.clipTagsInFile = true;
					if(c=='=') c = _streamGet(*file);
					clipLeft = _parseReadNumber(*file,c);
					c = _streamGet(*file);
					clipRight = _parseReadNumber(*file,c);
				}
				else if (current_tag == "count")
				{
					if(c=='=') c = _streamGet(*file);
					readCount = _parseReadNumber(*file,c);
				}
				else if (current_tag == "read")
				{
					clear(curr_read);
					readFound = true;
					if(c=='=') c = _streamGet(*file);
					while(c != ';' && !(c == '\n' || (c == '\r' && _streamPeek(*file) != '\n')))
					{
						appendValue(curr_read,c);
						c = _streamGet(*file);
					}
					if(mScore != 100)
					{
						editDist = (int)floor((length(curr_read) * ((100.0 - mScore + 0.001)/100.0)));
					}
				}
				else if (current_tag == "mutations")
				{
					if(!readFound)
					{
						temp_read = infix(readTemplate,0,gMatchLen);
						curr_read = temp_read;
						readFound = true;
					}
					while(c==',' || c=='=') // and add the mutated positions (misfragmentStore.alignedReadStore and insertions in read)
					{
						c = _streamGet(*file);
						pos = _parseReadNumber(*file,c);
						curr_read[pos-1]=(Dna5)c;
						c = _streamGet(*file);
						++editDist;
					}
				}
				else if (current_tag == "cigar")
				{
					pos = 0; pos2 = 0;
					int gPos = 0;
					readFound = true;
					while(c!=';'  && !(c == '\n' || (c == '\r' && _streamPeek(*file) != '\n')))
					{
						if(c=='=') c = _streamGet(*file);
						pos2 = _parseReadNumber(*file,c);
						if(c=='M')
						{
							unsigned k= 0;
							while(k<pos2)
							{
								appendValue(gAliPos,gPos,Generous());
								++gPos;
								++k;
							}
							appendValue(tmpCigarStr,TCigar('M',pos2));
							pos2 += pos;
							append(temp_read,infix(readTemplate,pos,pos2));
							pos = pos2;
							c = _streamGet(*file);
							continue;
						}
						if(c=='I')
						{ //insertion in the read
							//(*mIt).editDist += pos2; will be increased in mutations loop
							unsigned k= 0;
							while(k<pos2)
							{
								appendValue(gAliPos,-gPos,Generous());//no genome positions are used up
								++k;
							}
							appendValue(tmpCigarStr,TCigar('I',pos2));
							for(unsigned f = 0; f < pos2; ++f)
								append(temp_read,'A');  // will be replaced with correct base in "mutations" loop
							c = _streamGet(*file);
							rLen += pos2;
							options.hammingOnly = false;
							hasIndel = true;
							continue;
						}
						if(c=='D')
						{ //there is a deletion in the read
							unsigned k= 0;
							while(k<pos2)
							{
								++gPos;
								++k;
							}
							editDist += pos2;
							appendValue(tmpCigarStr,TCigar('D',pos2));
							pos += pos2;
							rLen -= pos2;
							options.hammingOnly = false;
							c = _streamGet(*file);
							hasIndel = true;
							continue;
						}
					}
					curr_read = temp_read;
				}
			}
			if(qualityFound ) {_parseSkipLine(*file,c); break;}
			
		}
		if (options._debugLevel>0&&(rSeq%1000000)==0)cout <<rSeq<<".."<<std::flush;
		if(mScore >= options.minMapQual && (!multi || options.keepMultiReads) && (!suboptimal || options.keepSuboptimalReads))// && (!((*mIt).hasIndel==1 && options.hammingOnly)))
		{
			if(!readFound)
			{	//neither quality nor read sequence found
				if(options._debugLevel>1)::std::cout << "neither quality nor read sequence found editDist = " << editDist <<"\n";			
				temp_read = infix(readTemplate,0,gMatchLen);
				curr_read = temp_read;
			}
			if(clipLeft + clipRight > (int)length(curr_read) - (int)options.minClippedLength)
			{
				if (options._debugLevel>1)cout <<"Discarding read "<<readName<<", too short after clipping.."<<std::endl;
				_parseSkipWhitespace(*file, c); continue;
			}
			if(options.realign && splitRead)
			{
				if(endPos-beginPos > (int)((float)length(curr_read)*1.5))
				{
					
					if (options._debugLevel>1)cout <<"Discarding split read "<<readName<<", deletion too large.."<<std::endl;
					_parseSkipWhitespace(*file, c); continue;
				}
			}
			TId readId = length(fragmentStore.readSeqStore);
			appendValue(fragmentStore.readSeqStore,curr_read,Generous());
			
#ifdef SNPSTORE_DEBUG
			if(clipLeft + clipRight > 76 )
				::std::cerr << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
#endif
			
			if(options._debugLevel > 1) 
				::std::cout<<fragmentStore.readSeqStore[rSeq]<<" with edit="<<editDist<<" at position "<< beginPos <<"\n";
			
			if(endPos - beginPos > (TContigPos)options.maxHitLength)
				options.maxHitLength = endPos - beginPos;
			
			// remember min and max positions seen
			if(beginPos < (TContigPos)options.minCoord) options.minCoord = (unsigned)beginPos;
			if(endPos > (TContigPos)options.maxCoord) options.maxCoord =  (unsigned)endPos;
			
			// create match m
			TMatch m;
			m.id = readId; //length(fragmentStore.alignedReadStore);
			if(orientation == 'F')
			{
				m.beginPos = beginPos;
				m.endPos = endPos;
			}
			else
			{
				m.beginPos = endPos;
				m.endPos = beginPos;
			}
			m.contigId = contigId;
			m.readId = m.id;
			//			std::cout << "readId=" << readId << ::std::endl;
			//			if(readId == 100) 
			//				std::cout << "wahnsinnnnn!!!!" << readId << ::std::endl;
			// corresponding match quality attributes are stored in q
			TMatchQuality q;
			q.errors = (char)editDist;
			q.score = (char) 0;
			if(options._debugLevel > 1)
			{
				if(splitRead && hasIndel)
					std::cout << "has indel!\n"; //TODO: how should split reads be treated in realignment?
			}
			if(!options.realign && splitRead && length(tmpCigarStr)<=3) hasIndel = false;
			if(splitRead) clipLeft = clipRight = 0;
			if(hasIndel)
				q.pairScore = 1;
			else
				q.pairScore = 0;
			
			typename Value<TReadStore>::Type r;
			r.matePairId = TReadStoreElement::INVALID_ID;
			if(readCount > 0) appendValue(readCounts, readCount, Generous());
			
			appendValue(fragmentStore.readStore, r, Generous());
			appendValue(fragmentStore.alignedReadStore, m, Generous());
			appendValue(fragmentStore.readNameStore, readName, Generous());
			appendValue(fragmentStore.alignQualityStore, q, Generous());
			appendValue(readClips,Pair<int,int>(clipLeft,clipRight));
			
			if(!splitRead)
				clear(tmpCigarStr); // split reads store their cigar string explicitly
			
			appendValue(readCigars,tmpCigarStr);
			++rSeq;
			if(options._debugLevel > 1)
			{
				::std::cout<<"Parsed: id= " <<m.readId<<" name="<<readName<<"="<<curr_read<<" with edit="<<editDist<<" at position "<< beginPos<<"\n";
				::std::cout << "mScore=" << mScore << " m.beginPos=" << m.beginPos << "m.endPos="<<m.endPos<<std::endl;
				if(q.pairScore==1) ::std::cout << "indel! pairScore=" << q.pairScore <<std::endl;
				if(q.pairScore==0) ::std::cout << "no indel! pairScore=" << q.pairScore <<std::endl;
				
			}
		}
		else 
		{
			if(options._debugLevel > 1 )
			{
				::std::cout<<"Discarded: "<<curr_read<<" with edit="<<editDist<<" at position "<< beginPos<<"\n";
				::std::cout << "mScore = " << mScore << std::endl;
			}
		}
		
		_parseSkipWhitespace(*file, c);
	}
	if(options._debugLevel > 0) 
		::std::cout << ::std::endl << "Parsed "<<length(fragmentStore.alignedReadStore)<<" matches of "<<length(fragmentStore.readSeqStore)<<" reads." << ::std::endl;
	
	
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int detectSNPs(
			   const char *genomeFileName,
			   String<CharString> & readFNames,
			   String<CharString> & qualityFNames,
			   SNPCallingOptions<TSpec> &options)
{
	
	typedef FragmentStore<SnpStoreSpec_>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeq		TReadSeq;				// TRead
	typedef typename TFragmentStore::TContigSeq		TContigSeq;				// TGenome
	//typedef typename Position<TReadSeq>::Type		TReadPos;				// TPos
	typedef typename TFragmentStore::TReadPos		TReadPos;				// TPos
	//typedef typename Position<TContigSeq>::Type		TContigPos;				// TContigPos
	typedef typename TFragmentStore::TContigPos 		TContigPos;
	typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;	 	// TMatches
	typedef typename TFragmentStore::TAlignQualityStore	TAlignQualityStore;	 	// TMatchQualities
	typedef typename TFragmentStore::TReadStore		TReadStore;				// TReadSet
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;				// TReadSet
	typedef typename TFragmentStore::TContigStore		TContigStore;			// TGenomeSet
	typedef typename Value<TContigStore>::Type		TContig;
	typedef              TContigSeq                    TGenome;
	typedef StringSet<TGenome>                         TGenomeSet;
	
	typedef String<SimplePosition<TContigPos > >		TPositions;
	typedef ::std::map<CharString,unsigned> 		TGenomeMap;
	typedef typename TGenomeMap::iterator 			TMapIter;
	typedef String<unsigned>				TReadCounts;
	typedef String<Pair<int,int> >				TReadClips;
	typedef StringSet<String<Pair<char,int> > >		TReadCigars;
	
	TGenomeSet				genomes;
	StringSet<CharString> 			genomeFileNameList; // filenamen
	StringSet<CharString> 			genomeNames;				// todo: raus
	TGenomeMap				gIdStringToIdNumMap;	// name to id
	
	
	// dump configuration in verbose mode
	if (options._debugLevel >= 1) 
	{
		::std::cerr << "___SETTINGS____________" << ::std::endl;
		::std::cerr << "Genome file:                             \t" << genomeFileName << ::std::endl;
		::std::cerr << "Read files:                              \t" << readFNames[0] << ::std::endl;
		for(unsigned i = 1; i < length(readFNames); ++i)
			::std::cerr << "                                         \t" << readFNames[i] << ::std::endl;
		if(options.inputFormat == 1) 
		{
			::std::cerr << "Quality files:                           \t" << qualityFNames[0] << ::std::endl;
			for(unsigned i = 1; i < length(qualityFNames); ++i)
				::std::cerr << "                                         \t" << qualityFNames[i] << ::std::endl;
		}::std::cerr << "MaxPile:                                 \t" << options.maxPile << ::std::endl;
		if(options.laneSpecificMaxPile)::std::cerr << "Lane specific:                           \tYES" << ::std::endl;
		else ::std::cerr << "Lane specific:                           \tNO" << ::std::endl;
		::std::cerr << "MinCoverage:                             \t" << options.minCoverage << ::std::endl;
		if(options.method == 0)
		{
			::std::cerr << "MinMutThreshold:                         \t" << options.minMutT << ::std::endl;
			::std::cerr << "MinPercentageThreshold:                  \t" << options.percentageT << ::std::endl;
			::std::cerr << "MinQualityThreshold:                     \t" << options.avgQualT << ::std::endl;
		}
		else
		{
			::std::cerr << "MinMappingQuality:                       \t" << options.minMapQual << ::std::endl;
		}
		if(options.doIndelCalling && *options.outputIndel != 0)
		{
			::std::cerr << "IndelCountThreshold:                     \t" << options.indelCountThreshold << ::std::endl;
			::std::cerr << "IndelPercentageThreshold:                \t" << options.indelPercentageT << ::std::endl;
			::std::cerr << "IndelWindow:                             \t" << options.indelWindow << ::std::endl;
		}
		::std::cerr << ::std::endl;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Determine genome file type and load genomes
	SEQAN_PROTIMESTART(load_time);
	
	int result = getGenomeFileNameList(genomeFileName, genomeFileNameList, options);
	if(result == CALLSNPS_GENOME_FAILED || !loadGenomes(genomes, genomeFileNameList,gIdStringToIdNumMap))
	{
		::std::cerr << "Failed to open genome file " << genomeFileName << ::std::endl;
		return result;
	}
	
	
	resize(genomeNames,gIdStringToIdNumMap.size()); //prepare genomeNames
	TMapIter gIt=gIdStringToIdNumMap.begin();
	for(unsigned i=0; i < length(genomeNames); ++i,++gIt)
		genomeNames[i] = gIt->first;
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Load fragmentStore.readSeqStore and fragmentStore.alignedReadStore
	// open read files and  store open file pointers
	String<int> highestChrId;
	resize(highestChrId,length(readFNames),0);
	vector< ::std::fstream* > readFileStreams;
	readFileStreams.resize(length(readFNames)); 
	for(unsigned i = 0; i < length(readFNames); ++i)
	{
		readFileStreams[i] = new fstream(toCString(readFNames[i]), ios_base::in | ios::binary);
		if(!(*(readFileStreams[i])).is_open())
		{
			::std::cerr << "Failed to open read file " << readFNames[i] << ::std::endl;
			return CALLSNPS_GFF_FAILED;
		}
	}
	
	/////////////////////////////////////////////////////////////////////
	// open out file streams and store open file pointers
	::std::ofstream snpFileStream; 
	if (*options.outputSNP != 0)
	{
		
		// prepare lookup tables for maq statistics
		if (options.method == 1 )
		{
			computeCnks(options.cnks,options.fks,options);
			options.priorHetQ = computeHetTable(options.hetTable,options);
		}
		
		snpFileStream.open(options.outputSNP,::std::ios_base::out);
		if(!snpFileStream.is_open())
			return CALLSNPS_OUT_FAILED;
		snpFileStream << "#" << (options.programCall).str() << "\n";
		if(options.outputFormat < 2)
		{
			if(options.orientationAware)
				snpFileStream << "#chr\tpos\tref\t[A+]\t[C+]\t[G+]\t[T+]\t[A-]\t[C-]\t[G-]\t[T-]\tcov\tcall";
			else
				snpFileStream << "#chr\tpos\tref\tA\tC\tG\tT\tcov\tcall";
			//if(options.method==1)
			snpFileStream << "\tquality\n";
			//else
			//	file <<"\n";
		}
	}
	::std::ofstream indelFileStream; 
	if (*options.outputIndel != 0)
	{
		indelFileStream.open(options.outputIndel,::std::ios_base::out);
		if(!indelFileStream.is_open())
			return CALLSNPS_OUT_FAILED;
	}
	//	::std::ofstream cnvFileStream; 
	//	if (*options.outputCNV != 0)
	//	{
	//		cnvFileStream.open(options.outputCNV,::std::ios_base::out);
	//		if(!cnvFileStream.is_open())
	//			return CALLSNPS_OUT_FAILED;
	//	}
	
	/////////////////////////////////////////////////////////////////////////////
	// helper variables
	Pair<int,int> zeroPair(0,0);
	int sumreads = 0;
	int sumwindows = 0;
	
	
	/////////////////////////////////////////////////////////////////////////////
	// Start scanning for SNPs/indels
	// for each chromosome
	for(unsigned i=0; i < length(genomeNames); ++i)
	{
		// parse matches batch by batch
		TContigPos currentWindowBegin = 0;
		if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << i << " ..." << ::std::endl;
		
		// containers for those matches that overlap two windows	
		TAlignedReadStore tmpMatches;
		TAlignQualityStore tmpQualities;
		TReadStore tmpRs;
		TReadSeqStore tmpReads;
		TReadCounts tmpReadCounts;
		TReadClips tmpReadClips;
		TReadCigars tmpReadCigars;
		options.minCoord = MaxValue<unsigned>::VALUE;
		options.maxCoord = 0;
		
		// snp calling is done for all positions between windowBegin and windowEnd
		while(currentWindowBegin < (TContigPos)length(genomes[i]))
		{
			TContigPos currentWindowEnd = currentWindowBegin + options.windowSize;
			if(currentWindowEnd > (TContigPos)length(genomes[i])) currentWindowEnd = (TContigPos)length(genomes[i]);
			
			if(options._debugLevel > 0)
				::std::cout << "Sequence number " << i << " window " << currentWindowBegin << ".." << currentWindowEnd << "\n";
			
			TFragmentStore fragmentStore;  
			TReadCounts readCounts;
			TReadClips readClips;
			TReadCigars readCigars;	// currently only stored for split-mapped reads
			
			// add the matches that were overlapping with this and the last window (copied in order to avoid 2 x makeGlobal)
			if(!empty(tmpMatches))
			{
				sumreads -=  length(tmpReads);	// count these reads only once
				resize(fragmentStore.alignQualityStore,length(tmpMatches));
				resize(fragmentStore.alignedReadStore,length(tmpMatches));
				resize(fragmentStore.readSeqStore,length(tmpMatches));
				resize(fragmentStore.readStore,length(tmpMatches));
				if(!empty(tmpReadClips))resize(readClips,length(tmpMatches));
				if(!empty(tmpReadCounts)) resize(readCounts,length(tmpMatches));
				if(!empty(tmpReadCigars))resize(readCigars,length(tmpMatches));
				
				arrayMoveForward(begin(tmpQualities,Standard()), end(tmpQualities,Standard()), begin(fragmentStore.alignQualityStore,Standard()));
				arrayMoveForward(begin(tmpMatches,Standard()), end(tmpMatches,Standard()), begin(fragmentStore.alignedReadStore,Standard()));
				arrayMoveForward(begin(tmpReads,Standard()), end(tmpReads,Standard()), begin(fragmentStore.readSeqStore,Standard()));
				arrayMoveForward(begin(tmpRs,Standard()), end(tmpRs,Standard()), begin(fragmentStore.readStore,Standard()));
				if(!empty(tmpReadCounts)) arrayMoveForward(begin(tmpReadCounts,Standard()), end(tmpReadCounts,Standard()), begin(readCounts,Standard()));
				if(!empty(tmpReadCigars)) arrayMoveForward(begin(tmpReadCigars,Standard()), end(tmpReadCigars,Standard()), begin(readCigars,Standard()));
				if(!empty(tmpReadClips)) arrayMoveForward(begin(tmpReadClips,Standard()), end(tmpReadClips,Standard()), begin(readClips,Standard()));
				
			}	
			
			// parse matches for current window
			if(options._debugLevel > 0)
				::std::cout << "Parsing reads up to position " << currentWindowEnd << "...\n";
			for(unsigned j = 0; j < length(readFNames); ++j)
			{
				unsigned sizeBefore = length(fragmentStore.alignedReadStore);
				
				// currently only gff supported
				if(options.inputFormat ==0)
					result = readMatchesFromGFF_Batch(readFileStreams[j], fragmentStore, readCounts, readClips,
													  readCigars, genomes[i], gIdStringToIdNumMap, 
													  i, currentWindowBegin, currentWindowEnd, highestChrId[j], options);
				if(result == CALLSNPS_GFF_FAILED)
				{
					::std::cerr << "Failed to open read file " << readFNames[j] << ::std::endl;
					::std::cerr << "or reads are not sorted correctly. " << ::std::endl;
					return result;
				}
				if(result > 0)
					return result;
				
				if(options._debugLevel > 0)
					::std::cout << "parsed reads of file " << j << "\n";
				
				// store average quality of each read
				addReadQualityToMatches(fragmentStore,sizeBefore,(unsigned)length(fragmentStore.alignedReadStore),options);
				
				// do pile up correction if lane-based
				if(options.maxPile != 0 && options.laneSpecificMaxPile) 
					applyPileupCorrection(fragmentStore,sizeBefore,(unsigned)length(fragmentStore.alignedReadStore),options);
				
			}
			if (options._debugLevel > 1)  // number of reads currently in memory
				::std::cerr << lengthSum(fragmentStore.readSeqStore) << " bps of " << length(fragmentStore.readSeqStore) << " reads in memory." << ::std::endl;
			sumreads +=  length(fragmentStore.readSeqStore);  // total count of reads
			
			// do merged pile up correction
			if(options.maxPile != 0 && !options.laneSpecificMaxPile)
				applyPileupCorrection(fragmentStore,(unsigned)0,(unsigned)length(fragmentStore.alignedReadStore),options);
			
			// these were set while parsing matches, first and last position of parsed matches
			//			TContigPos startCoord = options.minCoord;// can be < currentWindowBegin
			//			TContigPos endCoord = options.maxCoord; // can be > currentWindoEnd
			
			// check
			TContigPos startCoord = _max((int)options.minCoord-options.realignAddBorder,0);// can be < currentWindowBegin
			TContigPos endCoord = _min(options.maxCoord+options.realignAddBorder,length(genomes[i])); // can be > currentWindoEnd
			
			
			if(!empty(fragmentStore.alignedReadStore))
			{
				//initial values of min and max coords for next round are set here
				if(currentWindowEnd != (TContigPos)length(genomes[i]))
				{
					clear(tmpMatches);
					clear(tmpQualities);
					clear(tmpRs);
					clear(tmpReads);
					clear(tmpReadCounts);
					clear(tmpReadClips);
					clear(tmpReadCigars);
					copyNextWindowMatchesAndReads(fragmentStore,readCounts,readCigars,tmpReadCounts,tmpReads,tmpRs,tmpMatches,tmpQualities,tmpReadClips,tmpReadCigars,i,currentWindowEnd,options);
				}	
				
				//	::std::cout << "Min = " << options.minCoord << " Max = " << options.maxCoord << std::endl;
				//	::std::cout << "Min = " << startCoord << " Max = " << endCoord << std::endl;
				
				// coordinates are relative to current chromosomal window (segment)
				transformCoordinates(fragmentStore,startCoord,options);
				
				// set the current chromosomal segment as contig sequence 
				TContig conti;
				conti.seq = infix(genomes[i],startCoord,endCoord);
				appendValue(fragmentStore.contigStore, conti, Generous() );
				appendValue(fragmentStore.contigNameStore, genomeNames[i], Generous() );// internal id is always 0
				
				// clip Reads if clipping is switched on and there were clip tags in the gff file
				if(!options.dontClip && options.clipTagsInFile)
				{
					options.useBaseQuality = false;	// activate "average read quality"-mode for snp calling, low quality bases should be clipped anyway
					clipReads(fragmentStore,readClips,(unsigned)0,(unsigned)length(fragmentStore.alignedReadStore),options);
				}
				
				// check for indels
				if (*options.outputIndel != 0)
				{
					if(options._debugLevel > 1) ::std::cout << "Check for indels..." << std::endl;
					if(!options.realign) dumpShortIndelPolymorphismsBatch(fragmentStore, readCigars, fragmentStore.contigStore[0].seq, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, indelFileStream, options);
				}
				
				// // check for CNVs
				//				if (*options.outputCNV != 0)
				//					dumpCopyNumberPolymorphismsBatch(fragmentStore, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, cnvFileStream, options);
				
#ifdef SNPSTORE_DEBUG
				CharString strstr = "test";
				//				_dumpMatches(fragmentStore, strstr );
#endif
				if (*options.outputSNP != 0)
				{
					if(options._debugLevel > 1) ::std::cout << "Check for SNPs..." << std::endl;
					if(options.realign)
						dumpVariantsRealignBatchWrap(fragmentStore, readCigars, readCounts, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, snpFileStream,indelFileStream,options);
					else 
						dumpSNPsBatch(fragmentStore, readCigars, readCounts, genomeNames[i], startCoord, currentWindowBegin, currentWindowEnd, snpFileStream,options);
				}
			}
			currentWindowBegin = currentWindowEnd;
			++sumwindows;
		}
		
	}
	if (*options.outputSNP != 0)
		snpFileStream.close();
	
	if (*options.outputIndel != 0)
		indelFileStream.close();
	
	//	if (*options.outputCNV != 0)
	//		cnvFileStream.close();
	
	return 0;
}




// log file to keep track of happenings
template <typename TSpec>
int writeLogFile(
				 int argc, const char *argv[],
				 const char *genomeFileName,
				 String<CharString> & readFNames,
				 String<CharString> & ,
				 SNPCallingOptions<TSpec> &options)
{
	
	::std::ofstream logfile;
	logfile.open(options.outputLog, ::std::ios_base::out | ::std::ios_base::trunc);
	if (!logfile.is_open()) 
	{
		::std::cerr << "Failed to open log file" << ::std::endl;
		return 1;
	}
	logfile << "#call" << std::endl;
	for (int i = 0; i < argc; ++i)
		logfile << argv[i] << " ";
	logfile << std::endl;
	
	logfile << "#files" << ::std::endl;
	logfile << "Genome=\"" << genomeFileName << "\""<<::std::endl;
	logfile << "Reads=\"" << readFNames[0];
	for(unsigned i = 1; i < length(readFNames); ++i)
		logfile << " " << readFNames[i] << ::std::endl;
	logfile << "\"" << std::endl;
	if(*options.output != 0)
	{
		logfile << "OutputSnp=\"" << CharString(options.output) << "\"" << ::std::endl;
	}
	if(*options.outputIndel != 0)
	{
		logfile << "OutputIndel=\"" << CharString(options.outputIndel) << "\""<< ::std::endl;
	}
	logfile << "#settings" << ::std::endl;
	logfile << "MaxPile=" << options.maxPile << ::std::endl;
	logfile << "MinCov=" << options.minCoverage << ::std::endl;
	logfile << "Method=" << options.method << ::std::endl;
	if(options.method == 0)
	{
		logfile << "MinMutT=" << options.minMutT << ::std::endl;
		logfile << "MinPercT=" << options.percentageT << ::std::endl;
		logfile << "MinQualT=" << options.avgQualT << ::std::endl;
	}
	else
	{
		logfile << "MinMapQual=" << options.minMapQual << ::std::endl;
	}
	if(*options.outputIndel != 0)
	{
		logfile << "MinIndel=" << options.indelCountThreshold << ::std::endl;
		logfile << "MinPercIndelT=" << options.indelPercentageT << ::std::endl;
	}
	logfile.close();
	return 0;
}



template <typename TOptions>
void printHelp(int, const char *[], TOptions &options, bool longHelp = false)
{
	
	cerr << "Usage: snpStore [OPTION]... <GENOME FILE> <MAPPED READ FILE>" << endl;
	if (longHelp) {
		cerr << endl << "Options:" << endl;
		cerr << "  -o,  --output FILE               \t" << "change output filename (default <READ FILE>.snp)" << endl;
		cerr << "  -of, --output-format NUM         \t" << "output format:" << endl;
		cerr << "                                   \t" << "0 = output all candidate snps (default)" << endl;
		cerr << "                                   \t" << "1 = output succesful candidate snps only" << endl;
		cerr << "  -dc, --dont-clip                 \t" << "ignore clip tags in gff (off)" << endl;
		cerr << "  -mu, --multi                     \t" << "keep non-unique fragmentStore.alignedReadStore (off)" << endl;
		cerr << "  -hq, --hide-qualities            \t" << "only show coverage (no qualities) in output file (off)" << endl;
		cerr << "  -id, --indel-file FILE           \t" << "output file for called indels in gff format (off)" << endl;
		cerr << "  -m,  --method NUM                \t" << "set method used for SNP calling" << endl;
		cerr << "                                   \t" << "0 = threshold method" << endl;
		cerr << "                                   \t" << "1 = maq (default)" << endl;
		//		cerr << "                                   \t" << "(default = "<<options.method << ")" << endl;
		cerr << "  -mp, --max-pile NUM              \t" << "maximal number of matches allowed to pile up at the same genome position ("<<options.maxPile<<")" << endl;
		cerr << "  -mmp,--merged-max-pile           \t" << "do pile up correction on merged lanes (off)" << endl;
		cerr << "  -mc, --min-coverage NUM          \t" << "minimal number of fragmentStore.readSeqStore covering a candidate position ("<< options.minCoverage<<")" << endl;
		cerr << "  -fc, --force-call NUM            \t" << "always call base if count is >= fc, ignore other parameters (off)" << endl;
		cerr << "  -oa, --orientation-aware         \t" << "distinguish between forward and reverse fragmentStore.alignedReadStore (off)" << endl;
		cerr << "  -fl, --force-length NUM          \t" << "read length to be used (igonres suffix of read) (off)" << endl;
		cerr << endl;
		cerr << "SNP calling options: " << endl;
		cerr << "  Threshold method related: " << endl;
		cerr << "  -mm, --min-mutations NUM         \t" << "minimal number of observed mutations for mutation to be called ("<<options.minMutT<<")" << endl;
		cerr << "  -pt, --perc-threshold NUM        \t" << "minimal percentage of mutational base for mutation to be called (" << options.percentageT << ")" << endl;
		cerr << "  -mq, --min-quality NUM           \t" << "minimal average quality of mutational base for mutation to be called ("<<options.avgQualT <<")" << endl;
		cerr << "  Maq method related: " << endl;
		cerr << "  -th, --theta                     \t" << "dependency coefficient ("<< options.theta <<")" << endl;
		cerr << "  -hr, --hetero-rate               \t" << "heterozygote rate ("<< options.hetRate <<")" <<  endl;
		cerr << "  -mmq,--min-map-quality           \t" << "minimum mapping quality for a match to be considered ("<< options.minMapQual <<")" <<  endl;
		cerr << "  Indel calling related: " << endl;
		cerr << "  -it, --indel-threshold           \t" << "absolute count threshold for indel calling (" << options.indelCountThreshold<<")"<< endl;
		cerr << "  -ipt,--indel-perc-threshold NUM  \t" << "minimal ratio of indel-count/coverage for indel to be called (" << options.indelPercentageT<<")" << endl;
		cerr << "  -iw, --indel-window              \t" << "overlap window used for indel calling (" << options.indelWindow<<")"<< endl;
		
		cerr << endl<< "Other options: " << endl;
		cerr << "  -lf, --log-file FILE             \t" << "write log file to FILE" << endl;
		cerr << "  -v,  --verbose                   \t" << "verbose mode" << endl;
		cerr << "  -vv, --very-verbose              \t" << "very verbose mode" << endl;
		cerr << "  -h,  --help                      \t" << "print this help" << endl << endl;
	}
	else {
		cerr << "Try 'snpStore --help' for more information." << endl <<endl;
	}
}



int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line
	
	SNPCallingOptions<>		options;
	
	unsigned				fnameCount = 0;
	const char				*genomeFName = "";
	String<CharString> 		readFNames;
	String<CharString> 		qualityFNames;
	
	for(int arg = 0; arg < argc; ++arg) {
		options.programCall << argv[arg] << " ";
	}
	
	/*	std::cout << "lgamma(1) = " << lgamma(1) << std::endl;
	 std::cout << "lgamma(2) = " << lgamma(2) << std::endl;
	 std::cout << "lgamma(3) = " << lgamma(3) << std::endl;
	 std::cout << "lgamma(4) = " << lgamma(4) << std::endl;
	 std::cout << "lgamma(25) = " << lgamma(25) << std::endl;
	 std::cout << "lgamma(105) = " << lgamma(105) << std::endl;
	 std::cout << "lgamma(255) = " << lgamma(255) << std::endl;*/
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse options
			if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--method") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.method;
					if (!istr.fail())
					{
						if (options.method > 1)
							cerr << "Invalid method option." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--min-quality") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.avgQualT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mm") == 0 || strcmp(argv[arg], "--min-mutations") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.minMutT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-fl") == 0 || strcmp(argv[arg], "--force-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.forceReadLength;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-fc") == 0 || strcmp(argv[arg], "--force-call") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.forceCallCount;
					if (!istr.fail())
					{
						if(options.forceCallCount < 1)
							cerr << "--force-call expects a positive integer." << endl;	
						else continue;
					}
					
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-ipt") == 0 || strcmp(argv[arg], "--indel-perc-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.indelPercentageT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-pt") == 0 || strcmp(argv[arg], "--perc-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.percentageT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mp") == 0 || strcmp(argv[arg], "--max-pile") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.maxPile;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-dp") == 0 || strcmp(argv[arg], "--diff-pos") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.minDifferentReadPos;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-eb") == 0 || strcmp(argv[arg], "--exclude-border") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.excludeBorderPos;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mc") == 0 || strcmp(argv[arg], "--min-coverage") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.minCoverage;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-dc") == 0 || strcmp(argv[arg], "--dont-clip") == 0) {
				options.dontClip = true;
				continue;
			}
			if (strcmp(argv[arg], "-mu") == 0 || strcmp(argv[arg], "--multi") == 0) {
				options.keepMultiReads = true;
				continue;
			}
			if (strcmp(argv[arg], "-re") == 0 || strcmp(argv[arg], "--realign") == 0) {
				options.realign = true;
				continue;
			}
			if (strcmp(argv[arg], "-hq") == 0 || strcmp(argv[arg], "--hide-qualities") == 0) {
				options.showQualityStrings = false;
				continue;
			}
			if (strcmp(argv[arg], "-mmp") == 0 || strcmp(argv[arg], "--merged-max-pile") == 0) {
				options.laneSpecificMaxPile = false;
				continue;
			}
			if (strcmp(argv[arg], "-oa") == 0 || strcmp(argv[arg], "--orientation-aware") == 0) {
				options.orientationAware = true;
				continue;
			}
			if (strcmp(argv[arg], "-iw") == 0 || strcmp(argv[arg], "--indel-window") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.indelWindow;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-it") == 0 || strcmp(argv[arg], "--indel-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.indelCountThreshold;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mmq") == 0 || strcmp(argv[arg], "--min-map-quality") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.minMapQual;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-th") == 0 || strcmp(argv[arg], "--theta") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.theta;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-pws") == 0 || strcmp(argv[arg], "--parse-window-size") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.windowSize;
					if (!istr.fail() && options.windowSize > 0)
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-cws") == 0 || strcmp(argv[arg], "--cnv-window-size") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.cnvWindowSize;
					if (!istr.fail() && options.cnvWindowSize > 0)
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-hr") == 0 || strcmp(argv[arg], "--hetero-rate") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.hetRate;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			//if (strcmp(argv[arg], "-if") == 0 || strcmp(argv[arg], "--input-format") == 0) {
			//	if (arg + 1 < argc) {
			//		++arg;
			//		istringstream istr(argv[arg]);
			//		istr >> options.inputFormat;
			//		if (!istr.fail())
			//			continue;
			//	}
			//	printHelp(argc, argv, options, true);
			//	return 0;
			//}
			if (strcmp(argv[arg], "-of") == 0 || strcmp(argv[arg], "--output-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.outputFormat;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-lf") == 0 || strcmp(argv[arg], "--log-file") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.outputLog = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.outputSNP = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-id") == 0 || strcmp(argv[arg], "--indel-file") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.outputIndel = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-cnv") == 0 || strcmp(argv[arg], "--output-cnv") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.outputCNV = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--very-verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 2);
				continue;
			}
			cerr << "Unknown option: " << argv[arg] << endl << endl;
			printHelp(argc, argv, options);
			return 0;
		} else {
			// parse file name
			if (fnameCount == 0)
				genomeFName = argv[arg];
			if (fnameCount == 1)
			{
				if(argv[arg][0]=='[')
				{
					String<char> tempStr = argv[arg];
					appendValue(readFNames,suffix(tempStr,1),Generous());
					++arg;
					while(arg < argc && argv[arg][0] != '-')
					{
						tempStr = argv[arg];
						appendValue(readFNames,tempStr,Generous());
						++arg;
					}
					if(readFNames[length(readFNames)-1][length(readFNames[length(readFNames)-1])-1] != ']' )
						cerr << "Something wrong with read file list?" << endl;
					resize(readFNames[length(readFNames)-1],length(readFNames[length(readFNames)-1])-1);
					--arg;
				}
				else
				{
					//split by whitesapce and append each read file
					::std::string tempStr(argv[arg]);
					size_t firstPos = 0;
					size_t lastPos = tempStr.find(' ');
					::std::string tempFile = tempStr.substr(firstPos,lastPos);
					appendValue(readFNames,String<char>(tempFile),Generous());
					while (lastPos != 0 && lastPos != tempStr.npos)
					{
						while(tempStr[lastPos]==' ')
							++lastPos;
						firstPos = lastPos; 
						lastPos = tempStr.find(' ',firstPos);
						if (lastPos != tempStr.npos) tempFile = tempStr.substr(firstPos,lastPos-firstPos);
						else tempFile = tempStr.substr(firstPos,length(tempStr));
						appendValue(readFNames,String<char>(tempFile),Generous());
					}
				}
			}
			if (fnameCount == 2) {
				cerr << "More than 2 input files specified." <<endl;
				cerr << "If more than 2 mapped read files are to be parsed, use quotation marks directly before first file name and directly after last file name (e.g. \"lane1.gff lane2.gff\")." << endl << endl;
				printHelp(argc, argv, options);
				return 0;
			}
			++fnameCount;
		}
	}
	
	// some option checking
	if (fnameCount != 2) {
		if (argc > 1 && !options.printVersion)
			cerr << "Exactly 2 input files need to be specified." << endl << endl;
		printHelp(argc, argv, options);
		return 0;
	}
	if(options.inputFormat == 1 && (!options.qualityFile || (length(qualityFNames)!=length(readFNames))))
	{
		cerr << "If mapped read file is in Eland format, a .qual or .fastq file containing read qualities needs to be specified." << endl << endl;
		return 0;
	}
	if(*(options.positionFile) == 0 && *(options.outputPositionAnalysis) != 0)
	{
		cerr << "Position analysis output specified, but no position file given." << endl << endl;
		return 0;
	}
	
	if(options.realign || options.windowSize > 50000) 
		options.windowSize = 10000;
	
	if (*options.outputLog != 0)
		writeLogFile(argc, argv, genomeFName, readFNames, qualityFNames, options);
	
	// TODO: is forceReadLength still needed?
	if(options.inputFormat == 0 )  options.forceReadLength = 0; // for now this is safer
	
	if(options.runID == "")
	{
		::std::string tempStr = toCString(readFNames[0]);
		size_t lastPos = tempStr.find_last_of('/') + 1;
		if (lastPos == tempStr.npos) lastPos = tempStr.find_last_of('\\') + 1;
		if (lastPos == tempStr.npos) lastPos = 0;
		options.runID = tempStr.substr(lastPos);
	}
	
	
	
	//////////////////////////////////////////////////////////////////////////////
	// check for variants
	int result = detectSNPs(genomeFName, readFNames, qualityFNames, options);
	if (result > 0)
	{
		printHelp(argc, argv, options);
		return 0;
	}
	return result;
}
