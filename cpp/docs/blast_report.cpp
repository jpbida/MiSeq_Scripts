#include <iostream>
#include <seqan/blast.h>

using namespace seqan;


template <typename TFile>
void read_blast_report(TFile & strm)
{

	typedef BlastHsp<BlastN, FullInfo> TBlastHsp;

	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Iterator<TBlastReport,HitIterator>::Type THitIterator;
	typedef typename Iterator<TBlastHit,HspIterator>::Type THspIterator;

	
	TBlastReport blast;

	unsigned hspcount = 0;
	unsigned hitcount = 0;
	unsigned highsignif = 0;

	while(!atEnd(strm,blast)) 
	{
		read(strm,blast,Blast());
		::std::cout << "Query: "<<getQueryName(blast) <<"\n";
		::std::cout << "Database: "<<getDatabaseName(blast) <<"\n\n";

		THitIterator hit_it(blast); 
		for(; !atEnd(strm,hit_it); goNext(strm,hit_it)) 
		{
			++hitcount;
			TBlastHit hit = getValue(strm,hit_it);
			::std::cout << " Hit: " <<name(hit) <<"\n\n";

			/// iterate over alignments (HSPs)
			THspIterator hsp_it(hit);
			for(; !atEnd(strm,hsp_it); goNext(strm,hsp_it)) 
			{
				++hspcount;
 				TBlastHsp hsp = getValue(strm,hsp_it);
				
				::std::cout << "  Score  = " << getBitScore(hsp) << "\n";
				::std::cout << "  Length = " << length(hsp) << "\n\n";
				if(getEValue(hsp)<0.01)
					++highsignif;

			}
		}
	}
	::std::cout <<"Total number of Hits: "<< hitcount<<::std::endl;
	::std::cout <<"Total number of HSPs: "<< hspcount<<::std::endl;
	::std::cout <<"Number of highly significant HSPs: "<< highsignif<<::std::endl;

}


int main()
{
	::std::fstream strm;
	strm.open("ecoln.out", ::std::ios_base::in | ::std::ios_base::binary);
	read_blast_report(strm);
	return 0;
}
