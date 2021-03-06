#define SEQAN_PROFILE // enable time measurements
#include <seqan/misc/misc_cmdparser.h>
#include <tr1/unordered_map>
#include <seqan/file.h>
#include <iostream>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/basic.h>
#include <sstream>
#include <string>
#include <stdio.h>

using namespace seqan;
std::string IntToStr( int n )
{
    std::ostringstream result;
    result << n;
    return result.str();
}

//The known library of sequences is stored as a StringSet called THaystacks
//We will generate an index against this file to make the search faster
typedef StringSet<CharString> THaystacks;

//A hashmap is used to map the experimental ids to the descriptive labels
typedef std::tr1::unordered_map<std::string, std::string> hashmap;

//Versioning information
inline void
_addVersion(CommandLineParser& parser) {
    ::std::string rev = "$Revision: 0000 $";
    addVersionLine(parser, "Version 1.0 (9 August 2012) Revision: " + rev.substr(11, 4) + "");
}

int main(int argc, const char *argv[]) {
    int perfect=0;
//All command line arguments are parsed using SeqAn's command line parser
    CommandLineParser parser;
    _addVersion(parser);

    addTitleLine(parser, "                                                 ");
    addTitleLine(parser, "*************************************************");
    addTitleLine(parser, "* SHAPE-Seq Pairwise alignment - shapeAln       *");
    addTitleLine(parser, "* (c) Copyright 2012 by JP Bida                 *");
    addTitleLine(parser, "*************************************************");
    addTitleLine(parser, "                                                 ");

    addUsageLine(parser, "-c <constant sequence> -s <miseq output1> -e <miseq output2> -l <library> -b <experimetal barcodes> -n <sequence id length> -d <0/1 output type>");

    addSection(parser, "Main Options:");
    addOption(parser, addArgumentText(CommandLineOption("c", "cseq", "Constant sequence", OptionType::String), "<DNA sequence>"));
    addOption(parser, addArgumentText(CommandLineOption("s", "miseq1", "miseq output containing ids",OptionType::String), "<FASTAQ FILE>"));
    addOption(parser, addArgumentText(CommandLineOption("e", "miseq2", "miseq output containing 3' ends",OptionType::String), "<FASTAQ FILE>"));
    addOption(parser, addArgumentText(CommandLineOption("l", "library", "library of sequences to align against", OptionType::String),"<FASTA FILE>"));
    addOption(parser, addArgumentText(CommandLineOption("b", "barcodes", "fasta file containing experimental barcodes", OptionType::String), "<FASTA FILE>"));
    addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename", (int)OptionType::String, "out.fasta"), "<Filename>"));
    addOption(parser, addArgumentText(CommandLineOption("n", "sid.length", "sequence id length", OptionType::Int), "<Int>"));
    addOption(parser, addArgumentText(CommandLineOption("d", "debug", "full output =1 condensed output=0", OptionType::Int), "<Int>"));
    if (argc == 1)
    {
        shortHelp(parser, std::cerr);	// print short help and exit
        return 0;
    }

    if (!parse(parser, argc, argv, ::std::cerr)) return 1;
    if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit

//This isn't required but shows you how long the processing took
    SEQAN_PROTIMESTART(loadTime);
    std::string file1,file2,file3,file4,outfile;

//cseq is the constant region between the experimental id and the sequence id
//in the Das lab this is the tail2 sequence AAAGAAACAACAACAACAAC
    String<char> cseq="";
    getOptionValueLong(parser, "cseq",cseq);
    getOptionValueLong(parser, "miseq1",file1);
    getOptionValueLong(parser, "miseq2",file2);
    getOptionValueLong(parser, "library",file3);
    getOptionValueLong(parser, "barcodes",file4);
    getOptionValueLong(parser, "outfile",outfile);
    int seqid=-1;
    int debug=0;
    getOptionValueLong(parser,"sid.length",seqid);
    getOptionValueLong(parser,"debug",debug);

//Opening the output file and returning an error if it can't be opened
    FILE * oFile;
    oFile = fopen(outfile.c_str(),"w");

    MultiSeqFile multiSeqFile1;
    MultiSeqFile multiSeqFile2;
    MultiSeqFile multiSeqFile3;
    MultiSeqFile multiSeqFile4; //barcode patterns
    if (!open(multiSeqFile1.concat, file1.c_str(), OPEN_RDONLY) ||
            !open(multiSeqFile2.concat,file2.c_str(),OPEN_RDONLY) || !open(multiSeqFile3.concat,file3.c_str(),OPEN_RDONLY) ||
            !open(multiSeqFile4.concat, file4.c_str(), OPEN_RDONLY) || cseq=="" || seqid==-1 )
        return 1;

//The SeqAn library has a built in file parser that can guess the file format
//we use the AutoSeqFormat option for the MiSeq, Library, and barcode files

    //MiSeq File1
    AutoSeqFormat format1;
    guessFormat(multiSeqFile1.concat, format1);
    split(multiSeqFile1, format1);

    //MiSeq File2
    AutoSeqFormat format2;
    guessFormat(multiSeqFile2.concat, format2);
    split(multiSeqFile2, format2);

    //Library
    AutoSeqFormat format3;
    guessFormat(multiSeqFile3.concat, format3);
    split(multiSeqFile3, format3);

    //Barcodes
    AutoSeqFormat format4;
    guessFormat(multiSeqFile4.concat, format4);
    split(multiSeqFile4, format4);

    if(length(multiSeqFile1)!=length(multiSeqFile2)) {
        std::cout << "MiSeq input files contain different number of sequences" << std::endl;
        return 1;
    } else {
        std::cout << "Total Pairs: " << length(multiSeqFile1) << ":" << length(multiSeqFile2) << std::endl;
    }

//Getting the total number of sequences in each of the files
//The two MiSeq files should match in the number of sequences
//They should also be in the same order
    unsigned seqCount1 = length(multiSeqFile1);
    unsigned seqCount3 = length(multiSeqFile3);
    unsigned seqCount4 = length(multiSeqFile4);
    int unmatchedtail=0;
    int unmatchedlib=0;

// FRAGMENT(read_sequences)
    String<char> seq1;
    String<char> seq1id;
    String<char> seq2;
    String<char> seq3;
    String<char> seq3id;
    String<char> seq4;
    String<char> seq4id;
    CharString qual1;
    CharString qual2;
    CharString id1;
    CharString id2;
    THaystacks haystacks;

    std::cout << "Getting " << seqCount4 << " Barcodes.. " << std::endl;

    hashmap idmap;
    std::vector<int> eidlen;
    int max_eidlen=0;
    for(unsigned j=0; j< seqCount4; j++) {
        assignSeqId(seq4id, multiSeqFile4[j], format4);    // read sequence
        assignSeq(seq4, multiSeqFile4[j], format4);    // read sequence
        std::string desc_v=toCString(seq4id);
        std::string eid_v=toCString(seq4);
        std::cout << desc_v.c_str() << " " << eid_v.c_str() << std::endl;
        idmap[eid_v]=desc_v;
        int in=0;
        if(max_eidlen < length(seq4)) {
            max_eidlen=length(seq4);
        }
        for(int w=0; w<eidlen.size(); w++) {
            if(eidlen[w]==length(seq4)) {
                in=1;
                break;
            }
        }
        if(in==0) {
            eidlen.push_back(length(seq4));
        }
    }

//If barcodes are of different lengths print the lengths here

    std::cout << "Barcode Lengths(max=" << max_eidlen <<"):" << std::endl;
    for(int i=0; i<eidlen.size(); i++) {
        std::cout << "length: " << eidlen[i] << std::endl;
    }

// Index library sequences for alignment
    std::cout << "Indexing Sequences(N=" << seqCount3 << ")..";
    for(unsigned j=0; j< seqCount3; j++) {
        assignSeq(seq3, multiSeqFile3[j], format3);    // read sequence
        appendValue(haystacks, seq3);
    }
    Index<THaystacks> index(haystacks);
    Finder<Index<THaystacks> > finder(haystacks);
    std::cout << "completed" << std::endl;

    std::cout << "Running alignment, output should be appearing in " << outfile.c_str() << std::endl;
    for (unsigned i = 0; i < seqCount1; ++i)
    {
        assignSeq(seq1, multiSeqFile1[i], format1);    // read sequence
        assignQual(qual1, multiSeqFile1[i], format1);  // read ascii quality values
        assignSeqId(id1, multiSeqFile1[i], format1);   // read sequence id

        assignSeq(seq2, multiSeqFile2[i], format2);    // read sequence
        assignQual(qual2, multiSeqFile2[i], format2);  // read ascii quality values
        assignSeqId(id2, multiSeqFile2[i], format2);   // read sequence id
        reverseComplement(seq1);
//Get ID and MID
//		std::cout << seq1 << std::endl;
//There are multiple patterns to match however there are always 2 ids
        /*
        String<char> ndl = "[AGCT]{";
        append(ndl,IntToStr(seqid));
        append(ndl,"}");
        append(ndl,"AAAGAAACAACAACAACAAC");
        append(ndl,"[AGCT]{");
        append(ndl,IntToStr(eidlen[0]));
        append(ndl,"}");
        */
        String<char> ndl=cseq;
        Finder<String<char> > finder1(seq1);
//Set options for gap, mismatch,deletion
        Pattern<String<char>, DPSearch<SimpleScore> > pattern1(ndl,SimpleScore(0, -2, -1));
//Matches the Constant region and extracts the ids on either side
        if(find(finder1, pattern1,-2)) {
            int pos1=position(finder1);
            int score=getScore(pattern1);
            while(find(finder1,pattern1,-2)) {
                if(getScore(pattern1) > score) {
                    score=getScore(pattern1);
                    pos1=position(finder1);
                }
            }
//length is eidlen
            int max_pos=max_eidlen;
            if(length(seq1) < pos1+1+max_eidlen) {
                max_pos=length(seq1)-(pos1+1);
            }
            int min_pos=pos1-length(cseq)-seqid+1;
            if(min_pos < 0) {
                min_pos=0;
            }
            String<char> mid=infixWithLength(seq1,(pos1+1),max_pos);
//Length is seqid + constant
            String<char> uid=infixWithLength(seq1,min_pos,seqid);
            append(uid,cseq);
//std::cout << mid << ":" << uid << ::std::endl;
            Pattern<CharString> pattern(uid);
            if(find(finder, pattern)) {
                //std::cout << beginPosition(finder).i1 << std::endl;
                assignSeq(seq3, multiSeqFile3[beginPosition(finder).i1], format3);    // read sequence
                assignSeqId(seq3id,multiSeqFile3[beginPosition(finder).i1], format3);
                Finder<String<char> > finder3(seq3);
                Pattern<String<char>, MyersUkkonen> pattern3(seq2);
                setScoreLimit(pattern3, -10);//Edit Distance = -1. -2 would be worse
                int mpos=-1;
                int mscr=-100;
                while (find(finder3, pattern3)) {
                    int cscr=getScore(pattern3);
                    if(cscr > mscr) {
                        mscr=cscr;
                        mpos=position(finder3);
                    }
                }
                if(mpos >=0) {
                    mpos=mpos-length(seq2);
                    perfect++;
//map mid to description
                    std::string edescr=toCString(mid);
                    for(int j=0; j<eidlen.size(); j++) {
                        String<char> sub=infixWithLength(mid,0,eidlen[j]);
                        std::string str=toCString(sub);
                        hashmap::iterator it = idmap.find(str);
                        if(it != idmap.end()) {
                            edescr=it->second;
                        }
                    }
                    if(debug==1) {
                        fprintf(oFile,"%s,%s,%s,%s,%d,%d,%s,%s\n",toCString(id1),toCString(id2),edescr.c_str(),toCString(seq3id),(mpos+1),mscr,toCString(seq2),toCString(seq3));
                    }
                    else {
                        fprintf(oFile,"%s,%s,%d,%d\n",edescr.c_str(),toCString(seq3id),(mpos+1),mscr);
                    }
                }

            } else {
                unmatchedlib++;
//fprintf(uFile,"%s,%s,%d,%d,%s,%s\n",mid,uid,seq1id);
            }
            clear(finder);

        } else {
            unmatchedtail++;
        }
    }
    std::cout << "Loading " << seqCount1 << " sequences took " << SEQAN_PROTIMEDIFF(loadTime);
    std::cout << " seconds." << std::endl << std::endl;
    std::cout << "Total Sequence Pairs: " << length(multiSeqFile1) << std::endl;
    std::cout << "Perfect ID matches: " << perfect << std::endl;
    std::cout << "Unmatched Tail: " << unmatchedtail << std::endl;
    std::cout << "Unmatched Lib: " << unmatchedlib << std::endl;

    return 0;
}

