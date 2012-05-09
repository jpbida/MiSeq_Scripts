// ==========================================================================
//                          Mason - A Read Simulator
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
// Globally shared code for the Mason read simulator.
// ==========================================================================

#ifndef MASON_H_
#define MASON_H_

#include <numeric>

#include <seqan/random.h>
#include <seqan/sequence_journaled.h>

#include "store_config.h"
#include "util.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// Enum describing the read type to be simulated.
enum ReadsType
{
    READS_TYPE_ILLUMINA,
    READS_TYPE_454,
    READS_TYPE_SANGER
};

// Tag for global options.
typedef void Global;

// Class options.  We will use template inheritance for specializing
// the options class and C++ inheritance to prevent redundant code.
template <typename TTag>
struct Options;

// Program-wide options.
template <>
struct Options<Global>
{
    // true iff help is to be shown.
    bool showHelp;
    // true iff verbosity is enabled.
    bool verbose;
    // true if very verbose is enabled.
    bool veryVerbose;

    // The type of the reads to be simulated.
    ReadsType readsType;

    // Basic Read Simulation Parameters.

    // If set to true then reads are also sampled from regions in the
    // genome that have Ns in them.
    bool allowNFromGenome;
    // Seed to use for the random number generator.
    unsigned seed;
    // Number of reads (pairs) to simulate.
    unsigned numReads;
    // true iff a random reference sequence is to be used.
    bool useRandomSequence;
    // Probability for A, C, G in the random simulated sequence.
    double sourceProbabilityA;
    double sourceProbabilityC;
    double sourceProbabilityG;
    // Length of random sequence to be simulated.
    unsigned randomSourceLength;
    // true iff only the forward strand is to be simulated.
    bool onlyForward;
    // true iff only the reverse strand is to be simulated.
    bool onlyReverse;
    // The output file.  Defaults to REFERENCE-FILE.reads.fastq,
    // possibly with a ".1" or ".2" before the ".fastq" if mate pairs
    // are simulated.
    CharString outputFile;
    // Path to the Sam file to generate.  Defaults to fastq file name
    // with suffix ".sam"
    CharString samFile;
    // true iff qualities are to be simulated.
    bool simulateQualities;
    // true iff additional information is to be included in the reads file.
    bool includeReadInformation;
    // Path to file with sample counts from each contig.  The file is a
    // TSV file containing pairs mapping the contig id (non-whitespace
    // sequence in FASTA header) to a read count.  Overrides parameter
    // numReads if given.
    CharString sampleCountsFilename;

    // Mate-Pair Related Options.

    // true iff generating mate pairs is enabled.
    bool generateMatePairs;
    // true iff mate pair library sizes are to be uniformly distributed,
    // otherwise standard distribution is used.
    bool libraryLengthIsUniform;
    // Mate-pair library mean length.
    double libraryLengthMean;
    // Mate-pair library length error.  Standard deviation for normally
    // distributed library lengths, interval length around mean for uniform
    // distribution.
    double libraryLengthError;

    // Haplotype parameters.

    // Number of haplotypes to generated.  All are generated from the input genome.
    unsigned numHaplotypes;
    // SNP rate.
    double haplotypeSnpRate;
    // Indel rate.
    double haplotypeIndelRate;
    // Smallest number of indels.
    unsigned haplotypeIndelRangeMin;
    // Largest number of indels.
    unsigned haplotypeIndelRangeMax;
    // If true then no Ns are substituted or inserted into the haplotype.
    bool haplotypeNoN;

    Options()
            : showHelp(false),
              verbose(false),
              veryVerbose(false),
              allowNFromGenome(false),
              seed(0),
              numReads(1000),
              useRandomSequence(false),
              sourceProbabilityA(0.25),
              sourceProbabilityC(0.25),
              sourceProbabilityG(0.25),
              randomSourceLength(1000*1000),
              onlyForward(false),
              onlyReverse(false),
              outputFile(""),
              samFile(""),
              simulateQualities(false),
              includeReadInformation(false),
              generateMatePairs(false),
              libraryLengthIsUniform(false),
              libraryLengthMean(1000),
              libraryLengthError(100),
              numHaplotypes(1),
              haplotypeSnpRate(0.001),
              haplotypeIndelRate(0.001),
              haplotypeIndelRangeMin(1),
              haplotypeIndelRangeMax(6),
              haplotypeNoN(false)
    {}
};

// Use this container for model specific parameters generated before the actual
// simulation.  Setup in void simulateReadsSetupModelSpecificData(...).
template <typename TTag>
struct ModelParameters;

// Global model parameters.
template <>
struct ModelParameters<Global>
{
    // If non-empty, sampleCounts[i] gives the number of reads to
    // sample from contig i (the i-th sequence in the FASTA input
    // file).
    String<size_t> sampleCounts;
};

// Enum describing the type of an error.
enum ErrorType
{
    ERROR_TYPE_MATCH    = 0,
    ERROR_TYPE_MISMATCH = 1,
    ERROR_TYPE_INSERT   = 2,
    ERROR_TYPE_DELETE   = 3
};

// Class for storing the read simulation instructions.  Will be
// specialized for each technology to simulate.
template <typename TReadTypeTag>
struct ReadSimulationInstruction;

// Read simulation instructions used by all simulated technologies.
template <>
struct ReadSimulationInstruction<Global>
{
    // Index of the haplotype to sample the read from.
    unsigned haplotype;
    // Index of the contig to sample the read from.
    unsigned contigId;
    // Whether or not the read is to be sampled from the forward strand.
    bool isForward;
    // Begin and end position of the infix to sample the read from.
    size_t beginPos;
    size_t endPos;
    // Number of characters added/removed to the string by indels.
    unsigned delCount;
    unsigned insCount;
    // Edit string of the read.
    String<ErrorType> editString;
    // String of qualities for the bases written out.
    String<int> qualities;

    ReadSimulationInstruction() : delCount(0), insCount(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Prints the global options to stream.
template <typename TStream>
TStream & operator<<(TStream & stream, Options<Global> const & options) {
    stream << "global-options {" << std::endl
           << "  allowNFromGenome:       " << (options.allowNFromGenome ? "true" : "false") << std::endl
           << "  seed:                   " << options.seed << std::endl
           << "  numReads:               " << options.numReads << std::endl
           << "  useRandomSequence:      " << (options.useRandomSequence ? "true" : "false") << std::endl
           << "  sourceProbabilityA:     " << options.sourceProbabilityA << std::endl
           << "  sourceProbabilityC:     " << options.sourceProbabilityC << std::endl
           << "  sourceProbabilityG:     " << options.sourceProbabilityG << std::endl
           << "  useRandomSequence:      " << (options.useRandomSequence ? "true" : "false") << std::endl
           << "  randomSourceLength:     " << options.randomSourceLength << std::endl
           << "  onlyForward:            " << (options.onlyForward ? "true" : "false") << std::endl
           << "  onlyReverse:            " << (options.onlyReverse ? "true" : "false") << std::endl
           << "  outputFile:             \"" << options.outputFile << "\"" << std::endl
           << "  samFile:                \"" << options.samFile << "\"" << std::endl
           << "  simulateQualities:      " << (options.simulateQualities ? "true" : "false") << std::endl
           << "  includeReadInformation: " << options.includeReadInformation << std::endl
           << "  generateMatePairs:      " << (options.generateMatePairs ? "true" : "false") << std::endl
           << "  libraryLengthMean:      " << options.libraryLengthMean << std::endl
           << "  libraryLengthError:     " << options.libraryLengthError << std::endl
           << "  numHaplotypes:          " << options.numHaplotypes << std::endl
           << "  haplotypeSnpRate:       " << options.haplotypeSnpRate << std::endl
           << "  haplotypeIndelRate:     " << options.haplotypeIndelRate << std::endl
           << "  haplotypeIndelRangeMin: " << options.haplotypeIndelRangeMin << std::endl
           << "  haplotypeIndelRangeMax: " << options.haplotypeIndelRangeMax << std::endl
           << "  haplotypeNoN:           " << options.haplotypeNoN << std::endl
           << "  sampleCountsFilename:   " << options.sampleCountsFilename << std::endl
           << "}" << std::endl;
    return stream;
}

// Stream operator for simulation instructions.
template <typename TStream>
TStream & operator<<(TStream & stream, ReadSimulationInstruction<Global> const & inst) {
    stream << "(haplotype=" << inst.haplotype << ", contigId=" << inst.contigId << ", isForward=" << inst.isForward << ", beginPos=" << inst.beginPos << ", endPos=" << inst.endPos << ", insCount=" << inst.insCount << ", delCount=" << inst.delCount << ", editString=";
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        stream << "MEID"[inst.editString[i]];
    }
    stream << ", qualities=[";
    for (unsigned i = 0; i < length(inst.qualities); ++i) {
        if (i != 0)
            stream << ", ";
        stream << inst.qualities[i];
    }
    stream << "])";
    return stream;
}

// Initialize the command line parser for the global options.
void setUpCommandLineParser(CommandLineParser & parser)
{
    addVersionLine(parser, "0.1");

    addTitleLine(parser, "Mason - A Read Simulator");
    addUsageLine(parser, "illumina [OPTIONS] SEQUENCE");
    addLine(parser, "");
    addLine(parser, "Use 'random' for the SEQUENCE file name to generate it randomly.");

    addSection(parser, "Main Options");
    
    addOption(parser, CommandLineOption("aNg",  "allow-N-from-genome", "Allow N from genome.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("s",  "seed", "The seed for Rng.  Default: 0.", OptionType::Integer | OptionType::Label));
    addOption(parser, CommandLineOption("N",  "num-reads", "Number of reads (or mate pairs) to simulate.  Default: 1000.", OptionType::Integer));
    addOption(parser, CommandLineOption("sn", "source-length", "Length of random source sequence.  Default: 1,000,000.", OptionType::Integer));
    addOption(parser, CommandLineOption("spA", "source-probability-A", "Propabilibty for A in randomly generated sequence. Default: 0.25", OptionType::Double));
    addOption(parser, CommandLineOption("spC", "source-probability-C", "Propabilibty for C in randomly generated sequence. Default: 0.25", OptionType::Double));
    addOption(parser, CommandLineOption("spG", "source-probability-G", "Propabilibty for G in randomly generated sequence. Default: 0.25", OptionType::Double));
    addOption(parser, CommandLineOption("snN", "source-no-N", "If set then no Ns are generated in the random source sequence.", OptionType::Bool));
    addOption(parser, CommandLineOption("f",  "forward-only", "Simulate from forward strand only.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("r",  "reverse-only", "Simulate from reverse strand only.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("o",  "output-file", "Write results to PARAM.fasta file instead of SEQUENCE.reads.fasta.  Default: \"\".", OptionType::String));
    addOption(parser, CommandLineOption("sq", "simulate-qualities", "Simulate qualities, generate FASTQ instead of FASTA.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("i", "include-read-information", "Include additional read information in reads file.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("v", "verbose", "Verbosity mode.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "High verbosity mode, implies verbosity mode.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("scf", "sample-counts-file", "Path to TSV file that maps contig ids to read counts.", OptionType::String));

    addSection(parser, "Mate-Pair Options");

    addOption(parser, CommandLineOption("ll", "library-length-mean", "Mate-pair mean library length.  Default: 1000.", OptionType::Double));
    addOption(parser, CommandLineOption("le", "library-length-error", "Mate-pair library tolerance.  Default: 100.", OptionType::Double));
    addOption(parser, CommandLineOption("mp", "mate-pairs", "Enable mate pair simulation.  Default: false.", OptionType::Bool));

    addSection(parser, "Haplotype Options");

    addOption(parser, CommandLineOption("hn", "num-haplotypes", "Number of haplotypes to simulate.  Default: 1.", OptionType::Integer));
    addOption(parser, CommandLineOption("hs", "haplotype-snp-rate", "Haplotype SNP rate.  Default: 0.001.", OptionType::Double));
    addOption(parser, CommandLineOption("hi", "haplotype-indel-rate", "Haplotype indel rate.  Default: 0.001.", OptionType::Double));
    addOption(parser, CommandLineOption("hm", "haplotype-indel-range-min", "Haplotype indel size min.  Default: 1.", OptionType::Integer));
    addOption(parser, CommandLineOption("hM", "haplotype-indel-range-max", "Haplotype indel size max.  Default: 6.", OptionType::Integer));
    addOption(parser, CommandLineOption("hnN", "haplotype-no-N", "Do not allow Ns to be substituted or inserted into N.  Default: Is allowed.", OptionType::Bool));

    // Need reads type and {SEQUENCE.fasta, random}.
    requiredArguments(parser, 2);
}

// Parse command line for global options and perform some checks.
template <typename TOptions>
int parseCommandLineAndCheck(TOptions & options,
                             CharString & referenceFilename,
                             CommandLineParser & parser,
                             const int argc,
                             const char * argv[]) {
    if (!parse(parser, argc, argv)) {
        if (!isSetShort(parser, "h"))
            shortHelp(parser, std::cerr);
        return 1;
    }
    if (isSetShort(parser, "h")) {
        options.showHelp = true;
        return 0;
    }

    if (isSetLong(parser, "allow-N-from-genome"))
        options.allowNFromGenome = true;
    if (isSetLong(parser, "seed"))
        getOptionValueLong(parser, "seed", options.seed);
    if (isSetLong(parser, "num-reads"))
        getOptionValueLong(parser, "num-reads", options.numReads);
    if (isSetLong(parser, "source-length"))
        getOptionValueLong(parser, "source-length", options.randomSourceLength);
    if (isSetLong(parser, "source-probability-A"))
        getOptionValueLong(parser, "source-probability-A", options.sourceProbabilityA);
    if (isSetLong(parser, "source-probability-C"))
        getOptionValueLong(parser, "source-probability-C", options.sourceProbabilityC);
    if (isSetLong(parser, "source-probability-G"))
        getOptionValueLong(parser, "source-probability-G", options.sourceProbabilityG);
    if (isSetLong(parser, "forward-only"))
        options.onlyForward = true;
    if (isSetLong(parser, "reverse-only"))
        options.onlyReverse = true;
    if (isSetLong(parser, "output-file"))
        getOptionValueLong(parser, "output-file", options.outputFile);
    if (isSetLong(parser, "simulate-qualities"))
        options.simulateQualities = true;
    if (isSetLong(parser, "include-read-information"))
        options.includeReadInformation = true;
    if (isSetLong(parser, "verbose"))
        options.verbose = true;
    if (isSetLong(parser, "very-verbose")) {
        options.verbose = true;
        options.veryVerbose = true;
    }

    if (isSetLong(parser, "library-length-mean"))
        getOptionValueLong(parser, "library-length-mean", options.libraryLengthMean);
    if (isSetLong(parser, "library-error"))
        getOptionValueLong(parser, "library-length-error", options.libraryLengthError);
    if (isSetLong(parser, "mate-pairs"))
        options.generateMatePairs = true;

    if (isSetLong(parser, "num-haplotypes"))
        getOptionValueLong(parser, "num-haplotypes", options.numHaplotypes);
    if (isSetLong(parser, "haplotype-snp-rate"))
        getOptionValueLong(parser, "haplotype-snp-rate", options.haplotypeSnpRate);
    if (isSetLong(parser, "haplotype-indel-rate"))
        getOptionValueLong(parser, "haplotype-indel-rate", options.haplotypeIndelRate);
    if (isSetLong(parser, "haplotype-indel-range-min"))
        getOptionValueLong(parser, "haplotype-indel-range-min", options.haplotypeIndelRangeMin);
    if (isSetLong(parser, "haplotype-indel-range-max"))
        getOptionValueLong(parser, "haplotype-indel-range-max", options.haplotypeIndelRangeMax);
    if (isSetLong(parser, "haplotype-no-N"))
        options.haplotypeNoN = true;

    if (isSetLong(parser, "sample-counts-file"))
        getOptionValueLong(parser, "sample-counts-file", options.sampleCountsFilename);

    // First argument is "illumina", second one name of reference file.
    referenceFilename = getArgumentValue(parser, 1);

    if (referenceFilename == "random")
        options.useRandomSequence = true;

    return parseCommandLineAndCheckModelSpecific(options, parser);
}

// Write a random DNA sequence of the given length to the file with the given name.
template <typename TRNG, typename TOptions>
int writeRandomSequence(TRNG & rng, size_t length, CharString const & fileName, TOptions const & options)
{
    DnaString randomSequence;
    reserve(randomSequence, length);

    // Simulate the random sequence with the given background probabilities.
    double xA = options.sourceProbabilityA;
    double xC = xA + options.sourceProbabilityC;
    double xG = xC + options.sourceProbabilityG;

    Pdf<Uniform<double> > pdf(0, 1);
    for (size_t i = 0; i < length; ++i) {
        double x = pickRandomNumber(rng, pdf);
        if (x < xA)
            appendValue(randomSequence, Dna('A'));
        else if (x < xC)
            appendValue(randomSequence, Dna('C'));
        else if (x < xG)
            appendValue(randomSequence, Dna('G'));
        else
            appendValue(randomSequence, Dna('T'));
    }

    std::ofstream file;
    file.open(toCString(fileName), std::ios_base::out | std::ios_base::trunc);
    if (!file.is_open()) {
        std::cerr << "Failed to write random sequence to " << fileName << std::endl;
        return 1;
    }
    write(file, randomSequence, "random_sequence", Fasta());
    file.close();
    return 0;
}

// Load sample counts as integers from a file.
template <typename TFragmentStore, typename TSpec, typename TOptions>
int loadSampleCounts(ModelParameters<TSpec> & modelParameters, TFragmentStore /*const*/ & fragmentStore, TOptions const & options)
{
    std::ifstream file;
    file.open(toCString(options.sampleCountsFilename), std::ios_base::in);
    if (!file.is_open()) {
        std::cerr << "Failed to open sample counts file " << options.sampleCountsFilename << std::endl;
        return 1;
    }

    resize(modelParameters.sampleCounts, length(fragmentStore.contigNameStore), 0);

    char c = _streamGet(file);
    CharString contigName;
    size_t sampleCount;

    refresh(fragmentStore.contigNameStoreCache);

    while (!_streamEOF(file)) {
        clear(contigName);

        _parseReadSamIdentifier(file, contigName, c);
        _parseSkipUntilChar(file, '\t', c);
        sampleCount = _parseReadNumber(file, c);
        _parseSkipLine(file, c);

        size_t contigId = 0;
        if (!getIdByName(fragmentStore.contigNameStore, contigName, contigId, fragmentStore.contigNameStoreCache)) {
            std::cerr << "ERROR: Could not find contig with name \"" << contigName << "\" (from read counts file) in contigs." << std::endl;
            return 1;
        }

        if (options.veryVerbose)
            std::cout << "Sample count for contig " << contigName << " (id=" << contigId << ") is " << sampleCount << std::endl;
        modelParameters.sampleCounts[contigId] = sampleCount;
    }

    return 0;
}

// Top level dispatched function that mainly contains the input/output
// logic and dispatches to simulateReadsMain() for the actual
// simulation steps.
template <typename TOptions, typename TReadsTypeTag>
int simulateReads(TOptions options, CharString referenceFilename, TReadsTypeTag const &) {
    // Print options.
    std::cerr << options;
    std::cerr << "reference file: " << referenceFilename << std::endl;

    // Create a Rng with the given seed.
    Rng<MersenneTwister> rng(options.seed);

    // Load or randomly generate the reference sequence.
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    if (options.useRandomSequence) {
        referenceFilename = "random.fasta";
        std::cerr << "Generating random sequence of length " << options.randomSourceLength
                  << " to file \"" << referenceFilename << "\"." << std::endl;
        int ret = writeRandomSequence(rng, options.randomSourceLength, referenceFilename, options);
        if (ret != 0)
            return ret;
    }
    // Generate output file name if necessary.
    if (options.outputFile == "") {
        options.outputFile = referenceFilename;
        if (options.simulateQualities)
            append(options.outputFile, ".fastq");
        else
            append(options.outputFile, ".fasta");
    }
    if (options.samFile == "") {
        options.samFile = options.outputFile;
        append(options.samFile, ".sam");
    }
    std::cerr << "Loading reference sequence from \"" << referenceFilename << "\"" << std::endl;
    if (!loadContigs(fragmentStore, referenceFilename)) {
        std::cerr << "Could not load reference sequence from " << referenceFilename << std::endl;
        return 1;
    }

    // Load and/or build the model specific parameters for the simulation.
    ModelParameters<TReadsTypeTag> modelParameters;
    int ret = simulateReadsSetupModelSpecificData(modelParameters, options);
    if (ret != 0)
        return ret;

    // Load sample counts file if given.
    if (length(options.sampleCountsFilename) > 0) {
        ret = loadSampleCounts(modelParameters, fragmentStore, options);
        if (ret != 0)
            return ret;
    }
    
    // Kick off the read generation.
    ret = simulateReadsMain(fragmentStore, rng, options, modelParameters);
    if (ret != 0)
        return ret;

    // Write out results.
    if (options.generateMatePairs) {
        // Build filename with '.1.' infix.
        CharString mateFilename = options.outputFile;
        size_t dotPos = 0;
        for (size_t i = 0; i < length(mateFilename); ++i)
            if (mateFilename[i] == '.')
                dotPos = i;
        infix(mateFilename, dotPos, dotPos + 1) = "_1.";
        // Write out first mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/1" << std::endl;
        StringSet<String<Dna5Q>, Dependent<> > reads;
        StringSet<CharString, Dependent<> > readNames;
        for (size_t i = 0; i < length(fragmentStore.matePairStore); ++i) {
            size_t readId = fragmentStore.matePairStore[i].readId[0];
            appendValue(readNames, fragmentStore.readNameStore[readId]);
            appendValue(reads, fragmentStore.readSeqStore[readId]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (!fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            if (options.simulateQualities)
                write(fstrm, readNames, reads, Fastq());
            else
                write(fstrm, readNames, reads, Fasta());
        }
        // Build filename with '_2.' infix.
        infix(mateFilename, dotPos + 1, dotPos + 2) = "2";
        // Write out second mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/2" << std::endl;
        clear(reads);
        clear(readNames);
        for (size_t i = 0; i < length(fragmentStore.matePairStore); ++i) {
            size_t readId = fragmentStore.matePairStore[i].readId[1];
            appendValue(readNames, fragmentStore.readNameStore[readId]);
            appendValue(reads, fragmentStore.readSeqStore[readId]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (!fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            if (options.simulateQualities)
                write(fstrm, readNames, reads, Fastq());
            else
                write(fstrm, readNames, reads, Fasta());
        }
    } else {
        std::cerr << "Writing resulting reads to \"" << options.outputFile << "\"" << std::endl;
        std::fstream fstrm(toCString(options.outputFile), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << options.outputFile << "\"" << std::endl;
            return 1;
        }
        if (options.simulateQualities)
            write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fastq());
        else
            write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fasta());
    }
    std::cerr << "Writing Sam file \"" << options.samFile << "\"" << std::endl;
    {
        std::fstream fstrm(toCString(options.samFile), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open Sam file \"" << options.samFile << "\"" << std::endl;
            return 1;
        }
        write(fstrm, fragmentStore, Sam());
    }
    return 0;
}

// Build a haplotype, based on the contigs from the given fragment store.
template <typename TRNG>
void buildHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotype,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    TRNG & rng,
                    Options<Global> const & options) {
    resize(haplotype, length(fragmentStore.contigStore), Exact());
    String<Dna5> buffer;
    reserve(buffer, options.haplotypeIndelRangeMax);

    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
        std::cout << "    contig # " << i+1 << "/" << length(fragmentStore.contigStore) << std::endl;
        clear(haplotype[i]);
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        String<Dna5, Journaled<Alloc<> > > & haplotypeContig = haplotype[i];

        // Only generate Ns in the haplotype if allowed.
        int maxOrdValue = options.haplotypeNoN ? 3 : 4;

        // j is position in original sequence, k is position in haplotype
        for (size_t j = 0, k = 0; j < length(contig);) {
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x < options.haplotypeSnpRate) {
                // SNP
                Dna5 c = Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue - 1)));
                if (c == contig[j])
                    c = Dna5(ordValue(c) + 1);
                if (options.haplotypeNoN)
                    SEQAN_ASSERT_TRUE(c != Dna5('N'));
                assignValue(haplotypeContig, k, c);
                j += 1;
                k += 1;
            } else if (x < options.haplotypeSnpRate + options.haplotypeIndelRate) {
                // Indel of random length.
                unsigned rangeLen = options.haplotypeIndelRangeMax - options.haplotypeIndelRangeMin;
                unsigned indelLen = options.haplotypeIndelRangeMin + static_cast<unsigned>(pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) * rangeLen);
                if (pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) < 0.5) {
                    // Insertion.
                    clear(buffer);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                        appendValue(buffer, Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue))));
                    insert(haplotypeContig, k, buffer);
                    k += indelLen;
                } else {
                    // Deletion.
                    indelLen = _min(indelLen, length(haplotypeContig) - k);
                    erase(haplotypeContig, k, k + indelLen);
                    j += indelLen;
                }
            } else {
                // Match.
                j += 1;
                k += 1;
            }
        }
    }
}

// Build a read simulation instructions for a haplotype.
//
// pick a contig, probability is proportional to the length
// pick a start position, end position = start position + read length
// pick whether to match on the forward or reverse strand
// simulate edit string
// build quality values
// possibly adjust mate if left read has insert at the beginning or right read has insert at the right
template <typename TReadsTag, typename TRNG>
int buildReadSimulationInstruction(
        String<ReadSimulationInstruction<TReadsTag> > & instructions,
        TRNG & rng,
        unsigned const & haplotypeId,
        StringSet<String<Dna5, Journaled<Alloc<> > > > const & haplotype,
        String<double> const & relativeContigLengths,
        size_t const & contigId,
        bool fixedContigId,
        ModelParameters<TReadsTag> const & parameters,
        Options<TReadsTag> const & options)
{
    ReadSimulationInstruction<TReadsTag> inst;
    inst.haplotype = haplotypeId;

    // We have to retry simulation if the mate pair did not fit in.
    bool invalid = false;
    do {
        clear(instructions);
        invalid = false;  // By default, we do not want to repeat.
        if (fixedContigId) {
            // Use precomputed contig id.
            inst.contigId = contigId;
        } else {
            // Pick contig id, probability is proportional to the length.
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            for (unsigned i = 0; i < length(relativeContigLengths); ++i) {
                if (x < relativeContigLengths[i]) {
                    inst.contigId = i - 1;
                    break;
                }
            }
        }
        // Pick whether on forward or reverse strand.
        if (options.onlyForward)
            inst.isForward = true;
        else if (options.onlyReverse)
            inst.isForward = false;
        else
            inst.isForward = pickRandomNumber(rng, Pdf<Uniform<int> >(0, 1));
        // Pick the length in the haplotype infix the read comes from, possibly randomly.
        unsigned readLength = pickReadLength(rng, options);
        // This cannot work if the haplotype is shorter than the length of the read to simulate.
        if (length(haplotype[inst.contigId]) < readLength) {
            std::cerr << "ERROR: haplotype (== " << length(haplotype[inst.contigId]) << ") < read length!" << std::endl;
            return 1;
        }
        // Pick a start and end position.
        inst.beginPos = pickRandomNumber(rng, Pdf<Uniform<size_t> >(0, length(haplotype[inst.contigId]) - readLength - 1));
        inst.endPos = inst.beginPos + readLength;
        // Simulate the read with these parameters.
        buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
        // Append read to result list.
        appendValue(instructions, inst);

        ReadSimulationInstruction<TReadsTag> inst2(inst);

        // Maybe create a mate for this read.
        if (options.generateMatePairs) {
            // Pick a read length, possibly randomly.
            unsigned readLength = pickReadLength(rng, options);
            // Pick a library length, according to the options.
            size_t libraryLength = pickLibraryLength(rng, options);
            // Compute start and end position.
            inst.endPos = inst.beginPos + libraryLength;
            inst.beginPos = inst.endPos - readLength;
            // Verify that the mate fits right of the originally simulated read.
            size_t contigLength = length(haplotype[inst.contigId]);
            if ((inst.beginPos > contigLength) || (inst.endPos > contigLength)) {
                // Mate did not fit!  Remove previously added read and set
                // invalid to true so we repeat this simulation.
                SEQAN_ASSERT_GT(length(instructions), 0u);
                eraseBack(instructions);
                invalid = true;
                if (options.verbose) {
                    std::cerr << "INFO: Mate did not fit! Repeating..." << std::endl;
                    std::cerr << "      inst2.beginPos == " << inst2.beginPos << ", inst2.endPos == " << inst2.endPos << std::endl;
                    std::cerr << "       inst.beginPos == " << inst.beginPos << ",  inst.endPos == " << inst.endPos << std::endl;
                }
                continue;
            }
            // Simulate the read with these parameters.
            buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
            // Append read to result list.
            appendValue(instructions, inst);
        }

        // Check whether there are Ns in the selected areas.
        if (!options.allowNFromGenome) {
            for (unsigned i = 0; i < length(instructions); ++i) {
                String<Dna5, Journaled<Alloc<> > > const & contig = haplotype[instructions[i].contigId];
                typedef typename Position<Dna5String>::Type TPosition;
                TPosition beginPos = instructions[i].beginPos;
                TPosition endPos = instructions[i].endPos;
                if (beginPos > endPos)
                    std::swap(beginPos, endPos);
                for (unsigned i = beginPos; i != endPos; ++i) {
                    if (contig[i] == Dna5('N')) {
                        invalid = true;
                        break;
                    }
                }
            }
        }
    } while (invalid);
	
	if (options.generateMatePairs)
		SEQAN_ASSERT_EQ(length(instructions), 2u);
	else
		SEQAN_ASSERT_EQ(length(instructions), 1u);
	
    return 0;
}

// Pick library length, based on the configuration in options.
template <typename TRNG>
inline
unsigned pickLibraryLength(TRNG & rng, Options<Global> const & options)
{
    if (options.libraryLengthIsUniform) {
        // Pick uniformly.
        double minLen = options.libraryLengthMean - options.libraryLengthError;
        double maxLen = options.libraryLengthMean + options.libraryLengthError;
        double len = pickRandomNumber(rng, Pdf<Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(_max(0.0, round(len)));
    } else {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, Pdf<Normal>(options.libraryLengthMean, options.libraryLengthError));
        return static_cast<unsigned>(_max(0.0, round(len)));
    }
}

// Performs the actual read simulation.
template <typename TRNG, typename TReadsTag, typename TOptions>
int simulateReadsMain(FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                      TRNG & rng,
                      TOptions const & options,
                      ModelParameters<TReadsTag> const & parameters) {
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;
    typedef Value<TFragmentStore::TMatePairStore>::Type TMatePairStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type	TAlignedElement;
    typedef typename TAlignedElement::TGapAnchors								TReadGapAnchors;

    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

    typedef Position<CharString>::Type TPos;

    // Number of reads comes from command line by default.  If sample
    // counts are given, we compute it from this instead.
    size_t numReads = options.numReads;
    if (length(parameters.sampleCounts) > 0u) {
        numReads = 0;
        for (unsigned i = 0; i < length(parameters.sampleCounts); ++i)
            numReads += parameters.sampleCounts[i];
    }

    // First, we randomly pick the haplotype for each read to be
    // simulated or read it from the sample counts in parameters.
    String<unsigned> haplotypeIds;
    // Pick random haplotype origin.
    reserve(haplotypeIds, numReads);
    for (size_t i = 0; i < numReads; ++i)
        appendValue(haplotypeIds, pickRandomNumber(rng, Pdf<Uniform<unsigned> >(0, options.numHaplotypes - 1)));

    // Maybe pick contig ids to sample from.
    String<unsigned> contigIds;
    if (length(parameters.sampleCounts) > 0) {
        for (unsigned i = 0; i < length(fragmentStore.contigNameStore); ++i) {
            resize(contigIds, length(contigIds) + parameters.sampleCounts[i], i);
            if (options.veryVerbose)
                std::cerr << parameters.sampleCounts[i] << " reads from haplotype " << i << "..." << std::endl;
        }
        shuffle(contigIds, rng);
    }

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // TODO(holtgrew): Would only have to switch pointers to journals which is possible.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    reserve(fragmentStore.readSeqStore, numReads, Exact());
    reserve(fragmentStore.readNameStore, numReads, Exact());
    char readName[1024];
    char outFileName[151];
    snprintf(outFileName, 150, "%s", toCString(options.outputFile));
    String<bool> flipped;
    for (unsigned haplotypeId = 0; haplotypeId < options.numHaplotypes; ++haplotypeId) {
        std::cerr << "Simulating for haplotype #" << haplotypeId << "..." << std::endl;
        std::cout << "  Building haplotype..." << std::endl;
        StringSet<String<Dna5, Journaled<Alloc<> > > > haplotypeContigs;
        buildHaplotype(haplotypeContigs, fragmentStore, rng, options);

        // Build partial sums over relative contig lengths so we can pick the contigs later on.
        size_t totalLength = 0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
            totalLength += length(fragmentStore.contigStore[i].seq);
        String<double> relativeContigLengths;
        resize(relativeContigLengths, length(fragmentStore.contigStore) + 1, Exact());
        front(relativeContigLengths) = 0.0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            double l = static_cast<double>(length(fragmentStore.contigStore[i].seq));
            relativeContigLengths[i + 1] = l / totalLength;
        }
        std::partial_sum(begin(relativeContigLengths), end(relativeContigLengths), begin(relativeContigLengths));
        back(relativeContigLengths) = 1.0;

        // Simulate the reads...
        std::cerr << "  Simulating reads for haplotype #" << haplotypeId << "..." << std::endl;

//         std::cerr << "Journal: " << haplotypeContigs[0]._journalEntries << std::endl;

        for (unsigned j = 0; j < length(haplotypeIds); ++j) {
            if (haplotypeIds[j] != haplotypeId)
                continue;  // Guard against instructions on wrong haplotype.

            // Build simulation instructions.
            String<ReadSimulationInstruction<TReadsTag> > instructions;
            // TODO(holtgrew): Pick contig id outside of instructions.
            size_t contigId = 0;
            bool fixedContigId = false;
            if (length(parameters.sampleCounts) > 0) {
                fixedContigId = true;
                if (options.generateMatePairs)
                    contigId = contigIds[length(fragmentStore.readSeqStore) / 2];
                else
                    contigId = contigIds[length(fragmentStore.readSeqStore)];
            }
            int res = buildReadSimulationInstruction(instructions, rng, haplotypeId, haplotypeContigs, relativeContigLengths, contigId, fixedContigId, parameters, options);
            if (res != 0)
                return res;

            int previousMateNum = 0;
            for (unsigned k = 0; k < length(instructions); ++k) {
                ReadSimulationInstruction<TReadsTag> const & inst = instructions[k];
                // Apply simulation instructions.
                SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.readNameStore));
                // Cut out segment from haplotype.
                String<Dna5Q> read = infix(haplotypeContigs[inst.contigId], inst.beginPos, inst.endPos);
                String<Dna5Q> haplotypeInfix = read;  // Copy for printing later on.
                applySimulationInstructions(read, rng, inst, options);
                // Append read sequence to read seq store and mate pair to read
                // name store.  This also yields the read id.  We will generate
                // and append the read name below, depending on the read id.
                unsigned readId;
                if (options.generateMatePairs)
                    readId = appendRead(fragmentStore, read, length(fragmentStore.matePairStore));
                else
                    readId = appendRead(fragmentStore, read);

                // Get expected begin/end position in the original sequence.  If we decide to flip this read later on, we will modify the WIT store.
                TPos origBeginPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.beginPos);
                TPos origEndPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.endPos);

                // Generate read name.
                // TODO(holtgrew): Remove mateNum, not required?
                if (options.generateMatePairs) {
                    // Generate the mate num \in {1, 2}, randomly but consistent so two entries belonging together have different nums.  This also decides about the flipping.
                    int mateNum = 3 - previousMateNum;
                    if (readId % 2 == 0) {
                        mateNum = pickRandomNumber(rng, Pdf<Uniform<int> >(1, 2));
						SEQAN_ASSERT_GEQ(mateNum, 1);
						SEQAN_ASSERT_LEQ(mateNum, 2);
                        previousMateNum = mateNum;
                        appendValue(flipped, mateNum == 2);
                    } else {
						SEQAN_ASSERT_EQ(flipped[readId - 1], mateNum == 1);
                        appendValue(flipped, mateNum == 1);
                    }
                    if (options.includeReadInformation)
                        sprintf(readName, "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu haplotype_infix=%s edit_string=", outFileName, readId / 2, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos), toCString(CharString(haplotypeInfix)));
                    else
                        sprintf(readName, "%s.%09u", outFileName, readId / 2);
                } else {
                    if (options.includeReadInformation)
                        sprintf(readName, "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu haplotype_infix=%s edit_string=", outFileName, readId, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos), toCString(CharString(haplotypeInfix)));
                    else
                        sprintf(readName, "%s.%09u", outFileName, readId);
                }
                if (options.includeReadInformation) {
                    for (unsigned i = 0; i < length(inst.editString); ++i) {
                        char buffer[2] = "*";
                        buffer[0] = "MEID"[static_cast<int>(inst.editString[i])];
                        strcat(readName, buffer);
                    }
                }
                appendValue(fragmentStore.readNameStore, readName);

                // Print info about read and haplotype.
                if (options.veryVerbose) {
                    std::cout << ",-- Read #" << readId << std::endl
                              << "| inst.beginPos    " << inst.beginPos << std::endl
                              << "| inst.endPos      " << inst.endPos << std::endl
                              << "| origBeginPos     " << origBeginPos << std::endl
                              << "| origEndPos       " << origEndPos << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos-1) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos+1) << std::endl
                              << "| name:            " << readName << std::endl
                              << "| original infix:  " << infix(fragmentStore.contigStore[inst.contigId].seq, origBeginPos, origEndPos) << std::endl
                              << "| haplotype infix: " << infix(haplotypeContigs[inst.contigId], inst.beginPos, inst.endPos) << std::endl
                              << "| read:            " << read << std::endl
                              << "`-- " << std::endl;
                }

                // Tentatively add matches to aligned read store.  We will
                // maybe flip begin and end position below in the "flipping and
                // reordering" step and convert the matches to a global
                // alignment in the "convertMatchesToGlobalAlignment" call.
                // std::cout << "origBeginPos = " << origBeginPos << ", origEndPos = " << origEndPos << std::endl;
                if (options.generateMatePairs)
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos, length(fragmentStore.matePairStore));
                else
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos);

                // Perform flipping and reordering.
                if (options.generateMatePairs) {
                    if (readId % 2 == 1) {  // Only flip and append mate pair info after simulating second mate.
                        // Append mate pair element to fragment store's mate pair store.
                        TMatePairStoreElement matePair;
                        matePair.readId[0] = readId - 1 + flipped[readId];
                        matePair.readId[1] = readId - flipped[readId];
                        appendValue(fragmentStore.matePairStore, matePair);

                        // The first mate always comes from the forward strand.
                        if (options.includeReadInformation)
                            append(fragmentStore.readNameStore[readId - 1], " strand=forward");
                        // The second read always comes from the reverse strand.
                        reverseComplement(fragmentStore.readSeqStore[readId]);
                        if (options.includeReadInformation)
                            append(fragmentStore.readNameStore[readId], " strand=reverse");
                        // Note: readId is also last index of aligned read store because we only have one alignment per read!
                        std::swap(fragmentStore.alignedReadStore[readId].beginPos, fragmentStore.alignedReadStore[readId].endPos);
                    }
                } else {
                    if (pickRandomNumber(rng, Pdf<Uniform<double> >(0.0, 1.0)) < 0.5) {
                        reverseComplement(back(fragmentStore.readSeqStore));
                        if (options.includeReadInformation)
                            append(back(fragmentStore.readNameStore), " strand=reverse");
                        // Note: readId is also last index of aligned read store because we only have one alignment per read!
                        std::swap(fragmentStore.alignedReadStore[readId].beginPos, fragmentStore.alignedReadStore[readId].endPos);
                    } else {
                        if (options.includeReadInformation)
                            append(back(fragmentStore.readNameStore), " strand=forward");
                    }
                }
            }
			if (options.generateMatePairs) {
                // When generating mate pairs, an even number of reads is generated in each step.
 				SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readNameStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore) % 2, 0u);
			}
		}
    }

    // Last but not least, convert the matches collected before to a global alignment.
    convertMatchesToGlobalAlignment(fragmentStore, Score<int, EditDistance>(), True());
	
	// AlignedReadLayout layout;
	// layoutAlignment(layout, fragmentStore);
	// printAlignment(std::cout, Raw(), layout, fragmentStore, 0, 0, 300, 0, 100);
    
    if (options.verbose)
        std::cerr << "Simulated " << length(fragmentStore.readSeqStore) << " reads" << std::endl;

    return 0;
}

#endif  // MASON_H_
