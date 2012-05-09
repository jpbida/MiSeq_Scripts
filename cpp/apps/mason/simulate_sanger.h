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
// Code specific to the Sanger read simulation.
// ==========================================================================

#ifndef SIMULATE_SANGER_H_
#define SIMULATE_SANGER_H_

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct SangerReads_;
typedef Tag<SangerReads_> SangerReads;

template <>
struct Options<SangerReads> : public Options<Global>
{
    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthIsUniform;

    // Average read length.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation,
    // for uniform read length the interval around the average to use for
    // picking the read lengths.
    double readLengthError;

    // Base Calling Error Model Parameters.

    // Mismatch probability ramp.
    double probabilityMismatchBegin;
    double probabilityMismatchEnd;

    // Insert probability ramp.
    double probabilityInsertBegin;
    double probabilityInsertEnd;

    // Delete probability ramp.
    double probabilityDeleteBegin;
    double probabilityDeleteEnd;

    // Quality mean and standard deviation ramp specification for matches.
    double qualityMatchStartMean;
    double qualityMatchEndMean;
    double qualityMatchStartStdDev;
    double qualityMatchEndStdDev;

    // Quality mean and standard deviation ramp specification for errors.
    double qualityErrorStartMean;
    double qualityErrorEndMean;
    double qualityErrorStartStdDev;
    double qualityErrorEndStdDev;

    Options()
            : readLengthIsUniform(false),
              readLengthMean(400),
              readLengthError(40),
              probabilityMismatchBegin(0.005),
              probabilityMismatchEnd(0.01),
              probabilityInsertBegin(0.0025),
              probabilityInsertEnd(0.005),
              probabilityDeleteBegin(0.0025),
              probabilityDeleteEnd(0.005),
              qualityMatchStartMean(40),
              qualityMatchEndMean(39),
              qualityMatchStartStdDev(0.1),
              qualityMatchEndStdDev(2),
              qualityErrorStartMean(30),
              qualityErrorEndMean(20),
              qualityErrorStartStdDev(2),
              qualityErrorEndStdDev(5)
    {}
};

template <>
struct ReadSimulationInstruction<SangerReads> : public ReadSimulationInstruction<Global> {
};

template<>
struct ModelParameters<SangerReads> : public ModelParameters<Global>
{
};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<SangerReads> const & options) {
    stream << static_cast<Options<Global> >(options);
    stream << "sanger-options {" << std::endl
           << "  readLengthIsUniform:      " << options.readLengthIsUniform << std::endl
           << "  readLengthMean:           " << options.readLengthMean << std::endl
           << "  readLengthError:          " << options.readLengthError << std::endl
           << "  probabilityMismatchBegin: " << options.probabilityMismatchBegin << std::endl
           << "  probabilityMismatchEnd:   " << options.probabilityMismatchEnd << std::endl
           << "  probabilityInsertBegin:   " << options.probabilityInsertBegin << std::endl
           << "  probabilityInsertEnd:     " << options.probabilityInsertEnd << std::endl
           << "  probabilityDeleteBegin:   " << options.probabilityInsertEnd << std::endl
           << "  probabilityDeleteEnd:     " << options.probabilityDeleteEnd << std::endl
           << "  qualityMatchStartMean:    " << options.qualityMatchStartMean << std::endl
           << "  qualityMatchEndMean:      " << options.qualityMatchEndMean << std::endl
           << "  qualityMatchStartStdDev:  " << options.qualityMatchStartStdDev << std::endl
           << "  qualityMatchEndStdDev:    " << options.qualityMatchEndStdDev << std::endl
           << "  qualityErrorStartMean:    " << options.qualityErrorStartMean << std::endl
           << "  qualityErrorEndMean:      " << options.qualityErrorEndMean << std::endl
           << "  qualityErrorStartStdDev:  " << options.qualityErrorStartStdDev << std::endl
           << "  qualityErrorEndStdDev:    " << options.qualityErrorEndStdDev << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser,
                            SangerReads const &)
{
    setUpCommandLineParser(parser);

    addSection(parser, "Sanger Read Length Parameters");

    addOption(parser, CommandLineOption("nu",  "read-length-uniform", "If set, the read lengths are simulated with a uniform distribution, with standard distribution otherwise.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("nm",  "read-length-mean", "The mean of the read lengths.  Default: 400.", OptionType::Double));
    addOption(parser, CommandLineOption("ne",  "read-length-error", "The standard deviation (for standard distribution) and interval length (for uniform distribution) for the read length.  Default: 40.", OptionType::Double));

    addSection(parser, "Sanger Error Model Parameters");

    addOption(parser, CommandLineOption("pmb",  "probability-mismatch-begin", "Probability for a mismatch at begin of read.  Default: 0.005.", OptionType::Double));
    addOption(parser, CommandLineOption("pme",  "probability-mismatch-begin", "Probability for a mismatch at end of read.  Default: 0.01.", OptionType::Double));
    addOption(parser, CommandLineOption("pib",  "probability-insert-begin", "Probability for a insert at begin of read.  Default: 0.0025.", OptionType::Double));
    addOption(parser, CommandLineOption("pie",  "probability-insert-begin", "Probability for a insert at end of read.  Default: 0.005.", OptionType::Double));
    addOption(parser, CommandLineOption("pdb",  "probability-delete-begin", "Probability for a delete at begin of read.  Default: 0.0025.", OptionType::Double));
    addOption(parser, CommandLineOption("pde",  "probability-delete-begin", "Probability for a delete at end of read.  Default: 0.005.", OptionType::Double));

    addSection(parser, "Base Quality Model Parameters");
    addOption(parser, CommandLineOption("qmsm",  "quality-match-start-mean", "Mean quality for matches at start of read.  Default: 40.", OptionType::Double));
    addOption(parser, CommandLineOption("qmem",  "quality-match-end-mean", "Mean quality for matches at end of read.  Default: 39.", OptionType::Double));
    addOption(parser, CommandLineOption("qmss",  "quality-match-start-std-dev", "Quality standard deviation for matches at start of read.  Default: 0.1.", OptionType::Double));
    addOption(parser, CommandLineOption("qmes",  "quality-match-end-std-dev", "Quality standard deviation for matches at end of read.  Default: 2.", OptionType::Double));
    addOption(parser, CommandLineOption("qesm",  "quality-error-start-mean", "Mean quality for errors at start of read.  Default: 40.", OptionType::Double));
    addOption(parser, CommandLineOption("qeem",  "quality-error-end-mean", "Mean quality for errors at end of read.  Default: 39.", OptionType::Double));
    addOption(parser, CommandLineOption("qess",  "quality-error-start-std-dev", "Quality standard deviation for errors at start of read.  Default: 0.1.", OptionType::Double));
    addOption(parser, CommandLineOption("qees",  "quality-error-end-std-dev", "Quality standard deviation for errors at end of read.  Default: 2.", OptionType::Double));
}

int parseCommandLineAndCheckModelSpecific(Options<SangerReads> & options,
                                          CommandLineParser & parser)
{
    if (isSetLong(parser, "read-length-uniform"))
        options.readLengthIsUniform = true;
    if (isSetLong(parser, "read-length-mean"))
        getOptionValueLong(parser, "read-length-mean", options.readLengthMean);
    if (isSetLong(parser, "read-length-error"))
        getOptionValueLong(parser, "read-length-error", options.readLengthError);

    if (isSetLong(parser, "probability-mismatch-begin"))
        getOptionValueLong(parser, "probability-mismatch-begin", options.probabilityMismatchBegin);
    if (isSetLong(parser, "probability-mismatch-end"))
        getOptionValueLong(parser, "probability-mismatch-end", options.probabilityMismatchEnd);
    if (isSetLong(parser, "probability-insert-begin"))
        getOptionValueLong(parser, "probability-insert-begin", options.probabilityInsertBegin);
    if (isSetLong(parser, "probability-insert-end"))
        getOptionValueLong(parser, "probability-insert-end", options.probabilityInsertEnd);
    if (isSetLong(parser, "probability-delete-begin"))
        getOptionValueLong(parser, "probability-delete-begin", options.probabilityDeleteBegin);
    if (isSetLong(parser, "probability-delete-end"))
        getOptionValueLong(parser, "probability-delete-end", options.probabilityDeleteEnd);

    if (isSetLong(parser, "quality-match-start-mean"))
        getOptionValueLong(parser, "quality-match-start-mean", options.qualityMatchStartMean);
    if (isSetLong(parser, "quality-match-end-mean"))
        getOptionValueLong(parser, "quality-match-end-mean", options.qualityMatchEndMean);
    if (isSetLong(parser, "quality-match-start-std-dev"))
        getOptionValueLong(parser, "quality-match-start-std-dev", options.qualityMatchStartStdDev);
    if (isSetLong(parser, "quality-match-end-std-dev"))
        getOptionValueLong(parser, "quality-match-end-std-dev", options.qualityMatchEndStdDev);
    if (isSetLong(parser, "quality-Error-start-mean"))
        getOptionValueLong(parser, "quality-Error-start-mean", options.qualityErrorStartMean);
    if (isSetLong(parser, "quality-Error-end-mean"))
        getOptionValueLong(parser, "quality-Error-end-mean", options.qualityErrorEndMean);
    if (isSetLong(parser, "quality-Error-start-std-dev"))
        getOptionValueLong(parser, "quality-Error-start-std-dev", options.qualityErrorStartStdDev);
    if (isSetLong(parser, "quality-Error-end-std-dev"))
        getOptionValueLong(parser, "quality-Error-end-std-dev", options.qualityErrorEndStdDev);

    return 0;
}

// No model specific data for Sanger reads.
int simulateReadsSetupModelSpecificData(ModelParameters<SangerReads> & /*parameters*/,
                                        Options<SangerReads> const & /*options*/)
{
    return 0;
}

// TODO(holtgrew): Same as 454 reads! Remove redundancy.
template <typename TRNG>
inline
unsigned pickReadLength(TRNG & rng, Options<SangerReads> const & options)
{
    if (options.readLengthIsUniform) {
        // Pick uniformly.
        double minLen = options.readLengthMean - options.readLengthError;
        double maxLen = options.readLengthMean + options.readLengthError;
        double len = pickRandomNumber(rng, Pdf<Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(round(len));
    } else {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, Pdf<Normal>(options.readLengthMean, options.readLengthError));
        return static_cast<unsigned>(round(len));
    }
}


template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<SangerReads> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<SangerReads> const & /*parameters*/, Options<SangerReads> const & options)
{
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;
    if (options.simulateQualities) {
        reserve(inst.qualities, readLength, Generous());
        clear(inst.qualities);
    }

    //
    // Build Edit String.
    //
    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
        double pos = 1.0 * i / (readLength - 1);
        double pMismatch = options.probabilityMismatchBegin + pos * (options.probabilityMismatchEnd - options.probabilityMismatchBegin);
        double pInsert   = options.probabilityInsertBegin + pos * (options.probabilityInsertEnd - options.probabilityInsertBegin);
        double pDelete   = options.probabilityDeleteBegin + pos * (options.probabilityDeleteEnd - options.probabilityDeleteBegin);
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;
        if (x < pMatch) {
            // match
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        } else if (x < pMatch + pMismatch) {
            // mismatch
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
        } else if (x < pMatch + pMismatch + pInsert) {
            // insert
            if (length(inst.editString) > 0 && back(inst.editString == ERROR_TYPE_DELETE)) {
                inst.delCount -= 1;
                eraseBack(inst.editString);
            } else {
                i += 1;
                inst.insCount += 1;
                appendValue(inst.editString, ERROR_TYPE_INSERT);
            }
        } else {
            // Decrement string size, do not add a delete if string is
            // too short, possibly remove insert from edit string.
            if (length(inst.editString) > 0) {
                if (back(inst.editString == ERROR_TYPE_INSERT)) {
                    i -= 1;
                    inst.insCount -= 1;
                    eraseBack(inst.editString);
                } else {
                    inst.delCount += 1;
                    appendValue(inst.editString, ERROR_TYPE_DELETE);
                }
            }
        }
    }
    SEQAN_ASSERT_EQ(readLength, length(inst.editString) - inst.delCount);

    //
    // Adjust Positions.
    //

    // If the number of deletions does not equal the number of inserts
    // then we have to adjust the read positions.
    if (inst.delCount != inst.insCount) {
        int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
        inst.endPos += delta;
        if (inst.endPos > length(contig)) {
            delta = inst.endPos - length(contig);
            inst.endPos -= delta;
            inst.beginPos -= delta;
        }
        SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                        readLength);
    }

    //
    // Quality Simulation.
    //
    if (options.simulateQualities) {
        reserve(inst.qualities, readLength + inst.insCount - inst.delCount, Exact());
        clear(inst.qualities);

        for (unsigned i = 0; i < length(inst.editString); ++i) {
            double mean, stdDev;
            double pos = 1.0 * i / (readLength + inst.insCount - inst.delCount - 1);
            if (inst.editString[i] ==  ERROR_TYPE_DELETE) {
                continue;  // No quality to give out.
            } else if (inst.editString[i] ==  ERROR_TYPE_INSERT || inst.editString[i] ==  ERROR_TYPE_MISMATCH) {
                mean = options.qualityMatchStartMean + pos * (options.qualityMatchEndMean - options.qualityMatchStartMean);
                stdDev = options.qualityMatchStartStdDev + pos * (options.qualityMatchEndStdDev - options.qualityMatchStartStdDev);
            } else {
                mean = options.qualityErrorStartMean + pos * (options.qualityErrorEndMean - options.qualityErrorStartMean);
                stdDev = options.qualityErrorStartStdDev + pos * (options.qualityErrorEndStdDev - options.qualityErrorStartStdDev);
            }
            Pdf<Normal> pdf(mean, stdDev);
            appendValue(inst.qualities, static_cast<int>(pickRandomNumber(rng, pdf)));
        }
    }
}

template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & rng, ReadSimulationInstruction<SangerReads> const & inst, Options<SangerReads> const & options)
{
    typedef typename Value<TString>::Type TAlphabet;

    if (options.simulateQualities)
        SEQAN_ASSERT_EQ(length(inst.qualities) + inst.delCount, length(inst.editString));
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;  // Index in string.
    unsigned k = 0;  // Index in qualities.
    for (unsigned i = 0; i < length(inst.editString); ++i) {  // Index in edit string.
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[k]);
                k += 1;
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2)));  // -2, N allowed
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (ordValue(c) >= ordValue(read[j]))
                    c = TAlphabet(ordValue(c) + 1);
                appendValue(tmp, c);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[k]);
                k += 1;
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 1))));  // -1 == N allowed
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[k]);
                k += 1;
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }

    move(read, tmp);
}

#endif  // SIMULATE_SANGER_H_
