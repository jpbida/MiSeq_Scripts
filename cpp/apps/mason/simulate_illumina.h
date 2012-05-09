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
// Code specific to the simulation of Illumina reads.
// ==========================================================================

#ifndef SIMULATE_ILLUMINA_H_
#define SIMULATE_ILLUMINA_H_

#include <cmath>

#include <seqan/store.h>

#include "mason.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct IlluminaReads_;
typedef Tag<IlluminaReads_> IlluminaReads;

template<>
struct Options<IlluminaReads> : public Options<Global>
{
    // Length of the reads to simulate.
    unsigned readLength;

    // Base Calling Error Model Parameters.

    // Probability of an insertion.
    double probabilityInsert;
    // Probability of a deletion.
    double probabilityDelete;

    // Set to true if the probability distribution is to be loaded
    // from probabilityMismatchFile.
    bool probabilityMismatchFromFile;
    // Name of the file to load the probabilities from.
    CharString probabilityMismatchFile;
    // Scale factor to apply to mismatch probabilities,
    double probabilityMismatchScale;

    // Probability of a mismatch (single-base polymorphism).
    double probabilityMismatch;
    // Probability for a mismatch in the first base.
    double probabilityMismatchBegin;
    // Probability for a mismatch in the last base.
    double probabilityMismatchEnd;
    // Relative position in the read between 0 and 1 where the steeper curve begins.
    double positionRaise;

    // If set then no Ns will be introduced into the read.
    bool illuminaNoN;

    // Base Calling Quality Model Parameters.

    // Mean quality for non-mismatches at the first base.
    double meanQualityBegin;
    // Mean quality for non-mismatches at last base.
    double meanQualityEnd;
    // Standard deviation quality for non-mismatches at the first base.
    double stdDevQualityBegin;
    // Standard deviation quality for non-mismatches at the last base.
    double stdDevQualityEnd;

    // Mean quality for mismatches at the first base.
    double meanMismatchQualityBegin;
    // Mean quality for mismatches at last base.
    double meanMismatchQualityEnd;
    // Standard deviation quality for mismatches at the first base.
    double stdDevMismatchQualityBegin;
    // Standard deviation quality for mismatches at the last base.
    double stdDevMismatchQualityEnd;

    Options()
            : readLength(36),
              // Base Calling Error Model Parameters
              probabilityInsert(0.001),
              probabilityDelete(0.001),
              probabilityMismatchFromFile(false),
              probabilityMismatchScale(1.0),
              probabilityMismatch(0.004),
              probabilityMismatchBegin(0.002),
              probabilityMismatchEnd(0.012),
              positionRaise(0.66),
              illuminaNoN(false),
              // Base Calling Quality Model Parameters
              meanQualityBegin(40),
              meanQualityEnd(39.5),
              stdDevQualityBegin(0.05),
              stdDevQualityEnd(10),
              meanMismatchQualityBegin(39.5),
              meanMismatchQualityEnd(30),
              stdDevMismatchQualityBegin(3),
              stdDevMismatchQualityEnd(15)
    {}
};

template<>
struct ModelParameters<IlluminaReads> : public ModelParameters<Global>
{
    // Probabilities for a mismatch at a given position.
    String<double> mismatchProbabilities;

    // Standard deviations for the normal distributions of base
    // qualities for the mismatch case.
    String<double> mismatchQualityMeans;
    // Standard deviations for the normal distributions of base
    // qualities for the mismatch case.
    String<double> mismatchQualityStdDevs;

    // Standard deviations for the normal distributions of base
    // qualities for the non-mismatch case.
    String<double> qualityMeans;
    // Standard deviations for the normal distributions of base
    // qualities for the non-mismatch case.
    String<double> qualityStdDevs;
};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<IlluminaReads> const & options) {
    stream << static_cast<Options<Global> >(options);
    stream << "illumina-options {" << std::endl
           << "  readLength:                  " << options.readLength << std::endl
           << "  probabilityInsert:           " << options.probabilityInsert << std::endl
           << "  probabilityDelete:           " << options.probabilityDelete << std::endl
           << "  probabilityMismatchFromFile: " << options.probabilityMismatchFromFile << std::endl
           << "  probabilityMismatchScale:    " << options.probabilityMismatchScale << std::endl
           << "  probabilityMismatchFile:     " << options.probabilityMismatchFile << std::endl
           << "  probabilityMismatch:         " << options.probabilityMismatch << std::endl
           << "  probabilityMismatchBegin:    " << options.probabilityMismatchBegin << std::endl
           << "  probabilityMismatchEnd:      " << options.probabilityMismatchEnd << std::endl
           << "  positionRaise:               " << options.positionRaise << std::endl
           << "  illuminaNoN:                 " << options.illuminaNoN << std::endl
           << "  meanQualityBegin:            " << options.meanQualityBegin << std::endl
           << "  meanQualityEnd:              " << options.meanQualityEnd << std::endl
           << "  stdDevQualityBegin:          " << options.stdDevQualityBegin << std::endl
           << "  stdDevQualityEnd:            " << options.stdDevQualityEnd << std::endl
           << "  meanMismatchQualityBegin:    " << options.meanMismatchQualityBegin << std::endl
           << "  meanMismatchQualityEnd:      " << options.meanMismatchQualityEnd << std::endl
           << "  stdDevMismatchQualityBegin:  " << options.stdDevMismatchQualityBegin << std::endl
           << "  stdDevMismatchQualityEnd:    " << options.stdDevMismatchQualityEnd << std::endl
           << "}" << std::endl;
    return stream;
}

template <>
struct ReadSimulationInstruction<IlluminaReads> : public ReadSimulationInstruction<Global> {};

void setUpCommandLineParser(CommandLineParser & parser,
                            IlluminaReads const &)
{
    setUpCommandLineParser(parser);

    addSection(parser, "Illumina Read Lengths");

    addOption(parser, CommandLineOption("n",  "read-length", "The length of the reads to simulate.  Default: 36.", OptionType::Integer | OptionType::Label));
    addHelpLine(parser, "All resulting reads will have the same length.");

    addSection(parser, "Illumina Error Model");

    addOption(parser, CommandLineOption("pi", "prob-insert", "Probability of an insertion.  Default: 0.001.", OptionType::Double));
    addOption(parser, CommandLineOption("pd", "prob-delete", "Probability of a deletion.  Default: 0.001.", OptionType::Double));
    addOption(parser, CommandLineOption("pmmf", "prob-mismatch-file", "Mismatch probability path.  If set, probability distribution is loaded from argument.  Default: not set.", OptionType::String));
    addOption(parser, CommandLineOption("pmms", "prob-mismatch-scale", "Scale to apply for probability mismatch.  Default: 1.0", OptionType::Double));
    addOption(parser, CommandLineOption("pmm", "prob-mismatch", "Average mismatch probability.  Default: 0.004.", OptionType::Double));
    addOption(parser, CommandLineOption("pmmb", "prob-mismatch-begin", "Probability of a mismatch at the first base.  Default: 0.003.", OptionType::Double));
    addOption(parser, CommandLineOption("pmme", "prob-mismatch-end", "Probability of a mismatch at the last base.  Default: 0.012.", OptionType::Double));
    addOption(parser, CommandLineOption("pr", "position-raise", "Relative position of raise point.  Default: 0.66.", OptionType::Double));
    addOption(parser, CommandLineOption("nN", "no-N", "If set then no Ns will be introduced in the reads.  Default: Ns can be introduced.", OptionType::Bool));

    addOption(parser, CommandLineOption("qmb", "quality-mean-begin", "Quality mean at first base.  Default: 40.", OptionType::Double));
    addOption(parser, CommandLineOption("qme", "quality-mean-end", "Quality mean at last base.  Default: 39.5.", OptionType::Double));
    addOption(parser, CommandLineOption("qsdb", "quality-std-dev-begin", "Quality standard deviation at first base.  Default: 0.05.", OptionType::Double));
    addOption(parser, CommandLineOption("qsde", "quality-std-dev-end", "Quality standard deviation at last base.  Default: 10.", OptionType::Double));

    addOption(parser, CommandLineOption("mmqmb", "mismatch-quality-mean-begin", "Mismatch quality mean at first base.  Default: 39.4.", OptionType::Double));
    addOption(parser, CommandLineOption("mmqme", "mismatch-quality-mean-end", "Mismatch quality mean at last base.  Default: 30.", OptionType::Double));
    addOption(parser, CommandLineOption("mmqsdb", "mismatch-quality-std-dev-begin", "Mismatch quality standard deviation at first base.  Default: 3.", OptionType::Double));
    addOption(parser, CommandLineOption("mmqsde", "mismatch-quality-std-dev-end", "Mismatch quality standard deviation at last base.  Default: 15.", OptionType::Double));
}

int parseCommandLineAndCheckModelSpecific(Options<IlluminaReads> & options,
                                          CommandLineParser & parser)
{
    if (isSetLong(parser, "read-length"))
        getOptionValueLong(parser, "read-length", options.readLength);

    if (isSetLong(parser, "prob-insert"))
        getOptionValueLong(parser, "prob-insert", options.probabilityInsert);
    if (isSetLong(parser, "prob-delete"))
        getOptionValueLong(parser, "prob-delete", options.probabilityDelete);
    if (isSetLong(parser, "prob-mismatch-file")) {
        options.probabilityMismatchFromFile = true;
        getOptionValueLong(parser, "prob-mismatch-file", options.probabilityMismatchFile);
    }
    if (isSetLong(parser, "prob-mismatch-scale"))
        getOptionValueLong(parser, "prob-mismatch-scale", options.probabilityMismatchScale);
    if (isSetLong(parser, "prob-mismatch"))
        getOptionValueLong(parser, "prob-mismatch", options.probabilityMismatch);
    if (isSetLong(parser, "prob-mismatch-begin"))
        getOptionValueLong(parser, "prob-mismatch-begin", options.probabilityMismatchBegin);
    if (isSetLong(parser, "prob-mismatch-end"))
        getOptionValueLong(parser, "prob-mismatch-end", options.probabilityMismatchEnd);
    if (isSetLong(parser, "position-raise"))
        getOptionValueLong(parser, "positio-raise", options.positionRaise);
    if (isSetLong(parser, "no-N"))
        options.illuminaNoN = true;

    if (isSetLong(parser, "quality-mean-begin"))
        getOptionValueLong(parser, "quality-mean-begin", options.meanQualityBegin);
    if (isSetLong(parser, "quality-mean-end"))
        getOptionValueLong(parser, "quality-mean-end", options.meanQualityEnd);
    if (isSetLong(parser, "quality-std-dev-begin"))
        getOptionValueLong(parser, "quality-std-dev-begin", options.stdDevQualityBegin);
    if (isSetLong(parser, "quality-std-dev-end"))
        getOptionValueLong(parser, "quality-std-dev-end", options.stdDevQualityEnd);

    if (isSetLong(parser, "mismatch-quality-mean-begin"))
        getOptionValueLong(parser, "mismatch-quality-mean-begin", options.meanMismatchQualityBegin);
    if (isSetLong(parser, "mismatch-quality-mean-end"))
        getOptionValueLong(parser, "mismatch-quality-mean-end", options.meanMismatchQualityEnd);
    if (isSetLong(parser, "mismatch-quality-std-dev-begin"))
        getOptionValueLong(parser, "mismatch-quality-std-dev-begin", options.stdDevMismatchQualityBegin);
    if (isSetLong(parser, "mismatch-quality-std-dev-end"))
        getOptionValueLong(parser, "mismatch-quality-std-dev-end", options.stdDevMismatchQualityEnd);

    return 0;
}

// Called in simulateReads() to compute the model specific data from the options.
//
// For Illumina reads, this means computing probabilities, means and
// standard deviations for each position, depending on the options.
int simulateReadsSetupModelSpecificData(ModelParameters<IlluminaReads> & parameters,
                                        Options<IlluminaReads> const & options)
{
    // Compute mismatch probabilities, piecewise linear function.
    resize(parameters.mismatchProbabilities, options.readLength);
    // Compute probability at raise point.
    double y_r = 2 * options.probabilityMismatch - options.positionRaise * options.probabilityMismatchBegin - options.probabilityMismatchEnd + options.probabilityMismatchEnd * options.positionRaise;
    if (options.verbose) {
        std::cout << "Illumina error curve:" << std::endl
                  << "  (0, " << options.probabilityMismatchBegin << ") -- (" << options.positionRaise << ", " << y_r << ") -- (1, " << options.probabilityMismatchEnd << ")" << std::endl;
    }
    // std::cout << "y_r = " << y_r << std::endl;
    // Compute mismatch probability at each base.
    if (options.probabilityMismatchFromFile) {
        // Open file.
        std::fstream file;
        file.open(toCString(options.probabilityMismatchFile), std::ios_base::in);
        if (!file.is_open()) {
            std::cerr << "Failed to load mismatch probabilities from " << options.probabilityMismatchFile << std::endl;
            return 1;
        }
        // Load probabilities.
        double x;
        file >> x;
        unsigned i;
        for (i = 0; i < options.readLength && !file.eof(); ++i) {
            parameters.mismatchProbabilities[i] = x;
            file >> x;
        }
        if (i != options.readLength) {
            std::cerr << "Not enough mismatch probabilites in " << options.probabilityMismatchFile << " (" << i << " < " << options.readLength << ")!" << std::endl;
            return 1;
        }
    } else {
        // Use piecewise linear function for mismatch probability simulation.
        for (unsigned i = 0; i < options.readLength; ++i) {
            double x = static_cast<double>(i) / (options.readLength - 1);
            if (x < options.positionRaise) {
                double b = options.probabilityMismatchBegin;
                double m = (y_r - options.probabilityMismatchBegin) / options.positionRaise;
                parameters.mismatchProbabilities[i] = m * x + b;
                // std::cout << "parameters.mismatchProbabilities[" << i << "] = " << parameters.mismatchProbabilities[i] << std::endl;
            } else {
                double b = y_r;
                double m = (options.probabilityMismatchEnd - y_r) / (1 - options.positionRaise);
                x -= options.positionRaise;
                parameters.mismatchProbabilities[i] = m * x + b;
                // std::cout << "parameters.mismatchProbabilities[" << i << "] = " << parameters.mismatchProbabilities[i] << std::endl;
            }
        }
    }
    if (options.probabilityMismatchScale != 1.0) {
        for (unsigned i = 0; i < options.readLength; ++i)
          parameters.mismatchProbabilities[i] *= options.probabilityMismatchScale;
    }

    // Compute match/mismatch means and standard deviations.
    resize(parameters.mismatchQualityMeans, options.readLength);
    for (unsigned i = 0; i < options.readLength; ++i) {
        double b = options.meanMismatchQualityBegin;
        double x = static_cast<double>(i) / (options.readLength - 1);
        double m = (options.meanMismatchQualityEnd - options.meanMismatchQualityBegin);
        parameters.mismatchQualityMeans[i] = m * x + b;
        // std::cout << "parameters.mismatchQualityMeans[" << i << "] = " << parameters.mismatchQualityMeans[i] << std::endl;
    }
    resize(parameters.mismatchQualityStdDevs, options.readLength);
    for (unsigned i = 0; i < options.readLength; ++i) {
        double b = options.stdDevMismatchQualityBegin;
        double x = static_cast<double>(i) / (options.readLength - 1);
        double m = (options.stdDevMismatchQualityEnd - options.stdDevMismatchQualityBegin);
        parameters.mismatchQualityStdDevs[i] = m * x + b;
        // std::cout << "parameters.mismatchQualityStdDevs[" << i << "] = " << parameters.mismatchQualityStdDevs[i] << std::endl;
    }
    resize(parameters.qualityMeans, options.readLength);
    for (unsigned i = 0; i < options.readLength; ++i) {
        double b = options.meanQualityBegin;
        double x = static_cast<double>(i) / (options.readLength - 1);
        double m = (options.meanQualityEnd - options.meanQualityBegin);
        parameters.qualityMeans[i] = m * x + b;
        // std::cout << "parameters.qualityMeans[" << i << "] = " << parameters.qualityMeans[i] << std::endl;
    }
    resize(parameters.qualityStdDevs, options.readLength);
    for (unsigned i = 0; i < options.readLength; ++i) {
        double b = options.stdDevQualityBegin;
        double x = static_cast<double>(i) / (options.readLength - 1);
        double m = (options.stdDevQualityEnd - options.stdDevQualityBegin);
        parameters.qualityStdDevs[i] = m * x + b;
        // std::cout << "parameters.qualityStdDevs[" << i << "] = " << parameters.qualityStdDevs[i] << std::endl;
    }

    return 0;
}

template <typename TRNG>
unsigned pickReadLength(TRNG const &, Options<IlluminaReads> const & options)
{
    return options.readLength;
}

template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<IlluminaReads> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<IlluminaReads> const & parameters, Options<IlluminaReads> const & options) {
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    //
    // Build Edit String.
    //
    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
        double pMismatch = parameters.mismatchProbabilities[i];
        double pInsert   = options.probabilityInsert;
        double pDelete   = options.probabilityDelete;
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

    SEQAN_ASSERT_GT(length(inst.editString), 0u);
    if (options.simulateQualities) {
        clear(inst.qualities);
        resize(inst.qualities, length(inst.editString), 0, Exact());

        for (unsigned i = 0, j = 0; i < length(inst.editString); i++) {
            SEQAN_ASSERT_LEQ(j, inst.endPos - inst.beginPos + inst.delCount);
            if (inst.editString[i] == ERROR_TYPE_MISMATCH) {
                // std::cout << "i == " << i << ", j == " << j << ", parameters.mismatchQualityMeans[j] == " << parameters.mismatchQualityMeans[j] << ", parameters.mismatchQualityStdDevs[j] == " << parameters.mismatchQualityStdDevs[j] << std::endl;
                Pdf<Normal> pdf(parameters.mismatchQualityMeans[j], parameters.mismatchQualityStdDevs[j]);
                inst.qualities[i] = static_cast<int>(pickRandomNumber(rng, pdf));
            } else {
                Pdf<Normal> pdf(parameters.qualityMeans[j], parameters.qualityStdDevs[j]);
                inst.qualities[i] = static_cast<int>(pickRandomNumber(rng, pdf));
            }

            if (inst.editString[i] == ERROR_TYPE_MISMATCH || inst.editString[i] == ERROR_TYPE_MATCH)
                j += 1;
        }
    }
}


template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & rng, ReadSimulationInstruction<IlluminaReads> const & inst, Options<IlluminaReads> const & options)
{
    typedef typename Value<TString>::Type TAlphabet;

    if (options.simulateQualities)
        SEQAN_ASSERT_EQ(length(inst.qualities), length(inst.editString));
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        //int x, xold;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[i]);
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " match" << std::endl;
                //std::cout << back(tmp) << " " << read[j] << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                if (options.illuminaNoN) {
                    c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 3)));  // -3, N not allowed
                } else {
                    c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2)));  // -2, N allowed
                }
                //xold = ordValue(c);
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (ordValue(c) >= ordValue(read[j]))
                    c = TAlphabet(ordValue(c) + 1);
                if (options.illuminaNoN)
                    SEQAN_ASSERT_TRUE(c != TAlphabet('N'));
                //x = ordValue(c);
                appendValue(tmp, c);
                if (options.simulateQualities) {
                    if (options.illuminaNoN)  // Ns can be introduced through quality, too.
                        assignQualityValue(back(tmp), _max(1, inst.qualities[i]));
                    else
                        assignQualityValue(back(tmp), inst.qualities[i]);
                }
                // std::cout << i << " q(q_i)=" << getQualityValue(back(tmp)) << " q(i)=" << inst.qualities[i] << " char=" << convert<char>(back(tmp)) << " c_old=" << xold << " c=" << x << " r_j=" << ordValue(read[j]) << std::endl;
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " mismatch" << std::endl;
                //std::cout << "MM " << c << " " << back(tmp) << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                if (options.illuminaNoN)
                    appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2))));  // -2 == no N
                else
                    appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 1))));  // -1 == N allowed
                if (options.simulateQualities) {
                    if (options.illuminaNoN)  // Ns can be introduced through quality, too.
                        assignQualityValue(back(tmp), _max(1, inst.qualities[i]));
                    else
                        assignQualityValue(back(tmp), inst.qualities[i]);
                }
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " insertion" << std::endl;
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }
    SEQAN_ASSERT_EQ(j, length(read));
    SEQAN_ASSERT_GEQ(length(tmp), options.readLength);

    //std::cout << "tmp == " << tmp << std::endl;
    resize(tmp, options.readLength, Exact());
    move(read, tmp);
}

#endif  // SIMULATE_ILLUMINA_H_
