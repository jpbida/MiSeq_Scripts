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
// FragmentStore configuration used in the read simulator.
// ==========================================================================

#ifndef STORE_CONFIG_H_
#define STORE_CONFIG_H_

// Fragment Store Configuration.
//
// We change the default fragment store configuration to use normal
// string sets and not concat string sets for the reads since we need
// to be able to swap read sequences.

namespace seqan {
struct MyFragmentStoreConfig {};

template<>
struct FragmentStoreConfig<MyFragmentStoreConfig>
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void    TReadStoreElementSpec;
	typedef Owner<Default> TReadSeqStoreSpec;
	typedef void    TMatePairStoreElementSpec;
	typedef void    TLibraryStoreElementSpec;
	typedef void    TContigStoreElementSpec;
	typedef void    TContigFileSpec;
	typedef void    TAlignedReadStoreElementSpec;
	typedef Owner<Default>	TAlignedReadTagStoreSpec;
	typedef void    TAnnotationStoreElementSpec;
};
}  // namespace seqan

#endif  // #ifndef STORE_CONFIG_H_
