#ifndef MCMD_MAXIMUM_SF_CPP
#define MCMD_MAXIMUM_SF_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MaximumSf.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/crystal/PointGroup.h>
#include <util/crystal/PointSymmetry.h>
#include <util/archives/Serializable_includes.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   MaximumSf::MaximumSf(System& system) 
    : StructureFactor(system),
      hMax_(0),
      nStar_(0),
      lattice_(Triclinic),
      isInitialized_(false)
   {}

   /// Read parameters from file, and allocate data array.
   void MaximumSf::readParam(std::istream& in) 
   {
      nAtomType_ = system().simulation().nAtomType();

      // Read parameters
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "hMax", hMax_);
      read<Util::LatticeSystem>(in, "lattice", lattice_);

      // Allocate wavevectors arrays
      nWave_     = (2*hMax_ +1 )*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      structureFactors_.allocate(nWave_, nMode_);

      int i, j, h, k, l, m;
      IntVector g;

      // Cubic Symmetry
      if (lattice_ == Cubic) {

         nStar_ = (hMax_ +1 )*(hMax_ + 2)*(hMax_ + 3)/6;
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
                  
         // Create cubic point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,1) =  1;
         a.R(1,0) =  1;
         a.R(2,2) =  1;

         b.R(0,0) = -1;
         b.R(1,1) =  1;
         b.R(2,2) =  1;

         c.R(0,1) =  1;
         c.R(1,2) =  1;
         c.R(2,0) =  1;

         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();

         // Create grid of wavevectors
         FSArray<IntVector, 48> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= h; ++k) {
               g[1] = k;
               for (l = 0; l <= k; ++l) {
                  g[2] = l;
                  starIds_[i] = j;
                  group.makeStar(g, star);
                  starSizes_[i] = star.size();
                  for (m = 0; m < star.size(); ++m) {
                     waveIntVectors_[j] = star[m];
                     ++j;
                  }
                  ++i;
               }
            }
         }
         if (i != nStar_) {
            UTIL_THROW("Error");
         }
         if (j != nWave_) {
            UTIL_THROW("Error");
         }
      }

      // Clear accumulators
      for (i = 0; i < nWave_; ++i) {
         for (j = 0; j < nMode_; ++j) {
            structureFactors_(i, j) = 0.0;
         }
      }

      maximumValue_.allocate(Samples);
      maximumWaveIntVector_.allocate(Samples);
      maximumQ_.allocate(Samples);

      for (i=0; i < Samples; ++i) {
         maximumValue_[i] = 0.0;
      }
      nSample_ = 0;

      isInitialized_ = true;
   }

   void MaximumSf::setup() 
   {}

   void MaximumSf::output() 
   {
      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Outputs history of maximum structure factors
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      for (int i = 0; i < nSample_; ++i) {
         outputFile_ << maximumWaveIntVector_[i];
         outputFile_ << Dbl(maximumQ_[i], 20, 8);
         outputFile_ << Dbl(maximumValue_[i], 20, 8);
         outputFile_ << std::endl;
      }
      outputFile_.close();
   }
   
   /*
   * Save state to binary file archive.
   */
   void MaximumSf::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void MaximumSf::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}
#endif
