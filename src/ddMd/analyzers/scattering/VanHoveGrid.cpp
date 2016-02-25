/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VanHoveGrid.h"
#include <ddMd/simulation/Simulation.h>
#include <util/crystal/PointGroup.h>
#include <util/crystal/PointSymmetry.h>
#include <ddMd/simulation/SimulationAccess.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <util/boundary/Boundary.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <util/space/IntVector.h>
#include <util/misc/ioUtil.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   VanHoveGrid::VanHoveGrid(Simulation& simulation) 
    : StructureFactor(simulation),
      isInitialized_(false)
   {  setClassName("VanHoveGrid"); }

   /*
   * Destructor
   */
   VanHoveGrid::~VanHoveGrid() 
   {}

   /*
   * Read parameters from file, and allocate memory.
   */
   void VanHoveGrid::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      nAtomType_ = simulation().nAtomType();
      atomTypeCoeffs_.allocate(nAtomType_);
      readDArray<double>(in, "atomTypeCoeffs", atomTypeCoeffs_, nAtomType_);

      read<int>(in, "nBuffer", nBuffer_);
      read<int>(in, "hMax", hMax_);
      read<LatticeSystem>(in, "lattice", lattice_);

      // Allocate wavevectors arrays
      nWave_ = (2*hMax_ +1 )*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      totalFourierModes_.allocate(nWave_);

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
      } else if (lattice_ == Tetragonal) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 1)*(hMax_ + 2)/2;
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
         // Create tetragonal point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,0) =  1;
         a.R(1,2) =  1;
         a.R(2,1) =  1;

         b.R(0,0) =  -1;
         b.R(1,1) =  1;
         b.R(2,2) =  1;

         c.R(0,0) =  1;
         c.R(1,1) =  -1;
         c.R(2,2) =  1;

         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();

         // Create grid of wavevectors
         FSArray<IntVector, 16> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= hMax_; ++k) {
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
      } else if (lattice_ == Orthorhombic) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 1)*(hMax_ + 1);
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
         // Create tetragonal point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,0) =  -1;
         a.R(1,1) =  1;
         a.R(2,2) =  1;

         b.R(0,0) =  1;
         b.R(1,1) =  -1;
         b.R(2,2) =  1;

         c.R(0,0) =  1;
         c.R(1,1) =  1;
         c.R(2,2) =  -1;

         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();

         // Create grid of wavevectors
         FSArray<IntVector, 16> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= hMax_; ++k) {
               g[1] = k;
               for (l = 0; l <= hMax_; ++l) {
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

      if (simulation().domain().isMaster()) {
         accumulators_.allocate(nWave_);
         for (int i = 0; i < nWave_; ++i) {
            accumulators_[i].setParam(nBuffer_);
         }
      }

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void VanHoveGrid::loadParameters(Serializable::IArchive& ar)
   {  
      nAtomType_ = simulation().nAtomType();
      loadInterval(ar);
      loadOutputFileName(ar);
      atomTypeCoeffs_.allocate(nAtomType_);
      loadDArray<double>(ar, "atomTypeCoeffs", atomTypeCoeffs_, nAtomType_);
      loadParameter<int>(ar, "nBuffer", nBuffer_);
      loadParameter<int>(ar, "hMax", hMax_);
      loadParameter<LatticeSystem>(ar, "lattice", lattice_);
      waveIntVectors_.allocate(nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);

      // Load and broadcast nSample_
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nWave_);
      waveIntVectors_.allocate(nWave_);
      loader.load(waveIntVectors_, nWave_);
      loader.load(nStar_);
      starIds_.allocate(nStar_);
      loader.load(starIds_, nStar_);
      starSizes_.allocate(nStar_);
      loader.load(starSizes_, nStar_);
      
      if (simulation().domain().isMaster()) {
         accumulators_.allocate(nWave_);
         ar >> accumulators_;
      }
             
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_);
      totalFourierModes_.allocate(nWave_);

      isInitialized_ = true;
   }

   /*
   * Save state to an archive.
   */
   void VanHoveGrid::save(Serializable::OArchive& ar)
   { 
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << atomTypeCoeffs_;
      ar << nBuffer_;
      ar << hMax_;
      ar << lattice_;
      ar << nWave_;
      ar << waveIntVectors_;
      ar << nStar_;
      ar << starIds_;
      ar << starSizes_;

      ar << nSample_;

      if (simulation().domain().isMaster()) {
         ar << accumulators_;
      }
   }

   /*
   * Clear accumulators.
   */
   void VanHoveGrid::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      assert (nBuffer_ > 0);
      assert (nWave_ > 0);
      nSample_ = 0;
      if (simulation().domain().isMaster()) {
         for (int i = 0; i < nWave_; ++i) {
            accumulators_[i].clear();
         }
      }
   }
 
   /// Increment Structure Factor
   void VanHoveGrid::sample(long iStep) 
   {

      if (!isAtInterval(iStep))  {
         UTIL_THROW("Time step index is not a multiple of interval");
      }

      if (isAtInterval(iStep))  {

         StructureFactor::sample(iStep);
         // Log structure factors
         if (simulation().domain().isMaster()) {
            double volume = simulation().boundary().volume();
            double norm;
            for (int i = 0; i < nStar_; ++i) {
               int size = starSizes_[i];

               int k = starIds_[i];

               double average = 0.0;
               double value = 0.0;
               k = starIds_[i];
               for (int m = 0; m < size; ++m) {
                  norm = std::norm(totalFourierModes_[k]);
                  value = norm/volume;
                  average += value;
                  ++k;
               }
               average = average/double(size);
            }
         }

         
         if (simulation().domain().isMaster()) {
            // Add Fourier modes to autocorrelation accumulators
            double sqrtV = sqrt(simulation().boundary().volume());
            for (int i = 0; i < nWave_; ++i) {
               accumulators_[i].sample(totalFourierModes_[i]/sqrtV);
            }
         }
         ++nSample_;
      }

   }

   /**
   * Calculate floating point wavevectors.
   */
   void VanHoveGrid::makeWaveVectors() 
   {
      Vector    dWave;
      Boundary* boundaryPtr = &simulation().boundary();
      int       i, j;

      // Calculate wavevectors
      for (i = 0; i < nWave_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (j = 0; j < Dimension; ++j) {
            dWave  = boundaryPtr->reciprocalBasisVector(j);
            dWave *= waveIntVectors_[i][j];
            waveVectors_[i] += dWave;
         }
      }
   }

   void VanHoveGrid::output() 
   {
      if (simulation().domain().isMaster()) {
         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_.close();
        
         // Output structure factors to one *.dat file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
        
         int i, j;

         // Output autocorrelation functions to file
         for (i = 0; i < nWave_; ++i) {

            for (j = 0; j < Dimension; ++j) {
               outputFile_ << Int(waveIntVectors_[i][j], 5);
            }
            outputFile_ << Dbl(waveVectors_[i].abs(), 20, 8);
            outputFile_ << std::endl;
            accumulators_[i].output(outputFile_);
            outputFile_ << std::endl;
            outputFile_ << std::endl;
         }
         outputFile_.close();
      }

   }

}
