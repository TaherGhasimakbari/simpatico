#ifndef DDMD_TSTRUCTURE_FACTOR_GRID_CPP
#define DDMD_TSTRUCTURE_FACTOR_GRID_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TStructureFactorGrid.h"
#include <ddMd/simulation/Simulation.h>
#include <util/crystal/PointGroup.h>
#include <util/crystal/PointSymmetry.h>
#include <util/archives/Serializable_includes.h>
#include <util/mpi/MpiLoader.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>
namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   TStructureFactorGrid::TStructureFactorGrid(Simulation& simulation) 
    : StructureFactor(simulation),
      hMax_(0),
      nStar_(0),
      isInitialized_(false)
   {  setClassName("TStructureFactorGrid"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void TStructureFactorGrid::readParameters(std::istream& in) 
   {
      nAtomType_ = simulation().nAtomType();

      // Read parameters
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "hMax", hMax_);

      // Allocate wavevectors arrays
      nWave_     = (2*hMax_ +1 )*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      totalFourierModes_.allocate(nWave_, nMode_);

      int i, j, h, k, l, m;
      IntVector g;

      nStar_ = (2*hMax_ + 1 )*(2*hMax_ + 1)*(2*hMax_ + 1);
      starIds_.allocate(nStar_);
      starSizes_.allocate(nStar_);
      // Create tetragonal point group
      PointGroup group;
      PointSymmetry a, b, c;

      a.R(0,0) =  1;
      a.R(1,1) =  1;
      a.R(2,2) =  1;

      b.R(0,0) =  1;
      b.R(1,1) =  1;
      b.R(2,2) =  1;

      c.R(0,0) =  1;
      c.R(1,1) =  1;
      c.R(2,2) =  1;

      group.add(c);
      group.add(b);
      group.add(a);
      group.makeCompleteGroup();

      // Create grid of wavevectors
      FSArray<IntVector, 16> star;
      i = 0;
      j = 0;
      for (h = -hMax_; h <= hMax_; ++h) {
         g[0] = h;
         for (k = -hMax_; k <= hMax_; ++k) {
            g[1] = k;
            for (l = -hMax_; l <= hMax_; ++l) {
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

      if (simulation().domain().isMaster()) {

         structureFactors_.allocate(nWave_, nMode_);
         int i, j;
         for (i = 0; i < nWave_; ++i) {
            for (j = 0; j < nMode_; ++j) {
               structureFactors_(i, j) = 0.0;
            }
         }

      }
      nSample_ = 0;

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void TStructureFactorGrid::loadParameters(Serializable::IArchive &ar)
   {
      nAtomType_ = simulation().nAtomType();

      // Load parameter file parameters
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      loadDMatrix<double>(ar, "modes", modes_, nMode_, nAtomType_);
      loadParameter<int>(ar, "hMax", hMax_);

      // Load and broadcast other distributed members
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nWave_);
      waveIntVectors_.allocate(nWave_);
      loader.load(waveIntVectors_);
      loader.load(nStar_);
      starIds_.allocate(nStar_);
      loader.load(starIds_);
      starSizes_.allocate(nStar_);
      loader.load(starSizes_);
      loader.load(nSample_);

      if (simulation().domain().isMaster()) {
         structureFactors_.allocate(nWave_, nMode_);
         ar >> structureFactors_;
      }

      // Allocate work space
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      totalFourierModes_.allocate(nWave_, nMode_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void TStructureFactorGrid::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nMode_;
      ar << modes_;
      ar << hMax_;

      ar << nWave_;
      ar << waveIntVectors_;
      ar << nStar_;
      ar << starIds_;
      ar << starSizes_;
      ar << nSample_;

      ar << structureFactors_;
   }

   void TStructureFactorGrid::clear()
   {}

   void TStructureFactorGrid::sample(long iStep)
   {
      if (isAtInterval(iStep))  {

         if (simulation().domain().isMaster()) {
            simulation().fileMaster().openOutputFile(outputFileName(".dat"), logFile_, !isFirstStep_);
         }

         StructureFactor::sample(iStep);
         
         // Log structure factors
         if (simulation().domain().isMaster()) {
            double volume = simulation().boundary().volume();
            double norm;
            double arg;
            for (int i = 0; i < nStar_; ++i) {
               int size = starSizes_[i];

               int k = starIds_[i];
               for (int n = 0; n < Dimension; ++n) {
                  logFile_ << Int(waveIntVectors_[k][n], 5);
               }
               logFile_ << Dbl(waveVectors_[k].abs(), 20, 8);

               for (int j = 0; j < nMode_; ++j) {
                  double average_norm = 0.0;
                  double average_arg = 0.0;
                  double value = 0.0;
                  k = starIds_[i];
                  for (int m = 0; m < size; ++m) {
                     norm = std::norm(totalFourierModes_(k, j));
                     arg = std::arg(totalFourierModes_(k, j));
                     value = norm/volume;
                     average_norm += value;
                     average_arg += arg;
                     ++k;
                  }
                  average_norm = average_norm/double(size);
                  average_arg = average_arg/double(size);
                  logFile_ << Dbl(average_norm, 20, 8);
                  logFile_ << Dbl(average_arg, 20, 8);
               }
               logFile_ << std::endl;
            }
            logFile_ << std::endl;
            logFile_.close();
         }
      }
   }

   void TStructureFactorGrid::output()
   {
      if (simulation().domain().isMaster()) {
            
         double  value, average, size;
         int     i, j, k, m, n;

         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_.close();

         // Output all structure factors to one file
         simulation().fileMaster().openOutputFile(outputFileName("_avg.dat"), outputFile_);
         for (i = 0; i < nStar_; ++i) {
            size = starSizes_[i];
            k = starIds_[i];
            for (n = 0; n < Dimension; ++n) {
               outputFile_ << Int(waveIntVectors_[k][n], 5);
            }
            outputFile_ << Dbl(waveVectors_[k].abs(), 20, 8);
            for (j = 0; j < nMode_; ++j) {
               k = starIds_[i];
               average = 0.0;
               for (m = 0; m < size; ++m) {
                  value = structureFactors_(k, j)/double(nSample_);
                  average += value;
                  ++k;
               }
               average = average/double(size);
               outputFile_ << Dbl(average, 20, 8);
            }
            outputFile_ << std::endl;
         }
         outputFile_.close();
      
      }
   }

}
#endif
