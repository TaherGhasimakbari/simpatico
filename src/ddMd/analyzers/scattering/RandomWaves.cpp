/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RandomWaves.h"
#include <ddMd/simulation/Simulation.h>
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
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace DdMd
{

   using namespace Util;

   /// Constructor.
   RandomWaves::RandomWaves(Simulation& simulation) 
    : Analyzer(simulation),
      isInitialized_(false),
      isFirstStep_(true)
   {  setClassName("RandomWaves"); }

   RandomWaves::~RandomWaves() 
   {}

   /// Read parameters from file, and allocate data array.
   void RandomWaves::readParameters(std::istream& in) 
   {
      nAtomType_ = simulation().nAtomType();

      readInterval(in);
      readOutputFileName(in);
      mode_.allocate(nAtomType_);
      readDArray<double>(in, "mode", mode_, nAtomType_);
      read<int>(in, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWave_);
      phases_.allocate(nWave_);
      readDArray<double>(in, "phases", phases_, nWave_);
      grid_.allocate(Dimension);
      readDArray<int>(in, "grid", grid_, 3);

      waveVectors_.allocate(nWave_);
      randomWaves_.allocate(nWave_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void RandomWaves::loadParameters(Serializable::IArchive &ar)
   {
      nAtomType_ = simulation().nAtomType();

      // Load and broadcast parameter file parameters
      loadInterval(ar);
      loadOutputFileName(ar);
      mode_.allocate(nAtomType_);
      loadDArray<double>(ar, "mode", mode_, nAtomType_);
      loadParameter<int>(ar, "nWave", nWave_);
      waveIntVectors_.allocate(nWave_);
      loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);
      phases_.allocate(nWave_);
      loadDArray<double>(ar, "phases", phases_, nWave_);
      grid_.allocate(Dimension);
      loadDArray<int>(ar, "grid", grid_, 3);

      // Load and broadcast nSample_
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      randomWaves_.allocate(nWave_);
      ar >> randomWaves_;

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void RandomWaves::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << mode_;
      ar << nWave_;
      ar << waveIntVectors_;
      ar << phases_;
      ar << grid_;
      ar << nSample_;
      ar << accumulatorU_;
      ar << accumulatorN_;
      ar << accumulatorUavg_;
      ar << accumulatorNavg_;
   }
  
   /*
   * Clear accumulators.
   */
   void RandomWaves::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      assert (nWave_ > 0);
      accumulatorU_.clear();
      accumulatorN_.clear();
      accumulatorUavg_.clear();
      accumulatorNavg_.clear();

      nSample_ = 0;
   }

   /*
   * Increment structure factor.
   */
   void RandomWaves::sample(long iStep) 
   {
      if (!isAtInterval(iStep))  {
         UTIL_THROW("Time step index not a multiple of interval");
      }
      if (simulation().domain().isMaster()) {
         std::ios_base::openmode mode = std::ios_base::out;
         if (!isFirstStep_) {
            mode = std::ios_base::out | std::ios_base::app;
         }
         std::string name = outputFileName(".dat");
         simulation().fileMaster().openOutputFile(name, outputFile_, mode);
      }

      makeWaveVectors();

      for (int i = 0; i < nWave_; ++i) {
         randomWaves_[i] = std::complex<double>(gaussianDistribution_(randomGenerator_), gaussianDistribution_(randomGenerator_));
      }
 
      Boundary* boundaryPtr = &simulation().boundary();
      if(simulation().domain().isMaster()) {

         std::complex<double> CV;
         double CVmax = 0.0;
         double CVavg = 0.0;
         double CVreal = 0.0; 
         double abs = 0.0;

         int imax = 0;
         int jmax = 0;
         int kmax = 0;

         Vector R;
         R.zero();
         Vector dR0;
         Vector dR1;
         Vector dR2;

         for (int w = 0; w < nWave_; ++w) {
            abs += std::abs(randomWaves_[w]);
         }

         for (int i = 0; i < grid_[0]; ++i) {
            dR0 = boundaryPtr->bravaisBasisVector(0);
            dR0 *= i;
            dR0 /= grid_[0];
            for (int j = 0; j < grid_[1]; ++j) {
               dR1 = boundaryPtr->bravaisBasisVector(1);
               dR1 *= j;
               dR1 /= grid_[1];
               for (int k = 0; k < grid_[2]; ++k) {
                  dR2 = boundaryPtr->bravaisBasisVector(2);
                  dR2 *= k;
                  dR2 /= grid_[2];

                  R.zero();
                  R += dR0;
                  R += dR1;
                  R += dR2;
                  CV = 0.0;

                  for (int w = 0; w < nWave_; ++w) {
                     CV += exp(-Constants::Im*phases_[w])*exp(-Constants::Im*R.dot(waveVectors_[w]))*randomWaves_[w];
                  }
                  CVavg = CVavg + std::real(CV);

                  CVreal = std::real(CV);
                  if (CVmax < CVreal) {
                     imax = i;
                     jmax = j;
                     kmax = k;
                     CVmax = CVreal;
                  }
               }
            }
         }
         CVavg /= grid_[0]*grid_[1]*grid_[2];

         Vector Rmax;
         Rmax.zero();
         dR0 = boundaryPtr->bravaisBasisVector(0);
         dR0 *= imax;
         dR0 /= grid_[0];
         Rmax += dR0;
         dR1 = boundaryPtr->bravaisBasisVector(1);
         dR1 *= jmax;
         dR1 /= grid_[1];
         Rmax += dR1;
         dR2 = boundaryPtr->bravaisBasisVector(2);
         dR2 *= kmax;
         dR2 /= grid_[2];
         Rmax += dR2;
         

         for (int i = -grid_[0]; i < grid_[0]; ++i) {
            dR0 = boundaryPtr->bravaisBasisVector(0);
            dR0 *= i;
            dR0 /= grid_[0]*grid_[0];
            for (int j = -grid_[1]; j < grid_[1]; ++j) {
               dR1 = boundaryPtr->bravaisBasisVector(1);
               dR1 *= j;
               dR1 /= grid_[1]*grid_[1];
               for (int k = -grid_[2]; k < grid_[2]; ++k) {
                  dR2 = boundaryPtr->bravaisBasisVector(2);
                  dR2 *= k;
                  dR2 /= grid_[2]*grid_[2];

                  R.zero();
                  R += dR0;
                  R += dR1;
                  R += dR2;
                  R += Rmax;
                  CV = 0.0;

                  for (int w = 0; w < nWave_; ++w) {
                     CV += exp(-Constants::Im*phases_[w])*exp(-Constants::Im*R.dot(waveVectors_[w]))*randomWaves_[w];
                  }

                  CVreal = std::real(CV);
                  if (CVmax < CVreal) {
                     imax = i;
                     jmax = j;
                     kmax = k;
                     CVmax = CVreal;
                  }
               }
            }
         }

         dR0 = boundaryPtr->bravaisBasisVector(0);
         dR0 *= imax;
         dR0 /= grid_[0]*grid_[0];
         Rmax += dR0;
         dR1 = boundaryPtr->bravaisBasisVector(1);
         dR1 *= jmax;
         dR1 /= grid_[1]*grid_[1];
         Rmax += dR1;
         dR2 = boundaryPtr->bravaisBasisVector(2);
         dR2 *= kmax;
         dR2 /= grid_[2]*grid_[2];
         Rmax += dR2;
         
         for (int i = -grid_[0]; i < grid_[0]; ++i) {
            dR0 = boundaryPtr->bravaisBasisVector(0);
            dR0 *= i;
            dR0 /= grid_[0]*grid_[0]*grid_[0];
            for (int j = -grid_[1]; j < grid_[1]; ++j) {
               dR1 = boundaryPtr->bravaisBasisVector(1);
               dR1 *= j;
               dR1 /= grid_[1]*grid_[1]*grid_[1];
               for (int k = -grid_[2]; k < grid_[2]; ++k) {
                  dR2 = boundaryPtr->bravaisBasisVector(2);
                  dR2 *= k;
                  dR2 /= grid_[2]*grid_[2]*grid_[2];

                  R.zero();
                  R += dR0;
                  R += dR1;
                  R += dR2;
                  R += Rmax;
                  CV = 0.0;

                  for (int w = 0; w < nWave_; ++w) {
                     CV += exp(-Constants::Im*phases_[w])*exp(-Constants::Im*R.dot(waveVectors_[w]))*randomWaves_[w];
                  }

                  CVreal = std::real(CV);
                  if (CVmax < CVreal) {
                     imax = i;
                     jmax = j;
                     kmax = k;
                     CVmax = CVreal;
                  }
               }
            }
         }

         accumulatorU_.sample(CVmax);
         accumulatorN_.sample(CVmax/abs);
         accumulatorUavg_.sample(CVavg);
         accumulatorNavg_.sample(CVavg/abs);
         outputFile_ << Int(iStep, 10)
            << Dbl(CVmax, 20)
            << Dbl(CVmax/abs, 20)
            << Dbl(CVavg, 20)
            << Dbl(CVavg/abs, 20)
            << Dbl(abs, 20)
            << std::endl;
         outputFile_.close();
      }

      ++nSample_;

   }

   /*
   * Calculate floating point wavevectors.
   */
   void RandomWaves::makeWaveVectors() 
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

   /*
   * Write data to three output files.
   */
   void RandomWaves::output()
   {
      if (simulation().domain().isMaster()) {

         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), 
                                                  outputFile_);
         writeParam(outputFile_);
         outputFile_.close();

         // Write parameters to a *_avg.dat file
         simulation().fileMaster().openOutputFile(outputFileName("_Uavg.dat"), 
                                                  outputFile_);
         accumulatorU_.output(outputFile_);
         outputFile_.close();

         // Write parameters to a *_avg.dat file
         simulation().fileMaster().openOutputFile(outputFileName("_Navg.dat"), 
                                                  outputFile_);
         accumulatorN_.output(outputFile_);
         outputFile_.close();

      }
   }

}
