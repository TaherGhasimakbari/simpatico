/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CollectiveVariable.h"
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
   CollectiveVariable::CollectiveVariable(Simulation& simulation) 
    : Analyzer(simulation),
      isInitialized_(false),
      isFirstStep_(true)
   {  setClassName("CollectiveVariable"); }

   CollectiveVariable::~CollectiveVariable() 
   {}

   /// Read parameters from file, and allocate data array.
   void CollectiveVariable::readParameters(std::istream& in) 
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
      fourier_.allocate(nWave_);
      totalFourier_.allocate(nWave_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void CollectiveVariable::loadParameters(Serializable::IArchive &ar)
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

      // Allocate and load accumulators that exist only on master.
      if (simulation().domain().isMaster()) {
         CollectiveVariables_.allocate(nWave_);
         ar >> CollectiveVariables_;
      }

      // Allocate work space (all processors).
      waveVectors_.allocate(nWave_);
      fourier_.allocate(nWave_);
      totalFourier_.allocate(nWave_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void CollectiveVariable::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << mode_;
      ar << nWave_;
      ar << waveIntVectors_;
      ar << phases_;
      ar << grid_;
      ar << nSample_;
      ar << accumulator_;
   }
  
   /*
   * Clear accumulators.
   */
   void CollectiveVariable::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      assert (nWave_ > 0);
      accumulator_.clear();

      nSample_ = 0;
   }

   /*
   * Increment structure factor.
   */
   void CollectiveVariable::sample(long iStep) 
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

      isFirstStep_ = false;
      Vector position;
      std::complex<double> expFactor;
      double  product;
      AtomIterator  atomIter;
      int i, typeId;

      makeWaveVectors();

      // Set all Fourier modes to zero
      for (i = 0; i < nWave_; ++i) {
         fourier_[i] = std::complex<double>(0.0, 0.0);
      }

      simulation().atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         position = atomIter->position();
         typeId   = atomIter->typeId();
 
         // Loop over wavevectors
         for (i = 0; i < nWave_; ++i) {
            
            product = position.dot(waveVectors_[i]);
            expFactor = exp( product*Constants::Im );
            fourier_[i] += mode_[typeId]*expFactor;
         }
      }

      for (i = 0; i < nWave_; ++i) {
         totalFourier_[i] = std::complex<double>(0.0, 0.0);
      }
 
      #ifdef UTIL_MPI
      // Loop over wavevectors
      for (int i = 0; i < nWave_; ++i) {
         //Sum values from all processors.
         simulation().domain().communicator().
                      Reduce(&fourier_[i], &totalFourier_[i],
                             1, MPI::DOUBLE_COMPLEX, MPI::SUM, 0);
      }
      #else
      for (int i = 0; i < nWave_; ++i) {
         totalFourier_[i] = fourier_[i];
      }
      #endif

      Vector spacing;
      for (int i = 0; i < 3; ++i) {
         spacing[i] = 1.0/grid_[i];
      }
      
      Boundary* boundaryPtr = &simulation().boundary();

      std::complex<double> CV[grid_[0]][grid_[1]][grid_[2]];
      double maxCV = 0.0;
      double absCV = 0.0; 
      double abs = 0.0;

      int maxi = 0;
      int maxj = 0;
      int maxk = 0;

      Vector R;
      R.zero();
      Vector dR0;
      Vector dR1;
      Vector dR2;

      for (int w = 0; w < nWave_; ++w) {
         abs += std::norm(totalFourier_[w]);
      }
      abs = std::sqrt(abs);

      for (int i = 0; i < grid_[0]; ++i) {
         dR0 = boundaryPtr->bravaisBasisVector(0);
         dR0 *= i*spacing[0];
         for (int j = 0; j < grid_[1]; ++j) {
            dR1 = boundaryPtr->bravaisBasisVector(1);
            dR1 *= j*spacing[1];
            for (int k = 0; k < grid_[2]; ++k) {
               dR2 = boundaryPtr->bravaisBasisVector(2);
               dR2 *= k*spacing[2];

               R.zero();
               R += dR0;
               R += dR1;
               R += dR2;
               CV[i][j][k] = 0.0;

               for (int w = 0; w < nWave_; ++w) {
                  CV[i][j][k] += exp(-Constants::Im*phases_[w])*exp(-Constants::Im*R.dot(waveVectors_[w]))*totalFourier_[w];
               }

               absCV = std::abs(CV[i][j][k])/abs;
               if (maxCV < absCV) {
                  maxi = i;
                  maxj = j;
                  maxk = k;
                  maxCV = absCV;
               }
            }
         }
      }

      #if 0
      R.zero();
      double b0 = (std::norm(CollectiveVariable[(maxi-1)%grid_[0]][maxj][maxk])-std::norm(CollectiveVariable[(maxi-2)%grid_[0]][maxj][maxk]))/spacing[0];
      double a0 = std::norm(CollectiveVariable[(maxi-1)%grid_[0]][maxj][maxk])-b0*spacing[0]*(maxi-1);
      double d0 = (std::norm(CollectiveVariable[(maxi+1)%grid_[0]][maxj][maxk])-std::norm(CollectiveVariable[(maxi+0)%grid_[0]][maxj][maxk]))/spacing[0];
      double c0 = std::norm(CollectiveVariable[(maxi+1)%grid_[0]][maxj][maxk])-d0*spacing[0]*(maxi+1);
      dR0.zero();
      dR0 = boundaryPtr->bravaisBasisVector(0);
      dR0 *= (-(a0-c0)/(b0-d0))*spacing[0];

      double b1 = (std::norm(CollectiveVariable[maxi][(maxj-1)%grid_[1]][maxk])-std::norm(CollectiveVariable[maxi][(maxj-2)%grid_[1]][maxk]))/spacing[1];
      double a1 = std::norm(CollectiveVariable[maxi][(maxj-1)%grid_[1]][maxk])-b1*spacing[1]*(maxj-1);
      double d1 = (std::norm(CollectiveVariable[maxi][(maxj+1)%grid_[1]][maxk])-std::norm(CollectiveVariable[maxi][(maxj+0)%grid_[1]][maxk]))/spacing[1];
      double c1 = std::norm(CollectiveVariable[maxi][(maxj+1)%grid_[1]][maxk])-d1*spacing[1]*(maxj+1);
      dR1.zero();
      dR1 = boundaryPtr->bravaisBasisVector(1);
      dR1 *= (-(a1-c1)/(b1-d1))*spacing[1];

      double b2 = (std::norm(CollectiveVariable[maxi][maxj][(maxk-1)%grid_[2]])-std::norm(CollectiveVariable[maxi][maxj][(maxk-2)%grid_[2]]))/spacing[2];
      double a2 = std::norm(CollectiveVariable[maxi][maxj][(maxk-1)%grid_[2]])-b2*spacing[2]*(maxk-1);
      double d2 = (std::norm(CollectiveVariable[maxi][maxj][(maxk+1)%grid_[2]])-std::norm(CollectiveVariable[maxi][maxj][(maxk+0)%grid_[2]]))/spacing[2];
      double c2 = std::norm(CollectiveVariable[maxi][maxj][(maxk+1)%grid_[2]])-d2*spacing[2]*(maxk+1);
      dR2.zero();
      dR2 = boundaryPtr->bravaisBasisVector(0);
      dR2 *= (-(a2-c2)/(b2-d2))*spacing[0];
      R += dR0;
      R += dR1;
      R += dR2;
      std::complex<double> maxCV;
      for (int w = 0; w < nWave_; ++w) { 
         maxCV += exp(-phases_[w]*Constants::Im)*exp(-R.dot(waveVectors_[w])*Constants::Im)*totalFourier_[w];
      }
      #endif

      if(simulation().domain().isMaster()) {
         accumulator_.sample(maxCV);
         outputFile_ << Int(iStep, 10)
                     << Dbl(maxCV, 20)
                     << std::endl;
         outputFile_.close();
      }

      ++nSample_;

   }

   /*
   * Calculate floating point wavevectors.
   */
   void CollectiveVariable::makeWaveVectors() 
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
   void CollectiveVariable::output()
   {
      if (simulation().domain().isMaster()) {

         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), 
                                                  outputFile_);
         writeParam(outputFile_);
         outputFile_.close();

         // Write parameters to a *_avg.dat file
         simulation().fileMaster().openOutputFile(outputFileName("_avg.dat"), 
                                                  outputFile_);
         accumulator_.output(outputFile_);
         outputFile_.close();

      }
   }

}
