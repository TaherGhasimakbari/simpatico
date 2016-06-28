/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StrainModulator.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StrainModulator::StrainModulator(Simulation& simulation)
    : Modifier(simulation)
   {setClassName("StrainModulator");}

   /*
   * Destructor.
   */
   StrainModulator::~StrainModulator()
   {}

   /*
   * Read interval and outputFileName. 
   */
   void StrainModulator::readParameters(std::istream& in)
   {
      // Read interval value (inherited from Interval)
      readInterval(in);
      read<std::string>(in, "outputFileName", outputFileName_);
      read<double>(in,"mFactor", mFactor_);
      factor_ = mFactor_;
      read<int>(in,"mInterval", mInterval_);
      arrayFlag_ = true;
      array1_.allocate(static_cast<int>(mInterval_/interval()));
      array2_.allocate(static_cast<int>(mInterval_/interval()));

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void StrainModulator::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "outputFileName", outputFileName_);
      loadParameter(ar,"mFactor", mFactor_);
      loadParameter(ar,"mInterval", mInterval_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void StrainModulator::save(Serializable::OArchive &ar)
   {
      ar << outputFileName_;
   }

   /*
   * Get the outputFileName string with an added suffix
   */
   std::string StrainModulator::outputFileName(const std::string& suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

   /*
   * Set actual number of molecules and clear accumulator.
   */
   void StrainModulator::setup()
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized.");
      }

      if (simulation().domain().isMaster()) {
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
   }

   /*
    * Save internal state to an archive.
   */
   void StrainModulator::preIntegrate1(long iStep)
   {
      if (isAtInterval(iStep)) {
         simulation().computeVirialStress();
         simulation().computeKineticStress();
         
         factor_ = 1.0/factor_;
         if (iStep % mInterval_ == 0) {
            simulation().boundary().setTetragonal(factor_,std::sqrt(1.0/factor_));
            arrayFlag_ = ~arrayFlag_;
         }
            
         if (simulation().domain().isMaster()) {

            Tensor virial  = simulation().virialStress();
            Tensor kinetic = simulation().kineticStress();
            Tensor total = total.add(virial, kinetic);

            double factor = sqrt(simulation().boundary().volume()/10.0);
            for (int i = 0; i < Dimension; ++i) {
               for (int j = 0; j < Dimension; ++j) {
                  total(i,j) *= factor;
               }
            }
            
            if (arrayFlag_) {
               array1_[iStep % mInterval_] = array1_[iStep % mInterval_] + (total(1,1)+total(2,2))/2.0-total(0,0);
            } else {
               array2_[iStep % mInterval_] = array1_[iStep % mInterval_] + (total(1,1)+total(2,2))/2.0-total(0,0);
            }     

            outputFile_ << Int(iStep, 10) << std::endl;
            outputFile_.close();
         }
      }
   }

   /*
   * Output results.
   */
   void StrainModulator::output()
   {
      outputFile_.close();
      simulation().fileMaster().openOutputFile(outputFileName("_1.dat"), outputFile_);
      for (int i = 0; i < array1_.capacity(); i++) {
         outputFile_ << Int(i, 10) << array1_[i] << std::endl;
      }
      outputFile_.close();

      simulation().fileMaster().openOutputFile(outputFileName("_2.dat"), outputFile_);
      for (int i = 0; i < array2_.capacity(); i++) {
         outputFile_ << Int(i, 10) << array2_[i] << std::endl;
      }
      outputFile_.close();
   }

}
