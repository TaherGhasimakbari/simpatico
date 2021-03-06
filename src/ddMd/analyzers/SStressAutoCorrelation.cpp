#ifndef DDMD_SSTRESS_TENSOR_AUTO_CORRELATION_CPP
#define DDMD_SSTRESS_TENSOR_AUTO_CORRELATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SStressAutoCorrelation.h"
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SStressAutoCorrelation::SStressAutoCorrelation(Simulation& simulation) 
    : Analyzer(simulation),
      temperature_(1),
      accumulator_(),
      capacity_(-1),
      isInitialized_(false)
   {  setClassName("SStressAutoCorrelation"); }

   /*
   * Read interval and outputFileName. 
   */
   void SStressAutoCorrelation::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<double>(in,"temperature", temperature_);
      read<int>(in,"capacity", capacity_);
      accumulator_.setParam(capacity_);

      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void SStressAutoCorrelation::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(capacity_);

      if (simulation().domain().isMaster()) {
         accumulator_.loadParameters(ar);
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void SStressAutoCorrelation::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);

      if (simulation().domain().isMaster()) {
         ar << accumulator_;
      }

   }

  
   /*
   * Read interval and outputFileName. 
   */
   void SStressAutoCorrelation::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }

      if (simulation().domain().isMaster()) {
         accumulator_.clear();  
      }
   }


   /*
 *    * Set actual number of molecules and clear accumulator.
 *       */
   void SStressAutoCorrelation::setup()
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized.");
      }

       // Set number of molecules and clear accumulator
       accumulator_.clear();
   }

   /*
   * Sample the stress tensor.
   */
   void SStressAutoCorrelation::sample(long iStep) 
   {  
         double element;
         Simulation& sys = simulation();
         double volume = sys.boundary().volume();

      if (isAtInterval(iStep))  {
         sys.computeVirialStress();
         sys.computeKineticStress();
         sys.atomStorage().computeNAtomTotal(sys.domain().communicator());

         if (sys.domain().isMaster()) {
            Tensor virial  = sys.virialStress();
            Tensor kinetic = sys.kineticStress();
            Tensor total = total.add(virial, kinetic);

            element = (total(0,1) + total(1,0)) / 2.0 * sqrt(volume/temperature_);          
            accumulator_.sample(element);
         }
      }
   }

   /*
   * Dump StressTensor Measurment results.
   */
   void SStressAutoCorrelation::output() 
   {

      if (simulation().domain().isMaster()) {
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << "bufferCapacity  " << accumulator_.bufferCapacity() << std::endl;
         outputFile_ << "nSample         " << accumulator_.nSample() << std::endl;
         outputFile_ << std::endl;
         outputFile_ << "average   " << accumulator_.average() << std::endl;
         outputFile_ << std::endl;
         outputFile_ << "Format of *.dat file" << std::endl;
         outputFile_ << "[int time in samples]  [double autocorrelation function]"
                     << std::endl;
         outputFile_ << std::endl;
         outputFile_.close();

         // Write xy autocorrelation function to data file
         simulation().fileMaster().openOutputFile(outputFileName(".corr"), outputFile_);
         accumulator_.output(outputFile_);
         outputFile_.close();
      }     
   }

}
#endif  
