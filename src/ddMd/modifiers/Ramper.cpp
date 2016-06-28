/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Ramper.h"
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
   Ramper::Ramper(Simulation& simulation)
    : Modifier(simulation)
   {set(Flags::PreIntegrate1);
    set(Flags::Setup);
    setClassName("Ramper");}

   /*
   * Destructor.
   */
   Ramper::~Ramper()
   {}

   /*
   * Read interval and outputFileName. 
   */
   void Ramper::readParameters(std::istream& in)
   {
      // Read interval value (inherited from Interval)
      readInterval(in);
      read<std::string>(in, "outputFileName", outputFileName_);
      read<double>(in,"epsilonStart", epsilonStart_);
      read<double>(in,"epsilonSlope", epsilonSlope_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void Ramper::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "outputFileName", outputFileName_);
      loadParameter(ar,"epsilonStart", epsilonStart_);
      loadParameter(ar,"epsilonSlope", epsilonSlope_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void Ramper::save(Serializable::OArchive &ar)
   {
      ar << outputFileName_;
   }

   /*
   * Get the outputFileName string with an added suffix
   */
   std::string Ramper::outputFileName(const std::string& suffix) const
   {
      std::string filename = outputFileName_;
      filename += suffix;
      return filename;
   }

   /*
   * Set actual number of molecules and clear accumulator.
   */
   void Ramper::setup()
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
   void Ramper::preIntegrate1(long iStep)
   {
      if (isAtInterval(iStep)) {

         double epsilon = epsilonStart_ + epsilonSlope_*iStep;
         simulation().pairPotential().set("epsilon", 0, 1, epsilon);
         
         if (simulation().domain().isMaster()) {
            outputFile_ << Int(iStep, 15) << Dbl(epsilon, 20) << std::endl;
         }
      }
   }

}
