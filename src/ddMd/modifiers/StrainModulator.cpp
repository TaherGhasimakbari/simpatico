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
      read<int>(in,"measurementInterval", measurementInterval_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void StrainModulator::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "outputFileName", outputFileName_);
      loadParameter(ar,"measurementInterval", measurementInterval_);

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
   }

   /*
    * Save internal state to an archive.
   */
   void StrainModulator::preIntegrate1(long iStep)
   {
      if (isAtInterval(iStep))  {
         Simulation& sim = simulation();
         sim.computeVirialStress();
         sim.computeKineticStress();

         if (sim.domain().isMaster()) {

            Tensor virial  = sim.virialStress();
            Tensor kinetic = sim.kineticStress();
            Tensor total = total.add(virial, kinetic);

            double factor = sqrt(sim.boundary().volume()/10.0);
            for (int i = 0; i < Dimension; ++i) {
               for (int j = 0; j < Dimension; ++j) {
                  total(i,j) *= factor;
               }
            }

            simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
            outputFile_ << Int(iStep, 10) << std::endl;
            outputFile_.close();
         }
      }
   }

}
