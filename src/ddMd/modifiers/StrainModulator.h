#ifndef DDMD_STRAIN_MODULATOR_H
#define DDMD_STRAIN_MODULATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/modifiers/Modifier.h>        // base class
#include <util/misc/FileMaster.h>           // member variable
#include <util/containers/DArray.h>         // member variable
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{

   using namespace Util;

   /**
   * The design of the StrainModulator class was inspired by the Lammps "Fix" 
   * base class, which provides a very flexible framework for designing
   * algorithms that can modify the state of a system.
   *
   * \ingroup DdMd_StrainModulator_Module
   */
   class StrainModulator : public Modifier
   {

   public:

      /**
      * Default constructor (for unit testing)
      */
      StrainModulator();

      /**
      * Constructor (for use in simulation).
      */
      StrainModulator(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~StrainModulator();

      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Setup before entering the main loop.
      */ 
      virtual void setup();

      /// \name Integration action functions
      //@{ 
 
      /** 
      * Call just before the first step of velocity-Verlet algorithm. 
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void preIntegrate1(long iStep);

      /** 
      * Call just before the first step of velocity-Verlet algorithm. 
      *
      * Atom positions are Cartesian on entry and return.
      */
      virtual void output();

      protected:

      /**
      * Return outputFileName string with added suffix.
      */
      std::string outputFileName(const std::string& suffix) const;

      private:

      /// Output file stream
      std::ofstream  outputFile_;

      /// Base name of output file(s).
      std::string outputFileName_;

      /// Stress Measurement Interval
      double  mFactor_;

      /// Stress Measurement Interval
      double  factor_;

      /// Stress Measurement Interval
      int  mInterval_;

      /// Stress Measurement Interval
      bool  arrayFlag_;

      /// Stress Measurement Interval
      Util::DArray<double>  array1_;

      /// Stress Measurement Interval
      Util::DArray<double>  array2_;

      /// Has readParam been called?
      long  isInitialized_;

   };

}
#endif
