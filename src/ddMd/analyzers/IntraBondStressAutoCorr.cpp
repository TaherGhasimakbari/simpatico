#ifndef DDMD_INTRA_BOND_STRESS_AUTO_CORR_INC_H
#define DDMD_INTRA_BOND_STRESS_AUTO_CORR_INC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/stress/IntraBondStressAutoCorr.h>
#include <ddMd/analyzers/AutoCorrAnalyzer.tpp>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   IntraBondStressAutoCorr::IntraBondStressAutoCorr(Simulation& simulation) 
    : AutoCorrAnalyzer<Tensor, double>(simulation)
   { setClassName("IntraBondStressAutoCorr"); }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   Tensor IntraBondStressAutoCorr::data() 
   { 
      Tensor    intraBondStress;
      Vector    dr, f;
      double    rsq, pressure;
      int       i, j;

      // Confirm that nMolecule has remained constant
      if (nMolecule_ != system().nMolecule(speciesId_)) {
            UTIL_THROW("Number of molecules has changed.");
      }

      // Loop over molecules
      System::ConstMoleculeIterator molIter;
      Molecule::ConstBondIterator   bondIter;
      i = 0;
      system().begin(speciesId_, molIter); 
      for ( ; molIter.notEnd(); ++ molIter) {
         data_[i].zero();

         molIter->begin(bondIter);
         for ( ; bondIter.notEnd(); ++bondIter) {
            rsq = system().boundary().distanceSq(
                  bondIter->atom(0).position(), 
                  bondIter->atom(1).position(), dr);
            f = dr;
            f *= system().bondPotential()
                         .forceOverR(rsq, bondIter->typeId());
            incrementPairStress(f, dr, data_[i]);
         }

         // Remove trace
         pressure  = data_[i].trace() / double(Dimension);
         for (j = 0; j < Dimension; ++j) {
            data_[i](j, j) -= pressure;
         }

         ++i;
      }
      
      return intraBondStress;
   }

   /// Output results after simulation is completed.
   template <class SystemType>
   void IntraBondStressAutoCorr<SystemType>::output() 
   {  

      // Echo parameters to analyzer log file
      fileMaster().openOutputFile(outputFileName(), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
#endif 
