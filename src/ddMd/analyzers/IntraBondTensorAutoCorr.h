#ifndef DDMD_INTRA_BOND_TENSOR_AUTO_CORR_H
#define DDMD_INTRA_BOND_TENSOR_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/AutoCorrAnalyzer.h>   // base class template
#include <util/space/Tensor.h>                   // member template parameter

namespace DdMd
{

   using namespace Util;

   /**
   * Autocorrelation for bond stress of a molecule.
   *
   * The bond orientation tensor for each molecule is defined as a sum
   * \f[ 
   * 
   *    S_{ij} = \sum_{a}( u_{ai}u_{aj} - \delta_{ij} )
   *
   * \f]
   * where \f$u_{ai}\f$ is the ith Cartesian component (i=0,..,2) of a 
   * unit vector \f${\bf u}_{a}\f$ parallel to bond number a, and the 
   * sum is over all bonds in a molecule. This analyzer calculates the 
   * quantity:
   * \f[
   *
   *  C(t) = \sum_{i,j=0}^{2} \langle S_{ij}(t)S_{ij}(0) \rangle
   *
   * \f]
   * for molecules of a specified species. 
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */

   class IntraBondTensorAutoCorr : public AutoCorrAnalyzer<Tensor, double>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      IntraBondTensorAutoCorr(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~IntraBondTensorAutoCorr();

      using AutoCorrAnalyzer<Tensor, double>::readParameters;
      using AutoCorrAnalyzer<Tensor, double>::loadParameters;
      using AutoCorrAnalyzer<Tensor, double>::save;
      using AutoCorrAnalyzer<Tensor, double>::clear;
      using AutoCorrAnalyzer<Tensor, double>::setup;
      using AutoCorrAnalyzer<Tensor, double>::sample;
      using AutoCorrAnalyzer<Tensor, double>::output;

   protected:

      virtual void computeData();
      virtual Tensor data();

   };

}
#endif
