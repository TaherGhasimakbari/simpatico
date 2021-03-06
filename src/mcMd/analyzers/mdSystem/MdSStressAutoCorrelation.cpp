#ifndef MCMD_MD_SSTRESS_AUTOCORRELATION_CPP
#define MCMD_MD_SSTRESS_AUTOCORRELATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/


#include "MdSStressAutoCorrelation.h"        // class header

namespace McMd
{

   /* 
   * Constructor.
   */
   MdSStressAutoCorrelation::MdSStressAutoCorrelation(MdSystem& system)
    : SStressAutoCorrelation<MdSystem>(system)
   {  setClassName("MdSStressAutoCorrelation"); }

   /* 
   * Destructor.
   */
   MdSStressAutoCorrelation::~MdSStressAutoCorrelation()
   {}

}
#endif
