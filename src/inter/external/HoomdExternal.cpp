#ifndef INTER_PERIODIC_EXTERNAL_CPP
#define INTER_PERIODIC_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdExternal.h"

#include <iostream>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   HoomdExternal::HoomdExternal() 
    : a_(),
      boundaryPtr_(0),
      nAtomType_(0), 
      isInitialized_(false)
   { setClassName("HoomdExternal"); }
   
   /* 
   * Copy constructor.
   */
   HoomdExternal::HoomdExternal(const HoomdExternal& other)
    : a_(other.a_),
      p_(other.p_),
      w_(other.w_),
      boundaryPtr_(other.boundaryPtr_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {
      a_.allocate(nAtomType_);
      for (int i=0; i < nAtomType_; ++i) {
        a_[i] = other.a_[i];
      }
      for (int i=0; i < Dimension; ++i) {
        b_[i] = other.b_[i];
      }
   } 
     
   /* 
   * Assignment operator.
   */
   HoomdExternal& HoomdExternal::operator = (const HoomdExternal& other)
   {
      a_             = other.a_;
      p_             = other.p_;
      w_             = other.w_;
      boundaryPtr_   = other.boundaryPtr_;
      nAtomType_     = other.nAtomType_;
      isInitialized_ = other.isInitialized_;
      for (int i=0; i < nAtomType_; ++i) {
        a_[i] = other.a_[i];
      }
      for (int i=0; i < Dimension; ++i) {
        b_[i] = other.b_[i];
      }
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void HoomdExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > HoomdExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void HoomdExternal::setExternalParameter(DArray<double> amplitude) 
   {  
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("HoomdExternal potential is not initialized");
      }

      for (int i=0; i < nAtomType_; ++i) {
        a_[i] = amplitude[i];
      }
   }

   
   /* 
   * Set pointer to the Boundary.
   */
   void HoomdExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /* 
   * Read potential parameters from file.
   */
   void HoomdExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
  
      // Read parameters
      a_.allocate(nAtomType_);
      readDArray<double>(in, "A", a_, nAtomType_);
      read<IntVector>(in, "b", b_);
      read<double>(in, "w", w_);
      read<int>(in, "p", p_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void HoomdExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_;
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      a_.allocate(nAtomType_);
      loadDArray<double>(ar, "A", a_, nAtomType_);
      loadParameter<IntVector>(ar, "b", b_);
      loadParameter<int>(ar, "p", p_);
      loadParameter<double>(ar, "w", w_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void HoomdExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << a_;
      ar << b_;
      ar << p_;
      ar << w_;
   }

   /*
   * Return name string "HoomdExternal".
   */
   std::string HoomdExternal::className() const
   {  return std::string("HoomdExternal"); }
 
} 
#endif
