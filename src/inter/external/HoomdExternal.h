#ifndef INTER_HOOMD_EXTERNAL_H
#define INTER_HOOMD_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>
#include <cmath>

namespace Inter
{

   using namespace Util;

   /**
   *
   * Hoomd External Potential class is intended to be reproduced! 
   *
   * \ingroup Inter_External_Module
   */
   class HoomdExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      HoomdExternal();

      /**
      * Copy constructor.
      */
      HoomdExternal(const HoomdExternal& other);

      /**
      * Assignment.
      */
      HoomdExternal& operator = (const HoomdExternal& other);

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * Sets external parameter
      *
      * \param externalParameter external parameter of system
      */
      void setExternalParameter(DArray<double> amplitude);

      /**
      * Set pointer to Boundary.
      *
      * \param boundary Boundary object (used to calculate length along perpendicular direction).
      */
      void setBoundary(Boundary &boundary);

      /**
      * Read potential parameters, and initialize other variables.
      *
      * \pre nAtomType must have been set, by calling setNAtomType().
      * \pre Boundary must have been set, by calling setBoundary().
      *
      * \param in input stream 
      */
      void readParameters(std::istream &in);

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
      * Returns external parameter
      *
      * \return external parameter
      */
      double externalParameter() const;

      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return     external potential energy
      */
      double energy(const Vector& position, int i) const;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      void getForce(const Vector& position, int type, Vector& force) const;
 
      /**
      * Return name string "HoomdExternal".
      */
      std::string className() const;
 
   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 3;

      /// Prefactor array ofsize nAtomType.
      DArray<double> a_;

      /// Array of Miller index IntVectors for the reciprocal lattice vectors.
      IntVector  b_;

      /// Number of unit cells in box
      int p_;

      /// Interface width
      double w_;

      /// Pointer to associated Boundary object.
      Boundary *boundaryPtr_;
   
      /// Number of possible atom types.
      int    nAtomType_; 

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // inline methods 

   /* 
   * Calculate external potential energy for a single atom.
   */
   inline double HoomdExternal::energy(const Vector& position, int type) const
   {
      const Vector cellLengths = boundaryPtr_->lengths();
      double clipParameter = 1.0/(2.0*M_PI*p_*w_);
      
      double cosine = 0.0;
      Vector q;
      q[0] = 2.0*M_PI*p_*b_[0]/cellLengths[0];
      q[1] = 2.0*M_PI*p_*b_[1]/cellLengths[1]; 
      q[2] = 2.0*M_PI*p_*b_[2]/cellLengths[2];
      double arg = q.dot(position);
      cosine += cos(arg);
      cosine *= clipParameter;
      return a_[type]*tanh(cosine);
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void HoomdExternal::getForce(const Vector& position, int type, 
                                     Vector& force) const
   {
      const Vector cellLengths = boundaryPtr_->lengths();
      double clipParameter = 1.0/(2.0*M_PI*p_*w_);
 
      double cosine = 0.0;
      Vector deriv;
      deriv.zero();
      Vector q;
      q[0] = 2.0*M_PI*p_*b_[0]/cellLengths[0];
      q[1] = 2.0*M_PI*p_*b_[1]/cellLengths[1]; 
      q[2] = 2.0*M_PI*p_*b_[2]/cellLengths[2];
      double arg = q.dot(position);
      cosine += cos(arg);

      double sine = -1.0*sin(arg);
      q *= sine;
      deriv += q;
      cosine *= clipParameter;
      deriv *= clipParameter;
      double tanH = tanh(cosine);
      double sechSq = (1.0 - tanH*tanH);
      double f = a_[type]*sechSq;
      deriv *= -1.0*f;

      force = deriv;
   }
 
}
#endif
