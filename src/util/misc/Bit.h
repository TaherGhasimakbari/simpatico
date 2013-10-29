#ifndef UTIL_BIT_H
#define UTIL_BIT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Represents a specific bit location within an unsigned int.
   *
   * Provides methods to query, set or clear a particular bit.
   */
   class Bit 
   {
   public:
  
      /**
      * Constructor.
      *
      * \param shift location of the bit, 0 < shift <= 32.
      */ 
      Bit(unsigned int shift);
  
      /**
      * Set this bit in the flags parameter
      *
      * \param flags unsigned int to be modified
      */
      void set(unsigned int& flags) const;
   
      /**
      * Clear this bit in the flags parameter
      *
      * \param flags unsigned int to be modified
      */
      void clear(unsigned int& flags) const;

      /**
      * Is this bit set in the flags integer?
      *
      * \param flags unsigned int to be queried
      */
      bool isSet(unsigned int flags) const;
   
      /**
      * Return integer with only this bit set.
      */
      unsigned int value() const;

   private:

      /// Integer with this bit set, all others clear.
      unsigned int value_;
   
   };

   /*
   * Set this bit in the flags parameter.
   */
   inline void Bit::set(unsigned int& flags) const
   {  flags |= value_; }

   /*
   * Clear this bit in the flags parameter.
   */
   inline void Bit::clear(unsigned int& flags) const
   {  flags &= (~value_); }

   /*
   * Is this bit set in the flags integer?
   */
   inline bool Bit::isSet(unsigned int flags) const
   {  return bool(flags & value_); }

   /*
   * Return unsigned int with only this bit set.
   */
   inline unsigned int Bit::value() const
   {  return value_; }
   
}
#endif
