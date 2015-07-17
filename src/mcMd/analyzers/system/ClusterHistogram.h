#ifndef MCMD_CLUSTER_HISTOGRAM_H
#define MCMD_CLUSTER_HISTOGRAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>              // base class template
#include <mcMd/simulation/System.h>                     // base class template parameter
#include <mcMd/analyzers/system/ClusterIdentifier.h>    // member
#include <util/accumulators/IntDistribution.h>          // member

namespace McMd
{

   using namespace Util;

   /**
   * This class is intended to identify Clusters in polymeric systems.
   */
   class ClusterHistogram : public SystemAnalyzer<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      ClusterHistogram(System &system);
   
      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval         : sampling interval
      *   - string outputFileName   : base name for output file(s)
      *   - int    speciesId        : integer id for Species of interest
      *   - int    atomTypeId       : integer id for core atom type
      *   - double cutoff           : distance cutoff
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /**
      * Evaluate squared radii of gyration for all molecules, add to ensemble.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive. 
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// ClusterIdentifier
      ClusterIdentifier identifier_;

      /// Distribution of the Clusters.
      IntDistribution  hist_;
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Clusters species type
      int  speciesId_;

      /// Clusters species type
      int  atomTypeId_;

      /// Distance cutoff
      double cutoff_;

      /// Clusters species type
      Species*  speciesPtr_;

      /// Histogram Min bin.
      int  histMin_;

      /// Histogram Min bin.
      int  histMax_;

      /// Number of configurations dumped thus far(first dump is zero).
      long  nSample_;

      /// Has readParam been called?
      bool  isInitialized_;

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void ClusterHistogram::serialize(Archive& ar, const unsigned int version)
   {  
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & histMin_;
      ar & histMax_;
   }

}
#endif
