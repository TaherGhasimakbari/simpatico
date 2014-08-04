#ifndef MCMD_CLUSTERS_STATISTICS_CPP
#define MCMD_CLUSTERS_STATISTICS_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ClustersStatistics.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   ClustersStatistics::ClustersStatistics(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      speciesId_(),
      speciesPtr_(),
      coreId_(),
      cellList_(),
      clusters_(),
      clusterLengths_(),
      histMin_(1),
      histMax_(),
      hist_(),
      isInitialized_(false)
   {  setClassName("ClustersStatistics"); }

   /// Read parameters from file, and allocate arrays.
   void ClustersStatistics::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      speciesPtr_ = &system().simulation().species(speciesId_);

      read<int>(in, "coreId", coreId_);
      if (coreId_ < 0) {
         UTIL_THROW("Negative coreId");
      }
      
      read<double>(in, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }
      int nMolecule = speciesPtr_->capacity();

      clusters_.allocate(nMolecule);
      for(int i = 0; i < nMolecule; ++i) { 
          clusters_[i].self_ = 0;
          clusters_[i].clusterId_ = -1;
      }
      clusterLengths_.reserve(nMolecule);
 
      read<int>(in,"histMin", histMin_);
      read<int>(in,"histMax", histMax_);
      hist_.setParam(histMin_, histMax_);
      hist_.clear();

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void ClustersStatistics::loadParameters(Serializable::IArchive& ar)
   {
      // Load (everything but accumulators_)
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      speciesPtr_ = &system().simulation().species(speciesId_);

      loadParameter<int>(ar, "coreId", coreId_);
      if (coreId_ < 0) {
         UTIL_THROW("Negative coreId");
      }

      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }
      int nMolecule = speciesPtr_->capacity();

      clusters_.allocate(nMolecule);
      clusterLengths_.reserve(nMolecule);
 
      loadParameter<int>(ar,"histMin", histMin_);
      loadParameter<int>(ar,"histMax", histMax_);
      hist_.setParam(histMin_, histMax_);
      hist_.clear();

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void ClustersStatistics::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Clear accumulators.
   */
   void ClustersStatistics::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
   }

   /*
   * Clear accumulators.
   */
   void ClustersStatistics::findClusters(Molecule* molPtr, int clusterId) 
   { 
      GArray<int> mNeighbors;
      mNeighbors.reserve(25);
      
      CellList::NeighborArray aNeighbors;

      Molecule::AtomIterator atomIter;
      for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {                   // Identifying the neighbours! 
          if (atomIter->typeId() == coreId_) {                                         // Checks the atomType to make sure it has the right Type.
            cellList_.getNeighbors(atomIter->position(), aNeighbors);                  // Takes neighbors of molecules atom out of cellList.
            for (int i = 0; i < aNeighbors.size(); i++) {
                if (aNeighbors[i]->molecule().id() != atomIter->molecule().id()) {
                   if (clusters_[aNeighbors[i]->molecule().id()].clusterId_ == -1) {
                      clusters_[aNeighbors[i]->molecule().id()].clusterId_ = clusterId;
                      mNeighbors.append(aNeighbors[i]->molecule().id());               // Adds neighbor molecule to list of neighbors.
                   } else if (clusters_[aNeighbors[i]->molecule().id()].clusterId_ != clusterId) UTIL_THROW("Cluster Clash!"); 
                } 
            } 
          }
      }
      
      for (int i = 0; i < mNeighbors.size(); i++) {
          findClusters(clusters_[mNeighbors[i]].self_, clusterId);
      }

      if (!isInitialized_) UTIL_THROW("Object is not initialized");
   }

   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void ClustersStatistics::sample(long iStep) 
   { 
      if (isAtInterval(iStep)) {
         
         int nMolecule = system().nMolecule(speciesId_);
         int nAtom = nMolecule * speciesPtr_->nAtom();
         cellList_.allocate(nAtom, system().boundary(), cutoff_);

         cellList_.makeGrid(system().boundary(), cutoff_);
         cellList_.clear();
         
         int clusterId = 0;

         System::MoleculeIterator molIter;                                             // Loading cellList with atoms.  
         Molecule::AtomIterator atomIter;           
         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
             for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                 if (atomIter->typeId() == coreId_) {
                    system().boundary().shift(atomIter->position());
                    cellList_.addAtom(*atomIter);
                 }
             }              // Atom loop.
         }                  // Molecule loop.
         
         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
             if (clusters_[molIter->id()].clusterId_ == -1) {
                findClusters(molIter.get(), clusterId);
                clusterId++;
             }
         }
         
         for (int i = 0; i < nMolecule; i++) {
             clusterLengths_[clusters_[i].clusterId_]++;
         }
      }                     // If is at interval.
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ClustersStatistics::output() 
   {
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      outputFile_ << "Cluster Id" <<"\t"<< "Number of Molecules" <<"\n";
      for (int i = 0; i < clusterLengths_.size(); i++) {
          hist_.sample(clusterLengths_[i]);
          outputFile_ << i+1 <<"\t"<< clusterLengths_[i];
      }
      outputFile_.close();
  
      fileMaster().openOutputFile(outputFileName(".hist"), outputFile_);
      hist_.output(outputFile_);
      outputFile_.close();
   }

}
#endif 
