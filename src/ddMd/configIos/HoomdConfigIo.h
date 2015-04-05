#ifdef  DDMD_MOLECULES
#ifndef DDMD_HOOMD_CONFIG_IO_H
#define DDMD_HOOMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/configIos/ConfigIo.h>
#include <util/space/Vector.h>

#include <vector>

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Native / default format for configuration files.
   *
   * HoomdConfigIo is a ConfigIo that implements the default
   * configuration file format for ddSim, like DdMdConfigIo, but 
   * that outputs atoms in sequential order, sorted by atom id.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class HoomdConfigIo  : public ConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      HoomdConfigIo(bool hasMolecules = false);

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      HoomdConfigIo(Simulation& simulation, bool hasMolecules = false);

      /**
      * Read configuration file.
      *
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \pre  There are no atoms, ghosts, or groups.
      * \pre  AtomStorage is set for scaled / generalized coordinates
      *
      * \post atomic coordinates are scaled / generalized
      * \post there are no ghosts
      *
      * \param file input file stream
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master, collecting
      * atom and group data from all processors. Atoms and Groups are
      * output in order, sorted by global id.
      *
      * \param file output file stream
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      /*
      * Struct for atom data required in file format.  
      */
      struct IoAtom {

         Vector position;
         Vector velocity;
         int typeId;
         int id;
         AtomContext context;

         IoAtom()
          : position(0.0),
            velocity(0.0),
            typeId(-1),
            id(-1)
         {}

      };

      /*
      * Struct containing a Group<N> and a group id.
      */
      template <int N>
      struct IoGroup {
         int    id;
         Group<N> group;
      };

      /**
      * Array of atoms, ordered by global index.
      */
      std::vector<IoAtom> atoms_;

      /**
      * Flag to determine if molecule info is included in file format
      */
      bool hasMolecules_;

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      int readGroups(std::ifstream& file, 
                     const char* sectionLabel, const char* nGroupLabel,
                     GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ofstream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, GroupCollector<N>& collector);
   
   };

}
#endif
#endif
