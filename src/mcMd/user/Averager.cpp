/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/accumulators/Average.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

using namespace std;
using namespace Util;

#define Master 0
int main(int argc, char **argv)
{
   bool pFlag = false;
   bool iFlag = false;
   bool oFlag = false;
   bool nFlag = false;
   bool bFlag = false;
   bool cFlag = false;
   bool sFlag = false;

   std::string path;
   std::string iDir;
   std::string oDir;
   std::string name;
   unsigned int nBegin=0;
   unsigned int nColumns=0;
   unsigned int sColumn=0;

   int c;
   opterr = 0;
   while ((c = getopt(argc, argv, "p:i:o:n:b:c:s:")) != -1){
      switch (c) {
         case 'p':
            pFlag = true;
            path = std::string(optarg);
            break;
         case 'i':
            iFlag = true;
            iDir = std::string(optarg);
            break;
         case 'o': 
            oFlag = true;
            oDir = std::string(optarg);
            break;
         case 'n':
            nFlag = true;
            name = std::string(optarg);
            break;
         case 'b':
            bFlag = true;
            nBegin = static_cast<unsigned int>(atoi(optarg));
            break;
         case 'c':
            cFlag = true;
            nColumns = static_cast<unsigned int>(atoi(optarg));
            break;
         case 's':
            sFlag = true;
            sColumn = static_cast<unsigned int>(atoi(optarg));
            break;
         case '?':
           std::cout << "Unknown option -" << optopt << "\n";
         }
   }

   if(!(pFlag && iFlag && oFlag && nFlag && bFlag && cFlag && sFlag)) std::cout<<"Few input information!\n\n";

   MPI_Init(NULL, NULL);
   int rank;
   int size;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   
   ifstream inFile;
   Average accumulator;
   accumulator.clear();

   std::string path_i;
   std::string path_o;
   path_i=path+"/"+iDir+"/"+name+toString(rank);
   path_o=path+"/"+oDir+"/"+name+toString(rank)+".avg";

   std::cout<<"rank"<<setw(3)<<rank<<":\n";
   std::cout<<"input directory:"<<path_i<<"\n";
   std::cout<<"output directory:"<<path_o<<"\n\n";
   inFile.open(path_i.c_str(), ios_base::in);
   if (!inFile) cerr << "can't open Data!\n";
   MPI_Barrier(MPI_COMM_WORLD);

   std::cout<<"rank"<<setw(3)<<rank<<":\t lines ignored:\n";
   string lines_ignored;
   for(unsigned int l=0;l<nBegin;l++){
   getline(inFile, lines_ignored);
   std::cout<<lines_ignored<<"\n";}
   MPI_Barrier(MPI_COMM_WORLD);

   double data;
   for(unsigned int i=0; i<sColumn; i++){
   inFile >> data;
   std::cout<<data<<"\t";}

   while(!inFile.eof()){
      for(unsigned int j=0; j<nColumns; j++) inFile >> data;
      accumulator.sample(data);
   }

   inFile.close();

   ofstream outFile;
   outFile.open(path_o.c_str());
   std::cout<<setw(3)<<rank<<": writing auto-correlation results into file...";
   accumulator.output(outFile);
   std::cout<<"results were written to file: "<<path_o<<"\n";
   outFile.close();

   MPI_Finalize();
   return 0;
}
