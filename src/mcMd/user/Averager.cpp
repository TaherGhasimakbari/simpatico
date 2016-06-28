/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <mpi.h>
#include <util/accumulators/Average.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

using namespace std;
using namespace Util;
int main()
{
   MPI_Init(NULL, NULL);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   
   ifstream inFile;
   Average accumulator;
   accumulator.clear();

   std::string path;
   std::string in_dir;
   std::string out_dir;
   unsigned int N_Begin;
   unsigned int N_Finish;
   unsigned int N_Columns;
   unsigned int S_Column;

   

   if (rank==0) {
       std::cout<<"Path to Data: ";std::cin>>path;
       std::cout<<"Input Directory: ";std::cin>>in_dir;
       std::cout<<"Output Directory: ";std::cin>>out_dir;
       std::cout<<"Beginning Line of Sampling: ";std::cin>>N_Begin;
       std::cout<<"Finishing Line of Sampling: ";std::cin>>N_Finish;
       std::cout<<"Number of Columns: ";std::cin>>N_Columns;
       std::cout<<"The Column to Sample: ";std::cin>>S_Column;

   }

   std::string path_in;
   std::string path_out;
   path_in=path+in_dir;
   path_out=path+out_dir;
   inFile.open(path_in.c_str(), ios_base::in);

   string lines_ignored;
   for(unsigned int l=0;l<N_Begin;l++) getline(inFile, lines_ignored);
   MPI_Barrier(MPI_COMM_WORLD);

   std::cout<<setw(3)<<rank<<":"<<"\n";
   std::cout<<"input directory:"<<path_in<<"\n";
   std::cout<<"output directory:"<<path_out<<"\n";
   std::cout<<"Last Line Ignored:"<<"\n";
   std::cout<<lines_ignored<<"\n\n";
   MPI_Barrier(MPI_COMM_WORLD);

   if (!inFile) cerr << "can't open Data!\n";

   double data;
   for(unsigned int i=0; i<S_Column; i++) inFile >> data;
   for(unsigned int i = N_Begin; i < N_Finish; i++){
      for(unsigned int j=0; j<N_Columns; j++) inFile >> data;
      accumulator.sample(data);
   }

   inFile.close();

   ofstream outFile;
   outFile.open(path_out.c_str());
   std::cout<<setw(3)<<rank<<": writing auto-correlation results into file...";
   accumulator.output(outFile);
   std::cout<<"results were written to file: "<<path_out<<"\n";
   outFile.close();

   MPI_Finalize();
}
