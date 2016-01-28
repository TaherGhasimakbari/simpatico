/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <mpi.h>
#include <util/space/Tensor.h>
#include <util/accumulators/AutoCorrelation.tpp>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#define START      10000000
#define END        70000000
#define BOXLENGTH  22.04

using namespace std;
using namespace Util;
int main()
{
   MPI_Init(NULL, NULL);

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   
   ifstream inFile;
   AutoCorrelation<Tensor, double> accumulator;
   accumulator.setParam(64,10);
   accumulator.clear();

   double step;
   double junk;
   double lx;
   double ly;
   double lz;
   double pxx;
   double pyy;
   double pzz;
   double pxy;
   double pxz;
   double pyz;

   double factor=sqrt(BOXLENGTH*BOXLENGTH*BOXLENGTH/10.0);
   Tensor total;
   double pressure;
   
   std::string path ("/home/morsedc/tghasi/Simulation/Stress/StressRelaxation/S1_32_NVT/Box22.04_ODT/out/");
   std::string path_in;
   std::string path_out;

         
   path_in=path+toString(rank)+"/stress";
   std::cout<<path_in<<"\n";
   MPI_Barrier(MPI_COMM_WORLD);
   path_out=path+toString(rank)+"/stress.dat";

   inFile.open(path_in.c_str(), ios_base::in);

   if (!inFile)  {
      cerr << "can't open Data!\n";
   }

   for(int i = 0; i < END; i++)
   {

      inFile >> step;
      inFile >> junk;
      inFile >> lx;
      inFile >> ly;
      inFile >> lz;
      inFile >> pxx;
      inFile >> pyy;
      inFile >> pzz;
      inFile >> pxy;
      inFile >> pxz;
      inFile >> pyz;
      if (i == 0) {
         cout<<step<<"   "<<junk<<"   "<<lx<<"   "<<ly<<"   "<<lz<<"   "<<pxx<<"   "<<pyy<<"   "<<pzz<<"   "<<pxy<<"   "<<pxz<<"   "<<pyz<<"\n";
      }
            
      if (i == START) {
         cout<<"sampling starts!\n";
      }
      
      if (i>START) {
         if (i%10000000 == 0) {
            cout<<i<<" sampling\n";
         }
         pressure = (pxx+pyy+pzz)/3.0;

         total(0,0)=factor*(pxx-pressure);
         total(1,1)=factor*(pyy-pressure);
         total(2,2)=factor*(pzz-pressure);
         total(0,1)=factor*pxy;
         total(0,2)=factor*pxz;
         total(1,2)=factor*pyz;
         total(1,0)=factor*pxy;
         total(2,0)=factor*pxz;
         total(2,1)=factor*pyz;

         accumulator.sample(total);
      }

   }

   inFile.close();

   ofstream outFile;
   outFile.open(path_out.c_str());
   std::cout<<"writing auto-correlation results into file...";
   accumulator.output(outFile);
   std::cout<<"results were written to file: "<<path_out<<"\n";
   outFile.close();

   MPI_Finalize();
}
