#ifndef PGSM_H
#define PGSM_H

#include "stringtools.h"
#include "pTable.h"

class PGSM {
  
  private:
  
   int natoms;
   int* anumbers;
   string* anames;


  public:

   double gstring_se();
   double gstring_se(string filename);
   double gstring_lst();
   double gstring_lst(int run);
   double get_lst(int run);
   double gstring_dft();
   double gstring_dft(string filename);
   void gstring_dft_dnr(string filename);
   double get_energy(string filename);

   void alloc(int natoms);
   void init(int natoms, int* anumbers, string* anames, double* xyz0, double* xyz1);
   void reset(int natoms, int* anumbers, string* anames, double* xyz0, double* xyz1);
   void freemem();
   void xyz_read(string file);

   double energyts;

   double* xyz0;
   double* xyz1;
   double* xyzts;

};

#endif
