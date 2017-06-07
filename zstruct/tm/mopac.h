#ifndef MOPAC_H
#define MOPAC_H

#include "stringtools.h"
#include "pTable.h"
#include "icoord.h"

class Mopac {
  
  private:
  
   int natoms;
   int* anumbers;
   string* anames;
   int nfrz; //total frozen atoms
   int nfrz0; //total "moved" frozen atoms
   int* frzlist;
   int* frzlistb;
   int* solvent;
   int nsolvent;

   void opt_header(ofstream& inpfile);
   void write_ic_input(ofstream& inpfile, int anum, ICoord icoords);

  public:

   double energy_sp(string filename);
   double opt();
   double opt(string filename);
   double opt(string filename, ICoord icoords);
   void opt_write();
   void opt_write(string filename);
   void opt_write(string filename, ICoord icoords);

   double read_output(string filename);
   void xyz_read(string filename);
   void xyz_save(string filename);


   void alloc(int natoms);
   void init(int natoms, int* anumbers, string* anames, double* xyz);
   void reset(int natoms, int* anumbers, string* anames, double* xyz);
   void freemem();

   void freeze(int* frzlist, int nfrz, int nfrz0);
   void addSolvent(int* solvent0, int nsolvent0);

   double energy0;
   double energy;

   double* xyz0;
   double* xyz;

};

#endif