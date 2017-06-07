#ifndef ATOMS_H
#define ATOMS_H

#include "stringtools.h"
#include "pTable.h"
#include "icoord.h"
//#include "geombasis.h"

class Atom {

  private:

    
//    void list_read(string xyzlist);
    

  public:
    
    int coordn;
    ICoord lic; // local internals
    void alloc(int size);
    void init(int natoms, string* anam, int* anum, double* localxyz);
    void reset(int natoms, string* anam, int* anum, double* localxyz);

};



#endif

