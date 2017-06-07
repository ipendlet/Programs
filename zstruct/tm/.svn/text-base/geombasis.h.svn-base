#ifndef GEOMBASIS_H
#define GEOMBASIS_H

#include "stringtools.h"
#include "pTable.h"
#include "icoord.h"
#include "atom.h"

class GeomBasis {

  private:

    int nbasis;
    int natomt;
    string* xyzfile; //from xyzlist
    ICoord* icoord;
    
    void list_read(string xyzlist);
    int compare(Atom atom);
    int locate(Atom atom);
    int compare1(Atom atom1, Atom atom2);


  public:
    
    Atom* ats;
    int dummy;
    int init(string xyzlist);
    void freemem();
    int find(int i, ICoord icoord);
    int nplus(int i, ICoord icoord, int atype0);
    int nminus(int i, ICoord icoord, int atype0);
    int nplus2(int i, ICoord icoord, int atype0);
    int nminus2(int i, ICoord icoord, int atype0);

};



#endif

