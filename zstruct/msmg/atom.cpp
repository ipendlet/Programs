#include "atom.h"
#include "geombasis.h"
#include "icoord.h"
#include "utils.h"
using namespace std;


void Atom::init(int natoms, string* anam, int* anum, double* xyz){

 coordn = natoms - 1;
// printf("  Atom init \n");


 lic.init(natoms,anam,anum,xyz);

// printf("\n\n");

 return;
}


void Atom::reset(int natoms, string* anam, int* anum, double* xyz){

 coordn = natoms - 1;
// printf(" Atom reset \n");


 lic.reset(natoms,anam,anum,xyz);
 lic.ic_create();

// printf("\n\n");

 return;
}

void Atom::alloc(int size){

// printf(" Atom alloc \n");


 lic.alloc(size);

// printf(" NOT YET making atom centered types \n");



 return;
}


