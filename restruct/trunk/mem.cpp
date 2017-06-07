#include "icoord.h"

void ICoord::freemem(){

 // return;
   
#if 1
//looks okay
 for (int i=0;i<max_bonds;i++)
   delete [] bonds[i];
 delete [] bonds;
#endif

#if 0
//problem likely
 for (int i=0;i<max_angles;i++)
   delete [] angles[i];
 for (int i=0;i<max_angles;i++)
   delete [] angles0[i];
 //delete [] angles0; 
#endif

#if 0
//problem likely
 delete [] anglev;
 delete [] angled;
 delete [] angled0;
#endif

#if 1
//problem?
 for (int i=0;i<max_nonbond;i++)
   delete [] nonbond[i];
 delete [] nonbond;
//problem
 //delete [] nonbondd;
#endif

#if 0
//problem (are these allocated?)
 for (int i=0;i<max_torsions;i++)
   delete [] torsions[i];
 delete [] torsions;
 delete [] torv;
#endif

#if 1
//looks okay
 for (int i=0;i<max_imptor;i++)
   delete [] imptor[i];
 delete [] imptor;
 delete [] imptorv;
#endif

//bad!
// delete [] coordn;

#if 1
  for (int i=0;i<16;i++)
    delete [] planes[i];
  delete [] planes;

  for (int i=0;i<natoms0*2;i++)
    delete [] tmnoangles[i];
  delete [] tmnoangles;

  for (int i=0;i<natoms0*2;i++)
    delete [] tmlinangles[i];
  delete [] tmlinangles;

 for (int i=0;i<24;i++)
   delete [] tmtriangles[i];
  delete [] tmtriangles;
#endif

#if 1
  delete [] anumbers;
  delete [] amasses;
  delete [] anames;
  delete [] coords0;
  delete [] coordsts;
  delete [] coords; 
  delete [] ffR;
  delete [] ffeps;
  delete [] ffq;
#endif

#if 1
  delete [] bondstm;
  delete [] oatomtm;

  delete [] grad;

  delete [] solvent;
#endif

  return;

}


void ICoord::alloc_mem(){

 grad = new double[3*natoms];
 for (int i=0;i<natoms0*3;i++)
   grad[i]=0;

 nbonds = 0;
 max_bonds=natoms0*8;
 bonds = new int*[max_bonds];
 for (int i=0;i<max_bonds;i++)
   bonds[i]=new int[2];
 bondd = new double[max_bonds];

 nangles = 0;
 max_angles=natoms0*24;
 angles = new int*[max_angles];
 for (int i=0;i<max_angles;i++)
   angles[i]=new int[3];
 angles0 = new int*[max_angles];
 for (int i=0;i<max_angles;i++)
   angles0[i]=new int[3];
 anglev = new double[max_angles];
 angled = new double[max_angles];
 angled0 = new double[max_angles];

 ntor = 0;
 max_torsions=natoms0*60;
#if 1
 torsions = new int*[max_torsions];
 for (int i=0;i<max_torsions;i++)
   torsions[i] = new int[4];
#endif
 torv = new double[max_torsions];

 nimptor = 0;
 max_imptor=natoms0*10;
 imptor = new int*[max_imptor];
 for (int i=0;i<max_imptor;i++)
   imptor[i] = new int[4];
 imptorv = new double[max_imptor];


 n_nonbond = 0;
 max_nonbond = (2+natoms)*(2+natoms);
 nonbond = new int*[max_nonbond];
 for (int i=0;i<max_nonbond;i++)
   nonbond[i]=new int[2];
 nonbondd = new double[max_nonbond];

 coordn = new int[4+natoms];
 for (int i=0;i<4+natoms;i++)
    coordn[i]=0;

 planes = new int*[16];
 for (int i=0;i<16;i++)
   planes[i] = new int[6];

 tmnoangles = new int*[natoms0*2]; 
 for (int i=0;i<natoms0*2;i++)
   tmnoangles[i] = new int[3];

 tmlinangles = new int*[natoms0*2]; 
 for (int i=0;i<natoms0*2;i++)
   tmlinangles[i] = new int[3];

 tmtriangles = new int*[24];
 for (int i=0;i<24;i++)
   tmtriangles[i] = new int[10];

 ffR = new double[4+natoms];
 ffeps = new double[4+natoms];
 ffq = new double[4+natoms];

 bondstm = new int[4+natoms];
 oatomtm = new int[4+natoms];

 nsolvent = 0;
 solvent = new int[4+natoms];
 for (int i=0;i<4+natoms;i++)
   solvent[i]=0;
 //solvent_xyz = new double[3*(4+natoms)];

 return;
}

