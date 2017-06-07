#include <iostream>
#include <fstream>
#include <stdio.h>

#include "zstruct.h"
#include "icoord.h"
#include "geombasis.h"
#include "mopac.h"
#include "iso.h"
#include "pgsm.h"

void print_triple_xyz(int natoms, string* anames, int* anumbers, double* c0, double* c1, double* c2, double* e);
void print_double_xyz_save(string file, int natoms, string* anames, int* anumbers, double* c0, double* c1, double* e);
void print_triple_xyz_save(string file, int natoms, string* anames, int* anumbers, double* c0, double* c1, double* c2, double* e);

using namespace std;

int ZStruct::generate_paths_lst(int id, int nsave) {

  printf("\n connecting all low energy paths via LST (DFT) \n");

  double* e = new double[3]; 
  double** isos = new double*[nsave];
  for (int i=0;i<nsave;i++)
    isos[i] = new double[3*geoms[id].natoms];

  PGSM gsm;

  gsm.alloc(geoms[id].natoms);

  e[0] = e[1] = e[2] = -1;
  for (int i=1;i<niso;i++)
  {
//    e[0] = geoms[id].dftenergy; e[1] = geoms[i].dftenergy;
    e[0] = e[1] = -1;
    string nstr=StringTools::int2str(i,4,"0"); // 2,"0" is number of 0s in total
    string xyzfilename = "scratch/initial"+nstr+".xyz";
    print_double_xyz_save(xyzfilename,geoms[id].natoms,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,geoms[i].coords,e);
  }

  for (int i=1;i<nsave;i++)
  {
    if (savelist[i])
    {
      printf(" about to do: %i: ",i);
      gsm.reset(geoms[id].natoms,geoms[id].anumbers,geoms[id].anames,geoms[id].coords,geoms[i].coords);
      dftelistlst[i] = gsm.gstring_lst(i);

      string nstr=StringTools::int2str(i,4,"0"); 
//      string xyzfilename = "scratch/stringfile"+nstr+".xyz";
//      string cmd = "mv stringfile.xyz "+xyzfilename;
//      system(cmd.c_str());
      string cmd = "mv lst.log scratch/lstout"+nstr;
      system(cmd.c_str());

      for (int j=0;j<3*geoms[id].natoms;j++)
        isos[i][j] = gsm.xyzts[j];

    } // if savelist[i]
  }

  for (int i=1;i<nts;i++)
    dftelistlst[i] = (dftelistlst[i] - dftelist[0])*627.5;


  printf("\n printing pathways from LST (DFT), E[0]: %1.4f \n",dftelist[id]);
  for (int i=1;i<nsave;i++)
    printf(" connection: %i to %i, energies: %1.1f - %1.1f - %1.1f \n",id,i,0.0,dftelistlst[i],dftelist[i]);

  for (int i=1;i<nsave;i++)
  {
    e[0] = dftelist[id]; e[1] = dftelistlst[i]; e[2] = dftelistlst[i];
    string nstr=StringTools::int2str(i,4,"0"); 
    string xyzfilename = "scratch/tslst"+nstr+".xyz";
    print_triple_xyz_save(xyzfilename,geoms[id].natoms,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,isos[i],geoms[i].coords,e);
  }

  delete [] e;
  for (int i=0;i<nts;i++)
    delete [] isos[i];
  delete [] isos;  

  return nsave;
}

int ZStruct::generate_paths_se(int id, int nsave) {

  printf(" connecting all low energy paths via gsm (SE) \n");

  double* e = new double[3]; 
  double** isos = new double*[nsave];
  for (int i=0;i<nsave;i++)
    isos[i] = new double[3*geoms[id].natoms];

  PGSM gsm;

  gsm.alloc(geoms[id].natoms);

  write_initialxyz(id,nsave);
#if 0
  e[0] = e[1] = e[2] = -1;
  for (int i=1;i<niso;i++)
  {
    e[0] = geoms[id].energy; e[1] = geoms[i].energy;
    string nstr=StringTools::int2str(i,4,"0"); // 2,"0" is number of 0s in total
    string xyzfilename = "scratch/initial"+nstr+".xyz";
    print_double_xyz_save(xyzfilename,geoms[id].natoms,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,geoms[i].coords,e);
  }
#endif

  for (int i=1;i<nsave;i++)
  {
    printf(" about to do: %i: ",i);
    gsm.reset(geoms[id].natoms,geoms[id].anumbers,geoms[id].anames,geoms[id].coords,geoms[i].coords);
    elistts[i] = gsm.gstring_se();

    string nstr=StringTools::int2str(i,4,"0"); 
    string xyzfilename = "scratch/stringfile"+nstr+".xyz";
    string cmd = "mv stringfile.xyz "+xyzfilename;
    system(cmd.c_str());
    cmd = "mv initial.xyz.gsm scratch/gsmout"+nstr;
    system(cmd.c_str());

    for (int j=0;j<3*geoms[id].natoms;j++)
      isos[i][j] = gsm.xyzts[j];
  }

  printf("\n printing pathways from GSM/SE, absolute energies \n");
  for (int i=1;i<nsave;i++)
    printf(" connection: %i to %i, energies: %1.1f - %1.1f - %1.1f \n",id,i,elist[id],elistts[i],elist[i]);
//    printf(" connection: %i to %i, e[start]: %1.3f ets: %1.3f e[end]: %1.3f \n",id,i,elist[id],elistts[i],elist[i]);

  for (int i=1;i<nsave;i++)
  {
    e[0] = elist[id]; e[1] = elistts[i]; e[2] = elist[i];
    string nstr=StringTools::int2str(i,4,"0"); 
    string xyzfilename = "scratch/ts"+nstr+".xyz";
    print_triple_xyz_save(xyzfilename,geoms[id].natoms,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,isos[i],geoms[i].coords,e);
  }

  delete [] e;
  for (int i=0;i<nts;i++)
    delete [] isos[i];
  delete [] isos;  

  return nsave;
}




int ZStruct::generate_paths_dft() {

  printf("\n connecting all low energy paths via gsm (DFT) \n");

  int id = 0;
  double** isos = new double*[nts];
  for (int i=0;i<nts;i++)
    isos[i] = new double[3*geoms[id].natoms];
  double* e = new double[3]; 

  PGSM gsm;

  gsm.alloc(geoms[id].natoms);

  e[0] = e[1] = e[2] = -1;

  for (int i=1;i<nts;i++)
  {
    if (savelist[i])
    {
      printf(" about to do: %i: ",i);
      string nstr=StringTools::int2str(i,4,"0"); 
      string initialfilename = "scratch/initial"+nstr+".xyz";
      dftelistts[i] = gsm.gstring_dft(initialfilename);
      dftelistts[i] = (dftelistts[i]-dftelist[id])*627.5;

      string stringfilename = "scratch/stringfiledft"+nstr+".xyz";
      string cmd = "mv stringfile.xyz "+stringfilename;
      system(cmd.c_str());
      cmd = "mv "+initialfilename+".gsmq scratch/gsmqout"+nstr;
      system(cmd.c_str());

    for (int j=0;j<3*geoms[id].natoms;j++)
      isos[i][j] = gsm.xyzts[j];
    } // if savelist[i]
  } //loop i over nts

  for (int i=1;i<nts;i++)
    if (savelist[i])
      printf(" connection: %i to %i, energies: %1.3f - %1.3f - %1.3f \n",id,i,0.0,dftelistts[i],dftelist[i]);

#if 1
  for (int i=1;i<nts;i++)
  {
    if (savelist[i])
    {
      e[0] = 0.; e[1] = dftelistts[i]; e[2] = dftelist[i];
      string nstr=StringTools::int2str(i,4,"0"); 
      string xyzfilename = "scratch/tsq"+nstr+".xyz";
      print_triple_xyz_save(xyzfilename,geoms[id].natoms,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,isos[i],geoms[i].coords,e);
    }
  }
#endif


  delete [] e;
  for (int i=0;i<nts;i++)
    delete [] isos[i];
  delete [] isos;  

  return 0;
}
