#ifndef ZSTRUCT_H
#define ZSTRUCT_H

//standard library includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <string>
#include <cctype>
#include <ctime>
#include <sys/stat.h>

#include "stringtools.h"
#include "icoord.h"
#include "mopac.h"
#include "geombasis.h"


#define NMAXSTRUCTS 4999


class ZStruct {
 private:

//  int* natoms;			//number of atoms in each structure
  string comment;		//comment line of first structure read

  int niso;
  int nts;
  int nfail;
  int nmopacdiff;
  int nmopaclargediff;
 
  double oemax; //max energy for consideration

  int* savelist;
  int nrepeat;
  int* repeat;
  ICoord geom0; // initial structure
  ICoord geom1; // working geometry
  ICoord geom1r; // working geometry, reverse direction
  ICoord geomshadow; //shadow geom
  ICoord* geoms; // all isomers
  double* elist;
  double* elistts;
  double* dftelist;
  double* dftelistlst;
  double* dftelistgsm;
  double* dftelistts;

  int* frozen;
  int* solvent;
  void addSolventMopac(int natoms, Mopac mopt);
  void read_frozen(int natoms);
  int read_solvent(int natoms);

  int CHARGE;                   //charge 
  int SPIN;

  GeomBasis basis1; // database xyz

//  int gradJobCount;		//keeps track of the total number of gradient jobs performed

 
  int generate_isomers(int id, int nchoose);
  void write_initialxyz(int id, int nsave);
  int generate_paths_se(int id, int nchoose);
  int generate_paths_lst(int id, int nchoose);
  int generate_paths_lst_para(int id, int nsave);
  int generate_paths_dft();
  int dft_sp();
  int dft_sp_para();
  int dft_opt_para();
  int gsm_dft_para();
  int dft_ts_para();

//  int remove_duplicates(int nbranch, int tniso);
  int remove_duplicates(int nbranch);
  int diff_structure(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structureq(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structurec(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  int diff_structurecq(int nat, string* anames, int* anumbers, double* xyz1, double* xyz2);
  void save_unique(int type);

  void alloc_elists();
  void free_elists();

 public:

  int init(string xyzfile, string xyzlist, double ecutoff);
  int init(string xyzfile, string xyzlist);
  void begin();
  void go_mech_it(int nsteps);
  int do_one_step(int* final_ints);

};

#endif


 
 
