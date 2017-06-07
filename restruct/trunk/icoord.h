#ifndef ICOORD_H
#define ICOORD_H

#include "stringtools.h"
#include "pTable.h"

class ICoord {

  private:

    int data;
//    int* atypes;                  //array of MM atom types  // assign values directly?
    double* amasses;              //array of atomic masses
//    double* charge;               //array of MM atomic charges
    string comment;
    double* dxm1;

  int* oatomtm; 
  int nbondstm;
  int* bondstm;
  int ntmnoangles;
  int** tmnoangles;
  int ntmlinangles;
  int** tmlinangles;
  int ntri;
  int** tmtriangles;
 
    double* bondd;

    double* anglev;

    double* torv;

    double* imptorv;
    int** imptor;
 
    int nplanes;
    int** planes;

    int max_bonds;
    int max_angles;
    int max_torsions;
    int max_imptor;

    int max_nonbond;
    int n_nonbond;
    int** nonbond;
    double* nonbondd;

  int nsolvent; 
  int* solvent;
  //double* solvent_xyz;

  double tmRval;

  void structure_read(string xyzfile);
  void alloc_mem();
  void make_bonds();
  void coord_num();
  void make_angles();
  void make_angles_tm();
  void save_angles();
  void make_torsions();

  void make_imptor(); 
// ic creation function when bonds are already made
  void make_imptor_nobonds(); 
  int isImptor(int i, int a1, int a2, int a3);

  int make_nonbond();

//
  void update_bonds();
  void update_angles();
  void update_torsion();
  void update_imptor();
  void update_nonbond();

  void create_xyz();

  double getR(int i);

  // MM force field params
  double* ffR;
  double* ffeps;
  double* ffq;
  double ffbondd(int i, int j);
  double ffbonde(int i, int j);
  double* angled; //angle array for TM
  double* angled0; //angle array for TM
  double ffangled(int i, int j, int k); //a is angle #
  double ffanglee(int i, int j, int k);
  double fftord(int i, int j, int k, int l); 
  double fftore(int i, int j, int k, int l); 
  double fftorm(int i, int j, int k, int l); //multiplicity 
  double ffimptore(int i, int j, int k, int l); 
  double ffimptord(int i, int j, int k, int l); 
  // function to make arrays?

  //TM items
  int isTM(int anum);
  int isLinAngle(int i, int k);
  int isH2(int anum);
  void TM_planes();
  void TM_setup(); //after bonding change, wrapper function
  void TM_planes_2(); //reassign planes after bonding change
  void remove_TM_pairs();
  void TM_angles();

  // Gradient terms
  double* grad;
  double gradrms;
  double pgradrms;
  void print_grad();
  void solvent_grad();
  void bond_grad_all();
  void bond_grad_1(int i, int j);
  double bond_stretch(int i, int j);
  void angle_grad_all();
  void angle_grad_1(int i, int j, int k, int a); //i,j,k in bond, a is angle #
  void torsion_grad_all();
  void torsion_grad_1(int i, int j, int k, int l);
  void imptor_grad_all();
  void imptor_grad_1(int i, int j, int k, int l);
  void vdw_grad_all();
  void vdw_grad_1(int i, int j);
  void elec_grad_all();
  void elec_grad_1(int i, int j);

  //Optimizer
  void update_xyz_sd();
  void update_xyz_cg();

//  ofstream xyzfile;
  void print_xyzf(ofstream xyzfile); // print xyz coords to file

  public:

  int** bonds;
  int nbonds;
  int** angles;
  int nangles;
  int** angles0;
  int nangles0;
  int** torsions;
  int ntor;

    
  int id; //for geoms[id] in zstruct
  int pid; // previous structure id
  double seenergy;
  double segsmenergy;
  double dftenergy;
  double dftlstenergy;
  double dftgsmenergy;
  double dfttsenergy;

  int natoms;
  int natoms0; // before adding duplicate TM atoms
  int natomstm;
  int skipTMadd;
  double* coords;
  double* coordsts;
  double* coords0;
  string* anames;               //array of atomic symbols 
  int* anumbers;                //array of atomic indices 
  int* coordn;                  //coordination number
  int nimptor;

  int ic_create();
  int ic_create_nobonds();
  int mm_grad();
  int mm_grad(ICoord shadow);
  int opt();
  int opt(string xyzfile, ICoord shadow, int* solvent);
  int opt(string xyzfile, ICoord shadow);
  int opt(int* solvent);
  int opt(string xyzfile);
  void update_ic();
  void mm_init();
  void mm_init_tm();

// help functions for iso
  int bond_exists(int b1, int b2);
  int bond_exists_tri(int b1, int b2);
  int bond_num(int b1, int b2);
  int hpair(int a1, int a2);
  int isTMpair(int a1, int a2);
  int h2count();

  int same_struct(double* xyz);

  int init(string xyzfile);
  int init(int natoms, string* anames, int* anumbers, double* xyz);
  int alloc(int size); 
  int reset(int natoms, string* anames, int* anumbers, double* xyz);
  int reset(double* xyz);
  void print_ic();
  void print_bonds();
  void print_xyz();
  void print_xyz_save(string filename);
  void print_xyz_save(string xyzfile_string, double energy);


  double distance(int i, int j);
  double distance_v0(int i); //distance compared to 0.0,0.0,0.0
  double angle_val(int i, int j, int k);
  double torsion_val(int i, int j, int k, int l); // for imptor and torsion
  double angle_val_v0(int i, int j); //angle compared to xyz = 0.0,0.0,0.0
  double torsion_val_v01(int i, int j); // compared to 0.0,0.0,0.0; 0.0,1.0,0.0

  void freemem();


};



#endif


