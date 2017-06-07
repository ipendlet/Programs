#include "icoord.h"
#include "utils.h"

//for CG optimizer
#define SCALESD 25
#define SCALE1 1.4
#define SCALEA 3.0 
#define OPTTHRESH 0.01
#define MAXAD 0.03

//#define OPTSTEPS 240
#define OPTSTEPS 720
#define SHADOWSTEPS 120
 
int ICoord::opt(string xyzfile_string){

//  printf("  \n"); 
  dxm1 = new double[3*natoms];

  mm_init(); //vdw params

  //printf(" in opt, printing ic \n");
  //print_ic();

  if (natomstm>0)
    TM_setup();

  //printf(" in opt, first grad() \n");
  mm_grad();
  //print_grad();

  ofstream xyzfile;
//  string xyzfile_string = "xyzfile.txt";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  for (int n=0;n<OPTSTEPS;n++)
  {
//    printf(" Opt step: %i ",n+1);

    if (n%5==0)
    {
      xyzfile << " " << natoms << endl;
      xyzfile << " " << gradrms << endl;
      for (int i=0;i<natoms;i++) 
      {
        xyzfile << "  " << anames[i];
        xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
        xyzfile << endl;
      }
    }

    if (n==30) 
      pgradrms = 10000; // reset CG optimizer
    if (n==0) update_xyz_sd();
    else update_xyz_cg();
    mm_grad();
//    print_grad();
//    printf(" gradient RMS: %1.4f \n",gradrms);
    if (gradrms<OPTTHRESH) break;
  }

  printf("\n final MM grad RMS: %1.4f \n",gradrms);

//  printf("\n updating IC's \n");
  update_ic();
  for (int i=0;i<nbonds;i++)
  {
    int a1,a2;
    a1 = bonds[i][0];
    a2 = bonds[i][1];
    if (bond_stretch(a1,a2)/ffbondd(a1,a2)>1.5)
      printf(" warning: bond %i %i far from eq. dist %1.2f, actual: %1.2f \n",a1,a2,ffbondd(a1,a2),bond_stretch(a1,a2));
  }

  //print_ic();
//  printf(" new XYZ \n");
//  print_xyz();
//  print_bonds();
//  printf("\n");

  delete [] dxm1;

#if 1
  if (gradrms>1 || gradrms!=gradrms)
  {
    printf(" opt fails \n");
    return 0;
  }
  else 
#endif
    return 1;
}

//Shadow mixes in old IC into gradient
int ICoord::opt(string xyzfile_string, ICoord shadow){

  //return opt(xyzfile_string);

  //printf(" in opt/shadow \n");
  dxm1 = new double[3*natoms];

  mm_init();
  //printf(" opt shadow natomstm: %i \n",natomstm);
  if (natomstm>0)
  {
   // printf(" new structure TM_setup \n");
    TM_setup();
    printf(" shadow ic_create/TM_setup \n");
    shadow.ic_create();
    shadow.TM_setup();
    //printf(" after shadow TM_setup \n");
  }

  //printf(" before shadow grad \n");
  //shadow.print_ic();
  //printf(" in opt(shadow) IC's for new struct \n");
  //print_ic();
  mm_grad(shadow);
  //print_grad();
  //printf(" after shadow grad \n");

  ofstream xyzfile;
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

  for (int n=0;n<OPTSTEPS;n++)
  {
//    printf(" Opt step: %i ",n+1);

    if (n%5==0)
    {
       xyzfile << " " << natoms << endl;
      xyzfile << " " << gradrms << endl;
      for (int i=0;i<natoms;i++) 
      {
        xyzfile << "  " << anames[i];
        xyzfile << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2];
        xyzfile << endl;
      }
    }

    if (n==SHADOWSTEPS) 
      pgradrms = 10000; // reset CG optimizer
    if (n==0) update_xyz_sd();
    else update_xyz_cg();

    if (n<SHADOWSTEPS)
    {
      for (int i=0;i<3*natoms;i++)
        shadow.coords[i] = coords[i];
      mm_grad(shadow);
    }
    else
      mm_grad();
//    print_grad();
//    printf(" gradient RMS: %1.4f \n",gradrms);
    if (gradrms<OPTTHRESH) break;
  }

  printf(" final MM grad RMS: %1.4f \n",gradrms);

//  printf("\n updating IC's \n");
  update_ic();
  for (int i=0;i<nbonds;i++)
  {
    int a1,a2;
    a1 = bonds[i][0];
    a2 = bonds[i][1];
    if (bond_stretch(a1,a2)/ffbondd(a1,a2)>1.5)
      printf(" warning: bond %i %i far from eq. dist %1.2f, actual: %1.2f \n",a1,a2,ffbondd(a1,a2),bond_stretch(a1,a2));
  }

  //print_ic();
//  printf(" new XYZ \n");
//  print_xyz();
//  print_ic();
//  printf("\n");

  delete [] dxm1;

  if (gradrms>1 || gradrms!=gradrms)
  {
    printf(" opt fails \n");
    print_bonds();
    return 0;
  }
  else 
    return 1;
}

int ICoord::opt(string xyzfile_string, ICoord shadow, int* slist){

  nsolvent = 0;
  for (int i=0;i<natoms0;i++)
  {
    solvent[i] = slist[i];
    shadow.solvent[i] = slist[i];
    if (slist[i]) nsolvent++;
  }
  if (nsolvent>0)
    printf(" opting with solvent, natoms: %i nsolvent: %i \n",natoms,nsolvent);

  return opt(xyzfile_string,shadow);

}

int ICoord::opt(int* slist){

  string xyzfile_string = "xyzfile.xyz";
  printf(" opting with solvent, natoms: %i \n",natoms);
  nsolvent = 0;
  for (int i=0;i<natoms0;i++)
  {
    solvent[i] = slist[i];
    if (slist[i]) nsolvent++;
  }

  return opt(xyzfile_string);

}

// steepest descent
void ICoord::update_xyz_sd(){
  
  for (int i=0;i<3*natoms;i++)
    if(abs(grad[i])>MAXAD)
      grad[i]=sign(grad[i])*0.1;
  double SCALE = SCALESD;
  for (int i=0;i<3*natoms;i++)
  {
    dxm1[i] = grad[i]/SCALE;
    coords[i] += dxm1[i];
  }

  return;
} 

void ICoord::update_xyz_cg(){
  
  for (int i=0;i<3*natoms;i++)
    if(abs(grad[i])>MAXAD)
      grad[i]=sign(grad[i])*MAXAD;
  double SCALE = SCALE1; // 1.5 was fine
  double SCALE2 = gradrms*gradrms/pgradrms/pgradrms / SCALEA;
  if (SCALE2 > 1/SCALEA) SCALE2 = 1/SCALEA;
//  printf(" GRADRMS: %1.3f SCALE2: %1.3f \n",gradrms,SCALE2);
  for (int i=0;i<3*natoms;i++)
    coords[i] += grad[i]/SCALE+dxm1[i]*SCALE2;

  return;
} 

int ICoord::opt(){

  string xyzfile = "scratch/xyzfile.xyz";
  return opt(xyzfile);
}

#if 0
    printf(" adding bond between 0 and 2 \n");
    bonds[nbonds++][0] = 0;
    bonds[nbonds-1][1] = 2;
    bondd[nbonds-1] = distance(0,2);
#endif
#if 0
    printf(" moving bond to eq position, attaching to axial \n");
    bonds[3][0] = 3;
    bonds[3][1] = 0;
    bonds[nbonds++][0] = 8;
    bonds[nbonds-1][1] = 2;
    
    bondd[3] = distance(3,8);
    bondd[nbonds-1] = distance(0,2);
#endif
#if 0
    printf(" moving bond to axial position \n");
    bonds[3][0] = 5;
    bonds[3][1] = 4;
    bondd[3] = distance(5,4);
#endif
#if 0
    printf(" swapping bond between 4 and 11 \n");
    printf(" bond was: %i %i \n",bonds[9][0],bonds[9][1]);
    bonds[9][0] = 11;
    printf(" bond is: %i %i \n",bonds[9][0],bonds[9][1]);
    bondd[9] = distance(11,14);
#endif
#if 0
    printf(" swapping bond between 4 and 11 \n");
    printf(" bond was: %i %i \n",bonds[2][0],bonds[2][1]);
    bonds[2][0] = 2;
    printf(" bond is: %i %i \n",bonds[2][0],bonds[2][1]);
    bondd[2] = distance(0,2);
#endif
