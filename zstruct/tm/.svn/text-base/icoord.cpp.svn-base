#include "icoord.h"
#include "utils.h"
using namespace std;


int ICoord::init(string xyzfile){

// printf(" xyzfile: %s \n",xyzfile);
 printf("\n");
 cout << " xyzfile: " << xyzfile << endl;
 structure_read(xyzfile);
 
 natoms0 = natoms;
 natomstm = 0;
 skipTMadd = 0;
 //nimptor = 0;

 print_xyz();

 alloc_mem();
 // printf(" done allocating memory\n");

 int done = ic_create();

 //printf(" initializing MM parameters \n");
 mm_init();

 return 1;
}



// initialize by feeding in xyz coordinates
int ICoord::init(int nat, string* anam, int* anum, double* xyz){

// printf(" initializing icoord via xyz structure \n");
 skipTMadd = 0;
 natoms = natoms0 = nat;
 natomstm = 0;
 //nimptor = 0;

// printf(" natoms: %i \n",nat);
// for (int i=0;i<natoms;i++)
//    printf(" %1.3f %1.3f %1.3f \n",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);

//otherwise allocated in structure_read
 anumbers = new int[4+natoms];
 amasses = new double[4+natoms];
 anames = new string[4+natoms];
 coords = new double[(natoms+4)*3];
 coordsts = new double[(natoms+4)*3];
 coords0 = new double[(natoms+4)*3];

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
 //print_xyz();

 alloc_mem();

 int done = ic_create();
 //printf(" after ic_create \n");
 //print_ic(); 

 //printf(" initializing MM parameters \n");
 mm_init();
 //printf("\n\n");

 return 1;
}

// initialize memory only
int ICoord::alloc(int size){

 natoms = natoms0 = size;
 skipTMadd = 0;

//otherwise allocated in structure_read
 anumbers = new int[4+natoms];
 amasses = new double[4+natoms];
 anames = new string[4+natoms];
 coords = new double[(natoms+4)*3];
 coordsts = new double[(natoms+4)*3];
 coords0 = new double[(natoms+4)*3];

 alloc_mem();

 return 1;
}


// initialize by feeding in xyz coordinates
int ICoord::reset(double* xyz){

// printf(" resetting icoord via xyz structure \n");

// natoms0 = natoms;
// natomstm = 0;
 //nimptor = 0;

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// print_xyz();

// printf(" initializing MM parameters \n");
 mm_init();

 return 1;
}

// initialize by feeding in xyz coordinates
int ICoord::reset(int nat, string* anam, int* anum, double* xyz){

// printf(" resetting icoord via xyz structure \n");
 natoms = natoms0 = nat;
 natomstm = 0;
 //nimptor = 0;

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// print_xyz();

 mm_init();

 //printf("\n\n");

 return 1;
}

void ICoord::update_ic(){

  update_bonds();
  update_angles();
  //update_torsion();
  //printf(" nimptor: %i \n",nimptor);
  update_imptor();
  update_nonbond();

  return;
} 
 

// incomplete function to convert IC directly to XYZ
void ICoord::create_xyz()
{
  double* nxyz = new double[3*natoms];
  printf ("xyz_create not implemented\n");
  int* adone = new int[natoms];
  for (int i=0;i<natoms;i++) adone[i]=0;
  adone[0]=1;
  nxyz[0] = coords[0];
  nxyz[1] = coords[1];
  nxyz[2] = coords[2];
  
  double* v1 = new double[3];
  double* u1 = new double[3];
  v1[0] = coords[3] - coords[0];
  v1[1] = coords[4] - coords[1];
  v1[2] = coords[5] - coords[2];

  double norm1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  u1[0]=v1[0]/norm1;
  u1[1]=v1[1]/norm1;
  u1[2]=v1[2]/norm1;

  double R;
  for (int i=0;i<nbonds;i++)
    if ((bonds[i][0]==0 && bonds[i][1]==1) ||
        (bonds[i][1]==0 && bonds[i][0]==1) )
      R = bondd[i];

//alternatively, just put this along the x axis
  nxyz[3] = u1[0]*R;
  nxyz[4] = u1[1]*R;
  nxyz[5] = u1[2]*R;

  for (int i=0;i<nangles;i++)
  {

  }
 
  delete [] nxyz;
  delete [] u1;
  delete [] v1;
  delete [] adone;

  return;
}


int ICoord::ic_create()
{
  ntmnoangles = 0;
  ntmlinangles = 0;

  //printf(" Creating internals from xyz \n");
  make_bonds();

  //coord_num(); // counts # surrounding species
  //printf(" now making angles \n");
  nangles=0;
  make_angles();
  make_angles_tm();
  save_angles();

  //printf(" now making torsions \n");
  make_torsions();
  //printf(" now making improper torsions \n");
  make_imptor(); //for 3 coord species only

  //printf(" now counting nonbond\n");
  n_nonbond = make_nonbond(); //anything not connected by bond or angle

  //printf(" doing TM_planes \n");
  TM_planes();
  //printf(" after TM_planes \n");

}

int ICoord::ic_create_nobonds()
{
//  printf(" Creating internals from xyz, skipping bond making \n");

  coord_num(); // counts # surrounding species
  nangles=0;
  make_angles();
  make_torsions();
  make_imptor_nobonds();

  n_nonbond = make_nonbond(); //anything not connected by bond or angle

//  print_ic();
}


void ICoord::update_bonds(){  
  for (int i=0;i<nbonds;i++)
    bondd[i] = distance(bonds[i][0],bonds[i][1]);
  return;
}

void ICoord::update_angles(){
  for (int i=0;i<nangles;i++)
    anglev[i] = angle_val(angles[i][0],angles[i][1],angles[i][2]);
  return;
}

void ICoord::update_torsion(){
  for (int i=0;i<ntor;i++)
    torv[i]=torsion_val(torsions[i][0],torsions[i][1],torsions[i][2],torsions[i][3]);
  return;
}

void ICoord::update_imptor(){
  for (int i=0;i<nimptor;i++)
    imptorv[i]=torsion_val(imptor[i][0],imptor[i][1],imptor[i][2],imptor[i][3]);
  return;
}

void ICoord::update_nonbond(){
  for (int i=0;i<n_nonbond;i++)
    nonbondd[i] = distance(nonbond[i][0],nonbond[i][1]);
  return;
}

void ICoord::make_bonds()
{

  //printf("\n in make_bonds, natoms: %i natoms0: %i \n",natoms,natoms0);
  double MAX_BOND_DIST; 
  nbonds=0;
  for (int i=0;i<natoms0;i++)
    if (!isTM(i))
    for (int j=0;j<i;j++)
    {
      if (!isTM(j))
      {
        MAX_BOND_DIST = (getR(i) + getR(j))/2;
        if ((anumbers[i]==1 && anumbers[j]==8) ||
            (anumbers[i]==8 && anumbers[j]==1) )
          MAX_BOND_DIST = 1.2;
        double d = distance(i,j);
        if (d<MAX_BOND_DIST)
        {
         // printf(" found bond: %i %i dist: %f \n",i+1,j+1,d);
          bonds[nbonds][0]=i;
          bonds[nbonds][1]=j;
          bondd[nbonds]=d;
          nbonds++;
        }
      }
    } // loop j over !TM

  int* tmadded = new int[natoms+4];
  for (int i=0;i<natoms+4;i++)
    tmadded[i] = 0;
  for (int i=0;i<natoms+4;i++)
    oatomtm[i] = -1;
  nbondstm = 0;

  for (int i=0;i<natoms0;i++)
    for (int j=0;j<i;j++)
    if (isTM(i) || isTM(j))
    {
      MAX_BOND_DIST = (getR(i) + getR(j))/2;
      double d = distance(i,j);
      if (0.1<d && d<MAX_BOND_DIST)
      {
       // printf(" found TM bond: %i %i dist: %f \n",i+1,j+1,d);
        bonds[nbonds][0]=i;
        bonds[nbonds][1]=j;
        bondd[nbonds]=d;
        bondstm[nbondstm]=nbonds;
        nbonds++;
        nbondstm++;
      }
#if 1
      else if (d<0.1 && isTM(i))
      {
        printf(" assigned oatomtm[%i] \n",j);
        tmadded[i] = 1;
        oatomtm[i] = j;
      }
      else if (d<0.1 && isTM(j))
      {
        printf(" never called? \n");
        tmadded[j] = 1;
        oatomtm[i] = j;        
      }
#endif
    }

  //printf(" found nbonds: %i nbondstm: %i \n",nbonds,nbondstm);

  //CPMZ id triangles, disconnect C from C-H
  ntri = 0;
#if 1
  for (int i=0;i<nbonds;i++)
  for (int j=0;j<i;j++)
  {
    //printf(" i,j %i %i \n",i,j);
    if (bonds[i][0]==bonds[j][0])
    {
      for (int k=0;k<j;k++)
      {
        if (bonds[k][0]==bonds[i][1] && bonds[k][1]==bonds[j][1])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][0],bonds[j][1],bonds[k][0]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][0]; tmtriangles[ntri][2]=bonds[j][1]; tmtriangles[ntri][3]=bonds[k][0];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
        else if (bonds[k][1]==bonds[i][1] && bonds[k][0]==bonds[j][1])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][0],bonds[j][1],bonds[k][1]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][0]; tmtriangles[ntri][2]=bonds[j][1]; tmtriangles[ntri][3]=bonds[k][1];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
      }
    }
    else if (bonds[i][1]==bonds[j][1])
    {
      for (int k=0;k<j;k++)
      {
        if (bonds[k][0]==bonds[i][0] && bonds[k][1]==bonds[j][0])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][1],bonds[j][0],bonds[k][0]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][1]; tmtriangles[ntri][2]=bonds[j][0]; tmtriangles[ntri][3]=bonds[k][0];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
        else if (bonds[k][1]==bonds[i][0] && bonds[k][0]==bonds[j][0])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][1],bonds[j][0],bonds[k][1]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][1]; tmtriangles[ntri][2]=bonds[j][0]; tmtriangles[ntri][3]=bonds[k][1];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
      }
    }
    else if (bonds[i][0]==bonds[j][1])
    {
      for (int k=0;k<j;k++)
      {
        if (bonds[k][0]==bonds[i][1] && bonds[k][1]==bonds[j][0])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][0],bonds[j][0],bonds[k][0]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][0]; tmtriangles[ntri][2]=bonds[j][0]; tmtriangles[ntri][3]=bonds[k][0];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
        else if (bonds[k][1]==bonds[i][1] && bonds[k][0]==bonds[j][0])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][0],bonds[j][0],bonds[k][1]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][0]; tmtriangles[ntri][2]=bonds[j][0]; tmtriangles[ntri][3]=bonds[k][1];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
      }
    }
    else if (bonds[i][1]==bonds[j][0])
    {
      for (int k=0;k<j;k++)
      {
        if (bonds[k][0]==bonds[i][0] && bonds[k][1]==bonds[j][1])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][1],bonds[j][1],bonds[k][0]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][1]; tmtriangles[ntri][2]=bonds[j][1]; tmtriangles[ntri][3]=bonds[k][0];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
        else if (bonds[k][1]==bonds[i][0] && bonds[k][0]==bonds[j][1])
        {
          //printf(" found triangle: %i %i %i \n",bonds[i][1],bonds[j][1],bonds[k][1]);
          tmtriangles[ntri][0] = -1; 
          tmtriangles[ntri][1]=bonds[i][1]; tmtriangles[ntri][2]=bonds[j][1]; tmtriangles[ntri][3]=bonds[k][1];
          tmtriangles[ntri][4]=i; tmtriangles[ntri][5]=j; tmtriangles[ntri][6]=k;
          ntri++;
        }
      }
    }  //loop j<i
  } //loop i over nbonds
#endif

 // printf(" ntri: %i \n",ntri);
  for (int i=0;i<ntri;i++)
    printf(" triangle: %i %i %i \n",tmtriangles[i][1],tmtriangles[i][2],tmtriangles[i][3]);
  for (int i=0;i<ntri;i++)
  {
    if (!isTM(tmtriangles[i][1]) && !isTM(tmtriangles[i][2]) && !isTM(tmtriangles[i][3]))
    {
      printf(" triangle on non-TM, deleting \n");
      for (int j=i;j<ntri-1;j++)
        for (int k=0;k<7;k++)
          tmtriangles[j][k] = tmtriangles[j+1][k];
      ntri--; i--;
    }
  }
  for (int i=0;i<ntri;i++)
  {
    int t1,a1,a2;
    double d1,d2;
    if (isTM(tmtriangles[i][1]))
    {
      t1 = tmtriangles[i][1];
      a1 = tmtriangles[i][2];
      a2 = tmtriangles[i][3];
    }
    else if (isTM(tmtriangles[i][2]))
    {
      t1 = tmtriangles[i][2];
      a1 = tmtriangles[i][1];
      a2 = tmtriangles[i][3];
    }
    else if (isTM(tmtriangles[i][3]))
    {
      t1 = tmtriangles[i][3];
      a1 = tmtriangles[i][1];
      a2 = tmtriangles[i][2];
    }
    else
    {
      printf(" triangle on non-TM! \n");
      t1 = -1;
      ntri--;
    }
   // printf(" found TM %i for triangles, others: %i %i \n",t1,a1,a2);
    d1 = distance(t1,a1);
    d2 = distance(t1,a2);
   // printf(" distances: %1.2f %1.2f \n",d1,d2);

    int rem1, a12;
//    printf(" i: %i ntri: %i \n",i,ntri);

// first handle C-H interactions
    if ( (anumbers[a1]==1 && anumbers[a2]==6) ||
         (anumbers[a2]==1 && anumbers[a1]==6) )
    {
//CPMZ check, could just remove C
      if (d1>d2 && bond_num(t1,a1)!=-1) 
      {
        rem1 = bond_num(t1,a1);
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a2;
        tmtriangles[i][3] = a1;
        a12 = a1;
      }
      else if (bond_num(t1,a2)!=-1)
      {
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a1;
        tmtriangles[i][3] = a2;
        rem1 = bond_num(t1,a2);
        a12 = a2;
      }
      else rem1 = -1;
    } else rem1 = -1;
// then C-C interactions
    if ( (anumbers[a1]==6 && anumbers[a2]==6) ||
         (anumbers[a2]==6 && anumbers[a1]==6) )
    {
#if 0
//CPMZ need to check plane, keep atom in plane
      if (d1>d2 && bond_num(t1,a1)!=-1) 
      {
        rem1 = bond_num(t1,a1);
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a2;
        tmtriangles[i][3] = a1;
        a12 = a1;
      }
      else if (bond_num(t1,a2)!=-1)
      {
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a1;
        tmtriangles[i][3] = a2;
        rem1 = bond_num(t1,a2);
        a12 = a2;
      }
      else rem1 = -1;
#endif
    } else rem1 = -1;
// then handle B-H interactions
    if ( (anumbers[a1]==5 && anumbers[a2]==1) ||
         (anumbers[a2]==5 && anumbers[a1]==1))
    {
      if (anumbers[a1]==5 && bond_num(t1,a1)!=-1) 
      {
        rem1 = bond_num(t1,a1);
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a2;
        tmtriangles[i][3] = a1;
        a12 = a1;
      }
      else if (bond_num(t1,a2)!=-1)
      {
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a1;
        tmtriangles[i][3] = a2;
        rem1 = bond_num(t1,a2);
        a12 = a2;
      }
      else rem1 = -1;
    } else rem1 = -1;
// then handle B-N interactions
    if ( (anumbers[a1]==5 && anumbers[a2]==7) ||
         (anumbers[a2]==5 && anumbers[a1]==7))
    {
      if (anumbers[a1]==5 && bond_num(t1,a1)!=-1) 
      {
        rem1 = bond_num(t1,a1);
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a2;
        tmtriangles[i][3] = a1;
        a12 = a1;
      }
      else if (bond_num(t1,a2)!=-1)
      {
        tmtriangles[i][1] = t1;
        tmtriangles[i][2] = a1;
        tmtriangles[i][3] = a2;
        rem1 = bond_num(t1,a2);
        a12 = a2;
      }
      else rem1 = -1;
    } else rem1 = -1;

    if (rem1>-1)
    {
      //print_bonds();
      //print_xyz();
      printf(" found bond % to delete: %i %i \n",rem1,bonds[rem1][0],bonds[rem1][1]);
      coord_num();
      tmtriangles[i][0] = 1;
      tmtriangles[i][5] = coordn[t1]-1;
      tmtriangles[i][6] = coordn[a12]-1;
      tmtriangles[i][7] = t1;
      tmtriangles[i][8] = a12;
      for (int l=rem1;l<nbonds-1;l++)
      {
        bonds[l][0] = bonds[l+1][0];
        bonds[l][1] = bonds[l+1][1];
      }
      nbonds--;
    } //if found bond to delete
    else
      ntri--;

  } // loop i over ntri

  coord_num();
  if (!skipTMadd)
  for (int i=0;i<natoms0;i++)
  {
    if (isTM(i) && coordn[i]>0)
    {
    //  printf(" adding duplicate for TM %i with coordn: %i \n",i,coordn[i]);
      
      if (!tmadded[i])
      {
        oatomtm[natoms] = i;
        anumbers[natoms] = anumbers[i]+1000;
        anames[natoms] = "X"+anames[i];
        coords[3*natoms+0] = coords[3*i+0];
        coords[3*natoms+1] = coords[3*i+1];
        coords[3*natoms+2] = coords[3*i+2];
        tmRval = getR(i);

        natomstm++;
        natoms++;
      }
    }
  }
  //printf(" natoms: %i \n",natoms);
  //for (int i=0;i<natoms;i++)
  //  if (oatomtm[i]>-1)
  //    printf(" dummy atom %i of type %s is linked to atom %i \n",i,anames[i].c_str(),oatomtm[i]);

//  for (int i=0;i<natoms;i++)
//    printf(" coord num for atom %i is %i \n",i,coordn[i]);
  coord_num();

  delete [] tmadded;

  //printf(" \n done making bonds \n\n");
}

#if 0
//based on existing imptors
int ICoord::isImptor(int i, int a1, int a2, int a3) {

  //printf(" isImptor: %i, %i %i %i \n",i,a1,a2,a3);
  int t1 = 0;
  for (int j=0;j<nimptor;j++)
    if (i==imptor[j][2])
       if (a1==imptor[j][0] && a2==imptor[j][1] && a3==imptor[j][3] ||
           a1==imptor[j][0] && a3==imptor[j][1] && a2==imptor[j][3] ||
           a2==imptor[j][0] && a1==imptor[j][1] && a3==imptor[j][3] ||
           a2==imptor[j][0] && a3==imptor[j][1] && a1==imptor[j][3] ||
           a3==imptor[j][0] && a2==imptor[j][1] && a1==imptor[j][3] ||
           a3==imptor[j][0] && a1==imptor[j][1] && a2==imptor[j][3] )
       { t1=1; break; }

  return t1;
}
#endif

//based on computing imptor angles
int ICoord::isImptor(int i, int a1, int a2, int a3) {

  if (i==a1 || i==a2 || i==a3)
    return 0;
  if (a1==a2 || a1==a3 || a2==a3)
    return 0;
  if (a1<0 || a2<0 || a3<0)
    return 0;
  //printf(" isImptor: %i, %i %i %i \n",i,a1,a2,a3);
  double imptorvt = torsion_val(a1,a2,i,a3);
  //printf(" first try: %1.1f \n",imptorvt);
  if (abs(abs(imptorvt)-180.) < 25.0)
    return 1;

  imptorvt = torsion_val(a1,a3,i,a2);
  //printf(" second try: %1.1f \n",imptorvt);
  if (abs(abs(imptorvt)-180.) < 25.0)
    return -1;

  return 0;
}

void ICoord::coord_num()
{ 
  for (int i=0;i<natoms;i++)
    coordn[i] = 0;
  for (int i=0;i<nbonds;i++)
  {
    coordn[bonds[i][0]]++;
    coordn[bonds[i][1]]++;
  }
}

void ICoord::make_angles()
{
// makes all non_tm angles

  //include all consecutive connections 
  for (int i=0;i<nbonds;i++)
  {
     for (int j=0;j<i;j++)
     {
        int found = 0;
        if (bonds[i][0]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][1];
          nangles++; found = 1;
        }
        else if (bonds[i][0]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][0];
          nangles++; found = 1;
        }
        else if (bonds[i][1]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][1];
          nangles++; found = 1;
        }
        else if (bonds[i][1]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][0];
          nangles++; found = 1;
        }
        if (nangles>0 && found)
        {
          anglev[nangles-1]=angle_val(angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
          if (anglev[nangles-1]>170 || isTM(angles[nangles-1][1]))
          { nangles--; }
        }

     } //loop j
  } //loop i


  return;
}

void ICoord::make_angles_tm()
{
//called from TM_angles or ic_create

  if (natomstm==0) return;

  //include all consecutive connections 
  for (int i=0;i<nbonds;i++)
  {
     for (int j=0;j<i;j++)
     {
        int found=0;
        if (bonds[i][0]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][1];
          nangles++; found=1;
        }
        else if (bonds[i][0]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][0];
          nangles++; found=1;
        }
        else if (bonds[i][1]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][1];
          nangles++; found=1;
        }
        else if (bonds[i][1]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][0];
          nangles++; found=1;
        }
        if (found)
        {
#if 1
          if (isTM(angles[nangles-1][1]))
          {
           // printf(" potential angle: %i %i %i \n",angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
            for (int k=0;k<ntmnoangles;k++)     
            {    
              //printf(" check: %i %i to %i %i \n",tmnoangles[k][0],tmnoangles[k][2],angles[nangles-1][0],angles[nangles-1][2]);
              if (tmnoangles[k][0]==angles[nangles-1][0] && tmnoangles[k][2]==angles[nangles-1][2])
              { nangles--; found=0; break; }
              else if (tmnoangles[k][2]==angles[nangles-1][0] && tmnoangles[k][0]==angles[nangles-1][2])
              { nangles--; found=0; break; }
            }
          }
#endif
          if (found)
          {
            anglev[nangles-1]=angle_val(angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
//CPMZ change this later
//            if (anglev[nangles-1]<80 || anglev[nangles-1]>100 || !isTM(angles[nangles-1][1]))
//            if (anglev[nangles-1]>160 || !isTM(angles[nangles-1][1])
            if (!isTM(angles[nangles-1][1]))
            {
              //printf("stopped 3: %i %i %i \n",angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
              nangles--;
            }
          }
        }

     } //loop j
  } //loop i


  return;
}

void ICoord::save_angles(){

  for (int i=0;i<nangles;i++)
    for (int j=0;j<3;j++)
      angles0[i][j]=angles[i][j];
  for (int i=0;i<nangles;i++)
    angled0[i]=anglev[i];

  return;
}

void ICoord::make_torsions()
{
  int a1,b1,c1,a2,b2,c2;
  bool found;

  ntor = 0;

//  return;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       b1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       b2=angles[j][1];
       c2=angles[j][2];

      // printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,b1,c1,a2,b2,c2);

       if (a1==b2 && b1==a2)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=a2;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       else if (a1==b2 && b1==c2)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=c2;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
       else if (c1==b2 && b1==a2)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=c1;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       else if (c1==b2 && b1==c2)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=c1;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }

       if (found)
       {
         torv[ntor-1]=torsion_val(torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
        // printf(" torsion found: %i %i %i %i \n",torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
       }
//CPMZ maybe need to use abs val
       if (found && ((abs(torv[ntor-1])<165. && abs(torv[ntor-1])>15.0) || isTM(torsions[ntor-1][1]) || isTM(torsions[ntor-1][2])))
       {
         // printf(" not planar \n");
         ntor--;
       }
    }
  } 

  return;
}

void ICoord::make_imptor()
{
  int a1,m1,c1,a2,m2,c2;
  bool found;
  nimptor = 0;
  double imptorvt;

  // for 3 coord centers only
  for (int i=0;i<nangles;i++)
  {
//    for (int j=0;j<nangles;j++)
//    if (i!=j)
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       m1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       m2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,m1,c1,a2,m2,c2);

       if (m1==m2 && !isTM(m1) && coordn[m1]<4)
       {
         imptor[nimptor][2]=m1;
         if (a1==a2)
         {
             imptor[nimptor][0]=c1;
             imptor[nimptor][1]=a1;
             imptor[nimptor][3]=c2;
             nimptor++; found=true;
         }
         else if (a1==c2)
         {
             imptor[nimptor][0]=c1;
             imptor[nimptor][1]=a1;
             imptor[nimptor][3]=a2;
             nimptor++; found=true;
         }
         else if (c1==c2)
         {
             imptor[nimptor][0]=c1;
             imptor[nimptor][1]=a1;
             imptor[nimptor][3]=a2;
             nimptor++; found=true;
         }
         else if (c1==a2)
         {
             imptor[nimptor][0]=c1;
             imptor[nimptor][1]=a1;
             imptor[nimptor][3]=c2;
             nimptor++; found=true;
         }
       } // if m1==m2
       //if (nimptor>0)
       //  printf(" found? %i imptor: %i %i %i %i \n",found,imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
       if (found)
       {
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1 && coordn[imptor[k][2]]==3)
             { found = false; nimptor--; break; }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
         //printf(" imptorv[%i]: %1.4f is %i %i %i %i\n",nimptor,imptorvt,imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
         //if (abs(imptorvt) < 15.0 && abs(abs(imptorvt) - 180) < 15.0 ) { found = false; nimptor--; }
         if (abs(abs(imptorvt) - 180) > 15.0 ) { found = false; nimptor--; }
       }
       if (found) imptorv[nimptor-1] = imptorvt;
    }
  } 

  return;
}

void ICoord::make_imptor_nobonds()
{
  int a1,m1,c1,a2,m2,c2;
  bool found;
  double imptorvt;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       m1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       m2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,m1,c1,a2,m2,c2);

       if (m1==m2 && !isTM(m1) && coordn[m1]<4)
       {
         if (a1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
         else if (a1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
       } // if m1==m2
       if (found)
       {
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }
       }
       //printf(" anumber of imptor: %i \n",anumbers[m1]);
       if (found && anumbers[m1]==8) { found = false; nimptor--; }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
//       make all 3 centered atoms planar?
//         printf(" atom: %i has coordn %i \n",imptor[nimptor-1][2],coordn[imptor[nimptor-1][2]]);
       }
       if (found) imptorv[nimptor-1] = imptorvt;
    }
  } 

  return;
}

double ICoord::torsion_val(int i, int j, int k, int l)
{
  double tval = -999;

  double x1 = coords[3*j+0] - coords[3*i+0];
  double y1 = coords[3*j+1] - coords[3*i+1];
  double z1 = coords[3*j+2] - coords[3*i+2];
  double x2 = coords[3*k+0] - coords[3*j+0];
  double y2 = coords[3*k+1] - coords[3*j+1];
  double z2 = coords[3*k+2] - coords[3*j+2];
  
  double ux1 = y1*z2-z1*y2;
  double uy1 = z1*x2-x1*z2;
  double uz1 = x1*y2-y1*x2;

  double x3 = coords[3*l+0] - coords[3*k+0];
  double y3 = coords[3*l+1] - coords[3*k+1];
  double z3 = coords[3*l+2] - coords[3*k+2];

  double ux2 = z3*y2 - y3*z2;
  double uy2 = x3*z2 - z3*x2;
  double uz2 = y3*x2 - x3*y2;

  double u = (ux1*ux1+uy1*uy1+uz1*uz1)*(ux2*ux2+uy2*uy2+uz2*uz2);

  if (u!=0.0)
  {
     double a = (ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u);
     if (a>1) a=1; else if (a<-1) a=-1;
     tval = acos(a);
     if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
         uz1*(ux2*y2-uy2*x2) < 0.0) tval *=-1;
  }

  return tval * 180/3.14;
}

double ICoord::torsion_val_v01(int i, int j)
{
  double tval = -999;

  double x1 = coords[3*j+0] - coords[3*i+0];
  double y1 = coords[3*j+1] - coords[3*i+1];
  double z1 = coords[3*j+2] - coords[3*i+2];
  double x2 = 0.0 - coords[3*j+0];
  double y2 = 0.0 - coords[3*j+1];
  double z2 = 0.0 - coords[3*j+2];
  
  double ux1 = y1*z2-z1*y2;
  double uy1 = z1*x2-x1*z2;
  double uz1 = x1*y2-y1*x2;

  double x3 = 0.0 - 0.0;
  double y3 = 1.0 - 0.0;
  double z3 = 0.0 - 0.0;

  double ux2 = z3*y2 - y3*z2;
  double uy2 = x3*z2 - z3*x2;
  double uz2 = y3*x2 - x3*y2;

  double u = (ux1*ux1+uy1*uy1+uz1*uz1)*(ux2*ux2+uy2*uy2+uz2*uz2);

  if (u!=0.0)
  {
     double a = (ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u);
     if (a>1) a=1; else if (a<-1) a=-1;
     tval = acos(a);
     if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
         uz1*(ux2*y2-uy2*x2) < 0.0) tval *=-1;
  }

  return tval * 180/3.14;
}

double ICoord::angle_val(int i, int j, int k)
{
   double D1 = distance(i,j);
   double D2 = distance(j,k);
   double D3 = distance(i,k);
   
   double cos = ( D1*D1 + D2*D2 - D3*D3 ) / ( 2*D1*D2);
 
   if (cos > 1) cos = 1;
   if (cos < -1) cos = -1;

  // printf(" cos is: %f \n",cos);
 
   return acos(cos) * 180/3.14;
}

double ICoord::angle_val_v0(int i, int j)
{
   double D1 = distance(i,j);
   double D2 = distance_v0(j);
   double D3 = distance_v0(i);
   
   double cos = ( D1*D1 + D2*D2 - D3*D3 ) / ( 2*D1*D2);
 
   if (cos > 1) cos = 1;
   if (cos < -1) cos = -1;

  // printf(" cos is: %f \n",cos);
 
   return acos(cos) * 180/3.14;
}

int ICoord::make_nonbond(){

  int n = 0;
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    {
      if (oatomtm[i]==-1 && oatomtm[j]==-1)
      {
        bool found = false;
        for (int k=0;k<nbonds;k++)
        {
           if (found) break;
           if ((bonds[k][0]==i && bonds[k][1]==j) ||
               (bonds[k][0]==j && bonds[k][1]==i)) found = true;
        }
        //printf(" checking for pair: %i %i \n",i,j);
        for (int k=0;k<nangles;k++)
        {
          if (found) break;
          //printf(" angle %i bonds: %i %i %i \n",k,angles[k][0],angles[k][1],angles[k][2]);
          if (angles[k][0]==i)
          {
             if (angles[k][1]==j) found = true;
             else if (angles[k][2]==j) found = true;
          }
          else if (angles[k][1]==i)
          {
             if (angles[k][0]==j) found = true;
             else if (angles[k][2]==j) found = true;
          }
          else if (angles[k][2]==i)
          {
             if (angles[k][0]==j) found = true;
             else if (angles[k][1]==j) found = true;
          }
        } // loop k over angles
        if (!found)
        {
          //printf(" not found\n");
          nonbondd[n] = distance(i,j);
          nonbond[n][0] = i;
          nonbond[n][1] = j;
          n++;
        } // if !found
      } // if not identical
    }
  }
  //printf(" n_nonbond: %i \n",n);

  return n;
}

//returns whether atom a is bonded to another H
int ICoord::isH2(int a) {

  if (anumbers[a]!=1) return 0;
  for (int i=0;i<nbonds;i++)
    if (bonds[i][0]==a)
    {
      if (anumbers[bonds[i][1]]==1)
        return 1;
    }
    else if (bonds[i][1]==a)
    {
      if (anumbers[bonds[i][0]]==1)
        return 1;
    }
  return 0;
}

int ICoord::isTM(int a) {

//may later be extended to all 5+ coord types
  int anum;
  if (a>-1)
    anum = anumbers[a];
  else
    return 0;

  int TM = 0;
  if (anum > 1000)
    TM = 2;
  else if (anum > 20)
  {
    if (anum < 31)
      TM = 1;
    else if (38 < anum && anum < 49)
      TM = 1;
    else if (71 < anum && anum < 81)
      TM = 1;
  }

  return TM;
}


double ICoord::getR(int i){

  double value;
 
  if (i < 0) return 1.0;
  int anum = anumbers[i];
  if (anum > 1000) anum -= 1000;
  //printf(" getR: %i,%i \n",anum,i);
  if      (anum==1) value = 1.3;
  else if (anum==5) value = 1.7; //was 1.65, but changed due to N-B bond of AB
  else if (anum==6) value = 1.7;
  else if (anum==7) value = 1.7; //same as #5
  else if (anum==8) value = 1.65;
  else if (anum==9) value = 1.6;
  else if (anum==13) value = 2.6;
  else if (anum==14) value = 2.6;
  else if (anum==15) value = 2.5;
  else if (anum==16) value = 2.3;
  else if (anum==17) value = 2.1;
  else if (anum==27) value = 3.0;
  else if (anum==28) value = 3.0;
  else if (anum==46) value = 3.15;
  else if (anum==78) value = 3.35;
//  else if (anum>1000) value = tmRval;

  return value;
}

double ICoord::distance(int i, int j)
{
  //printf("in distance: %i %i\n",i+1,j+1);
  return sqrt((coords[3*i+0]-coords[3*j+0])*(coords[3*i+0]-coords[3*j+0])+
              (coords[3*i+1]-coords[3*j+1])*(coords[3*i+1]-coords[3*j+1])+
              (coords[3*i+2]-coords[3*j+2])*(coords[3*i+2]-coords[3*j+2])); 
}

double ICoord::distance_v0(int i)
{
  //printf("in distance: %i %i\n",i);
  return sqrt((coords[3*i+0])*(coords[3*i+0])+
              (coords[3*i+1])*(coords[3*i+1])+
              (coords[3*i+2])*(coords[3*i+2])); 
}

int ICoord::bond_exists(int b1, int b2) {

#if 0
   if (oatomtm[b1]>-1)
     if (b2==oatomtm[b1])
       printf(" ");
#endif

   int found = 0;
   if (bond_num(b1,b2)>-1)
     found = 1;

   return found;
}

int ICoord::bond_exists_tri(int b1, int b2) {

   int found = 0;
   if (bond_num(b1,b2)>-1)
     found = 1;

   // prevent triangle from reattaching 
   for (int i=0;i<ntri;i++)
     if ( (tmtriangles[i][7]==b1 && tmtriangles[i][8]==b2)
     ||   (tmtriangles[i][7]==b2 && tmtriangles[i][8]==b1) )
       found = 1;

   return found;
}

int ICoord::bond_num(int b1, int b2) {

   int found = -1;

   for (int k1=0;k1<nbonds;k1++)
     if ( (bonds[k1][0] == b1 && bonds[k1][1] == b2)
       || (bonds[k1][1] == b1 && bonds[k1][0] == b2))
     {
       found = k1;
       break;
     }

   return found;
}

int ICoord::hpair(int a1, int a2) {
  if (anumbers[a1]==1 && anumbers[a2]==1)
    return 1;
  else
    return 0;
}

int ICoord::isTMpair(int a1, int a2) {
  if (oatomtm[a1]==a2 || oatomtm[a2]==a1)
    return 1;
  else
    return 0;
}

int ICoord::h2count() {

  int count = 0;
  for (int i=0;i<nbonds;i++)
  {
    if (anumbers[bonds[i][0]]==1 && anumbers[bonds[i][1]]==1)
      count++;
  }

  return count;
}




void ICoord::structure_read(string xyzfile){ 
   
 // cout <<" Reading and initializing string coordinates" << endl;
 // cout <<"  -Opening structure file" << endl;
  
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    cout << "!!!!Error opening xyz file!!!!" << endl;
    exit(-1);
  } 
  
 // cout <<"  -reading file..." << endl;
  
  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
 // cout <<"  natoms: " << natoms << endl;
  
  success=getline(infile, line);
  if (success){  
    comment=line;
  }
  
  anumbers = new int[4+natoms];
  amasses = new double[4+natoms];
  anames = new string[4+natoms];
    
  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    anumbers[i]=PTable::atom_number(anames[i]);
    amasses[i]=PTable::atom_mass(anumbers[i]);
  }
  
  infile.close();
  
  coords = new double[(natoms+4)*3];
  coords0 = new double[(natoms+4)*3];
  coordsts = new double[(natoms+4)*3];
   
  //cout <<"  -Reading coordinates...";
 // cout << "Opening the xyz file" << endl;
  infile.open(xyzfile.c_str());
  fflush(stdout);
 // cout << "xyzfile opened" << endl;
  fflush(stdout);
  
  
//  for (int i=1;i<=2;i++){
    success=getline(infile, line);
    success=getline(infile, line);
    for (int j=0;j<natoms;j++){
      success=getline(infile, line);
      int length=StringTools::cleanstring(line);
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      coords[3*j+0]=atof(tok_line[1].c_str());
      coords[3*j+1]=atof(tok_line[2].c_str());
      coords[3*j+2]=atof(tok_line[3].c_str());
    
    }
//  }
  
  for (int i=0;i<3*natoms;i++)
     coords0[i] = coords[i];
   
 // cout << " done" << endl;
  infile.close();
  
 // cout << "Finished reading information from structure file" << endl;
}   




#if 0
        if (isTM(angles[nangles-1][1]) && tmnoangles>0)
        {
          for (int k=0;k<ntmnoangles;k++)
          {
            if (tmnoangles[k][0]==angles[nangles-1][0] && tmnoangles[k][2]==angles[nangles-1][2])
            { nangles--; break; }
            else if (tmnoangles[k][2]==angles[nangles-1][0] && tmnoangles[k][2]==angles[nangles-1][0])
            { nangles--; break; }
          }
        }
#endif
