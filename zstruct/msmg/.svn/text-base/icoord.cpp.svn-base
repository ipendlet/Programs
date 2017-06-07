#include "icoord.h"
#include "utils.h"
using namespace std;


int ICoord::init(string xyzfile){


// printf(" xyzfile: %s \n",xyzfile);
 printf("\n");
 cout << " xyzfile: " << xyzfile << endl;
 structure_read(xyzfile);
 
 print_xyz();

 alloc_mem();
 // printf(" done allocating memory\n");

 int done = ic_create();

 printf(" initializing MM parameters \n");
 mm_init();

 //printf("\n\n");

 return 1;
}



// initialize by feeding in xyz coordinates
int ICoord::init(int nat, string* anam, int* anum, double* xyz){

// printf(" initializing icoord via xyz structure \n");
 natoms = nat;
// printf(" natoms: %i \n",nat);
// for (int i=0;i<natoms;i++)
//    printf(" %1.3f %1.3f %1.3f \n",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);

//otherwise allocated in structure_read
 anumbers = new int[1+natoms];
 amasses = new double[1+natoms];
 anames = new string[1+natoms];
 coords = new double[natoms*3];
 coordsts = new double[natoms*3];
 coords0 = new double[natoms*3];

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

 alloc_mem();

 int done = ic_create();

// printf(" initializing MM parameters \n");
 mm_init();

 //printf("\n\n");

 return 1;
}

// initialize memory only
int ICoord::alloc(int size){

 natoms = size;

//otherwise allocated in structure_read
 anumbers = new int[1+natoms];
 amasses = new double[1+natoms];
 anames = new string[1+natoms];
 coords = new double[natoms*3];
 coordsts = new double[natoms*3];
 coords0 = new double[natoms*3];

 alloc_mem();

 return 1;
}


// initialize by feeding in xyz coordinates
int ICoord::reset(double* xyz){

// printf(" resetting icoord via xyz structure \n");

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

// int done = ic_create();

// printf(" initializing MM parameters \n");
// mm_init();

 //printf("\n\n");

 return 1;
}

// initialize by feeding in xyz coordinates
int ICoord::reset(int nat, string* anam, int* anum, double* xyz){

// printf(" resetting icoord via xyz structure \n");
 natoms = nat;

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

// int done = ic_create();

// printf(" initializing MM parameters \n");
// mm_init();

 //printf("\n\n");

 return 1;
}

void ICoord::update_ic(){

  update_bonds();
  update_angles();
  update_torsion();
  update_imptor();
  update_nonbond();

  return;
} 
 
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
 
  delete [] adone;

  return;
}


int ICoord::ic_create()
{
//  printf(" Creating internals from xyz \n");
  make_bonds();
  coord_num(); // counts # surrounding species
//  printf(" now making angles \n");
  make_angles();
//  printf(" now making torsions \n");
  make_torsions();
//  printf(" now making improper torsions \n");
  make_imptor();

//  printf(" now counting nonbond\n");
  n_nonbond = make_nonbond(); //anything not connected by bond or angle

//  print_ic();
}

int ICoord::ic_create_nobonds()
{
//  printf(" Creating internals from xyz, skipping bond making \n");
  coord_num(); // counts # surrounding species
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
  //printf(" in make_bonds, natoms: %i\n",natoms);
  double MAX_BOND_DIST; 
  nbonds=0;
  for (int i=0;i<natoms;i++)
    for (int j=0;j<i;j++)
    {
       MAX_BOND_DIST = (getR(i) + getR(j))/2;
       double d = distance(i,j);
       if (d<MAX_BOND_DIST)
       {
          //printf(" found bond: %i %i dist: %f \n",i+1,j+1,d);
          bonds[nbonds][0]=i;
          bonds[nbonds][1]=j;
          bondd[nbonds]=d;
          nbonds++;
       }
    }

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
  //include all consecutive connections 
  nangles=0;
  for (int i=0;i<nbonds;i++)
  {
     for (int j=0;j<i;j++)
     {
        if (bonds[i][0]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][1];
          nangles++;
        }
        else if (bonds[i][0]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][0];
          nangles++;
        }
        else if (bonds[i][1]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][1];
          nangles++;
        }
        else if (bonds[i][1]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][0];
          nangles++;
        }
        if (nangles>0)
          anglev[nangles-1]=angle_val(angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
     } //loop j
  } //loop i


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

       if (a1==a2 && b1==b2)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=a1;
          torsions[ntor][2]=b1;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       else if (b1==b2 && c1==c2)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=c1;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
       else if (a1==b2 && b1==c2)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=c2;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       else if (c1==b2 && b1==c2)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=c2;
          torsions[ntor][2]=c1;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
/*       else if (a1==b2 && b1==a2)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=a2;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }*/
       if (found) torv[ntor-1]=torsion_val(torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
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

       if (m1==m2)
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
//FIX ME
//the following only works when center is 3 coordinate
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
         if ((abs(imptorvt) < 90.0 || abs(imptorvt) > 120.0) && anumbers[imptor[nimptor-1][2]]==7) { found = false; nimptor--; }
         else if (abs(imptorvt) > 12.0 && abs(imptorvt - 180.) > 12.0) { found = false; nimptor--; }
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
  nimptor = 0;
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

       if (m1==m2)
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
//FIX ME
//the following only works when center is 3 coordinate
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
//       make all 3 centered atoms planar?
//         printf(" atom: %i has coordn %i \n",imptor[nimptor-1][2],coordn[imptor[nimptor-1][2]]);
         if (coordn[imptor[nimptor-1][2]]!=3)
         {
           found = false;
           nimptor--;
         }
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

int ICoord::make_nonbond(){

  int n = 0;
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
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
      }
    }
  }
  //printf(" n_nonbond: %i \n",n);

  return n;
}

double ICoord::getR(int i){

  double value;
 
  if	  (anumbers[i]==1) value = 1.3;
  else if (anumbers[i]==5) value = 1.7;	//was 1.65, but	changed	due to N-B bond	of AB
  else if (anumbers[i]==6) value = 1.65;
  else if (anumbers[i]==7) value = 1.7;	//same as #5
  else if (anumbers[i]==8) value = 1.65;
  else if (anumbers[i]==9) value = 1.6;
  else if (anumbers[i]==13) value = 2.6;
  else if (anumbers[i]==14) value = 2.6;
  else if (anumbers[i]==15) value = 2.5;
  else if (anumbers[i]==16) value = 2.3;
  else if (anumbers[i]==17) value = 2.1;


  return value;
}

double ICoord::distance(int i, int j)
{
  //printf("in distance: %i %i\n",i+1,j+1);
  return sqrt((coords[3*i+0]-coords[3*j+0])*(coords[3*i+0]-coords[3*j+0])+
              (coords[3*i+1]-coords[3*j+1])*(coords[3*i+1]-coords[3*j+1])+
              (coords[3*i+2]-coords[3*j+2])*(coords[3*i+2]-coords[3*j+2])); 
}

int ICoord::bond_exists(int b1, int b2) {

   int found = 0;
   if (bond_num(b1,b2)>-1)
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
  cout <<"  natoms: " << natoms << endl;
  
  success=getline(infile, line);
  if (success){  
    comment=line;
  }
  
  anumbers = new int[1+natoms];
  amasses = new double[1+natoms];
  anames = new string[1+natoms];
    
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
  
//  V_profile = new double[1+nnmax];
//  S = new double[1+nnmax];
  
  coords = new double[natoms*3];
  coords0 = new double[natoms*3];
   
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

