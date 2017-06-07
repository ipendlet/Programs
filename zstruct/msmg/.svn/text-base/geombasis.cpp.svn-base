#include "geombasis.h"
#include "icoord.h"
#include "utils.h"
using namespace std;


int GeomBasis::init(string xyzlist){

 printf(" test of GeomBasis, reading in xyz for ICoord \n");
 natomt = 0;

 printf("\n");
 cout << " xyzlist: " << xyzlist << endl;

 list_read(xyzlist);
 icoord = new ICoord[nbasis];

 printf(" files: \n");
 for (int i=0;i<nbasis;i++)
   cout << " " << xyzfile[i] << endl;
 for (int i=0;i<nbasis;i++)
   icoord[i].init(xyzfile[i]);

 printf("\n\n **now making atom centered types** \n");

 ats = new Atom[100];
 Atom att; att.alloc(50); // temporary atom type
// for (int i=0;i<icoord[0].natoms;i++)
//    printf(" %1.3f %1.3f %1.3f \n",icoord[0].coords[3*i+0],icoord[0].coords[3*i+1],icoord[0].coords[3*i+2]);

//temporary arrays 
 int natt;
 double* txyz = new double[150];
 string* tanam = new string[50];
 int* tanum = new int[50];
 int found;

// for (int i=0;i<1;i++)
 for (int i=0;i<nbasis;i++)
 {
    for (int j=0;j<icoord[i].natoms;j++)
    {
      //adding first atom
      txyz[3*0+0] = icoord[i].coords[3*j+0];
      txyz[3*0+1] = icoord[i].coords[3*j+1];
      txyz[3*0+2] = icoord[i].coords[3*j+2];
      tanam[0] = icoord[i].anames[j];
      tanum[0] = icoord[i].anumbers[j];
      int k = 1; // CPMZ CHECK THIS, was 1
      //search for all atoms attached by bonds
      for (int l=0;l<icoord[i].nbonds;l++)
      {
        int b1 = icoord[i].bonds[l][0];
        int b2 = icoord[i].bonds[l][1];
        int bn;
        if (b1==j || b2==j)
        {
//          printf(" atom j: %i has bond: %i %i \n",j,icoord[i].bonds[l][0],icoord[i].bonds[l][1]);
          if (b1==j) bn=b2; else bn=b1;
          txyz[3*k+0] = icoord[i].coords[3*bn+0];
          txyz[3*k+1] = icoord[i].coords[3*bn+1];
          txyz[3*k+2] = icoord[i].coords[3*bn+2];
          tanam[k] = icoord[i].anames[bn];
          tanum[k] = icoord[i].anumbers[bn];
          k++;
        }
      } // loop over nbonds
 
      printf("\n new atom coordination #: %i \n",k-1);

      natt = k;
      att.reset(natt,tanam,tanum,txyz);
      att.lic.ic_create();
      found = compare(att);
      if (found==0)
      {
        ats[natomt].init(natt,tanam,tanum,txyz);
        natomt++;
      }
    } // loop j over atoms in ic
 } // loop i over ic's

 printf("\n\n");

 printf(" found %i types of atoms \n",natomt);
 for (int i=0;i<natomt;i++)
   ats[i].lic.print_xyz();

 delete [] txyz;
 delete [] tanam;
 delete [] tanum;

 return 1;
}

int GeomBasis::find(int j, ICoord icoord){

  printf("\n locating atom %i \n",j);
//  icoord.print_xyz();

  int type = -1;

//temporary arrays 
  int natt;
  double* txyz = new double[150];
  string* tanam = new string[50];
  int* tanum = new int[50];
  int found;

      //adding first atom
      txyz[3*0+0] = icoord.coords[3*j+0];
      txyz[3*0+1] = icoord.coords[3*j+1];
      txyz[3*0+2] = icoord.coords[3*j+2];
      tanam[0] = icoord.anames[j];
      tanum[0] = icoord.anumbers[j];
      int k = 1;
      //search for all atoms attached by bonds
      for (int l=0;l<icoord.nbonds;l++)
      {
        int b1 = icoord.bonds[l][0];
        int b2 = icoord.bonds[l][1];
        int bn;
        if (b1==j || b2==j)
        {
//          printf(" atom j: %i has bond: %i %i \n",j,icoord.bonds[l][0],icoord.bonds[l][1]);
          if (b1==j) bn=b2; else bn=b1;
          txyz[3*k+0] = icoord.coords[3*bn+0];
          txyz[3*k+1] = icoord.coords[3*bn+1];
          txyz[3*k+2] = icoord.coords[3*bn+2];
          tanam[k] = icoord.anames[bn];
          tanum[k] = icoord.anumbers[bn];
          k++;
        }
      } // loop over nbonds
 
      natt = k;
      printf(" atom coordination #: %i \n",k-1);

      Atom att; att.alloc(50); // temporary atom type
      att.reset(natt,tanam,tanum,txyz);
      att.lic.ic_create();
      type = locate(att);
      if (type==-1) printf(" ERROR: couldn't find! \n");

  return type;
}

int GeomBasis::compare(Atom atom){

  printf(" comparing atom to rest of list \n");
 
  int found = 0;

  for (int i=0;i<natomt;i++)
  {
    found = compare1(atom,ats[i]);
    if (found) { printf(" atom is same as %i \n",i); break; }
  }

  return found;
}



int GeomBasis::nplus(int n, ICoord icoord, int atype0){

  int type = -1;

//  printf(" nplus for %i %i \n",n,atype0);

// currently allow hydrogen to be 2 coordinate only if attached to B
//  if (icoord.coordn[n]==1 && icoord.anumbers[n]==1)
//    return -1;

  int coordp1 = icoord.coordn[n]+1;
  for (int i=0;i<natomt;i++)
  {
    int anum = ats[i].lic.anumbers[0];
    int coordni = ats[i].lic.coordn[0];
//    printf(" in nplus loop atombasis %i, anum: %i coordni: %i \n",i,anum,coordni);
    if (anum==icoord.anumbers[n])
    {
      if (coordp1==coordni && icoord.anumbers[n]!=1)
      {
        //printf(" found nplus type: %i \n",i);
        type = i;
        break;
      }
#if 1
      else if (coordp1==coordni && icoord.anumbers[n]==1)
      {
//detecting H attached to B (extend later)
        int hb = -1;
        for (int j=0;j<icoord.nbonds;j++)
        {
          if (icoord.bonds[j][0]==n && icoord.anumbers[icoord.bonds[j][1]]==5)
            hb = icoord.bonds[j][1];
          else if (icoord.bonds[j][1]==n && icoord.anumbers[icoord.bonds[j][0]]==5)
            hb = icoord.bonds[j][0];
          if (hb>-1)
          {
            printf(" found H (atom %i), attached to: %i \n",n,hb);
            type = hb;
            break;
          }
        } //loop j over nbonds
      } // if hydrogen
#endif
    } //if atom found
  } // loop over natomtypes

#if 0
  if (type>=0)
  {
    printf(" printing atom type %i: \n",type);
    ats[type].lic.print_xyz();
  }
#endif

  return type;
}

int GeomBasis::nplus2(int n, ICoord icoord, int atype0){

  int type = -1;

// currently don't allow hydrogen to be 3 coordinate
  if (icoord.coordn[n]==1 && icoord.anumbers[n]==1)
    return -1;

  int coordp2 = icoord.coordn[n]+2;
  for (int i=0;i<natomt;i++)
  {
    int anum = ats[i].lic.anumbers[0];
    int coordni = ats[i].lic.coordn[0];
    if (anum==icoord.anumbers[n])
    {
      if (coordp2==coordni)
      {
//        printf(" found nplus2 type: %i \n",i);
        type = i;
        break;
      }
    }
  }

#if 0
  if (type>=0)
  {
    printf(" printing atom type %i: \n",type);
    ats[type].lic.print_xyz();
  }
#endif

  return type;
}

int GeomBasis::nminus(int n, ICoord icoord, int atype0){

  int type = -1;
  if (icoord.coordn[n]<2)
    return -1;

  int coordm1 = icoord.coordn[n]-1;
  for (int i=0;i<natomt;i++)
  {
    int anum = ats[i].lic.anumbers[0];
    int coordni = ats[i].lic.coordn[0];
    if (anum==icoord.anumbers[n])
    {
      if (coordm1==coordni)
      {
//        printf(" found nminus type: %i \n",i);
        type = i;
        break;
      }
    }
  }

#if 0
  if (type>=0)
  {
    printf(" printing atom type %i: \n",type);
    ats[type].lic.print_xyz();
  }
#endif

  return type;
}

int GeomBasis::nminus2(int n, ICoord icoord, int atype0){

  int type = -1;
  if (icoord.coordn[n]<3)
    return -1;

  int coordm2 = icoord.coordn[n]-2;
  for (int i=0;i<natomt;i++)
  {
    int anum = ats[i].lic.anumbers[0];
    int coordni = ats[i].lic.coordn[0];
    if (anum==icoord.anumbers[n])
    {
      if (coordm2==coordni)
      {
//        printf(" found nminus2 type: %i \n",i);
        type = i;
        break;
      }
    }
  }

#if 0
  if (type>=0)
  {
    printf(" printing atom type %i: \n",type);
    ats[type].lic.print_xyz();
  }
#endif

  return type;
}



int GeomBasis::locate(Atom atom){

//  printf(" locating atom type \n");
 
  int type = 0;

  for (int i=0;i<natomt;i++)
  {
    type = compare1(atom,ats[i]);
    if (type)
    {
     // printf(" atom is same as %i \n",i);
      type=i;
      break; 
    }
  }

  printf(" found atom type %i \n",type);
//  ats[type].lic.print_xyz();

  return type;
}

int GeomBasis::compare1(Atom atom1, Atom atom2){

  if (atom1.coordn != atom2.coordn)
    return 0;
  else if (atom1.lic.nimptor != atom2.lic.nimptor)
    return 0;
  else if (atom1.lic.anames[0] != atom2.lic.anames[0])
    return 0;

  return 1;
}



void GeomBasis::freemem(){

//  for (int i=0;i<nbasis;i++)
//    delete [] icoord[i];

  return;
}


void GeomBasis::list_read(string xyzlist){ 
   
  cout <<" Reading in list of basis structures" << endl;
  
  ifstream infile;
  infile.open(xyzlist.c_str());
  if (!infile){
    cout << "!!!!Error opening list!!!!" << endl;
    exit(-1);
  } 
  
  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    nbasis=atoi(line.c_str());
  }
  cout <<"  nbasis: " << nbasis << endl;

  xyzfile = new string[nbasis];
  
  //cout <<"  -Reading the file names";
  for (int i=0;i<nbasis;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    xyzfile[i]=tok_line[0];
  }
  
  infile.close();
  
  cout << "Finished reading list " << endl;
}   

