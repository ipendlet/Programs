#include "mopac.h"
//#include "utils.h"
using namespace std;
#include "constants.h"


void Mopac::write_ic_input(ofstream& inpfile, int anum, ICoord ic){

//  printf(" in write_ic_input() for atom %i \n",anum);

//  printf(" nbonds: %i \n",ic.nbonds);
 

  for (int i=0;i<nfrz0;i++)
  {
    if (anum==frzlist[i])
    {
      //printf(" atom %i was moved, doing xyz \n",anum);
      inpfile << " " << anames[anum] << " " << xyz[3*anum+0] << " 0 " << xyz[3*anum+1] << " 0 " << xyz[3*anum+2] << " 0 " << endl;
      return;
    }
  }

//  printf(" atom attached to moved atom, freezing bond \n"); 

  int b1,c1,d1;
  for (int i=0;i<ic.nbonds;i++)
  {
    b1 = -1;
    if (ic.bonds[i][0] == anum || ic.bonds[i][1] == anum)
    for (int j=0;j<nfrz0;j++)
    {
      //printf(" found %i, checking for bond partner %i \n",anum,frzlist[j]);
      if (ic.bonds[i][0]==frzlist[j]
        || ic.bonds[i][1]==frzlist[j])
      {
        b1 = frzlist[j];
        break;
      }
    }
// this could additionally check if b1 has another attachment,
//  but this case is handled by dummies anyway
//    if (ic.coordn[b1]<2)
//      b1 = -1;

    if (b1>-1)
    {
      //printf(" found attachment: %i %i \n",b1,anum);
      break;
    }
  }

//could change this to "xyz list" if statement
  c1 = -1; d1 = -1;
  for (int i=0;i<natoms;i++)
  {
    if (c1==-1)
    for (int j=0;j<i;j++)
    {
      if (!frzlistb[i] && !frzlistb[j])
      {
        c1 = i; d1 = j;
        //printf(" found angle/tor possibility: %i %i \n",c1,d1);
        break;
      }
    }
  }

  if (anum < 1)
  {
    //printf(" first atom %i, writing xyz \n",anum);
    inpfile << " " << anames[anum] << " " << xyz[3*anum+0] << " 0 " << xyz[3*anum+1] << " 0 " << xyz[3*anum+2] << " 0 " << endl;
  }
  else if (b1==-1 || c1==-1 || d1==-1 || (anum>b1 || anum>c1 || anum>d1))
  {
    //printf(" trying to find additional IC for atom %i \n",anum);
    if (ic.coordn[anum]==0)
    {
      printf(" failed to find any IC for atom %i, writing xyz \n",anum);
      inpfile << " " << anames[anum] << " " << xyz[3*anum+0] << " 0 " << xyz[3*anum+1] << " 0 " << xyz[3*anum+2] << " 0 " << endl;
    }
    else
    {
      //printf(" using dummies \n");
      double rv,anglev,torv;
      rv = ic.distance(anum,b1);
      anglev = ic.angle_val_v0(anum,b1);
      torv = ic.torsion_val_v01(anum,b1);
      //printf(" d r: %1.3f angle: %1.3f tor: %1.3f (bond %i %i) \n",rv,anglev,torv,anum,b1);
      inpfile << " " << anames[anum] << " " << rv << " 0 " << anglev << " 1 " << torv << " 1 " << b1+1+3 << " " << 1 << " " << 2 << endl;
    }

  }
  else
  {
    //printf(" found actual IC connection for atom %i \n",anum);
    double rv,anglev,torv;
    rv = ic.distance(anum,b1);
    anglev = ic.angle_val(anum,b1,c1);
    torv = ic.torsion_val(anum,b1,c1,d1);
    //printf(" n r: %1.3f angle: %1.3f tor: %1.3f \n",rv,anglev,torv);
    inpfile << " " << anames[anum] << " " << rv << " 0 " << anglev << " 1 " << torv << " 1 " << b1+1+3 << " " << c1+1+3 << " " << d1+1+3 << endl;
  }


  return;
}

double Mopac::opt() {

  string filename = "scratch/testmopac.mop";

  energy = opt(filename);
 
  return energy;
}

void Mopac::opt_write() {
      
  string filename = "scratch/testmopac.mop";
      
  opt_write(filename);
   
  return;
}


void Mopac::opt_header(ofstream& inpfile) {

//  inpfile << " PM7 C.I.=2 " << endl; 
//  inpfile << " PM7 NOSYM GNORM=0.7 " << endl;
#if NOMOOPT
  inpfile << " PM6 NOSYM 1SCF " << endl;
#else

#if !UNRESTRICTED
#if !CATION
  inpfile << " PM6 NOSYM GNORM=4.5 " << endl;
#else
  inpfile << " PM6 NOSYM GNORM=4.5 CHARGE=1 " << endl;
#endif
#else
#if !CATION 
#if !TRIPLET
  inpfile << " PM6 NOSYM GNORM=4.5 UHF " << endl;
#else
  inpfile << " PM6 NOSYM GNORM=4.5 UHF TRIPLET" << endl;
#endif
#else
  inpfile << " PM6 NOSYM GNORM=4.5 UHF CHARGE=1 " << endl;
#endif
#endif

#endif

#if 0

#if !UNRESTRICTED
#if !CATION
  inpfile << " PM6 NOSYM GNORM=4.5 " << endl;
#else
  inpfile << " PM6 NOSYM GNORM=4.5 CHARGE=1 " << endl;
#endif
#else
#if !CATION 
#if !TRIPLET
  inpfile << " PM6 NOSYM GNORM=4.5 UHF " << endl;
#else
  inpfile << " PM6 NOSYM GNORM=4.5 UHF TRIPLET" << endl;
#endif
#else
  inpfile << " PM6 NOSYM GNORM=4.5 UHF CHARGE=1 " << endl;
#endif
#endif

#endif
//  inpfile << " PM7 NOSYM GNORM=1.5 CHARGE=1 " << endl;
//  inpfile << " PM6-DH+ NOSYM GNORM=1.5" << endl;
  inpfile << "   MOPAC run " << endl;
  inpfile << "   ready to go? " << endl;

  return;
}

double Mopac::opt(string filename, ICoord icoords) {

#if 0
  printf(" WARNING: bypassing ICs! \n");
  return opt(filename);
#endif

//  printf(" in mopac/opt() w/IC freeze \n");
//  for (int i=0;i<natoms;i++)
//    printf(" %s %1.3f %1.3f %1.3f \n",anames[i].c_str(),xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);
//    printf(" %s %1.3f %1.3f %1.3f \n",anames[i].c_str(),ics.coords[3*i+0],ics.coords[3*i+1],ics.coords[3*i+2]);

  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
#if NOMOOPT
  xyzfile << " " << natoms-nsolvent << endl << endl;
#else
  xyzfile << " " << natoms << endl << endl;
#endif
  opt_header(inpfile);

  inpfile << " X 0.0 0 0.0 0 0.0 0 " << endl;
  inpfile << " X 0.0 0 1.0 0 0.0 0 " << endl;
  inpfile << " X 0.0 0 0.0 0 1.0 0 " << endl;

#if !NOMOOPT
  for (int i=0;i<natoms;i++)
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
    {
//CPMZ instead writing anames in write_ic_input function
//      inpfile << " " << anames[i] << " ";
      write_ic_input(inpfile,i,icoords);
    }
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
#else
  for (int i=0;i<natoms;i++)
  {
    if (!solvent[i])
    {
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
      xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
    }
  }
#endif

  string cmd = "/tmp/MOPAC2012.exe "+filename;
  system(cmd.c_str());
#endif

  energy = read_output(filename);
 
  printf(" initial energy: %1.2f final energy: %1.2f \n",energy0,energy); 

  // need to retrieve final geometry, write to xyz
#if !NOMOOPT
  xyz_read(inpfile_string);
  xyz_save(inpfile_string+".xyz");
#endif

  xyzfile.close();
  inpfile.close();

  if (abs(energy)<0.00001)
  {
    printf(" energy zero, mopac failed \n");
    return 10000;
  }

  return energy;
}


void Mopac::opt_write(string filename, ICoord icoords) {

#if 0
  printf(" WARNING: bypassing ICs! \n");
  return opt_write(filename);
#endif

#if 0
  printf(" nsolvent: %i \n",nsolvent);
  printf(" solvent: ");
  for (int i=0;i<natoms;i++)
    printf("%i ",solvent[i]);
  printf("\n");
#endif


  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
#if NOMOOPT
  xyzfile << " " << natoms-nsolvent << endl << endl;
#else
  xyzfile << " " << natoms << endl << endl;
#endif

  opt_header(inpfile);

  inpfile << " X 0.0 0 0.0 0 0.0 0 " << endl;
  inpfile << " X 0.0 0 1.0 0 0.0 0 " << endl;
  inpfile << " X 0.0 0 0.0 0 1.0 0 " << endl;

#if !NOMOOPT
  for (int i=0;i<natoms;i++)
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
    {
      write_ic_input(inpfile,i,icoords);
    }
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
#else
  for (int i=0;i<natoms;i++)
  {
    if (!solvent[i])
    {
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
      xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
    }
  }
#endif
#endif

  return;
}


// standard opt
double Mopac::opt(string filename) {

   //printf(" in mopac/opt() \n");

  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
#if NOMOOPT
  xyzfile << " " << natoms-nsolvent << endl << endl;
#else
  xyzfile << " " << natoms << endl << endl;
#endif

  opt_header(inpfile);

#if !NOMOOPT
  for (int i=0;i<natoms;i++)
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 0 " << xyz[3*i+1] << " 0 " << xyz[3*i+2] << " 0 " << endl;
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
#else
  for (int i=0;i<natoms;i++)
  {
    if (!solvent[i])
    {
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
      xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
    }
  }
#endif
//  xyzfile << " " << natoms << endl << endl;

  string cmd = "/tmp/MOPAC2012.exe "+filename;
  system(cmd.c_str());
#endif

  energy = read_output(filename);
 
  printf(" initial energy: %1.4f final energy: %1.4f \n",energy0,energy); 

  // need to retrieve final geometry, write to xyz
#if !NOMOOPT
  xyz_read(inpfile_string);
  xyz_save(inpfile_string+".xyz");
#endif

  xyzfile.close();
  inpfile.close();

  if (abs(energy)<0.00001)
  {
    printf(" energy zero, mopac failed \n");
    return 10000;
  }

  return energy;
}

// single point
double Mopac::energy_sp(string filename) {

   //printf(" in mopac/energy_sp() \n");

  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac! \n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
#if NOMOOPT
  xyzfile << " " << natoms-nsolvent << endl << endl;
#else
  xyzfile << " " << natoms << endl << endl;
#endif
  opt_header(inpfile);

#if !NOMOOPT
  for (int i=0;i<natoms;i++)
  {
    inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 0 " << xyz[3*i+1] << " 0 " << xyz[3*i+2] << " 0 " << endl;
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
#else
  for (int i=0;i<natoms;i++)
  {
    if (!solvent[i])
    {
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 0 " << xyz[3*i+1] << " 0 " << xyz[3*i+2] << " 0 " << endl;
      xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
    }
  }
#endif

  string cmd = "/tmp/MOPAC2012.exe "+filename;
  system(cmd.c_str());
#endif

  energy = read_output(filename);
 
  printf(" initial energy: %1.4f final energy: %1.4f \n",energy0,energy); 

  // need to retrieve final geometry, write to xyz
  xyz_save(inpfile_string+".xyz");

  xyzfile.close();
  inpfile.close();

  if (abs(energy)<0.00001)
  {
    printf(" energy zero, mopac failed \n");
    return 10000;
  }

  return energy;
}



void Mopac::opt_write(string filename) {
   //printf(" in mopac/opt() \n");

  energy0 = energy = 0;

#if SKIPMOPAC
  printf(" skipping mopac opt! \n");
#endif

#if 0
  printf(" nsolvent: %i \n",nsolvent);
  printf(" solvent: ");
  for (int i=0;i<natoms;i++)
    printf("%i ",solvent[i]);
  printf("\n");
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  ofstream xyzfile;
  string xyzfile_string = "scratch/testmopac.xyz";
#if !SKIPMOPAC
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);
#if NOMOOPT
  xyzfile << " " << natoms-nsolvent << endl << endl;
#else
  xyzfile << " " << natoms << endl << endl;
#endif

  opt_header(inpfile);

#if !NOMOOPT
  for (int i=0;i<natoms;i++)
  {
    if (!frzlistb[i])
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
    else
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 0 " << xyz[3*i+1] << " 0 " << xyz[3*i+2] << " 0 " << endl;
    xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
  }
#else
  for (int i=0;i<natoms;i++)
  {
    if (!solvent[i])
    {
      inpfile << " " << anames[i] << " " << xyz[3*i+0] << " 1 " << xyz[3*i+1] << " 1 " << xyz[3*i+2] << " 1 " << endl;
      xyzfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << endl;
    }
  }
#endif
#endif

  return;
}



double Mopac::read_output(string filename) {

  //double energy = -1;
  energy = -1;

  string oname = filename+".out";
  printf(" read_output in mopac: %s \n",oname.c_str());
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (line.find("CYCLE:     1")!=string::npos)
    {
      cout << "Initial: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
//      energy0=atof(tok_line[10].c_str());
    }
    if (line.find("FINAL")!=string::npos)
    {
      cout << "Final: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[5].c_str());
    }
  }

  return energy;
}

void Mopac::alloc(int natoms_i) {
 
//  printf(" in mopac/alloc() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  nfrz = 0;
  nsolvent = 0;

  frzlist = new int[natoms_i];
  frzlistb = new int[natoms_i];

  solvent = new int[natoms_i];

  return;
}

void Mopac::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in mopac/init() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  nfrz = 0;
  nsolvent = 0;

  frzlist = new int[natoms_i];
  frzlistb = new int[natoms_i];

  solvent = new int[natoms_i];

  return;
}

void Mopac::freemem() {

  delete [] xyz0;
  delete [] xyz;
  delete [] anumbers;
  delete [] anames;

  delete [] frzlist;
  delete [] frzlistb;

  return;
}

void Mopac::addSolvent(int* solvent0, int nsolvent0)
{
  nsolvent = nsolvent0;
  for (int i=0;i<natoms;i++)
    solvent[i] = solvent0[i];
#if 0
  printf(" solvent added to mopac: \n");
  for (int i=0;i<natoms;i++)
    printf("%i: %i \n",i,solvent[i]);
#endif

  return;
}

void Mopac::freeze(int* frzlist_new, int nfrz_new, int nfrz0_new) {

  nfrz0 = nfrz0_new;
  nfrz = nfrz_new;
  for (int i=0;i<natoms;i++)
    frzlistb[i] = 0;
  for (int i=0;i<nfrz;i++)
    frzlistb[frzlist_new[i]] = 1;
  for (int i=0;i<nfrz;i++)
    frzlist[i] = frzlist_new[i];

  printf(" freeze list: ");
  for (int i=0;i<nfrz;i++)
    printf("%i ",frzlist[i]);
  printf("\n");
  
  return;
}

void Mopac::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in mopac/reset() \n");

  if (natoms!=natoms_i)
  {
    printf(" mopac reset failed due to different # of atoms \n");
    return;
  }

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  nfrz = 0;

  return;
}

void Mopac::xyz_save(string filename){

  ofstream xyzfile;
  //string xyzfile_string = "xyzfile.txt";
  xyzfile.open(filename.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(6);

   xyzfile << " " << natoms << endl;
   xyzfile << " " << endl;
   for (int i=0;i<natoms;i++)
   {
     xyzfile << "  " << anames[i];
     xyzfile << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2];
     xyzfile << endl;
   }


}

void Mopac::xyz_read(string filename){ 
   
//  printf(" in xyz_read \n");
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
    //cout << " RR (count: " << count << ")" << line << endl;
    if (count == 2)
    {
      if (!StringTools::cleanstring(line)) break;
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      xyz[3*i+0]=atof(tok_line[2].c_str());
      xyz[3*i+1]=atof(tok_line[3].c_str());
      xyz[3*i+2]=atof(tok_line[4].c_str());
//      cout << tok_line[0] << " " << tok_line[1] << " " << tok_line[2] << " " << endl; 
      i++;
    }

    if (line.find("CARTESIAN COORDINATES")!=string::npos && count == 0)
    {
//      cout << "Initial: " << line << endl;
      count++;
    }
    else if (line.find("CARTESIAN COORDINATES")!=string::npos && count == 1)
    {
//      cout << "Final: " << line << endl;
      count++;
      getline(output,line);
      getline(output,line);
      getline(output,line);
    }
  }
#if 0
  printf(" xyz_read: \n");
  for (int i=0;i<natoms;i++)
    printf(" %1.3f %1.3f %1.3f \n",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);
  printf("\n");
#endif


}   


