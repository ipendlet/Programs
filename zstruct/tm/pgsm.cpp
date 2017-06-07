#include "pgsm.h"
using namespace std;

#define SKIPGSM 0

double PGSM::gstring_se() {

  string filename = "initial.xyz";
  energyts = gstring_se(filename);
  return energyts;
}

double PGSM::gstring_lst() {

  //string filename = "initial.xyz";
  energyts = gstring_lst(0);
  return energyts;
}


double PGSM::get_lst(int run) {

  string nstr=StringTools::int2str(run,4,"0");
//  string outfile = "scratch/lst"+nstr+".log";

  string oname = "scratch/tsl"+nstr+".xyz";
  ifstream output(oname.c_str(),ios::in);
  string line;
  getline(output,line);
  getline(output,line);
//  cout << " RR " << line << endl;

  energyts = atof(line.c_str());
  if (abs(energyts)<0.00001)
  { 
    printf(" GSM failed \n");
    energyts = 10000;
  }
  else
    printf(" final (absolute) TS energy is: %1.2f \n",energyts);
 
// need to retrieve final TS geometry, write to xyz (with ints)
  xyz_read(oname);


  return energyts;
}

double PGSM::gstring_lst(int run) {

  string nstr=StringTools::int2str(run,4,"0");
  string filename = "scratch/stringfile"+nstr+".xyzl";


  printf(" running gstring_lst  ");

  energyts = -9;

#if SKIPGSM
  printf(" skipping GSM! \n");
  return energy;
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);
 
  inpfile << " " << natoms << endl << endl;
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz0[3*i+0] << " " << xyz0[3*i+1] << " " << xyz0[3*i+2] << endl;
  inpfile << " " << natoms << endl << endl;
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz1[3*i+0] << " " << xyz1[3*i+1] << " " << xyz1[3*i+2] << endl;
  inpfile << endl;

//  string cmd0 = "mv -f "+filename+".gsm "+filename+".gsm_prev";
//  system(cmd0.c_str());
//  string outfile = filename+".gsm";
  string outfile = "lst.log";
  string cmd = "./gstringl.exe "+nstr+" >"+outfile;
//  cout << " running cmd " << cmd << endl;

  fflush(stdout);
  system(cmd.c_str());

//  printf(" getting final energy \n");
//  string oname = "scratch/tslst.xyz";
  string oname = "scratch/tsl"+nstr+".xyz";
  ifstream output(oname.c_str(),ios::in);
  string line;
  getline(output,line);
  getline(output,line);
//  cout << " RR " << line << endl;

//  energyts = atof(line.c_str()) * 627.5;
  energyts = atof(line.c_str());
  if (abs(energyts)<0.00001)
  { 
    printf(" GSM failed \n");
    energyts = 10000;
  }
  else
    printf(" final (absolute) TS energy is: %1.2f \n",energyts);
 
// need to retrieve final TS geometry, write to xyz (with ints)
  xyz_read(oname);


  inpfile.close();

  return energyts;
}

double PGSM::gstring_se(string filename) {

 //printf(" in pgsm/gstring() \n");
  printf(" running gstring_se  ");

  energyts = -9;

#if SKIPGSM
  printf(" skipping GSM! \n");
  return energy;
#endif

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);
 
  inpfile << " " << natoms << endl << endl;
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz0[3*i+0] << " " << xyz0[3*i+1] << " " << xyz0[3*i+2] << endl;
  inpfile << " " << natoms << endl << endl;
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz1[3*i+0] << " " << xyz1[3*i+1] << " " << xyz1[3*i+2] << endl;
  inpfile << endl;

  string outfile = filename+".gsm";
  string cmd = "./gstring.exe >"+outfile;

  fflush(stdout);
  system(cmd.c_str());

//  printf(" getting final energy \n");
  string oname = "ts.xyz";
  ifstream output(oname.c_str(),ios::in);
  string line;
  getline(output,line);
  getline(output,line);
//  cout << " RR " << line << endl;

  energyts = atof(line.c_str()) * 627.5;
  if (abs(energyts)<0.00001)
  { 
    printf(" GSM failed \n");
    energyts = 10000;
  }
  else
    printf(" final (absolute) TS energy is: %1.2f \n",energyts);
 
// need to retrieve final TS geometry, write to xyz (with ints)
  xyz_read(oname);


  inpfile.close();

  return energyts;
}




double PGSM::gstring_dft() {

  string filename = "initial.xyz";

  energyts = gstring_dft(filename);
 
  return energyts;
}

double PGSM::gstring_dft(string filename) {

  printf(" running gstring_dft  ");

  energyts = -9;

#if SKIPGSM
  printf(" skipping GSM! \n");
  return energy;
#endif

  cout << filename << endl;
  string outfile = filename+".gsmq";
  string cmd = "./gstringq.exe "+filename+" >"+outfile;
//  cout << " " << cmd << endl;

  system(cmd.c_str());

//  printf(" getting final energy \n");
  string oname = "ts.xyz";
  ifstream output(oname.c_str(),ios::in);
  string line;
  getline(output,line);
  getline(output,line);
//  cout << " RR " << line << endl;

//  energyts = (atof(line.c_str()) - dftelist[0]) * 627.5;
  energyts = atof(line.c_str()); // * 627.5;
  if (abs(energyts)<0.00001)
  { 
    printf(" GSM failed \n");
    energyts = 10000;
  }
  else
    printf(" final TS energy is: %1.4f \n",energyts);
 
  xyz_read(oname);

  return energyts;
}


#if 0
void PGSM::gstring_dft_dnr(string filename) {

//CPMZ this function may be unnecessary 

  printf(" setting up input for gstring_dft  ");

  energyts = -9;

  cout << filename << endl;
  string outfile = filename+".gsmq";
  string cmd = "./gstringq.exe "+filename+" >"+outfile;
  cout << " " << cmd << endl;

//  system(cmd.c_str());

}
#endif

double PGSM::get_energy(string filename){


  double eread = -1;
#if SKIPGSM
  printf(" skipping GSM! \n");
  return eread;
#endif

#if 1
//  printf(" getting final energy \n");
  string oname = "ts.xyz";
  ifstream output(oname.c_str(),ios::in);
  string line;
  getline(output,line);
  getline(output,line);
//  cout << " RR " << line << endl;

  energyts = atof(line.c_str());
  if (abs(energyts)<0.00001)
  { 
    printf(" GSM failed \n");
    energyts = 10000;
  }
  else
    printf(" final TS energy is: %1.2f \n",energyts);
 
// need to retrieve final TS geometry, write to xyz (with ints)
  xyz_read(oname);

#endif

  return energyts;
}

void PGSM::alloc(int natoms_i) {
 
//  printf(" in mopac/alloc() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz1 = new double[3*natoms];
  xyzts = new double[3*natoms];

  return;
}

void PGSM::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i, double* xyz_f) {
 
//  printf(" in mopac/init() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz1 = new double[3*natoms];
  xyzts = new double[3*natoms];

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz_i[i];  
  for (int i=0;i<3*natoms;i++)
    xyz1[i] = xyz_f[i];  

  return;
}

void PGSM::freemem() {

  delete [] xyz0;
  delete [] xyz1;
  delete [] xyzts;
  delete [] anumbers;
  delete [] anames;

  return;
}

void PGSM::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i, double* xyz_f) {
 
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
    xyz0[i] = xyz_i[i];  
  for (int i=0;i<3*natoms;i++)
    xyz1[i] = xyz_f[i];  

  return;
}

void PGSM::xyz_read(string filename){ 

//  printf(" in xyz_read \n");
  string oname = filename;
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int count = 0;
  int i = 0;
  int tnatoms = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
    //cout << " RR " << line << endl;
    if (count == 1)
    {
      if (!StringTools::cleanstring(line)) break;
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      xyzts[3*i+0]=atof(tok_line[1].c_str());
      xyzts[3*i+1]=atof(tok_line[2].c_str());
      xyzts[3*i+2]=atof(tok_line[3].c_str());
      i++;
    }

    if (count == 0)
    {
//      cout << "Initial: " << line << endl;
      tnatoms = atoi(line.c_str());
//      cout << " found " << tnatoms << " atoms " << endl;
      count++;
      getline(output,line);
    } 
  } // while reading file
 
}   




