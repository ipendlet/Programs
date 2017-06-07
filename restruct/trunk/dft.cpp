#include "dft.h"
#include "constants.h"
//#include "utils.h"
using namespace std;


// for para execution
void DFT::sp_dnr(string filename, string nstr) {

  //printf(" in DFT/sp_dnr() \n");

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

#if !QSTART
  inpfile << " $molecule " << endl;
#if !CATION
#if !TRIPLET
  inpfile << " 0 1 " << endl;
#else
  inpfile << " 0 3 " << endl;
#endif
#else
  inpfile << " 1 1 " << endl;
#endif
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  inpfile << " $end " << endl;
  inpfile << endl << endl;
  inpfile << " $rem " << endl;
  inpfile << " JOBTYPE SP " << endl;
#if DFTB3LYP
  inpfile << " EXCHANGE B3LYP " << endl;
  inpfile << " CORRELATION NONE " << endl;
#else
  inpfile << " EXCHANGE PBE " << endl;
  inpfile << " CORRELATION PBE " << endl;
#endif
#if UNRESTRICTED
  inpfile << " UNRESTRICTED TRUE " << endl;
#endif
  inpfile << " SCF_ALGORITHM diis " << endl;
//  inpfile << " SCF_ALGORITHM gdm " << endl;
  inpfile << " SCF_MAX_CYCLES 250 " << endl;
  inpfile << " SCF_CONVERGENCE 6 " << endl;
#if B631GSS
  inpfile << " BASIS 6-31G** " << endl;
#elif LANL2DZ
  inpfile << " BASIS LANL2DZ " << endl;
  inpfile << " ECP LANL2DZ " << endl;
#elif BASISMIX
  inpfile << " BASIS mixed " << endl;
  inpfile << " ECP gen " << endl;
#elif BASISGEN
  inpfile << " BASIS gen " << endl;
  inpfile << " ECP gen " << endl;
#else
  inpfile << " BASIS 6-31G " << endl;     
#endif
  inpfile << " WAVEFUNCTION_ANALYSIS FALSE " << endl;
  inpfile << " $end " << endl;
  inpfile << endl;

#if (BASISMIX || BASISGEN)
  inpfile << endl;
  string cmd = "cat qmix >> "+inpfile_string;
  system(cmd.c_str());
  inpfile << endl;
#endif

#else
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  string cmd2 = "mv "+filename+" scratch/molecule"+nstr;
  system(cmd2.c_str());
  cmd2 = "./gscreate0 "+nstr;
  system(cmd2.c_str());
#endif

//  string cmd = "qchem "+filename+" "+filename+".out";
//  system(cmd.c_str());

  return;

}

// for para execution
void DFT::opt_dnr(string filename, string nstr) {

  //printf(" in DFT/opt_dnr() \n");

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

#if !QSTART
  inpfile << " $molecule " << endl;
#if !CATION
#if !TRIPLET
  inpfile << " 0 1 " << endl;
#else
  inpfile << " 0 3 " << endl;
#endif
#else
  inpfile << " 1 1 " << endl;
#endif
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  inpfile << " $end " << endl;
  inpfile << endl << endl;
  inpfile << " $rem " << endl;
  inpfile << " JOBTYPE OPT " << endl;
#if DFTB3LYP
  inpfile << " EXCHANGE B3LYP " << endl;
  inpfile << " CORRELATION NONE " << endl;
#else
  inpfile << " EXCHANGE PBE " << endl;
  inpfile << " CORRELATION PBE " << endl;
#endif
#if UNRESTRICTED
  inpfile << " UNRESTRICTED TRUE " << endl;
#endif

  inpfile << " GEOM_OPT_MAX_CYCLES  " << DFT_OPT_CYCLES << endl;
  inpfile << " GEOM_OPT_TOL_DISPLACEMENT 2500 " << endl;
  inpfile << " GEOM_OPT_TOL_GRADIENT     800 " << endl;
  inpfile << " GEOM_OPT_TOL_ENERGY	5000 " << endl;

  inpfile << " SCF_ALGORITHM diis " << endl;
//  inpfile << " SCF_ALGORITHM gdm " << endl;
  inpfile << " SCF_MAX_CYCLES 250 " << endl;
  inpfile << " SCF_CONVERGENCE 6 " << endl;
#if B631GSS
  inpfile << " BASIS 6-31G** " << endl;
#elif LANL2DZ
  inpfile << " BASIS LANL2DZ " << endl;
  inpfile << " ECP LANL2DZ " << endl;   
#elif BASISMIX
  inpfile << " BASIS mixed " << endl;
  inpfile << " ECP gen " << endl;
#elif BASISGEN
  inpfile << " BASIS gen " << endl;
  inpfile << " ECP gen " << endl;
#else
  inpfile << " BASIS 6-31G " << endl;     
#endif
  inpfile << " WAVEFUNCTION_ANALYSIS FALSE " << endl;
  inpfile << " $end " << endl;
  inpfile << endl;

#if (BASISMIX || BASISGEN)
  inpfile << endl;
  string cmd = "cat qmix >> "+inpfile_string;
  system(cmd.c_str());
  inpfile << endl;
#endif
#else
  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  string cmd2 = "mv "+filename+" scratch/molecule"+nstr;
  system(cmd2.c_str());
  cmd2 = "./gscreate0 "+nstr;
  system(cmd2.c_str());
#endif

//  string cmd = "qchem "+filename+" "+filename+".out";
//  system(cmd.c_str());

  return;

}

// for para execution
void DFT::ts_dnr(string filename) {

  //printf(" in DFT/ts_dnr() \n");

  energy0 = energy = 0;

  ofstream inpfile;
  string inpfile_string = filename;
  inpfile.open(inpfile_string.c_str());
  inpfile.setf(ios::fixed);
  inpfile.setf(ios::left);
  inpfile << setprecision(6);

  inpfile << " $molecule " << endl;
#if !CATION
#if !TRIPLET
  inpfile << " 0 1 " << endl;
#else
  inpfile << " 0 3 " << endl;
#endif
#else
  inpfile << " 1 1 " << endl;
#endif

  for (int i=0;i<natoms;i++)
     inpfile << " " << anames[i] << " " << xyz[3*i+0] << " " << xyz[3*i+1] << " " << xyz[3*i+2] << " " << endl;
  inpfile << " $end " << endl;
  inpfile << endl << endl;
  inpfile << " $rem " << endl;
  inpfile << " JOBTYPE TS " << endl;
  inpfile << " GEOM_OPT_DMAX 80 " << endl;
  inpfile << " GEOM_OPT_MAX_CYCLES 150 " << endl;

  inpfile << " GEOM_OPT_TOL_DISPLACEMENT 2500 " << endl;
  inpfile << " GEOM_OPT_TOL_GRADIENT     800 " << endl;
  inpfile << " GEOM_OPT_TOL_ENERGY      5000 " << endl;

#if DFTB3LYP
  inpfile << " EXCHANGE B3LYP " << endl;
  inpfile << " CORRELATION NONE " << endl;
#else
  inpfile << " EXCHANGE PBE " << endl;
  inpfile << " CORRELATION PBE " << endl;
#endif
#if UNRESTRICTED
  inpfile << " UNRESTRICTED TRUE " << endl;
#endif
  inpfile << " SCF_ALGORITHM diis " << endl;
//  inpfile << " SCF_ALGORITHM gdm " << endl;
  inpfile << " SCF_MAX_CYCLES 250 " << endl;
  inpfile << " SCF_CONVERGENCE 6 " << endl;
#if B631GSS
  inpfile << " BASIS 6-31G** " << endl;
#elif LANL2DZ
  inpfile << " BASIS LANL2DZ " << endl;
  inpfile << " ECP LANL2DZ " << endl;   
#elif BASISMIX
  inpfile << " BASIS mixed " << endl;
  inpfile << " ECP gen " << endl;
#elif BASISGEN
  inpfile << " BASIS gen " << endl;
  inpfile << " ECP gen " << endl;
#else
  inpfile << " BASIS 6-31G " << endl;     
#endif
  inpfile << " WAVEFUNCTION_ANALYSIS FALSE " << endl;
  inpfile << " MOLDEN_FORMAT TRUE " << endl;
  inpfile << " $end " << endl;
  inpfile << endl;

#if (BASISMIX || BASISGEN)
  inpfile << endl;
  string cmd = "cat qmix >> "+inpfile_string;
  system(cmd.c_str());
  inpfile << endl;
#endif

  return;
}

void DFT::get_structure(string filename, double* xyzc) {

  xyz_read(filename);
  for (int i=0;i<natoms*3;i++)
    xyzc[i] = xyz[i];

  return;
}

double DFT::get_energy(string filename) {

  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening dft file: %s \n",oname.c_str()); return 0.0; }
  string line;
  vector<string> tok_line;
  energy = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (line.find("Total energy in the final basis set")!=string::npos)
    {
      cout << "  DFT out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[8].c_str());
      break;
    }
    else if (line.find("E_qmmm")!=string::npos)
    {
      cout << "  DFT out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[3].c_str());
      break;
    }
  }
 
//  printf(" DFT energy: %1.4f \n",energy); 

  if (abs(energy)<0.00001 || (energy != energy))
  {
    printf(" energy zero, DFT failed \n");
    return 10000;
  }

  return energy;
}


double DFT::get_opt_energy(string filename) {

  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening dft file: %s \n",oname.c_str()); return 0.0; }
  string line;
  vector<string> tok_line;
  energy = 0;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (line.find("Total energy in the final basis set")!=string::npos)
    {
//      cout << "  DFT out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[8].c_str());
    }
    else if (line.find("E_qmmm")!=string::npos)
    {
//      cout << "  DFT out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energy=atof(tok_line[3].c_str());
    }
  }
 
  printf(" DFT opt energy: %1.4f \n",energy); 

  if (abs(energy)<0.00001 || (energy != energy))
  {
    printf(" energy zero, DFT failed \n");
    return 10000;
  }

  return energy;
}

double DFT::get_energy_ts(string filename) {

  energyts = 0.0;
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  if (!output) { printf(" error opening ts file: %s \n",oname.c_str()); return 0.0; }
  string line;
  vector<string> tok_line;
  int complete = 0;
  int hline = 0;
  double hval1 = -1;
  double gradts = -1;
  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;

    if (line.find("Gradient   ")!=string::npos)
    {
      tok_line = StringTools::tokenize(line, " \t");
      gradts = atof(tok_line[1].c_str());
    }
    else if (line.find("Hessian Eigenvalues")!=string::npos)
    {
      hline++;
    }
    else if (hline)
    {
      tok_line = StringTools::tokenize(line, " \t");
      hval1 = atof(tok_line[0].c_str());
      hline--;
    }

    else if (line.find("Total energy in the final basis")!=string::npos)
    {
//      cout << "  DFT TS out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energyts=atof(tok_line[8].c_str());
    }
    else if (line.find("Final energy")!=string::npos)
    {
//      cout << "  DFT TS Final out: " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      energyts=atof(tok_line[3].c_str());
      complete = 1;
      break;
    }
  }
 
//  printf(" DFT TS energy: %1.4f \n",energyts); 

  if (abs(energyts)<0.00001 || (energyts != energyts))
  {
    printf(" energy zero, DFT failed \n");
    return 10000;
  }
  if (!complete)
  {
    printf(" WARNING: %s not converged, energy is: %1.4f grad is: %1.4f eigenval: %1.3f \n",filename.c_str(),energyts,gradts,hval1);
    converged = 0;
  }
  else
  {
    printf(" SUCCESS: %s converged,     energy is: %1.4f grad is: %1.4f eigenval: %1.3f \n",filename.c_str(),energyts,gradts,hval1);
    converged = 1;
  }

  return energyts;
}

void DFT::alloc(int natoms_i) {
 
//  printf(" in DFT/alloc() \n");

  natoms = natoms_i;
  anumbers = new int[natoms];
  anames = new string[natoms];

  xyz0 = new double[3*natoms];
  xyz = new double[3*natoms];

  return;
}

void DFT::init(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in DFT/init() \n");

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

  return;
}

void DFT::freemem() {

  delete [] xyz0;
  delete [] xyz;
  delete [] anumbers;
  delete [] anames;

  return;
}

void DFT::reset(int natoms_i, int* anumbers_i, string* anames_i, double* xyz_i) {
 
//  printf(" in DFT/reset() \n");

  if (natoms!=natoms_i)
  {
    printf(" DFT reset failed due to different # of atoms \n");
    return;
  }

  for (int i=0;i<natoms;i++)
    anumbers[i] = anumbers_i[i];
  for (int i=0;i<natoms;i++)
    anames[i] = anames_i[i];
  for (int i=0;i<3*natoms;i++)
    xyz0[i] = xyz[i] = xyz_i[i];  

  return;
}

void DFT::xyz_read(string filename){ 
   
//  printf(" in DFT - xyz_read \n");
  string oname = filename+".out";
  ifstream output(oname.c_str(),ios::in);
  string line;
  vector<string> tok_line;
  int count = 0;

  while(!output.eof()) 
  { 
    getline(output,line);
//    cout << " RR " << line << endl;
    if (count == 1)
    {
//      cout << " reading! " << endl;
      for (int i=0;i<natoms;i++)
      {
        getline(output,line);
        if (line.find("-----")!=string::npos) break;
        if (!StringTools::cleanstring(line)) break;
        vector<string> tok_line = StringTools::tokenize(line, " \t");
        xyz[3*i+0]=atof(tok_line[2].c_str());
        xyz[3*i+1]=atof(tok_line[3].c_str());
        xyz[3*i+2]=atof(tok_line[4].c_str());
//      cout << tok_line[0] << " " << tok_line[1] << " " << tok_line[2] << " " << endl; 
      }
      count = 0;
    }

    if (line.find("Standard Nuclear Orientation")!=string::npos)
    {
//      cout << "Final: " << line << endl;
      count++;
      getline(output,line);
    }
  }
 



}   


