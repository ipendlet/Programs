#include <iostream>
#include <fstream>
#include <stdio.h>

#include "zstruct.h"
#include "icoord.h"
#include "geombasis.h"

using namespace std;

void generate_isomers(string xyzfile, string xyzlist);

int main(int argc, char* argv[]){

  string inpfile;
  string xyzfile;
  double emax = 10000;
  int nsteps = 1;
  xyzfile="test.xyz";
  switch (argc){
  case 1:
//    printf(" case 1 \n");
    emax = atof(argv[0]);
    //inpfile="xyzfile";
    break;
  case 2:
//    printf(" case 2 \n");
    nsteps = atoi(argv[1]);
    break;
  case 3:
//    printf(" case 3 \n");
    emax = atof(argv[1]);
    nsteps = atoi(argv[2]);
    break;
  default:
    inpfile="xyzfile";
    break;
//    return -1;
  }

  int done=0;

#if 0
  printf(" Calling ICoord \n");
  ICoord ictest;
  done = ictest.init(xyzfile);
  printf(" doing opt \n");
  ictest.opt();
  ictest.freemem();
 
  printf(" done calling ICoord test\n");
#endif

  if (nsteps>4) 
  {
    printf(" nsteps>4 is a bad idea \n");
    nsteps = 4;
  }
  //printf("emax is: %1.4f \n",emax);
  //printf("nsteps is: %i \n",nsteps);

  string xyzlist = "basislist1";
  ZStruct zmain;
//  zmain.init(xyzfile,xyzlist);
  zmain.init(xyzfile,xyzlist,emax);

  //printf("\n skipping everything else \n");
  zmain.go_mech_it(nsteps);

  //int* final_ints = new int[4];
  //zmain.do_one_step(final_ints);

  return 0;
}

