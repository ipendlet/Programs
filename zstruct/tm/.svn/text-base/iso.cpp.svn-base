#include <iostream>
#include <fstream>
#include <stdio.h>

#include "zstruct.h"
#include "icoord.h"
#include "geombasis.h"
#include "mopac.h"
#include "iso.h"
#include "constants.h"
#include "utils.h"

using namespace std;

#define PPMM 0
#define OPTSHADOW 0

int ZStruct::generate_isomers(int id, int nsave) {

  printf("\n\n ------------------------------------ \n");
  printf(" ---- 3. choosing step ---- \n");

  int natoms = geoms[id].natoms;
  int natoms0 = geoms[id].natoms0;
  int N3 = natoms0*3;

  printf(" geom id: %i \n",id);
 // printf(" in generate_isomers, natoms: %i natoms0: %i \n\n",natoms,natoms0);
  geoms[id].print_xyz();
  geoms[id].print_ic();

  int done = 0;
  int tniso = 0;
  int same = 1;
  int nscratch = 0;

//  int* olist = new int[NMAXSTRUCTS];
//  for (int j=0;j<100;j++) olist[j]=-1;

 //isosr is the return geometry
  double** isos = new double*[NMAXSTRUCTS];
  for (int i=0;i<NMAXSTRUCTS;i++)
    isos[i] = new double[4*natoms];
//  double** isosr = new double*[NMAXSTRUCTS];
//  for (int i=0;i<NMAXSTRUCTS;i++)
//    isosr[i] = new double[3*natoms];
  double* energies = new double[NMAXSTRUCTS];
  //double* energiesr = new double[NMAXSTRUCTS];

  geom1r.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
//  geom1r.skipTMadd = 1;
  geom1r.ic_create();

  Mopac mopt;
  mopt.alloc(natoms0);
  int* frzlist = new int[natoms*2+100];
  int nfrz = 0;
  int nfrz0 = 0;

  mopt.reset(natoms0,geoms[id].anumbers,geoms[id].anames,geoms[id].coords);
  addSolventMopac(natoms0,mopt);

  //energies[0] = mopt.opt("scratch/xyzfile.xyzm_initial");
  energies[0] = mopt.energy_sp("scratch/xyzfile.xyzm_initial1");
  elist[0] = energies[0];
//  if (id==0) geoms[0].seenergy = energies[0];


  // Expect initial structure to already be optimized?
  for (int i=0;i<3*natoms0;i++)
    isos[0][i] = geoms[id].coords[i];

//  for (int i=0;i<3*natoms;i++)
//    isos[0][i] = mopt.xyz[i];

//  geoms[id].reset(natoms,geoms[id].anames,geoms[id].anumbers,isos[0]);
  
//  for (int i=0;i<natoms;i++)
//    printf(" %s %1.4f %1.4f %1.4f \n",geoms[id].anames[i].c_str(),isos[0][3*i],isos[0][3*i+1],isos[0][3*i+2]);

  printf("\n");
  for (int i=0;i<natoms;i++)
    printf(" atom %i is frozen? %i \n",i,frozen[i]);



  //determine allowed changes in coordination number

  // cmoves stores list of possible changes (currently max 3 per center)
  int** cmovesp = new int*[natoms];
  for (int i=0;i<natoms;i++)
  {
    cmovesp[i] = new int[3];
    cmovesp[i][0] = cmovesp[i][1] = cmovesp[i][2] = -1;
  }
  int* ncmovesp = new int[natoms];
  for (int i=0;i<natoms;i++) 
    ncmovesp[i]=0;

  int** cmovesm = new int*[natoms];
  for (int i=0;i<natoms;i++)
  {
    cmovesm[i] = new int[3];
    cmovesm[i][0] = cmovesm[i][1] = cmovesm[i][2] = -1;
  }
  int* ncmovesm = new int[natoms];
  for (int i=0;i<natoms;i++) 
    ncmovesm[i]=0;

  int** cmoves2p = new int*[natoms];
  for (int i=0;i<natoms;i++)
  {
    cmoves2p[i] = new int[3];
    cmoves2p[i][0] = cmoves2p[i][1] = cmoves2p[i][2] = -1;
  }
  int* ncmoves2p = new int[natoms];
  for (int i=0;i<natoms;i++) 
    ncmoves2p[i]=0;

  int** cmoves2m = new int*[natoms];
  for (int i=0;i<natoms;i++)
  {
    cmoves2m[i] = new int[3];
    cmoves2m[i][0] = cmoves2m[i][1] = cmoves2m[i][2] = -1;
  }
  int* ncmoves2m = new int[natoms];
  for (int i=0;i<natoms;i++) 
    ncmoves2m[i]=0;

  int anum;
  int atype0;
  int atype1;

  for (int i=0;i<natoms;i++)
  {
    //printf("\n\n looking for moves for atom %i \n",i);

    atype0 = basis1.find(i,geoms[id]);
    cmovesp[i][0] = atype0; //save original type
    cmovesm[i][0] = atype0; //save original type
    cmoves2p[i][0] = atype0; //save original type
    cmoves2m[i][0] = atype0; //save original type

    atype1 = basis1.nplus(i,geoms[id],atype0); 
    if (atype1 > -1 && !frozen[i])
    {
      cmovesp[i][ncmovesp[i]+1] = atype1;
      ncmovesp[i]++;
    }
    atype1 = basis1.nplus2(i,geoms[id],atype0); 
    if (atype1 > -1 && !frozen[i])
    {
      cmoves2p[i][ncmoves2p[i]+1] = atype1;
      ncmoves2p[i]++;
    }
    atype1 = basis1.nminus(i,geoms[id],atype0);
    if (atype1 > -1 && !frozen[i])
    {
      cmovesm[i][ncmovesm[i]+1] = atype1;
      ncmovesm[i]++;
    }
    atype1 = basis1.nminus2(i,geoms[id],atype0);
    if (atype1 > -1 && !frozen[i])
    {
      cmoves2m[i][ncmoves2m[i]+1] = atype1;
      ncmoves2m[i]++;
    }
    printf(" ncmovesp found for atom %i: %i \n",i,ncmovesp[i]);
    printf(" ncmovesm found: %i \n",ncmovesm[i]);
    printf(" ncmoves2p found: %i \n",ncmoves2p[i]);
    printf(" ncmoves2m found: %i \n",ncmoves2m[i]);

  } //loop i over natoms


  //look for move pairs

  printf("\n\n determining move pairs \n");

  int nasmm = 0; //number of moves 2x minus
  int naspp = 0; //number of moves 2x plus
  int nasmp = 0; //number of moves +/-
  int nasmp2 = 0; // number of moves 2x +/-

  int** addsubmm  = new int*[NMAXSTRUCTS];
  int** addsubpp  = new int*[NMAXSTRUCTS];
  int** addsubmp  = new int*[NMAXSTRUCTS];
  int** addsubmp2 = new int*[NMAXSTRUCTS];
  for (int i=0;i<NMAXSTRUCTS;i++)
    addsubmm[i] = new int[6];
  for (int i=0;i<NMAXSTRUCTS;i++)
    addsubpp[i] = new int[6];
  for (int i=0;i<NMAXSTRUCTS;i++)
    addsubmp[i] = new int[6];
  for (int i=0;i<NMAXSTRUCTS;i++)
    addsubmp2[i] = new int[6];



///////////////////////////////////////////////////
// determine allowed changes in coordn/atom type //
///////////////////////////////////////////////////


// +/- moves
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    {
      int k1,k2;
//      printf(" ncmovesp[i]: %i ncmovesm[j]: %i \n",ncmovesp[i],ncmovesm[j]);
//      for (int k1=1;k1<=ncmovesp[i];k1++)
//      for (int k2=1;k2<=ncmovesm[j];k2++)
      if (ncmovesp[i] && ncmovesm[j])
      {
        k1 = k2 = 1;
//        printf("i: atom %i can gain %i bonds \n",i,basis1.ats[cmovesp[i][k1]].coordn-basis1.ats[cmovesp[i][0]].coordn); 
//        printf("j: atom %i can gain %i bonds \n",j,basis1.ats[cmovesm[j][k2]].coordn-basis1.ats[cmovesm[j][0]].coordn); 

        printf(" found add+sub pair centers %i %i \n",i,j);
        addsubmp[nasmp][0] = i;
        addsubmp[nasmp][1] = cmovesp[i][0];
        addsubmp[nasmp][2] = cmovesp[i][k1];
        addsubmp[nasmp][3] = j;
        addsubmp[nasmp][4] = cmovesm[j][0];
        addsubmp[nasmp][5] = cmovesm[j][k2];
        nasmp++;
      } // loop pairs of (single) p/m combo moves

//      for (int k1=1;k1<=ncmovesm[i];k1++)
//      for (int k2=1;k2<=ncmovesp[j];k2++)
      if (ncmovesm[i] && ncmovesp[j])
      {
        k1 = k2 = 1;
//        printf("i: atom %i can gain %i bonds \n",i,basis1.ats[cmovesm[i][k1]].coordn-basis1.ats[cmovesm[i][0]].coordn); 
//        printf("j: atom %i can gain %i bonds \n",j,basis1.ats[cmovesp[j][k2]].coordn-basis1.ats[cmovesp[j][0]].coordn); 

          printf(" found add+sub pair centers %i %i\n",j,i);
          addsubmp[nasmp][3] = i;
          addsubmp[nasmp][4] = cmovesp[i][0];
          addsubmp[nasmp][5] = cmovesp[i][k1];
          addsubmp[nasmp][0] = j;
          addsubmp[nasmp][1] = cmovesm[j][0];
          addsubmp[nasmp][2] = cmovesm[j][k2];
          nasmp++;
      } // loop pairs of (single) p/m combo moves

//      for (int k1=1;k1<=ncmoves2p[i];k1++)
//      for (int k2=1;k2<=ncmoves2m[j];k2++)
      if (ncmoves2p[i] && ncmoves2m[j])
      {
        k1 = k2 = 1;
//        printf("i: atom %i can gain %i bonds \n",i,basis1.ats[cmoves2p[i][k1]].coordn-basis1.ats[cmoves2p[i][0]].coordn); 
//        printf("j: atom %i can gain %i bonds \n",j,basis1.ats[cmoves2m[j][k2]].coordn-basis1.ats[cmoves2m[j][0]].coordn); 

          printf(" found 2x add + 2x sub pair centers %i %i\n",i,j);
          addsubmp2[nasmp2][0] = i;
          addsubmp2[nasmp2][1] = cmoves2p[i][0];
          addsubmp2[nasmp2][2] = cmoves2p[i][k1];
          addsubmp2[nasmp2][3] = j;
          addsubmp2[nasmp2][4] = cmoves2m[j][0];
          addsubmp2[nasmp2][5] = cmoves2m[j][k2];
          nasmp2++;
      } // loop pairs of double p/m combo moves

//      for (int k1=1;k1<=ncmoves2m[i];k1++)
//      for (int k2=1;k2<=ncmoves2p[j];k2++)
      if (ncmoves2m[i] && ncmoves2p[j])
      {
        k1 = k2 = 1;
//        printf("i: atom %i can gain %i bonds \n",i,basis1.ats[cmoves2m[i][k1]].coordn-basis1.ats[cmoves2m[i][0]].coordn); 
//        printf("j: atom %i can gain %i bonds \n",j,basis1.ats[cmoves2p[j][k2]].coordn-basis1.ats[cmoves2p[j][0]].coordn); 

          printf(" found 2x add + 2x sub pair centers %i %i\n",j,i);
          addsubmp2[nasmp2][3] = i;
          addsubmp2[nasmp2][4] = cmoves2p[i][0];
          addsubmp2[nasmp2][5] = cmoves2p[i][k1];
          addsubmp2[nasmp2][0] = j;
          addsubmp2[nasmp2][1] = cmoves2m[j][0];
          addsubmp2[nasmp2][2] = cmoves2m[j][k2];
          nasmp2++;
      } // loop pairs of double m/p combo moves

    } // loop j<i
  } // loop over atoms i

// need to add 2x + two 1x

// 2x + moves
  printf(" locating 2x + moves \n");
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    for (int k=0;k<j;k++)
    for (int l=0;l<k;l++)
    {  
      if (ncmovesp[i] && ncmovesp[j] && ncmovesp[k] && ncmovesp[l])
      {
        int found;
        found = geoms[id].bond_exists_tri(i,j)+geoms[id].bond_exists_tri(k,l)
                + geoms[id].hpair(i,j) + geoms[id].hpair(k,l);
        found += geoms[id].isTMpair(i,j) + geoms[id].isTMpair(k,l);
        if (!found)
        {  
          printf(" found 2x add pair centers %i %i, %i %i \n",i,j,k,l);
          addsubpp[naspp][0] = i;
          addsubpp[naspp][1] = j;
          addsubpp[naspp][2] = k;
          addsubpp[naspp][3] = l;
          naspp++;
        } // if !found

        found = geoms[id].bond_exists_tri(i,k)+geoms[id].bond_exists_tri(j,l) 
                + geoms[id].hpair(i,k) + geoms[id].hpair(j,l);
        found += geoms[id].isTMpair(i,k) + geoms[id].isTMpair(j,l);
        if (!found)
        {  
          printf(" found 2x add pair centers %i %i, %i %i \n",i,k,j,l);
          addsubpp[naspp][0] = i;
          addsubpp[naspp][1] = k;
          addsubpp[naspp][2] = j;
          addsubpp[naspp][3] = l;
          naspp++;
        } // if !found

        found = geoms[id].bond_exists_tri(i,l)+geoms[id].bond_exists_tri(j,k) 
                + geoms[id].hpair(i,l) + geoms[id].hpair(j,k);
        found += geoms[id].isTMpair(i,l) + geoms[id].isTMpair(j,k);
        if (!found)
        {  
          printf(" found 2x add pair centers %i %i, %i %i \n",i,l,j,k);
          addsubpp[naspp][0] = i;
          addsubpp[naspp][1] = l;
          addsubpp[naspp][2] = j;
          addsubpp[naspp][3] = k;
          naspp++;
        } // if !found

      } // if all 4 atoms can make bonds 

    } //loop j<i, k<j, l<k
  } //loop i

// 2x - moves
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    for (int k=0;k<j;k++)
    for (int l=0;l<k;l++)
    {
      if (ncmovesm[i] && ncmovesm[j] && ncmovesm[k] && ncmovesm[l])
      {
        int found;
        found = geoms[id].bond_exists(i,j)+geoms[id].bond_exists(k,l);
        if (found>1)
        {  
          addsubmm[nasmm][0] = geoms[id].bond_num(i,j);
          addsubmm[nasmm][1] = geoms[id].bond_num(k,l);
          printf(" found 2x minus pair: %i %i, %i %i: bonds: %i %i \n",i,j,k,l,addsubmm[nasmm][0],addsubmm[nasmm][1]);
          nasmm++;
        } // if found

        found = geoms[id].bond_exists(i,k)+geoms[id].bond_exists(j,l);
        if (found>1)
        {  
          addsubmm[nasmm][0] = geoms[id].bond_num(i,k);
          addsubmm[nasmm][1] = geoms[id].bond_num(j,l);
          printf(" found 2x minus pair: %i %i, %i %i: bonds: %i %i \n",i,k,j,l,addsubmm[nasmm][0],addsubmm[nasmm][1]);
          nasmm++;
        } // if found

        found = geoms[id].bond_exists(i,l)+geoms[id].bond_exists(j,k);
        if (found>1)
        {  
          addsubmm[nasmm][0] = geoms[id].bond_num(i,l);
          addsubmm[nasmm][1] = geoms[id].bond_num(j,k);
          printf(" found 2x minus pair: %i %i, %i %i: bonds: %i %i \n",i,l,j,k,addsubmm[nasmm][0],addsubmm[nasmm][1]);
          nasmm++;
        } // if found

      } // if all 4 atoms can lose bonds

    } //loop j<i, k<j, l<k
  } //loop i

// need to add 2x + two 1x


  printf("\n nasmp: %i naspp: %i nasmm: %i \n",nasmp,naspp,nasmm);






//////////////////
// do the moves //
//////////////////

  int np = 0;
  int** movesp = new int*[natoms*8+1000];
  for (int i=0;i<natoms*8+1000;i++)
    movesp[i] = new int[2];
#if 1
  printf("\n\n add one bond moves \n");
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    if (!geoms[id].bond_exists_tri(i,j) && ncmovesp[i] && ncmovesp[j]
        && !geoms[id].hpair(i,j) && !geoms[id].isTMpair(i,j))
    {
      movesp[np][0] = i;
      movesp[np][1] = j;
      np++;
    }
  }

  printf("\n\n doing %i add one bond moves \n",np);
  for (int i=0;i<np;i++)
  {
    printf("\n adding bond: %i %i \n",movesp[i][0],movesp[i][1]);
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    printf(" natomstm: %i \n",geom1.natomstm);
    geom1.ic_create();
    printf(" natomstm: %i \n",geom1.natomstm);
    //geom1.print_bonds(); 
    geom1.bonds[geom1.nbonds][0]=movesp[i][0];
    geom1.bonds[geom1.nbonds][1]=movesp[i][1];
    geom1.nbonds = geoms[id].nbonds+1;
    geom1.ic_create_nobonds(); 
    geom1.mm_init();
    geom1.update_ic();
    //geom1.print_bonds(); 
    nfrz = 2;
    frzlist[0] = movesp[i][0];
    frzlist[1] = movesp[i][1];
    nfrz0 = nfrz;
    for (int j=0;j<geom1.nbonds;j++)
    {
      for (int k=0;k<nfrz0;k++)
      {
        if (frzlist[k]==geom1.bonds[j][0])
        { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
        if (frzlist[k]==geom1.bonds[j][1])
        { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
      }
    }

    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
    string optfiler = "scratch/xyzfile.xyz"+nstr+"r";
    //printf(" before opt \n");
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    done = geom1.opt(optfile,geomshadow,solvent);
#endif
    //printf(" after opt, printing ic \n");
    //geom1.print_ic();

// not doing this
//    geom1r.reset(geoms[id].coords);
//    done = geom1r.opt(optfiler);

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(natoms0,geom1.anumbers,geom1.anames,geom1.coords);
      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }

  } // loop i over np
#endif


  int nm = 0;
  int* movesm = new int[natoms*8];
#if 1
  printf("\n sub one bond moves \n");
  //printf(" natoms: %i \n",natoms);
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    if (geoms[id].bond_exists(i,j) && ncmovesm[i] && ncmovesm[j])
    {
      printf(" found bond to remove: %i %i \n",i,j);
      movesm[nm] = geoms[id].bond_num(i,j);
//      printf(" sub bond: %i \n",movesm[nm]);
      nm++;
    }
  }
#endif

  printf("\n doing %i sub one bond moves \n",nm);
  for (int i=0;i<nm;i++)
  {
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    printf("\n subtracting bond: %i atoms: %i %i \n",movesm[i],geom1.bonds[movesm[i]][0],geom1.bonds[movesm[i]][1]);
    //geom1.print_bonds(); 
    for (int l=movesm[i];l<geom1.nbonds-1;l++)
    {
      geom1.bonds[l][0] = geom1.bonds[l+1][0];
      geom1.bonds[l][1] = geom1.bonds[l+1][1];
    }
    geom1.nbonds = geoms[id].nbonds-1;
    geom1.ic_create_nobonds(); 
    geom1.mm_init();
    geom1.update_ic();
    //geom1.print_bonds(); 

    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow,solvent);
#endif

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
//      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }

  } // loop i over np



// +/- moves
#if 1
  printf("\n ******* creating mp moves +/- \n");
  for (int i=0;i<nasmp;i++)
  {
//    printf(" moving %i: %i type to %i type \n",addsubmp[i][0],addsubmp[i][1],addsubmp[i][2]);
//    printf(" moving %i: %i type to %i type \n",addsubmp[i][3],addsubmp[i][4],addsubmp[i][5]);
//    printf(" moving %i: %i coord to %i coord \n",addsubmp[i][0],basis1.ats[addsubmp[i][1]].coordn,basis1.ats[addsubmp[i][2]].coordn);
//    printf(" moving %i: %i coord to %i coord \n",addsubmp[i][3],basis1.ats[addsubmp[i][4]].coordn,basis1.ats[addsubmp[i][5]].coordn);
  
    int rem1,mov1;
    for (int j=0;j<geoms[id].nbonds;j++)
    {
      int skip = true;
       //must be bond connected to minus atom
//      if ( (addsubmp[i][3]==geoms[id].bonds[j][0] || addsubmp[i][3]==geoms[id].bonds[j][1]) )
      if ( (addsubmp[i][3]==geoms[id].bonds[j][0] || addsubmp[i][3]==geoms[id].bonds[j][1]) && 
         !((addsubmp[i][3]==geoms[id].bonds[j][0] && addsubmp[i][0]==geoms[id].bonds[j][1]) ||
           (addsubmp[i][3]==geoms[id].bonds[j][1] && addsubmp[i][0]==geoms[id].bonds[j][0])))
      {
         skip = false;
         //printf(" possible bond to destroy: %i %i (#%i) \n",geoms[id].bonds[j][0],geoms[id].bonds[j][1],j);
         rem1 = j;
         if (addsubmp[i][3]==geoms[id].bonds[j][0])
         {
           mov1=geoms[id].bonds[j][1];
         }
         else
         {
           mov1=geoms[id].bonds[j][0];
         }
//         printf(" with atom moving: %i \n",mov1);
         if (geoms[id].bond_exists_tri(mov1,addsubmp[i][0]) || 
             geoms[id].hpair(mov1,addsubmp[i][0]))
           skip = true;
         if(frozen[mov1]) skip = true;
         if(geoms[id].isTMpair(addsubmp[i][0],mov1))
         { 
           printf(" found tm pair, skipping this addition: %i %i \n",addsubmp[i][0],mov1);
           skip = true;
         }
      }  // if conditions on which bond
      if (!skip)
      {
       printf("\n bond to add: %i %i to remove: %i %i\n",addsubmp[i][0],mov1,
              geoms[id].bonds[rem1][0],geoms[id].bonds[rem1][1]);
     
       geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
       geom1.ic_create();
       geom1.nbonds = geoms[id].nbonds;
       geom1.bonds[rem1][0]=addsubmp[i][0];
       geom1.bonds[rem1][1]=mov1;
       geom1.ic_create_nobonds(); 
       geom1.mm_init();
       geom1.update_ic();
       //geom1.print_bonds(); 
       nfrz = 2;
       frzlist[0] = addsubmp[i][0];
       frzlist[1] = mov1;
       nfrz0 = nfrz;
       for (int j=0;j<geom1.nbonds;j++)
       {
         for (int k=0;k<nfrz0;k++)
         {
           if (frzlist[k]==geom1.bonds[j][0])
           { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
           if (frzlist[k]==geom1.bonds[j][1])
           { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
         }
       }

       string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
       string optfile = "scratch/xyzfile.xyz"+nstr;
       //printf(" before opt \n");
#if !OPTSHADOW
       done = geom1.opt(optfile);
#else
       geomshadow.reset(geoms[id].coords);
       done = geom1.opt(optfile,geomshadow,solvent);
#endif

       if (done)
       {
         //printf(" attempting mopac \n");
         optfile = "scratch/xyzfile.xyzm_"+nstr;
         mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
         mopt.freeze(frzlist, nfrz, nfrz0);
         printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
         mopt.opt_write(optfile);
         tniso++;
         for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
       }
       else { nfail++; nscratch--; }
      } // if !skip


    }
  } // loop i over nasmp
#endif

#if 1
  printf("\n creating mp moves +/- where + is not the detached atom \n");
  for (int i=0;i<nasmp;i++)
  {
//    printf(" moving %i: %i type to %i type \n",addsubmp[i][0],addsubmp[i][1],addsubmp[i][2]);
//    printf(" moving %i: %i type to %i type \n",addsubmp[i][3],addsubmp[i][4],addsubmp[i][5]);
//    printf(" moving %i: %i coord to %i coord \n",addsubmp[i][0],basis1.ats[addsubmp[i][1]].coordn,basis1.ats[addsubmp[i][2]].coordn);
//    printf(" moving %i: %i coord to %i coord \n",addsubmp[i][3],basis1.ats[addsubmp[i][4]].coordn,basis1.ats[addsubmp[i][5]].coordn);
  
    int rem1,mov1,mov2;
    for (int j=0;j<geoms[id].nbonds;j++)
    {
      int skip = true;
      //must be bond connected to minus atom
      if (addsubmp[i][3]==geoms[id].bonds[j][0] || addsubmp[i][3]==geoms[id].bonds[j][1])
      {
         skip = false;
         //printf(" possible bond to destroy: %i %i (#%i) \n",geoms[id].bonds[j][0],geoms[id].bonds[j][1],j);
         rem1 = j;
         if (addsubmp[i][3]==geoms[id].bonds[j][0])
           mov1=geoms[id].bonds[j][1];
         else
           mov1=geoms[id].bonds[j][0];
//         printf(" with atom detaching: %i \n",mov1);
         if (!ncmovesm[mov1]) skip = true;
         if(frozen[mov1]) skip = true;
         if(mov1>addsubmp[i][3]) skip = true;
//         if(geoms[id].isTMpair(addsubmp[i][0],mov1)) skip = true;
      }  // if conditions on which bond
      if (!skip)
      {
        printf(" in !skip, mov1,minus,plus: %i %i %i \n",mov1,addsubmp[i][3],addsubmp[i][0]);
        for (int k=0;k<addsubmp[i][0];k++)
        if(!geoms[id].isTMpair(addsubmp[i][0],k))
        if (!frozen[k] && mov1!=k && ncmovesp[k] && k<addsubmp[i][3])
        if (!geoms[id].bond_exists_tri(k,addsubmp[i][0]) && !geoms[id].hpair(k,addsubmp[i][0]))
        {
          mov2=k;
          printf("\n bond to add: %i %i to remove: %i %i\n",addsubmp[i][0],mov2,
                 geoms[id].bonds[rem1][0],geoms[id].bonds[rem1][1]);
     
          geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
          geom1.ic_create();
          geom1.nbonds = geoms[id].nbonds;
          geom1.bonds[rem1][0]=addsubmp[i][0];
          geom1.bonds[rem1][1]=mov2;
          geom1.ic_create_nobonds(); 
          //geom1.mm_init();
          geom1.update_ic();
          //geom1.print_bonds(); 
          nfrz = 2;
          frzlist[0] = addsubmp[i][0];
          frzlist[1] = mov2;
          nfrz0 = nfrz;
          for (int j=0;j<geom1.nbonds;j++)
          {
            for (int k=0;k<nfrz0;k++)
            {
              if (frzlist[k]==geom1.bonds[j][0])
              { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
              if (frzlist[k]==geom1.bonds[j][1])
              { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
            }
          }

          string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
          string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
          done = geom1.opt(optfile);
#else
          geomshadow.reset(geoms[id].coords);
          done = geom1.opt(optfile,geomshadow);
#endif

          if (done)
          {
            //printf(" attempting mopac \n");
           optfile = "scratch/xyzfile.xyzm_"+nstr;
            mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
            mopt.freeze(frzlist, nfrz, nfrz0);
            printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
            mopt.opt_write(optfile);
            tniso++;
            for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
          }
          else { nfail++; nscratch--; }
        } // loop over second + atom k
      } // if !skip


    }
  } // loop i over nasmp
#endif

  // do 2x +/- moves
  int nmoves2pm = 0;
  int** moves2pm = new int*[nasmp*nasmp*10+1000];
  for (int i=0;i<nasmp*nasmp*10+1000;i++)
    moves2pm[i] = new int[8];

#if PPMM
  printf("\n finding 2x +/- moves \n");
  for (int i=0;i<nasmp;i++)
  {
    for (int j=0;j<i;j++)
    {  
      int skip = true;
      int mov1,mov2,rem1,rem2;
      int a1,a2,c1,c2;
      a1=0;a2=0;c1=0;c2=0;mov1=0;mov2=0;rem1=0;rem2=0;
      mov1 = 0; //avoid mem problem

      // condition 1. all 4 unique
      // condition 2. 2 atoms repeated (so do exchange)
      // condition 3. 1 atom repeated
      // condition 4. 1 atom 2x +
      // condition 5. 1 atom 2x -
      if ( addsubmp[i][0] != addsubmp[j][0] && addsubmp[i][3] != addsubmp[j][0]
        && addsubmp[i][0] != addsubmp[j][3] && addsubmp[i][3] != addsubmp[j][3] )
      {
        printf(" not repeated: %i %i %i %i \n",addsubmp[i][0],addsubmp[i][3],addsubmp[j][0],addsubmp[j][3]);
        a1 = addsubmp[i][0]; //gaining bond to atom  
        a2 = addsubmp[j][0]; //gaining bond to atom 
        c1 = addsubmp[i][3]; //losing bonded atom
        c2 = addsubmp[j][3]; //losing bonded atom

        //printf(" a1,a2: %i %i c1,c2: %i %i \n",a1,a2,c1,c2);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
//        for (int k2=0;k2<k1;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1) 
            && (geoms[id].bonds[k2][0]==c2 || geoms[id].bonds[k2][1]==c2) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==c1) mov1=geoms[id].bonds[k1][1];
            else mov1=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==c2) mov2=geoms[id].bonds[k2][1];
            else mov2=geoms[id].bonds[k2][0];
          }
          if (rem1==rem2) skip = true;
          if (geoms[id].bond_exists_tri(mov1,a1) || geoms[id].bond_exists_tri(mov2,a2)) skip = true;
          if (a1==mov1 || a2==mov2) skip = true;
//          printf(" a1,mov1 %i %i a2,mov2 %i %i \n",a1,mov1,a2,mov2);
          //CPMZ allowing a single partly-attached H2 to form
          //if (geoms[id].hpair(a1,mov1) || geoms[id].hpair(a2,mov2)) skip = true;
          //if (geoms[id].hpair(a1,mov1)) skip = true;
          if (geoms[id].hpair(a1,mov1) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          if (geoms[id].hpair(a2,mov2) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          //if (geoms[id].hpair(a1,mov1) && geoms[id].hpair(a2,mov2)) skip = true;
          if(geoms[id].isTMpair(a1,mov1) || geoms[id].isTMpair(a2,mov2)) skip = true;

          if (frozen[mov1] || frozen[mov2]) skip = true;
          if (frozen[a1] || frozen[a2]) skip = true;
          if (!skip)
          {
            printf(" found rem1 rem2: %i %i for atoms %i %i, %i %i \n",rem1,rem2,a1,mov1,a2,mov2);
            moves2pm[nmoves2pm][0] = rem1;
            moves2pm[nmoves2pm][1] = rem2;
            moves2pm[nmoves2pm][2] = a1;
            moves2pm[nmoves2pm][3] = mov1;
            moves2pm[nmoves2pm][4] = a2;
            moves2pm[nmoves2pm][5] = mov2;
            nmoves2pm++;
          } // if !skip
        } // loop k1,k2 over all bond pairs

      } // atom not repeated
      else if ( addsubmp[i][0] == addsubmp[j][3] && addsubmp[i][3] == addsubmp[j][0])
      {
        //this is a 2x center exchange reaction
        printf(" double repeated: %i %i %i %i \n",addsubmp[i][0],addsubmp[i][3],addsubmp[j][0],addsubmp[j][3]);

        a1 = addsubmp[i][0]; //gaining/losing bonded atom (central) 
        a2 = addsubmp[j][0]; //gaining/losing bonded atom (central2)

        //printf(" a1,a2: %i %i \n",a1,a2);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
//        for (int k2=0;k2<k1;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==a1 || geoms[id].bonds[k1][1]==a1) 
            && (geoms[id].bonds[k2][0]==a2 || geoms[id].bonds[k2][1]==a2) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==a1) mov1=geoms[id].bonds[k1][1];
            else mov1=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==a2) mov2=geoms[id].bonds[k2][1];
            else mov2=geoms[id].bonds[k2][0];
          }
          if (rem1==rem2) skip = true;
          if (geoms[id].bond_exists_tri(mov1,a2) || geoms[id].bond_exists_tri(mov2,a1)) skip = true;
          if (a1==mov2 || a2==mov1) skip = true;
          if(geoms[id].isTMpair(a1,mov1) || geoms[id].isTMpair(a2,mov2)) skip = true;

          if (frozen[mov1] || frozen[mov2]) skip = true;
          if (frozen[a1] || frozen[a2]) skip = true;
          if (!skip)
          {
            //printf(" found rem1 rem2: %i %i for atoms %i %i, %i %i \n",rem1,rem2,a1,mov2,a2,mov1);
            moves2pm[nmoves2pm][0] = rem1;
            moves2pm[nmoves2pm][1] = rem2;
            moves2pm[nmoves2pm][2] = a1;
            moves2pm[nmoves2pm][3] = mov2;
            moves2pm[nmoves2pm][4] = a2;
            moves2pm[nmoves2pm][5] = mov1;
            nmoves2pm++;
          } // if !skip
        } // loop k1,k2 over all bond pairs

      }
      else if ( addsubmp[i][0] == addsubmp[j][3] || addsubmp[i][3] == addsubmp[j][0])
      {
        printf("\n repeated: %i %i %i %i \n",addsubmp[i][0],addsubmp[i][3],addsubmp[j][0],addsubmp[j][3]);
        if (addsubmp[i][0]==addsubmp[j][3])
        {
          a1 = addsubmp[i][0]; //gaining/losing bonded atom (central)
          a2 = addsubmp[j][0]; //gaining bond to atom
          c1 = addsubmp[i][3]; //losing bonded atom
        }
        else if (addsubmp[i][3]==addsubmp[j][0])
        {
          a1 = addsubmp[i][3]; //gaining/losing bonded atom (central)
          a2 = addsubmp[i][0]; //gaining bond to atom
          c1 = addsubmp[j][3]; //losing bonded atom
        }
        //printf(" a1,a2,c1: %i %i %i \n",a1,a2,c1);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
//        for (int k2=0;k2<k1;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1) 
            && (geoms[id].bonds[k2][0]==a1 || geoms[id].bonds[k2][1]==a1) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==c1) mov1=geoms[id].bonds[k1][1];
            else mov1=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==a1) mov2=geoms[id].bonds[k2][1];
            else mov2=geoms[id].bonds[k2][0];
          }
          if (rem1==rem2) skip = true;
          if (geoms[id].bond_exists_tri(mov1,a1) || geoms[id].bond_exists_tri(mov2,a2)) skip = true;
          if (a2==mov2 || a1==mov1) skip = true;
          //CPMZ allowing a single partly-attached H2 to form
          //if (geoms[id].hpair(a1,mov1)) skip = true;
          if (geoms[id].hpair(a1,mov1) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          if (geoms[id].hpair(a2,mov2) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          //if (geoms[id].hpair(a1,mov1) && geoms[id].hpair(a2,mov2)) skip = true;
          if(geoms[id].isTMpair(a1,mov1) || geoms[id].isTMpair(a2,mov2)) skip = true;

          if (frozen[mov1] || frozen[mov2]) skip = true;
          if (frozen[a1] || frozen[a2]) skip = true;
          if (!skip)
          {
            //printf(" found rem1 rem2: %i %i for atoms %i %i, %i %i \n",rem1,rem2,a1,mov1,a2,mov2);
            moves2pm[nmoves2pm][0] = rem1;
            moves2pm[nmoves2pm][1] = rem2;
            moves2pm[nmoves2pm][2] = a1;
            moves2pm[nmoves2pm][3] = mov1;
            moves2pm[nmoves2pm][4] = a2;
            moves2pm[nmoves2pm][5] = mov2;
            nmoves2pm++;
          } // if !skip

        } // loops k1,k2 over bonds
      } //case 1 repeated
      else if (addsubmp[i][0] == addsubmp[j][0] && ncmoves2p[addsubmp[i][0]] && addsubmp[i][3] != addsubmp[j][3])
      {
        printf("\n repeated 2x +: %i %i (rem centers: %i %i) \n",addsubmp[i][0],addsubmp[j][0],addsubmp[i][3],addsubmp[j][3]);
        a1 = addsubmp[i][0]; //gaining bond to atom (central)
        a2 = addsubmp[j][0]; //gaining bond to atom (central)
        c1 = addsubmp[i][3]; //losing bonded atom
        c2 = addsubmp[j][3]; //losing bonded atom
        //printf(" a1,a2,c1,c2: %i %i %i %i \n",a1,a2,c1,c2);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
//        for (int k2=0;k2<k1;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1) 
            && (geoms[id].bonds[k2][0]==c2 || geoms[id].bonds[k2][1]==c2) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==c1) mov1=geoms[id].bonds[k1][1];
            else mov1=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==c2) mov2=geoms[id].bonds[k2][1];
            else mov2=geoms[id].bonds[k2][0];
          }
          if (rem1==rem2) skip = true;
          if (geoms[id].bond_exists_tri(mov1,a1) || geoms[id].bond_exists_tri(mov2,a2)) skip = true;
          if (a2==mov2 || a1==mov1) skip = true;
          //CPMZ allowing a single partly-attached H2 to form
          //if (geoms[id].hpair(a1,mov1) || geoms[id].hpair(a2,mov2)) skip = true;
          //if (geoms[id].hpair(a1,mov1)) skip = true;
          if (geoms[id].hpair(a1,mov1) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          if (geoms[id].hpair(a2,mov2) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          //if (geoms[id].hpair(a1,mov1) && geoms[id].hpair(a2,mov2)) skip = true;
          if(geoms[id].isTMpair(a1,mov1) || geoms[id].isTMpair(a2,mov2)) skip = true;

          if (frozen[mov1] || frozen[mov2]) skip = true;
          if (frozen[a1] || frozen[a2]) skip = true;
          if (!skip)
          {
            //printf(" found rem1 rem2: %i %i for atoms %i %i, %i %i \n",rem1,rem2,a1,mov1,a2,mov2);
            moves2pm[nmoves2pm][0] = rem1;
            moves2pm[nmoves2pm][1] = rem2;
            moves2pm[nmoves2pm][2] = a1;
            moves2pm[nmoves2pm][3] = mov1;
            moves2pm[nmoves2pm][4] = a2;
            moves2pm[nmoves2pm][5] = mov2;
            nmoves2pm++;
          } // if !skip

        } // loops k1,k2 over bonds
      } //case 2x + on one center repeated
      else if (addsubmp[i][3] == addsubmp[j][3] && ncmoves2m[addsubmp[i][3]] && addsubmp[i][0] != addsubmp[j][0])
      {
        printf("\n repeated 2x +: %i %i (rem centers: %i %i) \n",addsubmp[i][0],addsubmp[j][0],addsubmp[i][3],addsubmp[j][3]);
        a1 = addsubmp[i][0]; //gaining bond to atom 
        a2 = addsubmp[j][0]; //gaining bond to atom 
        c1 = addsubmp[i][3]; //losing bonded atom (central)
        c2 = addsubmp[j][3]; //losing bonded atom (central)
        //printf(" a1,a2,c1,c2: %i %i %i %i \n",a1,a2,c1,c2);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
//        for (int k2=0;k2<k1;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1) 
            && (geoms[id].bonds[k2][0]==c2 || geoms[id].bonds[k2][1]==c2) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==c1) mov1=geoms[id].bonds[k1][1];
            else mov1=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==c2) mov2=geoms[id].bonds[k2][1];
            else mov2=geoms[id].bonds[k2][0];
          }
          if (rem1==rem2) skip = true;
          if (geoms[id].bond_exists_tri(mov1,a1) || geoms[id].bond_exists_tri(mov2,a2)) skip = true;
          if (a2==mov2 || a1==mov1) skip = true;
          //CPMZ allowing a single partly-attached H2 to form
          //if (geoms[id].hpair(a1,mov1) || geoms[id].hpair(a2,mov2)) skip = true;
          //if (geoms[id].hpair(a1,mov1)) skip = true;
          if (geoms[id].hpair(a1,mov1) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          if (geoms[id].hpair(a2,mov2) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
          //if (geoms[id].hpair(a1,mov1) && geoms[id].hpair(a2,mov2)) skip = true;
          if(geoms[id].isTMpair(a1,mov1) || geoms[id].isTMpair(a2,mov2)) skip = true;

          if (frozen[mov1] || frozen[mov2]) skip = true;
          if (frozen[a1] || frozen[a2]) skip = true;
          if (!skip)
          {
            //printf(" found rem1 rem2: %i %i for atoms %i %i, %i %i \n",rem1,rem2,a1,mov1,a2,mov2);
            moves2pm[nmoves2pm][0] = rem1;
            moves2pm[nmoves2pm][1] = rem2;
            moves2pm[nmoves2pm][2] = a1;
            moves2pm[nmoves2pm][3] = mov1;
            moves2pm[nmoves2pm][4] = a2;
            moves2pm[nmoves2pm][5] = mov2;
            nmoves2pm++;
          } // if !skip

        } // loops k1,k2 over bonds
      } //case 2x - on one center repeated
      //not doing 2x + on one center and 2x - on one center

    } //loop j over unique pairs of nasmp
  } //loop i over nasmp

  printf(" removing duplicates from %i nmoves2pm \n",nmoves2pm);
  for (int i=0;i<nmoves2pm;i++)
  for (int j=0;j<i;j++)
  {
    //printf(" checking duplicate %i %i: %i %i %i %i %i %i \n",i,j,moves2pm[j][0],moves2pm[j][1],moves2pm[j][2],moves2pm[j][3],moves2pm[j][4],moves2pm[j][5]);
    //printf(" against %i %i: %i %i %i %i %i %i \n",i,j,moves2pm[i][0],moves2pm[i][1],moves2pm[i][2],moves2pm[i][3],moves2pm[i][4],moves2pm[i][5]);
    if (moves2pm[i][0]==moves2pm[j][0] 
     && moves2pm[i][1]==moves2pm[j][1]
     && moves2pm[i][2]==moves2pm[j][4]
     && moves2pm[i][3]==moves2pm[j][5]
     && moves2pm[i][4]==moves2pm[j][2]
     && moves2pm[i][5]==moves2pm[j][3])
    {
      printf(" removing duplicate %i %i: %i %i %i %i %i %i \n",i,j,moves2pm[j][0],moves2pm[j][1],moves2pm[j][2],moves2pm[j][3],moves2pm[j][4],moves2pm[j][5]);
      nmoves2pm--;
      for (int k1=j;k1<nmoves2pm;k1++)
      {
        moves2pm[k1][0]=moves2pm[k1+1][0];
        moves2pm[k1][1]=moves2pm[k1+1][1];
        moves2pm[k1][2]=moves2pm[k1+1][2];
        moves2pm[k1][3]=moves2pm[k1+1][3];
        moves2pm[k1][4]=moves2pm[k1+1][4];
        moves2pm[k1][5]=moves2pm[k1+1][5];
      }
      i--; j--;
    }

  }

#endif

#if PPMM
  printf("\n finding 2x +/- moves which have one unsat. or sat. center +/- \n");
  for (int i=0;i<nasmp;i++)
  {
    for (int j=0;j<natoms;j++)
    {  
      int skip = true;
      int rem1,rem2;
      int a1,a2,a3,a4,c1,c2,c3,c4;
      a3 = a4 = c3 = c4 = 0;
      a1 = addsubmp[i][0];
      a2 = j;
      c1 = addsubmp[i][3];
      c2 = j;

     // a2/c2 must be min or max coord, non hydrogen
      if ( (!ncmovesm[a2] || !ncmovesp[a2]) && (geoms[id].coordn[a2]>1 || geoms[id].anumbers[a2]>1)
          && !frozen[a2] )
      {
        //printf(" a1,a2: %i %i c1,c2: %i %i \n",a1,a2,c1,c2);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1) 
            && (geoms[id].bonds[k2][0]==c2 || geoms[id].bonds[k2][1]==c2) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==c1) c3=geoms[id].bonds[k1][1];
            else c3=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==c2) c4=geoms[id].bonds[k2][1];
            else c4=geoms[id].bonds[k2][0];
          }

         //case one: a1 to c4, a2/c2 to c3
          if (a1!=c4 && a2!=c3 && !skip)
          {
            int skip1 = false;
            if (geoms[id].bond_exists_tri(a1,c4) || geoms[id].bond_exists_tri(a2,c3)) skip1 = true;
            if (a1==c4 || a2==c3) skip1 = true;
            //CPMZ allowing a single partly-attached H2 to form
            //if ( geoms[id].hpair(a1,c4) ) skip1 = true;
            //if ( geoms[id].hpair(a2,c3) ) skip1 = true;
            if (geoms[id].hpair(a1,c4) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
            if (geoms[id].hpair(a2,c3) && !geoms[id].bond_exists_tri(a1,a2)) skip = true;
            //if (geoms[id].hpair(a1,c4) && geoms[id].hpair(a2,c3)) skip = true;
            if(geoms[id].isTMpair(a1,c4) || geoms[id].isTMpair(a2,c3)) skip = true;

            if ( frozen[a1] || frozen[a2] ) skip1 = true;
            if ( frozen[c3] || frozen[c4] ) skip1 = true;

            if (!skip && !skip1)
            {
              //if ( geoms[id].hpair(a2,c3) || geoms[id].hpair(a1,c4) ) printf(" H-H bond formed \n");
              //printf(" c1,c2,c3,c4: %i %i %i %i ",c1,c2,c3,c4);
              //printf(" found rem1 rem2: %i %i a1,c4 %i %i, a2,c3 %i %i \n",rem1,rem2,a1,c4,a2,c3);
              moves2pm[nmoves2pm][0] = rem1;
              moves2pm[nmoves2pm][1] = rem2;
              moves2pm[nmoves2pm][2] = a1;
              moves2pm[nmoves2pm][3] = c4;
              moves2pm[nmoves2pm][4] = a2;
              moves2pm[nmoves2pm][5] = c3;
              nmoves2pm++;
            } // if !skip
          } // if a1!=c4
#if BADMOVE
       //May not use this case, currently breaks for B2H6
       //case two a1 to c3, a2/c2 to a3, c4 disconnects (doesn't follow disconnect to connect scheme)
          if (a1!=c3 && ncmovesm[c4] && !skip)
          {
            for (int k=0;k<natoms;k++)
            {
              a3 = k;
              int skip1 = false;
              if (geoms[id].bond_exists_tri(a1,c3) || geoms[id].bond_exists_tri(a2,a3)) skip1 = true;
              if (a2==a3) skip = true;
//            if ( geoms[id].hpair(a2,a3) ) skip1 = true;
//              if (!ncmovesp[a2]) skip1 = true;

              if (geoms[id].hpair(a1,c3))
              {
                if (c2==a1 || c4==a1)
                {
                  if (geoms[id].coordn[a1]!=1 || geoms[id].coordn[c3]!=1)
                    skip1 = true;
                }
                else skip1 = true;
              } //H2 checking
              if (!skip && !skip1)
              {
                printf(" c1,c2,c3,c4: %i %i %i %i ",c1,c2,c3,c4);
                printf(" found rem1 rem2: %i %i a1,c3 %i %i, a2,a3 %i %i \n",rem1,rem2,a1,c3,a2,a3);
                moves2pm[nmoves2pm][0] = rem1;
                moves2pm[nmoves2pm][1] = rem2;
                moves2pm[nmoves2pm][2] = a1;
                moves2pm[nmoves2pm][3] = c3;
                moves2pm[nmoves2pm][4] = a2;
                moves2pm[nmoves2pm][5] = a3;
                nmoves2pm++;
              } // if !skip
            } //loop k over natoms
          } // if a1!=c3 
#endif
        } // loop k1,k2 over all bond pairs

      } // atom not repeated

    } //loop j over unique pairs of nasmp
  } //loop i over nasmp
#endif


#if PPMM
  printf("\n\n finding 2x +/- moves which have two unsat. or sat. centers +/- \n");
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    {  
      int skip = true;
      int rem1,rem2;
      int a1,a2,a3,a4,c1,c2,c3,c4;
      a3 = a4 = c3 = c4 = 0;
      a1 = i;
      a2 = j;
      c1 = i;
      c2 = j;

     // a1/c1 and a2/c2 must be min or max coord, non hydrogen
      if ( (!ncmovesm[a1] || !ncmovesp[a1]) && (geoms[id].coordn[a1]>1 || geoms[id].anumbers[a1]>1) 
            && !frozen[a1] )
      if ( (!ncmovesm[a2] || !ncmovesp[a2]) && (geoms[id].coordn[a2]>1 || geoms[id].anumbers[a2]>1) 
            && !frozen[a1] )
      {
        printf(" a1,a2: %i %i c1,c2: %i %i \n",a1,a2,c1,c2);
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        for (int k2=0;k2<geoms[id].nbonds;k2++)
        {
          skip = true;
          if ( (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1) 
            && (geoms[id].bonds[k2][0]==c2 || geoms[id].bonds[k2][1]==c2) )
          {
            skip = false;
            rem1=k1; rem2=k2;
            if (geoms[id].bonds[k1][0]==c1) c3=geoms[id].bonds[k1][1];
            else c3=geoms[id].bonds[k1][0];
            if (geoms[id].bonds[k2][0]==c2) c4=geoms[id].bonds[k2][1];
            else c4=geoms[id].bonds[k2][0];
          }
          if (geoms[id].bond_exists_tri(a1,c4) || geoms[id].bond_exists_tri(a2,c3)) skip = true;
          if (a1==c4 || a2==c3) skip = true;
          if (rem1==rem2) skip = true;
         // conditions for H-H bonding (never)
          if (geoms[id].hpair(a1,c4) || geoms[id].hpair(a2,c3)) skip = true;
          //printf(" a1,mov1 %i %i a2,mov2 %i %i \n",a1,c3,a2,c4);
          if(geoms[id].isTMpair(a1,c4) || geoms[id].isTMpair(a2,c3)) skip = true;

          if (frozen[a1] || frozen[a2]) skip = true;
          if (frozen[c3] || frozen[c4]) skip = true;
          if (!skip)
          {
            printf(" found rem1 rem2: %i %i for atoms %i %i, %i %i \n",rem1,rem2,a1,c4,a2,c3);
            moves2pm[nmoves2pm][0] = rem1;
            moves2pm[nmoves2pm][1] = rem2;
            moves2pm[nmoves2pm][2] = a1;
            moves2pm[nmoves2pm][3] = c4;
            moves2pm[nmoves2pm][4] = a2;
            moves2pm[nmoves2pm][5] = c3;
            nmoves2pm++;
          } // if !skip
        } // loop k1,k2 over all bond pairs

      } // atom not repeated

    } //loop j over unique pairs of nasmp
  } //loop i over nasmp
#endif


  printf("\n doing %i nmoves2pm \n",nmoves2pm);
  for (int i=0;i<nmoves2pm;i++)
  {
      printf("\n rem1,rem2: %i %i new bonds: %i %i, %i %i\n",moves2pm[i][0],moves2pm[i][1],
      moves2pm[i][2],moves2pm[i][3],moves2pm[i][4],moves2pm[i][5]);
     
      geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
      geom1.ic_create();
      //geom1.print_bonds(); 
      geom1.nbonds = geoms[id].nbonds;
      geom1.bonds[moves2pm[i][0]][0]=moves2pm[i][2];
      geom1.bonds[moves2pm[i][0]][1]=moves2pm[i][3];
      geom1.bonds[moves2pm[i][1]][0]=moves2pm[i][4];
      geom1.bonds[moves2pm[i][1]][1]=moves2pm[i][5];
      geom1.ic_create_nobonds(); 
      //geom1.mm_init();
      geom1.update_ic();
      //geom1.print_bonds(); 
      nfrz = 4;
      frzlist[0] = moves2pm[i][2];
      frzlist[1] = moves2pm[i][3];
      frzlist[2] = moves2pm[i][4];
      frzlist[3] = moves2pm[i][5];
      nfrz0 = nfrz;
      for (int j=0;j<geom1.nbonds;j++)
      {
        for (int k=0;k<nfrz0;k++)
        {
          if (frzlist[k]==geom1.bonds[j][0])
          { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
          if (frzlist[k]==geom1.bonds[j][1])
          { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
        }
      }

      string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
      string optfile = "scratch/xyzfile.xyz"+nstr;
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
#if !OPTSHADOW
      done = geom1.opt(optfile);
#else
      geomshadow.reset(geoms[id].coords);
      done = geom1.opt(optfile,geomshadow);
#endif

      if (done)
      {
        //printf(" attempting mopac \n");
        optfile = "scratch/xyzfile.xyzm_"+nstr;
        mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
        mopt.freeze(frzlist, nfrz, nfrz0);
        printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
        mopt.opt_write(optfile);
        tniso++;
        for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
      }
      else { nfail++; nscratch--; }
    } // loop i over allowed moves





#if 1
// for 2x + moves
  printf("\n doing %i naspp moves: \n",naspp);
  for (int i=0;i<naspp;i++)
  {
    int tbonds = geoms[id].nbonds;
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    geom1.nbonds = geom1.nbonds+2;
    geom1.bonds[tbonds][0]=addsubpp[i][0];
    geom1.bonds[tbonds][1]=addsubpp[i][1];
    geom1.bonds[tbonds+1][0]=addsubpp[i][2];
    geom1.bonds[tbonds+1][1]=addsubpp[i][3];

    geom1.ic_create_nobonds(); 
    //geom1.mm_init();
    geom1.update_ic();
//    geom1.print_bonds();
    nfrz = 4;
    frzlist[0] = addsubpp[i][0];
    frzlist[1] = addsubpp[i][1];
    frzlist[2] = addsubpp[i][2];
    frzlist[3] = addsubpp[i][3];
    nfrz0 = nfrz;
    for (int j=0;j<geom1.nbonds;j++)
    {
      for (int k=0;k<nfrz0;k++)
      {
       	if (frzlist[k]==geom1.bonds[j][0])
        { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
        if (frzlist[k]==geom1.bonds[j][1])
        { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
      }
    }
	
    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow);
#endif

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }
  }
#endif

#if 1
// for 2x - moves
  printf("\n doing %i nasmm moves: \n",nasmm);
  for (int i=0;i<nasmm;i++)
  {
    int tbonds = geoms[id].nbonds;
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    int rem1 = addsubmm[i][0];
    int rem2 = addsubmm[i][1];

//    geom1.print_bonds();
    printf(" removing bonds: %i %i \n",rem1,rem2);
    for (int l=rem1;l<geom1.nbonds-1;l++)
    {
//      printf(" l: %i b2: %i nbonds: %i \n",l,rem1,geom1.nbonds);
      geom1.bonds[l][0] = geom1.bonds[l+1][0];
      geom1.bonds[l][1] = geom1.bonds[l+1][1];
    }
    for (int l=rem2;l<geom1.nbonds-1;l++)
    {
//      printf(" l: %i b2: %i nbonds: %i \n",l,rem2,geom1.nbonds);
      geom1.bonds[l][0] = geom1.bonds[l+1][0];
      geom1.bonds[l][1] = geom1.bonds[l+1][1];
    }
    geom1.nbonds = geom1.nbonds-2;

    geom1.ic_create_nobonds(); 
    //geom1.mm_init();
    geom1.update_ic();
//    geom1.print_bonds();
//    geom1.print_ic();


    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
    printf(" creating %s ",optfile.c_str());
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow);
#endif

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
//      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*natoms0;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }
  }
  printf("\n");
#endif




// for 2x - + moves 
  int n2m1p = 0;
  int** moves2m1p = new int*[nasmp*natoms*4+1000];
  for (int i=0;i<nasmp*natoms*4+1000;i++)
    moves2m1p[i] = new int[4];
#if 1
  printf("\n finding 2x - + moves \n");
  for (int i=0;i<nasmp;i++)
  {
    int a1,a2,a3,a13,c1,c2,c3,c4;
    int rem1,rem2;
    a1 = addsubmp[i][0];
    c1 = addsubmp[i][3];
//    for (int j=0;j<=a1;j++)
    for (int j=0;j<=c1;j++)
    {
      c2 = j;

    //would overcount 
   // for (int k2=0;k2<geoms[id].nbonds;k2++) 

      if ( ncmovesm[c2] || (c1==c2 && ncmoves2m[c2]) )
      for (int k1=0;k1<geoms[id].nbonds;k1++)
      for (int k2=0;k2<k1;k2++)
      {
        c3 = c4 = -1;
        rem1 = k1; rem2 = k2;
        if (geoms[id].bonds[k1][0]==c1)
          c3 = geoms[id].bonds[k1][1];
        else if (geoms[id].bonds[k1][1]==c1)
          c3 = geoms[id].bonds[k1][0];
        if (geoms[id].bonds[k2][0]==c2)
          c4 = geoms[id].bonds[k2][1];
        else if (geoms[id].bonds[k2][1]==c2)
          c4 = geoms[id].bonds[k2][0];
        
        int skip = true;
        a2 = 0;
        a3 = -1;
        //printf(" c1,c2,c3,c4 %i %i %i %i ",c1,c2,c3,c4);
        if (c3>-1 && c4>-1 && ncmovesm[c3] && ncmovesm[c4])
        {
          skip = false;
          a2 = c3;
        }
        else if (c3>-1 && c4>-1 && ncmovesm[c3] && !ncmovesm[c4]) 
        {
          skip = false;
          a2 = c4;
          if (!ncmovesp[a2] && geoms[id].coordn[a2]>1) skip = true;
        }
        else if (c3>-1 && c4>-1 && ncmovesm[c4] && !ncmovesm[c3]) 
        {
          skip = false;
          a2 = c3;
        }
        else if (c3>-1 && c4>-1 && !ncmovesm[c3] && !ncmovesm[c3]) 
        {
          skip = false;
          a2 = c3;
          a3 = c4;
        }

//not needed        if (rem1==rem2) skip = true;
        if (c1==c4 && !ncmoves2m[c1]) skip = true;
        if (c2==c3 && !ncmoves2m[c2]) skip = true;
        if (c3==c4 && !ncmoves2m[c3]) skip = true;
        if (a3>-1) a13 = a3; else a13 = a1;
       //make sure H2 only forms if two R-H bonds are broken
        if (geoms[id].hpair(a13,a2))
        { 
//          printf(" found hpair \n");
          if (c1==a13 || c3==a13)
          {
            if (c2==a2 || c4==a2)
            {
              if (geoms[id].coordn[a13]!=1 || geoms[id].coordn[a2]!=1)
                skip = true;
            }
            else skip = true;
          }
          else if (c2==a13 || c4==a13)
          {
            if (c1==a2 || c3==a2)
            {
              if (geoms[id].coordn[a13]!=1 || geoms[id].coordn[a2]!=1)
                skip = true;
            }
            else skip = true;
          }
          else skip = true;
        } //H2 checking

        if (!ncmovesp[a2] && geoms[id].coordn[a2]>1) skip = true;
// apparently doesn't matter
//        if (a2==c2) skip = true;
        if (a13==a2) skip = true;
        if (geoms[id].bond_exists_tri(a13,a2)) skip = true;
        if (geoms[id].isTMpair(a13,a2)) skip = true;

        if (frozen[a2] || frozen[a13]) skip = true;
        if (!skip)
        {
          printf(" c1,c2,c3,c4 %i %i %i %i ",c1,c2,c3,c4);
          printf(" saving rem1,rem2: %i %i a13,a2: %i %i \n",rem1,rem2,a13,a2);
          moves2m1p[n2m1p][0]=rem1;
          moves2m1p[n2m1p][1]=rem2;
          moves2m1p[n2m1p][2]=a13;
          moves2m1p[n2m1p][3]=a2;
          n2m1p++;
        } //if !skip
      } //loop k1,k2 over all unique bond pairs
    } // loop j over natoms
  } //loop i over nasmp

// removes duplicates in n2m1p (double counting over a1 when a3 is used)
  for (int i=0;i<n2m1p;i++)
  for (int j=0;j<i;j++)
  {
    if (moves2m1p[i][0]==moves2m1p[j][0] 
     && moves2m1p[i][1]==moves2m1p[j][1]
     && moves2m1p[i][2]==moves2m1p[j][2]
     && moves2m1p[i][3]==moves2m1p[j][3])
    {
      printf(" removing duplicate %i %i: %i %i %i %i \n",i,j,moves2m1p[j][0],moves2m1p[j][1],moves2m1p[j][2],moves2m1p[j][3]);
      n2m1p--;
      for (int k1=j;k1<n2m1p;k1++)
      {
        moves2m1p[k1][0]=moves2m1p[k1+1][0];
        moves2m1p[k1][1]=moves2m1p[k1+1][1];
        moves2m1p[k1][2]=moves2m1p[k1+1][2];
        moves2m1p[k1][3]=moves2m1p[k1+1][3];
      }
      j--; i--;
    }
  }
#endif


#if 1
// 2x - + second case for coord min saturated +/- centers
  printf("\n finding min sat centers 2x - + \n");
  for (int i=0;i<natoms;i++)
  {
    int a1,a2,a3,a13,c1,c2,c3,c4;
    int rem1,rem2;
    a1 = i;
    c1 = a1;
    if ( !ncmovesm[a1] && (geoms[id].anumbers[a1]>1 || geoms[id].coordn[a1]>1) 
         && !frozen[a1] )
    for (int j=0;j<natoms;j++)
    {
      c2 = j;

     //only do lowest sat a1==c1 center
      if (ncmovesm[c2])
      for (int k1=0;k1<geoms[id].nbonds;k1++)
      for (int k2=0;k2<geoms[id].nbonds;k2++)
      {
        int skip = false;

        c3 = c4 = -1;
        rem1 = k1; rem2 = k2;
        if (geoms[id].bonds[k1][0]==c1)
          c3 = geoms[id].bonds[k1][1];
        else if (geoms[id].bonds[k1][1]==c1)
          c3 = geoms[id].bonds[k1][0];
        if (geoms[id].bonds[k2][0]==c2)
          c4 = geoms[id].bonds[k2][1];
        else if (geoms[id].bonds[k2][1]==c2)
          c4 = geoms[id].bonds[k2][0];

        a2 = c4;
        if (rem1==rem2) skip = true;
        if (c1==c2 || c3==c4) skip = true;
        if (a1==a2) skip = true;
        if (!ncmovesm[c3]) skip = true;
        if (geoms[id].isTMpair(a1,a2)) skip = true;
        if (frozen[c1] || frozen[c2] || frozen[c3] || frozen[c4] || frozen[a2]) skip = true;

//        if (!ncmovesm[c4]) skip = true;
//        if (c3>-1 && c4>-1 && skip)
//          printf(" now c1,c2,c3,c4: %i %i %i %i, rem1,rem2 %i %i a1,a2 %i %i \n",c1,c2,c3,c4,rem1,rem2,a1,a2);
        if (c3>-1 && c4>-1 && !skip)
        {
          printf(" found c1,c2,c3,c4: %i %i %i %i, rem1,rem2 %i %i a1,a2 %i %i \n",c1,c2,c3,c4,rem1,rem2,a1,a2);
          moves2m1p[n2m1p][0]=rem1;
          moves2m1p[n2m1p][1]=rem2;
          moves2m1p[n2m1p][2]=a1;
          moves2m1p[n2m1p][3]=a2;
          n2m1p++;
        }
      } // loops k1,k2 over bonds
    } // loop j over natoms
  } // loop i over natoms
#endif



//  for (int i=0;i<n2m1p;i++)
//    printf(" found moves2m1p[%i]: rem1,rem2: %i %i a1,mov2: %i %i \n",i,moves2m1p[i][0],moves2m1p[i][1],moves2m1p[i][2],moves2m1p[i][3]);


  printf("\n\n performing %i moves2m1p \n",n2m1p);
  for (int i=0;i<n2m1p;i++)
  {
    printf("\n found moves2m1p[%i]: rem1,rem2: %i %i a1,mov2: %i %i \n",i,moves2m1p[i][0],moves2m1p[i][1],moves2m1p[i][2],moves2m1p[i][3]);
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    geom1.bonds[moves2m1p[i][0]][0]=moves2m1p[i][2];
    geom1.bonds[moves2m1p[i][0]][1]=moves2m1p[i][3];
    geom1.nbonds--;
    for (int l=moves2m1p[i][1];l<geom1.nbonds;l++)
    {
//     printf(" l: %i b2: %i nbonds: %i \n",l,b2,geom1.nbonds);
      geom1.bonds[l][0] = geom1.bonds[l+1][0];
      geom1.bonds[l][1] = geom1.bonds[l+1][1];
    }
    geom1.ic_create_nobonds(); 
    //geom1.mm_init();
    geom1.update_ic();
    //geom1.print_bonds();
    nfrz = 2;
    frzlist[0] = moves2m1p[i][2];
    frzlist[1] = moves2m1p[i][3];
    nfrz0 = nfrz;
    for (int j=0;j<geom1.nbonds;j++)
    {
      for (int k=0;k<nfrz0;k++)
      {
       	if (frzlist[k]==geom1.bonds[j][0])
        { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
        if (frzlist[k]==geom1.bonds[j][1])
        { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
      }
    }

    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow);
#endif

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*geom1.natoms;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }
  } // loop i over n2m1p




// third case for +/- at one center with one minus move
  int n_nmp2 = 0;
  int** movesnmp2 = new int*[natoms*100+1000];
  for (int i=0;i<natoms*100+1000;i++) 
    movesnmp2[i] = new int[4];
#if 1
  printf("\n +/- - at one saturated center \n");
  for (int i=0;i<natoms;i++)
  {
    if (geoms[id].coordn[i]>1 && !ncmovesp[i] && !frozen[i])
    for (int j=0;j<natoms;j++)
    {
      if (i!=j && !geoms[id].bond_exists_tri(i,j) && !frozen[j])
      {
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        {
        //labels are wrong!
          int a1 = -1; // breaking this atom from bond k1 (was attached to i)
          int b2 = -1; // breaks from j (j attaches to i)
          int c1,c2,c3,c4; //complements to subtracted atoms 
          if (geoms[id].bonds[k1][0]==i)
            a1 = geoms[id].bonds[k1][1];
          else if(geoms[id].bonds[k1][1]==i)
            a1 = geoms[id].bonds[k1][0];
          c1 = a1;

//old          if (a1>-1 && j!=a1 && ncmovesm[a1])
          if ( a1>-1 && ((ncmovesm[a1] && j!=a1) || (ncmoves2m[a1])) )
          {
            //printf(" in i==j: %i %i \n",i,j);
            for (int k2=0;k2<geoms[id].nbonds;k2++)
            {
              b2 = -1;
              if      (geoms[id].bonds[k2][0] == j && ncmovesm[geoms[id].bonds[k2][1]])
              {
                b2 = k2;
                c2 = geoms[id].bonds[k2][1];
              }
              else if (geoms[id].bonds[k2][1] == j && ncmovesm[geoms[id].bonds[k2][0]])
              {
                b2 = k2;  
                c2 = geoms[id].bonds[k2][0];
              }
              if (geoms[id].hpair(i,j)) b2 = -1;
              if (geoms[id].isTMpair(i,j)) b2 = -1;

//double check this (should be fine now)
              if ( b2>-1 && ncmovesm[c2] && (c1!=c2 || ncmoves2m[c1]) )
              {
                //printf("\n c1,c2 %i %i ",c1,c2);
                printf(" found +/- - combo: %i %i %i %i %i \n",i,j,a1,b2,k1);
//removing b2, k2 bonds, make bond between i,j
                movesnmp2[n_nmp2][0]=k1;
                movesnmp2[n_nmp2][1]=i;
                movesnmp2[n_nmp2][2]=j;
                movesnmp2[n_nmp2][3]=b2;
                n_nmp2++;
              } // if b2

            } // loop m over nbonds
          } // if b1>-1 and j!=b1
        } //loop over k1
      } // if ncmovesm[j]
    } // loop j over natoms
  } //loop i over natoms
#endif


// 2x - + when nasmp == 0 (when it cannot find + site above)
#if 1
  if (nasmp==0)
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  if (!ncmovesm[i] && !ncmovesm[j] && !frozen[i] && !frozen[j])
  {
    int a1,a2,c1,c2,c3,c4;
    int rem1,rem2;
    a1 = c1 = i;
    a2 = c2 = j;

    for (int k1=0;k1<geoms[id].nbonds;k1++)
    for (int k2=0;k2<geoms[id].nbonds;k2++)
    {
      int skip = false;

      c3 = c4 = -1;
      rem1 = k1; rem2 = k2;
      if (geoms[id].bonds[k1][0]==c1)
        c3 = geoms[id].bonds[k1][1];
      else if (geoms[id].bonds[k1][1]==c1)
        c3 = geoms[id].bonds[k1][0];
      if (geoms[id].bonds[k2][0]==c2)
        c4 = geoms[id].bonds[k2][1];
      else if (geoms[id].bonds[k2][1]==c2)
        c4 = geoms[id].bonds[k2][0];

      if (!ncmovesm[c3] || !ncmovesm[c4]) skip = true;
      if (c3==c4 && !ncmoves2m[c3]) skip = true;
      if (rem1==rem2) skip = true;
      if (geoms[id].isTMpair(a1,a2)) skip = true;

      if (c3>-1 && c4>-1 && !skip)
      {
        printf(" found 2x - + rem1,rem2: %i %i a1/c1,a2/c2: %i %i c3,c4: %i %i  \n",rem1,rem2,c1,c2,c3,c4);
        movesnmp2[n_nmp2][0]=rem1;
        movesnmp2[n_nmp2][1]=a1;
        movesnmp2[n_nmp2][2]=a2;
        movesnmp2[n_nmp2][3]=rem2;
        n_nmp2++;
      }
    } // loops k1,k2 over bonds

  } //loop i,j over unique pairs of atoms
#endif

  printf(" found %i nmp2: +/- - at one sat. center or 2x min. coord - \n",n_nmp2);
  for (int i=0;i<n_nmp2;i++)
  {
    printf("\n found +/- - combo: %i %i %i %i  \n",movesnmp2[i][1],movesnmp2[i][2],movesnmp2[i][3],movesnmp2[i][0]);
//    do_nmp(geoms[id],geom1,i,j,b1,b2,k1);
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    geom1.bonds[movesnmp2[i][0]][0]=movesnmp2[i][1];
    geom1.bonds[movesnmp2[i][0]][1]=movesnmp2[i][2];
    geom1.nbonds--;
    for (int l=movesnmp2[i][3];l<geom1.nbonds;l++)
    {
//    printf(" l: %i b2: %i nbonds: %i \n",l,b2,geom1.nbonds);
      geom1.bonds[l][0] = geom1.bonds[l+1][0];
      geom1.bonds[l][1] = geom1.bonds[l+1][1];
    }
    geom1.ic_create_nobonds(); 
    //geom1.mm_init();
    geom1.update_ic();
  //geom1.print_bonds();
    nfrz = 2;
    frzlist[0] = movesnmp2[i][1];
    frzlist[1] = movesnmp2[i][2];
    nfrz0 = nfrz;
    for (int j=0;j<geom1.nbonds;j++)
    {
      for (int k=0;k<nfrz0;k++)
      {
       	if (frzlist[k]==geom1.bonds[j][0])
        { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
        if (frzlist[k]==geom1.bonds[j][1])
        { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
      }
    }

    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow);
#endif

    if (done)
    {
     //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*geom1.natoms;l++) isos[tniso][l] = mopt.xyz[l]; 
    }
    else { nfail++; nscratch--; }
  } // loop i over n_nmp2








  int n2p1m = 0;
  int** moves2p1m = new int*[nasmp*natoms*20+1000];
  for (int i=0;i<nasmp*natoms*20+1000;i++)
    moves2p1m[i] = new int[5];

#if 1
// for 2x + - moves 
  printf("\n finding 2x + - moves \n");
  for (int i=0;i<nasmp;i++)
  {
//    printf(" in mp pair: %i %i \n",addsubmp[i][0],addsubmp[i][3]);
    int skip;
    for (int j=0;j<natoms;j++)
    for (int k=0;k<j;k++)
    { 
      int a1,mov1,a2,a3,a4;
      int c1; //atom will detach from c1
      int rem1;
      //adding two bonds, destroying one
      c1 = addsubmp[i][3];
      a1 = addsubmp[i][0]; //first attach position

      skip = false;
      a3 = j; a4 = k;
      if (geoms[id].bond_exists_tri(a3,a4)) skip = true;
      if (a3==a1 && !ncmoves2p[a1]) skip = true;
      if (a4==a1 && !ncmoves2p[a1]) skip = true;
//  line below definitely wrong
//      if ( (a3==c1 || a4==c1)) skip = true;
      if (!ncmovesp[a3]) skip = true;
      if (!ncmovesp[a4]) skip = true;

//      printf(" skip? %i a1,mov1,a3,a4: %i n/a %i %i c1: %i \n",skip,a1,a3,a4,c1);
//      printf(" ncmovesp[a3]: %i ncmovesp[a4] %i \n",ncmovesp[a3],ncmovesp[a4]);

      if (!skip)
      for (int k1=0;k1<geoms[id].nbonds;k1++)
      {
        skip = false;
        if (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1)
        {
          rem1 = k1;
          if (geoms[id].bonds[k1][0]==c1) 
            mov1=geoms[id].bonds[k1][1];
          else 
            mov1=geoms[id].bonds[k1][0];
     
          if (geoms[id].bond_exists_tri(a1,mov1)) skip = true;
          if (mov1==a1) skip = true;
          if (mov1==a3) skip = true;
          if (mov1==a4) skip = true;
          if (geoms[id].hpair(a1,mov1) || geoms[id].hpair(a3,a4)) skip = true;
          if (geoms[id].isTMpair(a1,mov1) || geoms[id].isTMpair(a3,a4)) skip = true;

          if (!skip)
          {
//            printf(" saving rem1: %i a1,mov1,a3,a4: %i %i %i %i c1: %i \n",rem1,a1,mov1,a3,a4,c1);
            moves2p1m[n2p1m][0]=rem1;
            moves2p1m[n2p1m][1]=a1;
            moves2p1m[n2p1m][2]=mov1;
            moves2p1m[n2p1m][3]=a3;
            moves2p1m[n2p1m][4]=a4;
            n2p1m++;
          } // if !skip

        } //if c1 is in bond k1
      } // loop k1 over nbonds
    } // loops j,k over unique atom pairs

  } //loop i over nasmp
#endif

#if 1
  printf(" finding +/- min coord centers with another + \n");
  for (int i=0;i<natoms;i++)
  {
//    printf(" central atom for +/-: %i \n",i);
    int skip;

    if ( !ncmovesm[i] && (geoms[id].coordn[i]>1 || geoms[id].anumbers[i]>1) 
        && !frozen[i])
    for (int j=0;j<i;j++)
    for (int k=0;k<j;k++)
    { 
      int a1,mov1,a2,a3,a4;
      int c1; //atom will detach from c1
      int rem1;
      //adding two bonds, destroying one
      c1 = i;
      a1 = i; //first attach position

      skip = false;
      a3 = j; a4 = k;
      if (geoms[id].bond_exists_tri(a1,a3)) skip = true;
      if (!ncmovesp[a3]) skip = true;
      if (!ncmovesp[a4]) skip = true;

//      printf(" skip? %i a1,mov1,a3,a4: %i n/a %i %i c1: %i \n",skip,a1,a3,a4,c1);
//      printf(" ncmovesp[a3]: %i ncmovesp[a4] %i \n",ncmovesp[a3],ncmovesp[a4]);

      if (!skip)
      for (int k1=0;k1<geoms[id].nbonds;k1++)
      {
        skip = false;
        if (geoms[id].bonds[k1][0]==c1 || geoms[id].bonds[k1][1]==c1)
        {
          rem1 = k1;
          if (geoms[id].bonds[k1][0]==c1) 
            mov1=geoms[id].bonds[k1][1];
          else 
            mov1=geoms[id].bonds[k1][0];
     
          if (geoms[id].bond_exists_tri(mov1,a4)) skip = true;
          if (mov1==a1) skip = true;
          if (mov1==a3) skip = true;
          if (mov1==a4) skip = true;
          if (geoms[id].hpair(a1,a3) || geoms[id].hpair(mov1,a4)) skip = true;
          if (geoms[id].isTMpair(a1,a3) || geoms[id].isTMpair(mov1,a4)) skip = true;

          if (!skip)
          {
//            printf(" saving rem1: %i a1,mov1,a3,a4: %i %i %i %i c1: %i \n",rem1,a1,mov1,a3,a4,c1);
            moves2p1m[n2p1m][0]=rem1;
            moves2p1m[n2p1m][1]=a1;
            moves2p1m[n2p1m][2]=a3;
            moves2p1m[n2p1m][3]=mov1;
            moves2p1m[n2p1m][4]=a4;
            n2p1m++;
          } // if !skip

        } //if c1 is in bond k1
      } // loop k1 over nbonds
    } // loops j,k over unique atom pairs

  } //loop i over natoms
#endif


  for (int i=0;i<n2p1m;i++)
    printf(" found moves2p1m[%i]: rem1: %i a1,mov1,a3,a4: %i %i %i %i \n",i,moves2p1m[i][0],moves2p1m[i][1],moves2p1m[i][2],moves2p1m[i][3],moves2p1m[i][4]);

  printf("\n performing %i moves2p1m \n",n2p1m);
  for (int i=0;i<n2p1m;i++)
  {
//  printf(" found moves2p1m[%i]: rem1,rem2: %i %i a1,mov2: %i %i \n",i,moves2m1p[i][0],moves2m1p[i][1],moves2m1p[i][2],moves2m1p[i][3]);
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    geom1.bonds[moves2p1m[i][0]][0]=moves2p1m[i][1];
    geom1.bonds[moves2p1m[i][0]][1]=moves2p1m[i][2];
    geom1.bonds[geom1.nbonds][0]=moves2p1m[i][3];
    geom1.bonds[geom1.nbonds][1]=moves2p1m[i][4];
    geom1.nbonds++;

    geom1.ic_create_nobonds(); 
    //geom1.mm_init();
    geom1.update_ic();
    //geom1.print_bonds();
    nfrz = 4;
    frzlist[0] = moves2p1m[i][1];
    frzlist[1] = moves2p1m[i][2];
    frzlist[2] = moves2p1m[i][3];
    frzlist[3] = moves2p1m[i][4];
    nfrz0 = nfrz;
    for (int j=0;j<geom1.nbonds;j++)
    {
      for (int k=0;k<nfrz0;k++)
      {
       	if (frzlist[k]==geom1.bonds[j][0])
        { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
        if (frzlist[k]==geom1.bonds[j][1])
        { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
      }
    }

    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow);
#endif

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*geom1.natoms;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }
  } // loop i over n2p1m




// for +/- at one sat. center + plus moves
  int n_nmp1 = 0;
  int** movesnmp1 = new int*[natoms*100+1000];
  for (int i=0;i<natoms*100+1000;i++) 
    movesnmp1[i] = new int[5];
#if 1
  printf("\n finding +/- + moves at sat. center \n");
  for (int i=0;i<natoms;i++)
  {
    if (geoms[id].coordn[i]>1 && !ncmovesp[i] && !frozen[i])
    for (int j=0;j<natoms;j++)
    {
      if (ncmovesp[j] && i!=j)
      {
        for (int k1=0;k1<geoms[id].nbonds;k1++)
        {
          int b1 = -1; // attaches to j by breaking k1 (was attached to i)
          int b2 = -1; // attaches to i
          if (geoms[id].bonds[k1][0]==i) b1 = geoms[id].bonds[k1][1];
          else if(geoms[id].bonds[k1][1]==i) b1 = geoms[id].bonds[k1][0];

          for (int m=0;m<geoms[id].nbonds;m++)
            if ((geoms[id].bonds[m][0] == b1 && geoms[id].bonds[m][1] == j)
             || (geoms[id].bonds[m][1] == b1 && geoms[id].bonds[m][0] == j))
              b1 = -1;

//need to prevent same bond formation
          if (b1>-1 && j!=b1)
          {
            for (int l=0;l<natoms;l++)
            {
              if (ncmovesp[l] && i!=l && j!=l)
              {
                b2 = l;
               //make sure b2 is not bonded to i
                for (int m=0;m<geoms[id].nbonds;m++)
                  if ((geoms[id].bonds[m][0] == l && geoms[id].bonds[m][1] == i)
                   || (geoms[id].bonds[m][1] == l && geoms[id].bonds[m][0] == i))
                    b2 = -1;
                if (geoms[id].hpair(i,b2) || geoms[id].hpair(j,b1)) b2 = -1;
                if (geoms[id].isTMpair(i,b2) || geoms[id].isTMpair(j,b1)) b2 = -1;
                if (b2>-1)
                {
                 printf(" found +/- + combo: %i %i %i %i %i \n",i,j,b1,b2,k1);
                 movesnmp1[n_nmp1][0] = k1;
                 movesnmp1[n_nmp1][1] = i;
                 movesnmp1[n_nmp1][2] = b2;
                 movesnmp1[n_nmp1][3] = j;
                 movesnmp1[n_nmp1][4] = b1;
                 n_nmp1++;

                } // if b2

              } // if ncmovesm[l]
            } // loop l over natoms
          }
        }
      } // if ncmovesm[j]
    } // loop j over natoms
  } //loop i over natoms
#endif


// for 2+ at one/two center with - at two bonded min coord centers
// picks up slack where nasmp didn't find minus possibility
#if 1
  printf(" finding 2x + - where - connects two min coord centers \n");
  for (int i=0;i<geoms[id].nbonds;i++)
  {
    int a1,a2,c1,c2;
    int rem1;
    c1 = geoms[id].bonds[i][0];
    c2 = geoms[id].bonds[i][1];
    if (!ncmovesm[c1] && !ncmovesm[c2] && !frozen[c1] && !frozen[c2])
    {
      rem1 = i;
      for (int j=0;j<natoms;j++)
      for (int k=0;k<natoms;k++)
      {
        int skip = true;
        a1 = j;
        a2 = k;

        if (ncmovesp[a1] && ncmovesp[a2]) skip = false;
        if (geoms[id].hpair(c1,a1) || geoms[id].hpair(c2,a2)) skip = true;
        if (geoms[id].isTMpair(c1,a1) || geoms[id].isTMpair(c2,a2)) skip = true;
        if (a1==a2 && !ncmoves2p[a1]) skip = true;
//        if (skip)
//          printf(" skip found rem1: %i connecting c1,c2: %i %i a1,a2: %i %i \n",rem1,c1,c2,a1,a2);
        if (!skip)
        {
          printf(" found rem1: %i connecting c1,c2: %i %i a1,a2: %i %i \n",rem1,c1,c2,a1,a2);
          movesnmp1[n_nmp1][0] = rem1;
          movesnmp1[n_nmp1][1] = c1;
          movesnmp1[n_nmp1][2] = a1;
          movesnmp1[n_nmp1][3] = c2;
          movesnmp1[n_nmp1][4] = a2;
          n_nmp1++;
        }
      } //loop j,k over all atom pairs
    } // if bond connects two min coord. centers
  } // loop i over bonds
#endif

  printf("\n performing %i +/- + at sat. center or 2x min. coord moves \n",n_nmp1);
  for (int i=0;i<n_nmp1;i++) 
  {
    printf("\n found +/- + combo: %i %i %i %i rem: %i  \n",movesnmp1[i][1],movesnmp1[i][2],movesnmp1[i][3],movesnmp1[i][4],movesnmp1[i][0]);
//  do_nmp(geoms[id],geom1,i,j,b1,b2,k1);
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords);
    geom1.ic_create();
    geom1.bonds[movesnmp1[i][0]][0]=movesnmp1[i][1];
    geom1.bonds[movesnmp1[i][0]][1]=movesnmp1[i][2];
    geom1.bonds[geom1.nbonds][0]=movesnmp1[i][3];
    geom1.bonds[geom1.nbonds][1]=movesnmp1[i][4];
    geom1.nbonds++;
    geom1.ic_create_nobonds(); 
    //geom1.mm_init();
    geom1.update_ic();
  //geom1.print_ic();
    nfrz = 4;
    frzlist[0] = movesnmp1[i][1];
    frzlist[1] = movesnmp1[i][2];
    frzlist[2] = movesnmp1[i][3];
    frzlist[3] = movesnmp1[i][4];
    nfrz0 = nfrz;
    for (int j=0;j<geom1.nbonds;j++)
    {
      for (int k=0;k<nfrz0;k++)
      {
       	if (frzlist[k]==geom1.bonds[j][0])
        { frzlist[nfrz] = geom1.bonds[j][1]; nfrz++; }
        if (frzlist[k]==geom1.bonds[j][1])
        { frzlist[nfrz] = geom1.bonds[j][0]; nfrz++; }
      }
    }

    string nstr=StringTools::int2str(nscratch++ +1,4,"0"); // 2,"0" is number of 0s in total
    string optfile = "scratch/xyzfile.xyz"+nstr;
#if !OPTSHADOW
    done = geom1.opt(optfile);
#else
    geomshadow.reset(geoms[id].coords);
    done = geom1.opt(optfile,geomshadow);
#endif

    if (done)
    {
      //printf(" attempting mopac \n");
      optfile = "scratch/xyzfile.xyzm_"+nstr;
      mopt.reset(geom1.natoms0,geom1.anumbers,geom1.anames,geom1.coords);
      mopt.freeze(frzlist, nfrz, nfrz0);
      printf(" file name: %s nscratch: %i tniso+1: %i \n",optfile.c_str(),nscratch,tniso+1);
      mopt.opt_write(optfile);
      tniso++;
      for (int l=0;l<3*geom1.natoms;l++) isos[tniso][l] = mopt.xyz[l];
    }
    else { nfail++; nscratch--; }
  }





  // parallel run mopac jobs

  ofstream mparafile;
  string mparafile_string = "scratch/go_mopac";

  int npara = tniso / MAXPARA + 1;
  int nmdone = 0;
#if !SKIPMOPAC
  printf("\n\n Starting MOPAC parallel runs \n");
  printf(" Loops: %i tniso: %i \n",npara,tniso);
  fflush(stdout);
  for (int i=0;i<npara;i++)
  {
    printf(" in MOPAC loop %i \n",i);
    string nstr1=StringTools::int2str(i,4,"0");
    string mparafile_string = "scratch/go_mopac"+nstr1;
    mparafile.open(mparafile_string.c_str());
    mparafile << "#!/bin/bash" << endl << endl;
    for (int j=0;j<MAXPARA;j++)
    {
      string nstr2=StringTools::int2str(nmdone++ +1,4,"0");
      if (nmdone<=tniso)
        mparafile << " /tmp/MOPAC2012.exe scratch/xyzfile.xyzm_" << nstr2 << " & " << endl;
    }
    mparafile << " wait " << endl;
    mparafile.close();
    string cmd = "chmod ug+x "+mparafile_string;
    system(cmd.c_str());
    system(mparafile_string.c_str());
  }
#endif

  printf(" done running MOPAC, reading output \n"); 
  double* txyz = new double[N3];
  int ttniso = 0;
  for (int i=1;i<=tniso;i++)
  {
    string nstr=StringTools::int2str(i,4,"0");
    string oname = "scratch/xyzfile.xyzm_" + nstr;
    printf("\n reading in: %s \n",oname.c_str());

    energies[ttniso+1] = mopt.read_output(oname);
#if NOMOOPT
    mopt.reset(natoms0,geoms[id].anumbers,geoms[id].anames,isos[i]);
#else
    mopt.xyz_read(oname);
#endif
    if (energies[ttniso+1]==-1)
    {
      printf(" Mopac failed \n");
      energies[ttniso+1] = 10000;
    }

    //check duplicate structure
    if (diff_structurec(geom1.natoms,geom1.anames,geom1.anumbers,isos[i],mopt.xyz)>MMVSMOPAC)
    {
      printf(" warning: significantly different structure deleted \n");
      nmopaclargediff++;
    }
    else if (diff_structure(geom1.natoms,geom1.anames,geom1.anumbers,mopt.xyz,geoms[id].coords)==0
        && close_value(energies[id],energies[ttniso+1],EVSORIGTOL))
    {
      printf(" warning: mopac returned to original structure, deleting \n");
      nmopaclargediff++;
    }
    else
    {
      ttniso++;
      for (int l=0;l<N3;l++) isos[ttniso][l] = mopt.xyz[l];
    }
  }
  tniso = ttniso;

  delete txyz;



//////////////////////////////////////
// post processing after generation //
//////////////////////////////////////

  printf("\n\n POST-PROCESSING AFTER ISO \n");
/*
  printf(" printing iso energies and geometries: \n");
  for (int i=0;i<=tniso;i++)
    printf(" %i: %1.4f \n",i,energies[i]);
  printf("\n");
*/

  if (nsave>tniso) nsave=tniso;


#if 0
  printf(" counting H2 and correcting energy \n");
  for (int i=0;i<niso+tniso;i++)
  {
    geom1.reset(natoms0,geoms[id].anames,geoms[id].anumbers,isos[i]);
    geom1.ic_create();
    int numh2 = geom1.h2count();
    energies[i] += numh2*H2PENALTY;
    printf(" xyzfile %i: counted %i H2 formed \n",i,numh2);
  }
#endif
  printf("\n now selecting %i lowest energy isomers among tniso: %i \n",nsave,tniso);
 
  double* telist = new double[tniso+1];
  int* order = new int[tniso+1];
  for (int i=0;i<=tniso;i++)
  {
    telist[i] = 50000;
    order[i] = i;
  }
  telist[0] = energies[0];
  order[0] = id;

  for (int i=1;i<=tniso;i++)
  {
//    printf(" order: ");
//    for (int j=0;j<=tniso;j++)
//      printf(" %i %1.3f ",order[j],telist[j]);
//    printf("\n");
    for (int j=1;j<=i;j++)
    {
      if (energies[i]<telist[j] && i!=j)
      {
        for (int k=tniso-1;k>=j;k--) 
        {
          order[k+1] = order[k]; 
          telist[k+1] = telist[k];
        }
        order[j] = i;
        telist[j] = energies[i];
        break;
      }
      else if (j==i)
      { 
        telist[j] = energies[i];
        order[j] = i;
      }
    } // loop j
  }
  printf("\n");

  for (int i=0;i<=tniso;i++)
    printf(" ordered %i: %1.4f was # %i \n",i,telist[i],order[i]);

  printf("\n creating geoms from optimized structures \n");
 
  int nmofail = 0;
  for (int i=0;i<=tniso;i++)
    if (telist[i]>9999)
      nmofail++;
  printf(" found %i MOPAC failures \n",nmofail);
  nsave -= nmofail;

  printf(" niso: %i \n",niso);
  printf(" nmopacdiff: %i nmopaclargediff: %i \n",nmopacdiff,nmopaclargediff);

// only saving lowest energy structures
  for (int i=0;i<nsave;i++)
  {
    printf(" origin: %i saving structure %i as geoms[%i] with E: %1.1f \n",id,order[i+1],i+niso,telist[i+1]);
//    geoms[i+niso].init(natoms,geoms[id].anames,geoms[id].anumbers,isos[order[i+1]]);
    geoms[i+niso].reset(natoms0,geoms[id].anames,geoms[id].anumbers,isos[order[i+1]]);
    geoms[i+niso].id = i+niso;
    geoms[i+niso].pid = id;
    geoms[i+niso].seenergy = telist[i+1];
    elist[i+niso] = telist[i+1];
//    geoms[1+i+niso].print_xyz();
//    geoms[1+i+niso].print_ic();
  }

  printf(" putting isos in same order as geoms \n");
  for (int i=1;i<=nsave;i++)
    for (int j=0;j<3*natoms0;j++)
      isos[i][j] = geoms[i].coords[j];

#define REMOVEDUPLICATES 0
#if REMOVEDUPLICATES
  int nunique = remove_duplicates(nsave);
  printf("\n found %i unique structures \n",nunique);
  for (int i=0;i<nsave;i++)
  {
    if (repeat[order[i+1]]==-1)
    {
      printf(" origin: %i saving structure %i as geoms[%i] with E: %1.3f \n",id,order[i+1],niso,telist[i+1]);
      geoms[niso].reset(isos[i+1]);
      geoms[niso].ic_create();
      //geoms[niso].mm_init();
      geoms[niso].id = niso;
      geoms[niso].pid = id;
      geoms[niso].seenergy = telist[i+1];
      elist[niso] = telist[i+1];
      niso++;
    }
  }
  nsave = nunique;
#else
  niso += nsave;
#endif
  printf("\n\n saved %i structures, %i failed \n",nsave,nfail);
  printf(" total structures: %i \n",niso);

  printf("\n\n");

//  niso += nunique;

//  basis1.freemem();

#if 0
  for (int i=0;i<natoms;i++)
    delete [] cmovesp[i];
  delete [] cmovesp;
  for (int i=0;i<natoms;i++)
    delete [] cmovesm[i];
  delete [] cmovesm;
  delete [] ncmovesp;
  delete [] ncmovesm;

  for (int i=0;i<natoms;i++)
    delete [] cmoves2p[i];
  delete [] cmoves2p;
  for (int i=0;i<natoms;i++)
    delete [] cmoves2m[i];
  delete [] cmoves2m;
  delete [] ncmoves2p;
  delete [] ncmoves2m;

  for (int i=0;i<NMAXSTRUCTS;i++)
    delete [] addsubmm[i];
  delete [] addsubmm;
  for (int i=0;i<NMAXSTRUCTS;i++)
    delete [] addsubpp[i];
  delete [] addsubpp;
  for (int i=0;i<NMAXSTRUCTS;i++)
    delete [] addsubmp[i];
  delete [] addsubmp;
  for (int i=0;i<NMAXSTRUCTS;i++)
    delete [] addsubmp2[i];
  delete [] addsubmp2;

  for (int i=0;i<natoms*8+1000;i++)
    delete [] movesp[i];
  delete [] movesp;
  delete [] movesm;
  for (int i=0;i<nasmp*nasmp*10+1000;i++)
    delete [] moves2pm[i];
  delete [] moves2pm;
  for (int i=0;i<nasmp*natoms*4+1000;i++)
    delete [] moves2m1p[i];
  delete [] moves2m1p;
  for (int i=0;i<natoms*100+1000;i++)
    delete [] movesnmp2[i];
  delete [] movesnmp2;
  for (int i=0;i<nasmp*natoms*20+1000;i++)
    delete [] moves2p1m[i];
  delete [] moves2p1m;
  for (int i=0;i<natoms*100+1000;i++)
    delete [] movesnmp1[i];
  delete [] movesnmp1;
#endif

#if 0
  for (int i=0;i<NMAXSTRUCTS;i++)
    delete [] isos[i];
  delete [] isos;
#endif

//  delete [] elist;
//  delete [] olist; 
  delete [] order;
  delete [] telist;

#if 0
//problem deallocating
  delete [] frzlist;
#endif

  delete [] energies;
  //delete [] energiesr;

  return nsave;

}



