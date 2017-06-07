#include "zstruct.h"
#include "utils.h"
#include "dft.h"
#include "pgsm.h"
#include "constants.h"
#include "print.h"
using namespace std;

#define NUM_INTS 2500

void ZStruct::go_mech_it(int nsteps) {

  if (nsteps>0)
  {
    geomsts = new ICoord[NMAXSTRUCTS];
    dfterelts = new double[NMAXSTRUCTS];
    for (int i=1;i<NMAXSTRUCTS;i++)
      geomsts[i].alloc(geom0.natoms);
    for (int i=1;i<NMAXSTRUCTS;i++)
      geomsts[i].reset(geom0.natoms,geoms[0].anames,geoms[0].anumbers,geoms[0].coords);
  }

  int* kept_ints = new int[NUM_INTS];
  int* nkept = new int[NUM_INTS];
  cMap = new int[NUM_INTS*NUM_INTS];
  dfterel = new double[NUM_INTS];
  savelist = new int[NUM_INTS]; 
  savebranch = new int[NUM_INTS]; 
  for (int i=0;i<NUM_INTS*NUM_INTS;i++) cMap[i] = -1;
  for (int i=0;i<NUM_INTS;i++) cMap[i*NUM_INTS+i] = 2; //diagonal means same struct
  for (int i=0;i<NUM_INTS;i++) kept_ints[i] = -1;
  for (int i=0;i<NUM_INTS;i++) nkept[i] = -1;
  for (int i=0;i<NUM_INTS;i++) dfterel[i] = 0.;
  for (int i=0;i<NUM_INTS;i++) savelist[i] = -1;
  for (int i=0;i<NUM_INTS;i++) savebranch[i] = -1;
  savelist[0] = 1;
  savebranch[0] = 1;

  if (nsteps==1) nsteps = 1;
  if (nsteps==-1) nsteps = 1;

  if (nsteps>1)
  {
    printf("\n\n*****************************\n");
    printf("  Getting a mechanism! \n");
    printf("*****************************\n");

    printf(" oemax is: %1.4f \n",oemax);
    printf(" oemin is: %1.4f \n",oemin);
    printf(" Steps remaining: %i \n",nsteps);
  }
  else
    nkept[0] = do_one_step(kept_ints);
    
  string bDir, cmd, nDirNum0, nDirNum1, bDirNum;
  int cont = 1; if(nsteps==1) cont = 0;
  int nstepsdone = 0;
  int success = 0;
  nbranches = 1;
  cbranch = 0;
  string emax, emin;

  if (nbranches==1 && nsteps>1)
  {
    int N3 = geoms[0].natoms*3;
    double** isos = new double*[NMAXSTRUCTS];
    for (int i=0;i<NMAXSTRUCTS;i++)
      isos[i] = new double[N3];
    alloc_elists();

   //do the first step
    printf("\n creating first step \n");
    bDirNum =  StringTools::int2str(0,4,"0");
    emax = StringTools::double2str1(oemax);
    emin = StringTools::double2str1(oemin);
#if !SKIPSUBMIT
#if !SKIPFIRST
    cmd = "./dircreate 0000 "+emax+" "+emin+" 1; cp test.xyz 0000/";
  //  printf("cmd: %s \n",cmd.c_str());
    system(cmd.c_str());
    cmd = "cd 0000/; qsub zstruct.qsh; cd ..";
    system(cmd.c_str());
#endif
#endif
    printf("  submitted %s \n",cmd.c_str());
  }

  while(cont)
  {
    printf("\n  Start of continue loop! \n\n");
    wait_for_branches();
    if (nbranches==1) get_branch_ints(0,0,nkept,kept_ints);
 
    printf(" #ints in each branch: \n");
    for (int i=0;i<nbranches;i++)
      printf("  nkept[%i]: %i \n",i,nkept[i]);
    printf("\n");

    int newbranches = 0;
    int nextInt = 1; 
    for (int i=cbranch;i<nbranches;i++)
    {
      printf("\n working on branch: %i \n",i);
     // read final intermediates from branch directories
      nextInt = 1;
      for (int k=0;k<nbranches;k++)
      if (nkept[k]>0)
        nextInt += nkept[k];
      printf(" nextInt: %i \n",nextInt);

      printf(" before get_branch_ints, savelist: %i \n",savelist[i]);
      if (savelist[i])
        get_branch_ints(i,nextInt,nkept,kept_ints);
      else nkept[i] = 0;

      int testInt = 1;
      for (int k=0;k<nbranches;k++)
      if (nkept[k]>0)
        testInt += nkept[k];
      printf(" testInt: %i \n",testInt);
      save_unique_multi(nextInt,testInt);

// now screen based on TS energies
// NOT IMPLEMENTED


      for (int j=0;j<nkept[i];j++)
      {
        if (savelist[nextInt+j])
          printf("\n  branch[%i] continues on %i with E: %1.1f E(GSM): %1.1f \n",i,kept_ints[j],dfterel[nextInt+j],(dftelistgsm[nextInt+j]-dftelist[0])*627.5); 
        else
          printf("\n  branch[%i] not continuing on %i with E: %1.1f E(GSM): %1.1f \n",i,kept_ints[j],dfterel[nextInt+j],(dftelistgsm[nextInt+j]-dftelist[0])*627.5);

       //Creates new branch folder, submits
       // also copies stringfile to current directory
        bDirNum =  StringTools::int2str(i,4,"0");
        bDir = bDirNum+"/";
        nDirNum0 =  StringTools::int2str(kept_ints[j],4,"0");
        nDirNum1 =  StringTools::int2str(nbranches+newbranches,4,"0");

        cmd = "cp "+bDir+"stringfile.xyz"+nDirNum0+" stringfile.xyz"+bDirNum+"to"+nDirNum1;
        printf("cmd: %s \n",cmd.c_str());
        system(cmd.c_str());
        cmd = "cp "+bDir+"scratch/paragsm"+nDirNum0+" scratch/paragsm"+bDirNum+"to"+nDirNum1;
        printf("cmd: %s \n",cmd.c_str());
        system(cmd.c_str());

        geoms[nextInt+j].print_xyz_save("next.xyz");
        cmd = "cd "+bDir+"; mv ../next.xyz .; cd ..";
      //  printf("cmd: %s \n",cmd.c_str());
        system(cmd.c_str());
        emax = StringTools::double2str1(oemax-dfterel[nextInt+j]);
        emin = StringTools::double2str1(oemin-dfterel[nextInt+j]);
#if !SKIPSUBMIT
        cmd = "cd "+bDir+"; ./dircreate "+nDirNum0+" "+emax+" "+emin+" 1; cd ..";
      //  cmd = "cd "+bDir+"; ./dircreate "+nDirNum0+" "+"10000"+" "+"1; cd ..";
      //  printf("cmd: %s \n",cmd.c_str());
        system(cmd.c_str());
        cmd = "mv "+bDir+nDirNum0+" "+nDirNum1;
        //printf("cmd: %s \n",cmd.c_str());
        system(cmd.c_str());
#endif
        printf(" %s%s is now %s \n",bDir.c_str(),nDirNum0.c_str(),nDirNum1.c_str());
        if (nsteps>1 && savelist[nextInt+j])
        {
          cmd = "cd "+nDirNum1+"; qsub zstruct.qsh; cd ..";
#if !SKIPSUBMIT
          system(cmd.c_str());
#endif
          printf("  submitted %s \n",cmd.c_str());
          fflush(stdout);
        }
        cMap[i*NUM_INTS+nbranches+newbranches] = 1;

        newbranches++;
      } // nkept submit loop
    } // current branch loop
    cbranch = nbranches;
    nbranches += newbranches;

    nstepsdone++;
    if (--nsteps<1) cont = 0;

  } // while continuing


  if (nstepsdone>0)
  {
    printf("\n saving all structures (%i) \n",nbranches);
    for (int i=0;i<nbranches;i++)
    {
      string nstr=StringTools::int2str(i,4,"0");
      string finfile = "final.xyz_"+nstr;
      geoms[i].print_xyz_save(finfile,(dftelist[i]-dftelist[0])*627.5);
    }
    system("cat final.xyz_* > final.xyz");
    for (int i=0;i<nbranches;i++)
    if (savelist[i])
    {
      string nstr=StringTools::int2str(i,4,"0");
      string finfile = "final.xyzu_"+nstr;
      geoms[i].print_xyz_save(finfile,(dftelist[i]-dftelist[0])*627.5);
    }
    system("cat final.xyzu_* > final.xyzu");
  }

  print_cmap();
  get_low_paths();

  system("rm tmp.xyz");
  system("rm final.xyzu_*");
  system("rm final.xyz_*");
  if (nstepsdone==0)
  system("rm final.xyzm_*");
  if (nstepsdone==0)
  system("rm molecule*");
 
  //Alert manager 
  if (nstepsdone==0) system("touch zdone");
  printf("\n Done moving forward! \n");

  delete [] savelist;
  delete [] nkept;
  delete [] kept_ints;
  free_elists();

  return;
}



void ZStruct::get_low_paths() {
  
  printf("\n getting low energy pathways \n");

  int* paths = new int[nbranches*nbranches];
  int* paths0 = new int[nbranches*nbranches]; //without reverting to duplicates
  for (int i=0;i<nbranches*nbranches;i++) paths[i] = -1;
  for (int i=0;i<nbranches*nbranches;i++) paths0[i] = -1;

  double* epath = new double[nbranches*nbranches];
  double* epathts = new double[nbranches*nbranches];
  for (int i=0;i<nbranches*nbranches;i++) epath[i] = -1.;
  for (int i=0;i<nbranches*nbranches;i++) epathts[i] = -1.;

  int npaths = 0;
  int s1 = 0;
  int s2 = 0;

  int* tpath = new int[nbranches];
  for (int i=0;i<1;i++)
  {
    for (int j=0;j<nbranches;j++) tpath[j] = -1;
    tpath[0] = i;
    int ncurr = 1;
    int done = 0;
    int r1 = tpath[ncurr-1];
    int r2 = tpath[ncurr-1];

    int bu = 0;

    if (cMap[i])
    while(!done)
    {
   //   printf(" path so far (ncurr: %i): \n",ncurr);
   //   for (int j=0;j<ncurr;j++)
   //     printf(" %i",tpath[j]);
   //   printf("\n");

   //   printf(" r2, r1: %i %i \n",r2,r1);
      int c=0;

   //   printf(" for this r1, row is: \n");
   //   for (int j=r2+1;j<nbranches;j++)
   //     printf(" %i",cMap[r1*NUM_INTS+j]);
   //   printf("\n");
      for (int j=r2+1;j<nbranches;j++)
      if (cMap[r1*NUM_INTS+j]==1)
      {
        tpath[ncurr++] = j;
        c++;
        bu = 1;
        break;
      }
     //check for duplicate
      if (c==0)
      for (int j=0;j<r2;j++)
      if (cMap[r1*NUM_INTS+j]==3)
      {
        tpath[ncurr++] = -j;
        bu = -1;
      //  printf(" this path reverts to %i \n",j);
        break;
      }
      if (c>0)
      {
        r1 = tpath[ncurr-1];
        r2 = tpath[ncurr-1];
      }
      else if (c==0)
      {
        if (bu)
        { 
          printf(" found path %i:",npaths);
          int c = 0;
          for (int j=0;j<ncurr;j++)
          {
            if (bu==-1 && j==ncurr-1)
              printf("-%i",abs(tpath[j]));
            else
              printf(" %i",tpath[j]);
            if ( bu>0 || (j!=ncurr-1) )
            {
              paths[npaths*nbranches+c] = tpath[j];
              paths0[npaths*nbranches+c] = tpath[j];
              epathts[npaths*nbranches+c] = (dftelistgsm[tpath[j]]-dftelist[0])*627.5;
              epath[npaths*nbranches+c] = dfterel[tpath[j]];
              c++;
            }
#if 1
            else if (bu<0 && j==ncurr-1)
            { 
             //assign previous ID instead of duplicate
              paths[npaths*nbranches+c-1] = abs(tpath[j]);
              epath[npaths*nbranches+c-1] = dfterel[abs(tpath[j])];
            }
#endif
          } //loop j over ncurr in path
          npaths++;
          printf("\n");
        }
       //go back 1 step
        //r2++;
        bu = 0;
        ncurr--;
        r1 = tpath[ncurr-1];
        r2 = tpath[ncurr];
      } //if c==0
      if (ncurr<1)
        done = 1;
    } //while !done
  }
  delete [] tpath;

  int* nmaxts = new int[npaths];
  double* emaxts = new double[npaths];
  for (int i=0;i<npaths;i++) nmaxts[i] = 0;
  for (int i=0;i<npaths;i++) emaxts[i] = -0.1;
  for (int i=0;i<npaths;i++)
  {
    for (int j=0;j<nbranches;j++)
    if (emaxts[i]<epathts[i*nbranches+j])
    {
      nmaxts[i] = j;
      emaxts[i] = epathts[i*nbranches+j];
    }
  }
  int nminall = 0;
  double eminall = emaxts[0];
  for (int i=1;i<npaths;i++)
  {
    if (emaxts[i]<eminall)
    {
      nminall = i;
      eminall = emaxts[i];
    }
  }


  double* telist = new double[npaths];
  int* porder = new int[npaths];
  for (int i=0;i<npaths;i++)
  {
    telist[i] = 50000;
    porder[i] = i;
  }

  for (int i=0;i<npaths;i++)
  {
    for (int j=0;j<=i;j++)
    {
      if (emaxts[i]<telist[j] && i!=j)
      {
        for (int k=npaths-2;k>=j;k--)
        {
          porder[k+1] = porder[k];
          telist[k+1] = telist[k];
        }
        porder[j] = i;
        telist[j] = emaxts[i];
        break;
      }
      else if (j==i)
      {
       	telist[j] = emaxts[i];
        porder[j] = i;
      }
    } // loop j
  }



  printf("\n **************************************** \n");
  printf(" ******* Printing Final Pathways! ******* \n");
  printf(" **************************************** \n");

  //removing duplicate paths
  int* uniquepaths = new int[npaths];
  for (int i=0;i<npaths;i++) uniquepaths[i] = 1;
  for (int i=0;i<npaths;i++)
  for (int j=0;j<i;j++)
  {
    int same = 1;
    for (int k=0;k<nbranches;k++)
    {
      if (paths[i*nbranches+k]!=paths[j*nbranches+k])
      {
        same = 0;
        break;
      }
     // if (paths[i*nbranches+k]==-1) break;
    }
    if (same)
    {
    //  printf(" paths %i %i are duplicate \n",i,j);
      uniquepaths[i] = 0;
    }
  }

  for (int i=0;i<npaths;i++)
  if (uniquepaths[i])
  {
    string nstr=StringTools::int2str(i,4,"0");
    string cmd;
    string cmd2;
    //system(cmd.c_str());
    cmd = "rm -f path.xyz"+nstr;
    cmd2 = "rm -f path.xyzs"+nstr;
    system(cmd.c_str());
    system(cmd2.c_str());
    cmd = "cat tmp.xyz >> path.xyz"+nstr;
    int t0 = 0;
    int first = 1;
    for (int j=0;j<nbranches;j++)
    if (paths[i*nbranches+j]>-1)
    {
      int t1 = paths[i*nbranches+j];
      int t2 = paths0[i*nbranches+j];
//      double e1 = dfterelts[s1];
//      double e2 = dfterel[s1];
      double e1 = epathts[i*nbranches+j];
      double e2 = epath[i*nbranches+j];
      if (!first)
      {
        geomsts[t1].print_xyz_save("tmp.xyz",e1);
        system(cmd.c_str());
        string nstr1=StringTools::int2str(t0,4,"0");
        string nstr2=StringTools::int2str(t2,4,"0");
        cmd2 = "cat stringfile.xyz"+nstr1+"to"+nstr2+" >> path.xyzs"+nstr;
        system(cmd2.c_str());
      }
      geoms[t1].print_xyz_save("tmp.xyz",e2);
      system(cmd.c_str());
      t0 = t1;
      first = 0;
    }

    printf("\n path %3i:",i);
    for (int j=0;j<nbranches;j++)
    if (paths[i*nbranches+j]>-1)
    {
      if (j!=0)
        printf("       TS %8i",paths[i*nbranches+j]);
      else
        printf("    %i",paths[i*nbranches+j]);
    }
    printf("\n");
//    printf(" epath: 0.");
    printf(" epath:     ");
    for (int j=0;j<nbranches;j++)
    if (paths[i*nbranches+j]>-1)
    {
      if (j!=0)
        printf(" to  %4.1f  to %4.1f",epathts[i*nbranches+j],epath[i*nbranches+j]);
      else
        printf(" %1.1f",epath[i*nbranches+j]);
    }
    printf("\n");
  }

  printf("\n   Lowest barrier path: %i at %1.1f kcal/mol \n",nminall,eminall);
  printf("\n   Ordered paths (unique): \n");
  for (int i=0;i<npaths;i++)
  if (uniquepaths[porder[i]])
    printf(" path %i with Ea(max): %3.1f \n",porder[i],telist[i]);

#if 0
  printf("\n dumping paths matrix \n");
  for (int i=0;i<npaths;i++)
  {
    for (int j=0;j<nbranches;j++)
       printf(" %i",paths[i*nbranches+j]);
    printf("\n");
  }
#endif

#if 1
  printf("\n **************************************** \n");
  printf(" ******* Printing Lowest Routes ******* \n");
  printf(" **************************************** \n");

  int* endint = new int[npaths];
  for (int i=0;i<npaths;i++) endint[i] = -1;
  for (int i=0;i<npaths;i++)
  {
    int found = 0;
    for (int j=0;j<nbranches;j++)
    if (paths[i*nbranches+j]==-1)
    {
     // printf(" path %i ends with %i \n",i,paths[i*nbranches+j-1]);
      endint[i] = paths[i*nbranches+j-1];
      break;
    }
  } //loop i over nbranches

  int nunique = 0;
  int* uniquelist = new int[npaths];
  for (int i=0;i<npaths;i++) uniquelist[i] = -2;
  for (int i=0;i<npaths;i++)
  { 
    int found = 0;
    for (int j=0;j<nunique;j++)
    if (endint[i]==uniquelist[j])
    {
      found = 1;
      break;
    }
    if (!found)
      uniquelist[nunique++] = endint[i];
  }

 // for (int i=0;i<nunique;i++)
 //   printf(" unique ending intermediate: %i \n",uniquelist[i]);

  for (int i=0;i<nbranches;i++)
  {
    int found = 0;
    for (int j=0;j<nunique;j++)
    if (uniquelist[j]==i)
      found = 1;

    if (found)
      printf("\n Paths to intermediate %2i (dH: %2.1f): \n",i,dfterel[i]);
    if (found)
    for (int j=0;j<npaths;j++)
    if (endint[j]==i && uniquepaths[j])
    {
      printf("  #%2i (Ea: %3.1f):",j,emaxts[j]);
      for (int k=0;k<nbranches;k++)
      {
        if (paths[j*nbranches+k]>-1)
          printf(" %2i",paths[j*nbranches+k]);
        else
          break;
      }
//      printf(" Ea: %1.1f \n",emaxts[j]);
      printf("\n");
    }
  }

#endif

  delete [] uniquelist;
  delete [] uniquepaths;
  delete [] endint;
 
  delete [] telist;
  delete [] porder;

  delete [] nmaxts;
  delete [] emaxts;

  delete [] paths0;
  delete [] paths;
  delete [] epath;
  delete [] epathts;

  return;
}


void ZStruct::print_cmap() {

  printf("\n printing cMap \n");
  for (int i=0;i<nbranches;i++)
  {
    printf(" %2i: ",i);
    for (int j=0;j<nbranches;j++)
    {
//      printf(" %2i",cMap[i*NUM_INTS+j]);
      if(cMap[i*NUM_INTS+j]==-1)
        printf(" .");
      else if(cMap[i*NUM_INTS+j]==2)
        printf(" x");
      else if(cMap[i*NUM_INTS+j]==1)
        printf(" c");
      else if(cMap[i*NUM_INTS+j]==3)
        printf(" d");
    }
    printf("\n");
  }

  // now creating graph
  printf("\n creating graph_n and graph_e.csv \n");

  double shift = 0.;
  double emax = -100.;
  double tsemax = 0.;
  for (int i=0;i<nbranches;i++)
  if (emax<dfterel[i])
    emax = dfterel[i];
  shift = emax + 1.;

  for (int i=0;i<nbranches;i++)
  if (tsemax<(dftelistgsm[i]-dftelist[0]))
    tsemax = dftelistgsm[i]-dftelist[0];
  tsemax *= 627.5;
  tsemax += 1.;

  ofstream grpfile;
  string grpfile_string = "graph_n.csv";
  grpfile.open(grpfile_string.c_str());

  grpfile << "Id;Label;dH" << endl;
  for (int i=0;i<nbranches;i++)
    grpfile << i << ";" << i << ";" << -dfterel[i]+shift << endl;
//    grpfile << i << ";S-" << i << ";" << -dfterel[i]+shift << endl;
  grpfile.close();

  grpfile_string = "graph_e.csv";
  grpfile.open(grpfile_string.c_str());

  grpfile << "Source;Target;Id;Weight" << endl;
  for (int i=0;i<nbranches;i++)
  {
    for (int j=i;j<nbranches;j++)
    {
      if(cMap[i*NUM_INTS+j]==1)
      {
        double gsme = (dftelistgsm[j]-dftelist[0])*627.5;
        grpfile << i << ";" << j << ";" << i << ";" << 1.-gsme/tsemax << endl;
      }
    }
  }
  grpfile.close();

  grpfile_string = "graph_e2.csv";
  grpfile.open(grpfile_string.c_str());

  grpfile << "Source;Target;Id;Weight" << endl;
  for (int i=0;i<nbranches;i++)
  {
    for (int j=i;j<nbranches;j++)
    {
      if(cMap[i*NUM_INTS+j]==1)
      {
        int label = j;
        for (int k=0;k<j;k++)
        if (cMap[label*NUM_INTS+k]==3)
          label = k;
        double gsme = (dftelistgsm[j]-dftelist[0])*627.5;
        grpfile << i << ";" << label << ";" << i << ";" << 1.-gsme/tsemax << endl;
      }
    }
  }
  grpfile.close();

  return;
}


void ZStruct::get_branch_ints(int bn, int nextInt, int* nkept, int* kept_ints) {

  printf(" in get_branch_ints \n");

  string bDirNum =  StringTools::int2str(bn,4,"0");
  string bDir = bDirNum+"/";
  string contfile; 
  if (bn==0 && nextInt==0) contfile = "geom0.xyz";
  else contfile = "continue.xyz";
  string cname = bDir+contfile;
  printf(" opening: %s \n",cname.c_str());

  ifstream output(cname.c_str(),ios::in);
  if (!output) { printf(" error opening file: %s \n",cname.c_str()); return;}
  string line;
  vector<string> tok_line;
  int c = 0;
  getline(output,line);
  while(!output.eof())
  {
    cout << " RR " << line << endl;
    tok_line = StringTools::tokenize(line, " \t");
    kept_ints[c] = atoi(tok_line[0].c_str());
    dftelistgsm[nextInt+c] = atof(tok_line[1].c_str());
    dftelist[nextInt+c] = atof(tok_line[2].c_str());
    dfterel[nextInt+c] = (dftelist[nextInt+c]-dftelist[0])*627.5;
    dfterelts[nextInt+c] = (dftelistgsm[nextInt+c]-dftelist[0])*627.5;
    c++;
    getline(output,line);
  }
  nkept[bn] = c;
  if (bn==0 && nextInt==0)
  {
    nkept[bn] = -1;
    elist[0] = dftelistgsm[0] = dftelist[0];
    dfterel[0] = 0.;
    printf(" E(DFT): %1.4f \n",dftelist[0]);
  }

//read in XYZ
  int natoms = geoms[0].natoms;
  if (bn>0 || nextInt>0)
  printf("\n reading geometries from initial.xyz for ");
  if (bn>0 || nextInt>0)
  for (int i=0;i<c;i++)
  {
    savelist[nextInt+i] = 1;
    printf(" %i",kept_ints[i]);
    string nstr = StringTools::int2str(kept_ints[i],4,"0");
    string xyzname = bDir+"scratch/initial"+nstr+".xyz"; 
    ifstream output2(xyzname.c_str(),ios::in);
    if (!output2) { printf(" error opening xyz file: %s \n",xyzname.c_str());}

    for (int k=0;k<natoms+4;k++)
      getline(output2,line);
    int j = 0;
    while(!output2.eof())
    {
      getline(output2,line);
      //cout << " RR " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      geoms[nextInt+i].coords[3*j+0] = atof(tok_line[1].c_str());
      geoms[nextInt+i].coords[3*j+1] = atof(tok_line[2].c_str());
      geoms[nextInt+i].coords[3*j+2] = atof(tok_line[3].c_str());
      j++;
    }
    output2.close();

    xyzname = bDir+"scratch/tsq"+nstr+".xyz"; 
    ifstream output3(xyzname.c_str(),ios::in);
    if (!output3) { printf(" error opening ts xyz file: %s \n",xyzname.c_str());}

    getline(output3,line);
    getline(output3,line);
    j = 0;
    while(!output3.eof())
    {
      getline(output3,line);
      //cout << " RR " << line << endl;
      tok_line = StringTools::tokenize(line, " \t");
      geomsts[nextInt+i].coords[3*j+0] = atof(tok_line[1].c_str());
      geomsts[nextInt+i].coords[3*j+1] = atof(tok_line[2].c_str());
      geomsts[nextInt+i].coords[3*j+2] = atof(tok_line[3].c_str());
      j++;
    }
    output3.close();
  }
  printf("\n");

#if 0
  for (int i=0;i<c;i++)
  {
    printf(" printing structure for: %i \n",kept_ints[i]);
    geoms[nextInt+i].print_xyz();
//    for (int j=0;j<natoms;j++)
//      printf(" %s %1.4f %1.4f %1.4f \n",geoms[nextInt+i].anames[j].c_str(),geoms[nextInt+i].coords[3*j+0],geoms[nextInt+i].coords[3*j+1],geoms[nextInt+i].coords[3*j+2]);
  }
#endif

  return;
}






void ZStruct::wait_for_branches() {

  printf(" in wait_for_branches (cbranch: %i nbranches: %i) \n",cbranch,nbranches);

//  for (int i=0;i<nbranches;i++)
//    nkept[i] = 2;
  int qnotdone = true;
  int* branchdone = new int[nbranches];
  for (int i=0;i<nbranches;i++) branchdone[i]=0;
  for (int i=0;i<nbranches;i++) 
  if (!savelist[i])
    branchdone[i]=1;

  int max_wait = MAX_BRANCH_WAIT*10;
  int tc = 0;
  do {
    tc++; if (tc>max_wait) { printf(" enough waiting! \n"); break; }
#if !SKIPSUBMIT
    sleep(15);
#endif
    qnotdone = false;
    for (int i=cbranch;i<nbranches;i++)
    {
      if (tc%400==0) { printf(" checking branches "); fflush(stdout); }
      if (!branchdone[i])
      {
        string nstr=StringTools::int2str(i,4,"0");
        string bstring = nstr+"/zdone";
        struct stat sts;
        if (tc%1600==0) { printf(" checking %s \n",bstring.c_str()); fflush(stdout); }
        if (stat(bstring.c_str(), &sts) != -1)
        {
          branchdone[i]=1;
          printf("\n done with Branch %i",i);
        }
        else
          qnotdone = true;
      } // if not done
    } //loop i over nts
  } while (qnotdone);
  printf("\n");

  return;
}




int ZStruct::do_one_step(int* final_ints) {
  
  //printf("\n one elementary step \n");

  int nbranch = NMAXSTRUCTS;
  int nunique;
  niso = 1;
  int nkept = 0;

  int* olist = new int[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) olist[j]=-1;
  alloc_elists();

  read_frozen(geom0.natoms);

  int ngen;
  ngen = generate_isomers(0,nbranch); 
  nunique = ngen;

  fflush(stdout);


  printf("\n saving all low-energy structures \n");
  system("rm final.xyzm*");
  for (int i=0;i<niso;i++)
  {
    string nstr=StringTools::int2str(i,4,"0"); 
    string finfile = "final.xyzm_"+nstr;
    geoms[i].print_xyz_save(finfile);
  }
  system("cat final.xyzm_* > final.xyzm");
  system("rm final.xyzm_*");

  double emin = elist[1];
  double max_check = emin + PM6_EMAX;
  if (max_check < elist[0]) max_check = elist[0] + 1.;

  nts = niso;
  for (int i=1;i<niso;i++)
    if (elist[i]>max_check)
    {
      nts = i;
      break;
    }
  printf(" found nts: %i from emax: %1.1f \n",nts,max_check);


  int nkeep;
  savelist = new int[niso];
  //criteria for saving
  for (int i=0;i<niso;i++) savelist[i] = 0;
  for (int i=0;i<nts;i++) savelist[i] = 1;

  fflush(stdout);
#if 1
//  dft_sp_para();
  dft_opt_para();
#else
  printf(" skipping dft opt \n");
#endif
  fflush(stdout);
  system("rm scratch/xyzfile*");


  printf(" about to order DFT structs, nts: %i \n",nts);
  int tniso = nts - 1;
  double* telist = new double[tniso+2];
  int* order = new int[tniso+2];
  for (int i=0;i<=tniso;i++)
  {
    telist[i] = 50000;
    order[i] = i;
  }
  telist[0] = dftelist[0];
  order[0] = 0;

  for (int i=1;i<=tniso;i++)
  {
    for (int j=1;j<=i;j++)
    {
      if (dftelist[i]<telist[j] && i!=j)
      {
        for (int k=tniso-1;k>=j;k--)
        {
          order[k+1] = order[k];
          telist[k+1] = telist[k];
        }
        order[j] = i;
        telist[j] = dftelist[i];
        break;
      }
      else if (j==i)
      {
       	telist[j] = dftelist[i];
        order[j] = i;
      }
    } // loop j
  }
  int natoms = geoms[0].natoms;
  double** isos = new double*[tniso+1];
  for (int i=0;i<tniso+1;i++)
    isos[i] = new double[3*natoms];
  for (int i=1;i<=tniso;i++)
    for (int j=0;j<3*natoms;j++)
      isos[i][j] = geoms[i].coords[j];

  for (int i=0;i<tniso;i++)
  {
    printf(" ordering dft structure %i as geoms[%i] with E: %1.1f \n",order[i+1],i+1,telist[i+1]);
    geoms[i+1].reset(natoms,geoms[0].anames,geoms[0].anumbers,isos[order[i+1]]);
    geoms[i+1].id = i;
    geoms[i+1].pid = 0;
    geoms[i+1].seenergy = telist[i+1];
    dftelist[i+1] = telist[i+1];
  }
  for (int i=0;i<tniso+1;i++)
    delete [] isos[i];
  delete [] isos; 
  delete [] telist;
  delete [] order;


  printf("\n saving all DFT opt'd structures \n");
  for (int i=0;i<nts;i++)
  {
    string nstr=StringTools::int2str(i,4,"0"); 
    string finfile = "final.xyz_"+nstr;
    if (i!=0)
      geoms[i].print_xyz_save(finfile,dftelist[i]);
    else
      geoms[i].print_xyz_save(finfile,0.0);
  }
  system("cat final.xyz_* > final.xyz");
//  system("rm final.xyz_*");


  emin = 25000;
  for (int i=1;i<nts;i++)
    if (dftelist[i]<emin) emin = dftelist[i];
//  max_check = emin + DFT_EMAX;
//  if (max_check<0) max_check=DFT_EMAX; // avoid exothermicity problems
  max_check = DFT_EMAX; //new criterion

  double min_check = -1000;
  if (oemax>0. || oemin<0.)
  {
    if (max_check>oemax)
      max_check = oemax;
    min_check = oemin;
  }

  printf("\n min/max E for saving: %1.1f/%1.1f \n",min_check,max_check);

  for (int i=1;i<niso;i++) savelist[i] = 0;
  for (int i=1;i<nts;i++)
  {
    if (min_check < dftelist[i] && dftelist[i]<max_check)
      savelist[i] = 1;
  }

  save_unique(1);

  printf(" saved structure %3i with E(DFT) (%1.4f) E(PM6) %1.1f \n",0,dftelist[0],elist[0]-elist[0]);
  for (int i=1;i<nts;i++) 
    if (savelist[i])
      printf(" saved structure %3i with E(DFT) %1.2f E(PM6) %1.1f \n",i,dftelist[i],elist[i]-elist[0]);

  printf("\n saving all DFT opt'd structures \n");
  for (int i=0;i<nts;i++)
  {
    if (savelist[i])
    {
      string nstr=StringTools::int2str(i,4,"0");
      string finfile = "final.xyzu_"+nstr;
      if (i!=0)
        geoms[i].print_xyz_save(finfile,dftelist[i]);
      else
        geoms[i].print_xyz_save(finfile,0.0);
    }
  }
  system("cat final.xyzu_* > final.xyzu");
  system("rm final.xyzu_*");

  tniso = 0;
  for (int i=0;i<nts;i++)
  if (savelist[i]) tniso=i;
  nts = tniso+1;

// for LST or GSM runs
  write_initialxyz(0,nts);


#if 0
  printf("\n Stopping after DFT opt! \n");
  return 0;
#endif

  fflush(stdout);
//  nkeep = generate_paths_lst(0,nts);
#if 0
  nkeep = generate_paths_lst_para(0,nts);  
#else
  printf(" skipping generate_paths_lst_para \n");
#endif
  fflush(stdout);


  printf(" \n recommending structures for DFT GSM \n");
  double max_check_ts;
  double dfteminlst = 25000;
  for (int i=1;i<nts;i++)
    if (dftelistlst[i]<dfteminlst && dftelistlst[i]>MINTSE) dfteminlst = dftelistlst[i];

  savelist[0] = 1;
#if 0
  max_check_ts = dfteminlst + LST_EMAX_TS;
  for (int i=1;i<nts;i++)
  {
//    printf(" checking: %i emin emints: %1.2f %1.2f elist elistts %1.2f %1.2f \n",i,emin,emints,dftelist[i],dftelistlst[i]);
    if (dftelistlst[i]<max_check_ts && savelist[i])
    {
      savelist[i] = 1;
     // if LST barrier is small and final structure is similar energy
      if ( close_value(dftelistlst[i],0.0,EVSORIGTOL) && close_value(dftelist[i],0.0,EVSORIGTOL) )
      {
        printf(" WARNING: probably duplicate structure %i \n",i);
        savelist[i] = 0;
      }
    }
    else
      savelist[i] = 0;
  }
  printf("\n max lst ts energy criteria: %2.1f \n",max_check_ts);
  for (int i=0;i<nts;i++) 
    if (savelist[i])
      printf(" saved structure %3i with E %2.1f E(LST/TS) %2.1f \n",i,dftelist[i],dftelistlst[i]);
#else
  printf("\n skipping lst energy check \n");
#endif





  fflush(stdout);
#if !SKIPGSM
//  nkeep = generate_paths_dft();
  nkeep = gsm_dft_para();
#else
  printf(" skipping gsm dft \n");
  printf(" kludging GSM energy \n");
  for (int i=1;i<nts;i++) dftelistgsm[i] = dftelist[i] + 10.;
#endif
  fflush(stdout);


  printf(" \n recommending structures after GSM TS search \n");
  double dftemingsm = 25000;
  for (int i=1;i<nts;i++)
    if (dftelistgsm[i]<dftemingsm && dftelistgsm[i]>MINTSE) dftemingsm = dftelistgsm[i];
//  max_check_ts = dftemingsm + GSM_EMAX_TS;
  max_check_ts = GSM_EMAX_TS + oemax; //new criterion
  if (oemax>max_check)
    max_check_ts = GSM_EMAX_TS + max_check;

#if 0
  if (oemax>0. && oemin<0.)
  {
    if (max_check>oemax)
      max_check = oemax;
  }
#endif

  for (int i=1;i<nts;i++)
  {
//    printf(" checking: %i emin emints: %1.2f %1.2f elist elistts %1.2f %1.2f \n",i,emin,emints,dftelist[i],dftelistgsm[i]);
    if (dftelistgsm[i]<max_check_ts && savelist[i])
    {
      savelist[i] = 1;
#if 0
      if (close_value(dftelistgsm[i],0.0,EVSORIGTOL) || dftelistgsm[i]<-EVSORIGTOL)
      {
        printf(" WARNING: probably duplicate structure (or GSM failed) %i \n",i);
        savelist[i] = 0;
      }
#endif
    }
    else
      savelist[i] = 0;
  }
  nkeep = 0;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      nkeep++;
  printf("\n max ts energy criteria after GSM: %1.1f \n",max_check_ts);
  printf(" saved %2i structures from TS search \n",nkeep);
  for (int i=0;i<nts;i++) 
    if (savelist[i])
      printf(" saved structure %2i with E %2.1f ETS(GSM) %2.1f \n",i,dftelist[i],dftelistgsm[i]);

  fflush(stdout);
#if 0
  nkeep = dft_ts_para();
  double dftemints = 25000;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      if (dftelistts[i]<dftemints && dftelistts[i]>MINTSE) dftemints = dftelistts[i];

//  max_check_ts = dftemints + GSM_EMAX_TS;
  max_check_ts = oemax + GSM_EMAX_TS;
  if (oemax>max_check)
    max_check_ts = GSM_EMAX_TS + max_check;
  printf(" lowest TS energy is: %1.1f max for saving: %1.1f \n",dftemints,max_check_ts);


// another round of checking
  for (int i=1;i<nts;i++)
  {
//    printf(" checking: %i emin emints: %1.2f %1.2f elist elistts %1.2f %1.2f \n",i,emin,emints,dftelist[i],dftelistts[i]);
    if (dftelistts[i]<max_check_ts && savelist[i])
      savelist[i] = 1;
    else
      savelist[i] = 0;
  }
#else
  printf("\n skipping local TS search \n");
  for (int i=1;i<nts;i++)
    dftelistts[i] = dftelistgsm[i];
#endif
  fflush(stdout);


// removes duplicates for subsequent steps
  save_unique(0);

#if 1
  printf("\n determing ints to continue \n");
  for (int i=1;i<nts;i++)
    if (savelist[i])
      final_ints[nkept++]=i;
  printf(" saved:");
  for (int i=0;i<nkept;i++)
    printf(" %i",final_ints[i]);

  printf("\n continue list generated, writing to file\n");

  ofstream nextfile;
  string nextfilestr = "continue.xyz";
  nextfile.open(nextfilestr.c_str());
  nextfile << setprecision(6);

  for (int i=0;i<nkept;i++)
  {
    int c = final_ints[i];
    string nextNum =  StringTools::int2str(c,4,"0");
    nextfile << " " << nextNum << " " << dftelistgsm[c]/627.5+dftelist[0] << " " << dftelist[c]/627.5+dftelist[0] << " " << endl;
//    nextfile << " " << nextNum << " " << elist[c]+10. << " " << elist[c] << " " << endl;
  }
  nextfile.close();

// write initial geom energy
  nextfilestr = "geom0.xyz";
  nextfile.open(nextfilestr.c_str());
  nextfile << " 0000 -1.0 " << dftelist[0] << " " << endl;
  nextfile.close();

#endif

  delete [] olist;

  printf("\n done with single elementary step search! \n");

  return nkept;
}



int ZStruct::generate_paths_lst_para(int id, int nsave) {

  printf("\n in generate_paths_lst_para \n");

  int maxprocs = MAXPARA;
  int gsmprocs = DFTPARA;
  int nkeep = nsave;

  string* rclist = new string[nts+100];
  ofstream cmdfile;
  string cmdfile_string = "scratch/go_lst_dft";

  int totalc = 0;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      totalc++;
  printf("\n found %i structures for DFT LST \n",totalc);
 
  printf("\n printing job array numbers for LST \n");
#if 0
  cmdfile.open(cmdfile_string.c_str());
  cmdfile << "#PBS -t ";
  int next = 1;
  for (int i=1;i<nts;i++)
    if (savelist[i])
    {
      cmdfile << i;
      next = i+1; break;
    }
  for (int i=next;i<nts;i++)
    if (savelist[i])
      cmdfile << "," << i;
//      printf("%i,",i);
  cmdfile.close();
#else
  int continue1 = 1;
  int start = 1;
  int next = 1;
  do 
  {
    string cmd = "mv "+cmdfile_string+" "+cmdfile_string+"_prev";
    //system(cmd.c_str());
    cmdfile.open(cmdfile_string.c_str());
    cmdfile << "#PBS -t ";
    int max = min(start + 100,nts);
    next = max;
    //printf(" start: %i max: %i next: %i nts: %i \n",start,max,next,nts);
    fflush(stdout);
    for (int i=start;i<max;i++)
      if (savelist[i])
      {
        cmdfile << i;
        next = i+1; break;
      }
    for (int i=next;i<max;i++)
      if (savelist[i])
        cmdfile << "," << i;
    cmdfile << endl;
    cmdfile.close();
    start = max;
    cmd = "./qmakel";
    system(cmd.c_str());
    cmd = "qsub "+cmdfile_string+".qsh";
#if !SKIPDFT
    system(cmd.c_str());
#endif
    if (max==nts) continue1 = 0;
  } while (continue1);
#endif

  printf("\n");

  int qnotdone = true;
  int* lstdone = new int[nts];
  for (int i=0;i<nts;i++) lstdone[i]=0;

#if !SKIPDFT
  int max_wait = MAX_TIME_WAIT;
  int tc = 0;
  do {
    tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
    sleep(10);
    qnotdone = false;
    for (int i=1;i<nts;i++)
    {
      if (savelist[i] && !lstdone[i])
      {
        string nstr=StringTools::int2str(i,1,"0");
        string lstfile_string = "scratch/lstdone"+nstr;
        struct stat sts;
        if (stat(lstfile_string.c_str(), &sts) != -1)
        {
          lstdone[i]=1;
          printf(" done with LST %i",i);
        }
        else
          qnotdone = true;
      } // if not done
    } //loop i over nts
  } while (qnotdone);
  printf("\n");
#endif


  printf(" done with LST calcs \n\n");

  fflush(stdout);

  for (int i=1;i<nts;i++)
  {
    if (savelist[i])
    {
      string nstr = StringTools::int2str(i,4,"0");
      string oname = "scratch/tsl"+nstr+".xyz";
      ifstream output(oname.c_str(),ios::in);
      string line;
      getline(output,line);
      getline(output,line);
//    cout << " RR " << line << endl;
 
      double energyts = (atof(line.c_str())-dftelist[0]) * 627.5;
      if (abs(energyts)<0.00001)
      {
        printf(" LST failed \n");
        energyts = 25000;
      }
      else
        printf(" final LST energy for %i is: %1.2f \n",i,energyts);

      dftelistlst[i] = energyts;

    } //if savelist[i]
  } //loop i over all nts

  delete [] rclist;

  return nkeep;
}

#if 0
int ZStruct::dft_sp() {

  printf("\n doing dft sp \n");

  DFT dft1;
  dft1.init(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);
  for (int i=0;i<nts;i++)
  {
    if (savelist[i])
    {
      printf(" DFT working on structure %i \n",i);
      dft1.reset(geoms[i].natoms,geoms[i].anumbers,geoms[i].anames,geoms[i].coords);
      dftelist[i] = dft1.sp();
    }
    else
    {
//      printf(" not working on structure %i \n",i);
      dftelist[i] = -1;
    }
  }
  for (int i=1;i<nts;i++)
    dftelist[i] = (dftelist[i] - dftelist[0])*627.5;
  dftelist[0] = 0.;
  printf("\n");
  for (int i=0;i<nts;i++)
    if (savelist[i])
      printf(" E[%i]: %1.3f",i,dftelist[i]);
  printf("\n");

  return 0;
}
#endif


int ZStruct::dft_sp_para() {

  int maxprocs = MAXPARA;
  int numq = DFTPARA;

  fflush(stdout);
  printf("\n doing dft sp \n");

  DFT dft1;
  dft1.init(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);
//  for (int i=0;i<maxprocs;i++)
//    dft1[i].init(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

//  string* rlist = new string[maxprocs];
  string* rclist = new string[nts]; 
  printf(" DFT preparing structure ");
  for (int i=0;i<nts;i++)
  {
    if (savelist[i])
    {
      printf(" %i",i);
      string nstr=StringTools::int2str(i,4,"0"); 
      string filename = "scratch/paradft"+nstr;
      //rlist[i%maxprocs] = filename;
      rclist[i] = filename;
      dft1.reset(geoms[i].natoms,geoms[i].anumbers,geoms[i].anames,geoms[i].coords);
      dft1.sp_dnr(filename);
    }
    else
    {
//      printf(" not working on structure %i \n",i);
      dftelist[i] = -1;
    }
  }
  printf("\n");


  ofstream cmdfile;
  string cmdfile_string = "scratch/go_dft";

  cmdfile.open(cmdfile_string.c_str());
  printf("\n printing job array numbers for DFT SP \n");
  cmdfile << "#PBS -t ";
#if 0
  int next = 0;
  for (int i=0;i<nts;i++)
    if (savelist[i])
    {
      cmdfile << i;
      next = i+1; break;
    }
  for (int i=next;i<nts;i++)
    if (savelist[i])
      cmdfile << "," << i;
//      printf("%i,",i);
#else
  cmdfile << "0-" << nts-1;
#endif
  cmdfile.close();

  printf("\n");
#if !SKIPDFT
  string cmd = "./qmaked";
  system(cmd.c_str());
  cmd = "qsub "+cmdfile_string+".qsh";
  system(cmd.c_str());
#endif

  int qnotdone = true;
  int* dftdone = new int[nts];
  for (int i=0;i<nts;i++) dftdone[i]=0;

#if !SKIPDFT
  int max_wait = MAX_TIME_WAIT;
  int tc = 0;
  do {
    tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
    sleep(10);
    qnotdone = false;
    for (int i=1;i<nts;i++)
    {
      if (savelist[i] && !dftdone[i])
      {
        string nstr=StringTools::int2str(i,1,"0");
        string dftfile_string = "scratch/dftdone"+nstr;
        struct stat sts;
        if (stat(dftfile_string.c_str(), &sts) != -1)
        {
          dftdone[i]=1;
          printf(" done with DFT SP %i",i);
        }
        else
          qnotdone = true;
      } // if not done
    } //loop i over nts
  } while (qnotdone);
  printf("\n");
#endif


  dftelist[0] = dft1.get_energy(rclist[0]);
  for (int i=1;i<nts;i++)
    if (savelist[i])
      dftelist[i] = dft1.get_energy(rclist[i]);

  for (int i=1;i<nts;i++)
    dftelist[i] = (dftelist[i] - dftelist[0])*627.5;
//checking for severe DFT failures
  for (int i=1;i<nts;i++)
    if (dftelist[i]<dftelist[0]-500)
      dftelist[i] = 25000;
//  dftelist[0] = 0.;
  printf("\n");
  for (int i=1;i<nts;i++)
    if (savelist[i])
      printf(" E[%i]: %1.3f",i,dftelist[i]);
  printf("\n");

  delete [] rclist;

  return 0;
}



int ZStruct::dft_opt_para() {

  int maxprocs = MAXPARA;
  int numq = DFTPARA;

  fflush(stdout);
  printf("\n doing dft opt \n");

  DFT dft1;
  dft1.init(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

  string* rclist = new string[nts]; 
  printf(" DFT preparing structure for opt ");
  for (int i=0;i<nts;i++)
  {
    if (savelist[i])
    {
      printf(" %i",i);
      string nstr=StringTools::int2str(i,4,"0"); 
      string filename = "scratch/paradft"+nstr;
      //rlist[i%maxprocs] = filename;
      rclist[i] = filename;
      dft1.reset(geoms[i].natoms,geoms[i].anumbers,geoms[i].anames,geoms[i].coords);
      if (i!=0)
        dft1.opt_dnr(filename);
      else
        dft1.sp_dnr(filename);
    }
    else
    {
//      printf(" not working on structure %i \n",i);
      dftelist[i] = -1;
    }
  }
  printf("\n");


  ofstream cmdfile;
  string cmdfile_string = "scratch/go_dft";

  cmdfile.open(cmdfile_string.c_str());
  printf("\n printing job array numbers for DFT OPT \n");
  cmdfile << "#PBS -t ";
#if 0
  int next = 0;
  for (int i=0;i<nts;i++)
    if (savelist[i])
    {
      cmdfile << i;
      next = i+1; break;
    }
  for (int i=next;i<nts;i++)
    if (savelist[i])
      cmdfile << "," << i;
//      printf("%i,",i);
#else
  cmdfile << "0-" << nts-1 << endl;
#endif
  cmdfile.close();
  printf("\n");
#if !SKIPDFT
  string cmd = "./qmaked";
  system(cmd.c_str());
  cmd = "qsub "+cmdfile_string+".qsh";
  system(cmd.c_str());
#endif

  int qnotdone = true;
  int* dftdone = new int[nts];
  for (int i=0;i<nts;i++) dftdone[i]=0;

#if !SKIPDFT
  int max_wait = MAX_TIME_WAIT;
  int tc = 0;
  do {
    tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
    sleep(10);
    qnotdone = false;
    for (int i=0;i<nts;i++)
    {
      if (savelist[i] && !dftdone[i])
      {
        string nstr=StringTools::int2str(i,1,"0");
        string dftfile_string = "scratch/dftdone"+nstr;
        struct stat sts;
        if (stat(dftfile_string.c_str(), &sts) != -1)
        {
          dftdone[i]=1;
          printf(" done with DFT OPT %i",i);
        }
        else
          qnotdone = true;
      } // if not done
    } //loop i over nts
  } while (qnotdone);
  printf("\n");
  sleep(5);
#endif

  dftelist[0] = dft1.get_energy(rclist[0]);
  for (int i=1;i<nts;i++)
    if (savelist[i])
    {
      dftelist[i] = dft1.get_opt_energy(rclist[i]);
      dft1.get_structure(rclist[i],geoms[i].coords);
      cout << " xyz read in from DFT " << i << " E is " << dftelist[i] << endl;
      //print_xyz_gen(geoms[i].natoms,geoms[i].anames,geoms[i].coords);
    }

  for (int i=1;i<nts;i++)
    dftelist[i] = (dftelist[i] - dftelist[0])*627.5;
//checking for severe DFT failures
  for (int i=1;i<nts;i++)
    if (dftelist[i]<dftelist[0]-500)
      dftelist[i] = 25000;
//  dftelist[0] = 0.;
  printf("\n");
  for (int i=1;i<nts;i++)
    if (savelist[i])
      printf(" E[%i]: %1.3f",i,dftelist[i]);
  printf("\n");

  delete [] rclist;

  return 0;
}


int ZStruct::dft_ts_para() {

//  printf("\n dft_ts_para needs ts initial geometries \n");
//  return 0;

  int maxprocs = MAXPARA;
  int numq = DFTPARA;

  fflush(stdout);
  printf("\n doing dft ts \n");

  DFT dft1;
  dft1.init(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

//  string* rlist = new string[maxprocs];
  string* rclist = new string[nts]; 
  printf(" DFT TS preparing structure ");
  for (int i=1;i<nts;i++)
  {
    if (savelist[i])
    {
      printf(" %i",i);
      string nstr=StringTools::int2str(i,4,"0"); 
      string filename = "scratch/paradftts"+nstr;
      //rlist[i%maxprocs] = filename;
      rclist[i] = filename;
      dft1.reset(geoms[i].natoms,geoms[i].anumbers,geoms[i].anames,geoms[i].coordsts);
      dft1.ts_dnr(filename);
    }
    else
    {
//      printf(" not working on structure %i \n",i);
      dftelistts[i] = -1;
    }
  }
  printf("\n");


  ofstream cmdfile;
  string cmdfile_string = "scratch/go_dft_ts";

  printf("\n printing job array numbers for DFT TS \n");
#if 0
  cmdfile.open(cmdfile_string.c_str());
  cmdfile << "#PBS -t ";
  int next = 0;
  for (int i=0;i<nts;i++)
    if (savelist[i])
    {
      cmdfile << i;
      next = i+1; break;
    }
  for (int i=next;i<nts;i++)
    if (savelist[i])
      cmdfile << "," << i;
//      printf("%i,",i);
  cmdfile.close();
#else
  int continue1 = 1;
  int start = 1;
  int next = 1;
  do 
  {
    cmdfile.open(cmdfile_string.c_str());
    cmdfile << "#PBS -t ";
    int max = min(start + 100,nts);
    next = max;
    for (int i=start;i<max;i++)
      if (savelist[i])
      {
        cmdfile << i;
        next = i+1; break;
      }
    for (int i=next;i<max;i++)
      if (savelist[i])
        cmdfile << "," << i;
    cmdfile << endl;
    cmdfile.close();
    string cmd = "./qmakedts";
    system(cmd.c_str());
    cmd = "qsub "+cmdfile_string+".qsh";
#if !SKIPDFT
    system(cmd.c_str());
#endif
    start = max;
    if (max==nts) continue1 = 0;
  } while (continue1);
#endif


  printf("\n");

  int qnotdone = true;
  int* dfttsdone = new int[nts];
  for (int i=0;i<nts;i++) dfttsdone[i]=0;

#if !SKIPDFT
  int max_wait = MAX_TIME_WAIT;
  int tc = 0;
  do {
    tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
    sleep(10);
    qnotdone = false;
    for (int i=1;i<nts;i++)
    {
      if (savelist[i] && !dfttsdone[i])
      {
        string nstr=StringTools::int2str(i,1,"0");
        string dfttsfile_string = "scratch/dfttsdone"+nstr;
        struct stat sts;
        if (stat(dfttsfile_string.c_str(), &sts) != -1)
        {
          dfttsdone[i]=1;
          printf(" done with DFT TS %i",i);
        }
        else
          qnotdone = true;
      } // if not done
    } //loop i over nts
  } while (qnotdone);
  printf("\n");
#endif

  int* tsavelist = new int[nts];
  for (int i=1;i<nts;i++)
    if (savelist[i])
    {
      dftelistts[i] = dft1.get_energy_ts(rclist[i]);
      dft1.get_structure(rclist[i],geoms[i].coordsts);
      tsavelist[i] = dft1.converged;
    }


  printf(" done getting energy ts \n");

  for (int i=1;i<nts;i++)
    if (savelist[i])
      dftelistts[i] = (dftelistts[i] - dftelist[0])*627.5;
  for (int i=1;i<nts;i++)
    if (dftelistts[i]>2000)
      dftelistts[i] = -1.;
  
  printf("\n");
  for (int i=1;i<nts;i++)
    if (savelist[i])
    {
      if (tsavelist[i])
        printf(" E(TS/DFT)[%i]: %1.1f",i,dftelistts[i]);
      else
        printf(" E*(TS/DFT)[%i]: %1.1f",i,dftelistts[i]);
      savelist[i]=tsavelist[i];
    }
  printf("\n");
  printf(" savelist updated for convergence\n");

  delete [] tsavelist;
  delete [] rclist;

  return 0;
}



int ZStruct::gsm_dft_para(){

  int maxprocs = MAXPARA;
  int gsmprocs = DFTPARA;
  int nkeep = nts;

  string* rclist = new string[nts+100];
  ofstream cmdfile;
  string cmdfile_string = "scratch/go_gsm_dft";

  int totalc = 0;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      totalc++;
  printf("\n found %i structures for DFT GSM \n",totalc);

  printf("\n printing job array numbers for GSM \n");
#if 0
  cmdfile.open(cmdfile_string.c_str());
  cmdfile << "#PBS -t ";
  int next = 1;
  for (int i=1;i<nts;i++)
    if (savelist[i])
    {
      cmdfile << i;
      next = i+1; break;
    }
  for (int i=next;i<nts;i++)
    if (savelist[i])
      cmdfile << "," << i;
//      printf("%i,",i);
#else
  int continue1 = 1;
  int start = 1;
  int next = 1;
  do 
  {
    cmdfile.open(cmdfile_string.c_str());
    cmdfile << "#PBS -t ";
    int max = min(start + 100,nts);
    next = max;
    for (int i=start;i<max;i++)
      if (savelist[i])
      {
        cmdfile << i;
        next = i+1; break;
      }
    for (int i=next;i<max;i++)
      if (savelist[i])
        cmdfile << "," << i;
    cmdfile << endl;
    cmdfile.close();
    string cmd = "./qmakeg";
    system(cmd.c_str());
    cmd = "qsub "+cmdfile_string+".qsh";
#if !SKIPDFT
    system(cmd.c_str());
#endif
    start = max;
    if (max==nts) continue1 = 0;
  } while (continue1);
#endif

  printf("\n");

  int qnotdone = true;
  int* gsmdone = new int[nts];
  for (int i=0;i<nts;i++) gsmdone[i]=0;

  int max_wait = MAX_TIME_WAIT*10;
  int tc = 0;
  do {
    tc++; if (tc>max_wait) { printf(" done waiting! \n"); break; }
    sleep(10);
    qnotdone = false;
    for (int i=1;i<nts;i++)
    {
      if (savelist[i] && !gsmdone[i])
      {
        string nstr=StringTools::int2str(i,1,"0");
        string gsmfile_string = "scratch/gsmdone"+nstr;
        struct stat sts;
        if (stat(gsmfile_string.c_str(), &sts) != -1)
        {
          gsmdone[i]=1;
          printf(" done with GSM %i",i);
        }
        else
          qnotdone = true;
      } // if not done
    } //loop i over nts
  } while (qnotdone);
  printf("\n");

#if 0
  //Old parallel way of running suitable for 1 node
  int maxcycles = nts/maxprocs + 1;
  int done = 1;
  for (int i=0;i<maxcycles;i++)
  { 
    string nstr;
    string file1;
    string file2;
    int c = 0;
    if (done<nts)
    {
      fflush(stdout);
      cmdfile.open(cmdfile_string.c_str());
      for (int j=0;j<maxprocs;)
      {
        nstr=StringTools::int2str(done+j+c,4,"0");
        file2 = "scratch/paragsmdft"+nstr+".out";
        rclist[done+j+c] = file2;
        if (savelist[done+j+c])
        {
          cmdfile << " ./gstringq.exe " << done+j+c << " " << gsmprocs << " > " << file2 << " &" << endl;
          j++;
        }
        else c++;
        if (done+c+j>nts) break;
      }
//      printf(" cmd file almost done \n");
//      cout << " cmd file almost done " << endl;
      fflush(stdout);
      cmdfile << " wait " << endl;
      cmdfile.close();

      printf(" executing up to %i GSM runs (DFT) \n",maxprocs);
      string cmd = "chmod ug+rx scratch/go_gsm_dft";
      system(cmd.c_str());
#if !SKIPDFT
      cmd = "scratch/go_gsm_dft";
      system(cmd.c_str());
#endif
    }
    done += maxprocs+c;

  }
  fflush(stdout);
#endif

  printf(" done with gstring calcs \n\n");

  fflush(stdout);

  PGSM gsm;
  gsm.alloc(geoms[0].natoms);

  for (int i=1;i<nts;i++)
  {
    if (savelist[i])
    {
      string nstr = StringTools::int2str(i,4,"0");
      string oname = "scratch/tsq"+nstr+".xyz";
      ifstream output(oname.c_str(),ios::in);
      string line;
      getline(output,line);
      getline(output,line);
//    cout << " RR " << line << endl;

      double energyts = (atof(line.c_str())-dftelist[0]) * 627.5;
      if (abs(energyts)<0.00001)
      {
        printf(" GSM failed \n");
        energyts = 25000;
      }
      else
        printf(" final GSM energy for %i is: %3.1f \n",i,energyts);

// need to retrieve final TS geometry, write to xyz (with ints)
      gsm.xyz_read(oname);
      geoms[i].dftgsmenergy = energyts;
      dftelistgsm[i] = energyts;
      for (int j=0;j<3*geoms[i].natoms;j++)
        geoms[i].coordsts[j] = gsm.xyzts[j];

    } //if savelist[i]
  } //loop i over all nts

#if 0
  DFT dft1;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      dftelistts[i] = dft1.get_energy(rclist[i]);
#endif

  delete [] rclist;
  return nkeep;
}





int ZStruct::remove_duplicates(int nbranch) {

  printf(" removing duplicate structures \n");
  printf("  currently only based on coordn \n");

  int n = nbranch; // number of unique structures
//  return nbranch;

  for (int i=1;i<nbranch;i++)
    repeat[i] = -1;

  for (int i=niso;i<niso+nbranch;i++) //new structure to check
  {
//    printf(" checking structure # %i \n",i);
    int found = false;
    for (int j=0;j<i;j++) //existing structures
    {
      int echeck = 0;
      if (abs(elist[i]-elist[j])*10<10) echeck = 1;
      int nsamecoordn = 0;
      for (int k=0;k<geoms[i].natoms;k++)
        if (geoms[i].coordn[k]==geoms[j].coordn[k]) nsamecoordn++;

//      if (echeck) // based on energy for now
      if (nsamecoordn == geoms[i].natoms && echeck) 
      {
        printf(" e[%i] - e[%i]: %1.2f ",i,j,elist[i]-elist[j]);
        printf(" structure %i same as structure %i \n",i,j);
        found = true;
        n--;
        repeat[i] = j;
        break;
      }
    } // loop j over prev geoms
  } //loop i over nbranch  
 
  return n; 
}




void ZStruct::save_unique_multi(int int0, int int1) {

  if (int0==int1) return;
  printf(" in save_unique_multi() int0: %i int1: %i \n",int0,int1);

  int nremain = 0;
  for (int i=int0;i<int1;i++)
    if (savelist[i])
      nremain++;
#if 0
  for (int i=0;i<int1;i++)
  if (savelist[i])
    printf(" savelist[%i] 1 \n",i);
#endif

  printf("\n comparing %i new structs \n",nremain);

  for (int i=int0;i<int1;i++)
    for (int j=0;j<i;j++)
      if (savelist[i] && savelist[j])
      {
        int diff;
//        double edft1 = (dftelist[i] - dftelist[0])*627.5;
//        double edft2 = (dftelist[j] - dftelist[0])*627.5;
        double edft1 = dfterel[i];
        double edft2 = dfterel[j];
        diff = diff_structureiq(geoms[0].natoms,geoms[0].anames,geoms[0].anumbers,geoms[i].coords,geoms[j].coords);
        if (diff==0 && close_value(edft1,edft2,EVSORIGTOL))
        {
          printf(" pair %i %i is not unique (%i), E: %1.1f %1.1f EGSM: %1.1f %1.1f \n",i,j,diff,edft1,edft2,(dftelistgsm[i]-dftelist[0])*627.5,(dftelistgsm[j]-dftelist[0])*627.5);
          cMap[i*NUM_INTS+j] = 3;
          nremain--;
          savelist[i] = 0;
        }
        else
          printf(" pair %i %i is unique (%i), E: %1.1f %1.1f EGSM: %1.1f %1.1f \n",i,j,diff,edft1,edft2,(dftelistgsm[i]-dftelist[0])*627.5,(dftelistgsm[j]-dftelist[0])*627.5);
      }

  printf("\n nremain after removing pairs: %i \n",nremain);
#if 0
  for (int i=0;i<nts;i++)
    if (savelist[i])
      printf(" %i remains \n",i);
#endif

  return;
}



void ZStruct::save_unique(int type) {

  //printf(" in save_unique() \n");

  int nremain = 0;
  for (int i=0;i<nts;i++)
    if (savelist[i])
      nremain++;

  printf("\n comparing %i structs \n",nremain);

  for (int i=0;i<nts;i++)
    for (int j=0;j<i;j++)
      if (savelist[i] && savelist[j])
      {
        int samec;
        if (type==0) samec = diff_structurecq(geoms[0].natoms,geoms[0].anames,geoms[0].anumbers,geoms[i].coords,geoms[j].coords);
        else samec = diff_structureq(geoms[0].natoms,geoms[0].anames,geoms[0].anumbers,geoms[i].coords,geoms[j].coords);
        if (samec==0 && close_value(dftelist[i],dftelist[j],EVSORIGTOL/2))
        {
//          printf(" pair %i %i is not unique, E's: %1.1f %1.1f \n",i,j,elist[i],elist[j]);
          printf(" pair %i %i is not unique, E: %1.1f %1.1f EGSM: %1.1f %1.1f \n",i,j,dftelist[i],dftelist[j],dftelistgsm[i],dftelistgsm[j]);
          nremain--;
          if (dftelistgsm[i]<dftelistgsm[j])
            savelist[j] = 0;
          else
            savelist[i] = 0;
        }
      }

  printf(" nremain after removing pairs: %i \n",nremain);
#if 0
  for (int i=0;i<nts;i++)
    if (savelist[i])
      printf(" %i remains \n",i);
#endif

  return;
}



int ZStruct::diff_structure(int natoms, string* anames, int* anumbers, double* xyz1, double* xyz2) {
 
  int diff = 0;

  ICoord test1;
  ICoord test2;
  test1.init(natoms,anames,anumbers,xyz1);
  test2.init(natoms,anames,anumbers,xyz2);

  //printf(" testing diff_structure \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms) 
  {
    //printf(" same coordination #s, checking connectivity \n");
    int found = 0;
    for (int i=0;i<test1.nbonds;i++)
    {
      found = 0;
      if (test2.bond_exists(test1.bonds[i][0],test1.bonds[i][1]))
        found = 1;
      else
        break;
    }
    if (found) diff = 0;
    else diff = 2;
  }
  else
  {
    diff = natoms-nsamecoordn;
    printf(" different coordination# (by %i) \n",diff);
    nmopacdiff++;
//    test1.print_xyz();
//    test2.print_xyz();
//    test1.print_bonds();
//    test2.print_bonds();
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}


int ZStruct::diff_structureq(int natoms, string* anames, int* anumbers, double* xyz1, double* xyz2) {
 
  int diff = 0;
  ICoord test1;
  ICoord test2;
  test1.init(natoms,anames,anumbers,xyz1);
  test2.init(natoms,anames,anumbers,xyz2);

  //printf(" testing diff_structure \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms) 
  {
    //printf(" same coordination #s, checking connectivity \n");
    int found = 0;
    for (int i=0;i<test1.nbonds;i++)
    {
      found = 0;
      if (test2.bond_exists(test1.bonds[i][0],test1.bonds[i][1]))
        found = 1;
      else
        break;
    }
    if (found) diff = 0;
    else diff = 2;
  }
  else
  {
    diff = natoms-nsamecoordn;
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}



int ZStruct::diff_structureiq(int natoms, string* anames, int* anumbers, double* xyz1, double* xyz2) {
 
  printf("\n\n diff_structureiq test \n");

  int diff = 0;
  double* xyz3 = new double[3*natoms];
  for (int i=0;i<3*natoms;i++)
    xyz3[i] = xyz2[i];

  ICoord test1;
  ICoord test2;
  test1.init(natoms,anames,anumbers,xyz1);
  test2.init(natoms,anames,anumbers,xyz2);

  int cont = 1;

// diff struct without rotation
  if (diff_structureq(natoms,anames,anumbers,xyz1,xyz2))
    printf(" first test no info \n");
  else
  {
    printf(" same found \n");
    diff = 0;
    cont = 0;
  }

// make element list
  int* atype = new int[natoms];
  for (int i=0;i<natoms;i++) 
    atype[i] = -1;
  atype[0] = anumbers[0];
  int ntype = 1;
  for (int i=0;i<natoms;i++) 
  {
    int found = 0;
    for (int j=0;j<ntype;j++) 
    if (anumbers[i]==atype[j])
      found = 1;
    if (!found)
      atype[ntype++] = anumbers[i];
  }  
  printf(" printing element types:");
  for (int i=0;i<ntype;i++)
    printf(" %i",atype[i]);
  printf("\n");

// count number of each coordn for each element
  int* coordn1 = new int[6];
  int* coordn2 = new int[6];
  for (int i=0;i<ntype;i++)
  {
    for (int j=0;j<6;j++) coordn1[j] = 0;
    for (int j=0;j<6;j++) coordn2[j] = 0;
    for (int j=0;j<natoms;j++)
    if (anumbers[j] == atype[i])
    {
      coordn1[test1.coordn[j]]++;
      coordn2[test2.coordn[j]]++;
    }
//    for (int j=0;j<6;j++)
//      printf(" AtNum: %i coord#: %i coordn1: %i coordn2: %i \n",atype[i],j,coordn2[j],coordn2[j]);
    for (int j=0;j<6;j++)
    if (coordn1[j]!=coordn2[j])
      diff++;
  }
  printf(" after checking elements for same coordn, diff: %i \n",diff);
  if (diff>0) cont = 0;

  if (cont)
  printf(" diff_structureiq element swapping \n");

// rotate heavy atoms by coordn
  if (cont)
  for (int i=0;i<ntype;i++)
  if (atype[i]!=1)
  {
    for (int j=0;j<natoms;j++)
    if (anumbers[j]==atype[i])
    if (test1.coordn[j]!=test2.coordn[j])
    {
      printf(" need a swap partner for %i \n",j);
      int s1 = -1;
      int cn = test1.coordn[j];
      for (int k=0;k<natoms;k++)
      if (j!=k && cn==test2.coordn[k] && anumbers[k]==atype[i])
      {
        s1 = k;
        printf(" found: %i \n",s1);
        break;
      }
      swap_atoms(test2.coords,j,s1);
      test2.ic_create();
    }
  }

  int* clist1 = new int[6*natoms];
  int* clist2 = new int[6*natoms];
  for (int i=0;i<6*natoms;i++) clist1[i] = -1;
  for (int i=0;i<6*natoms;i++) clist2[i] = -1;
  if (cont)
  for (int i=0;i<natoms;i++)
  {
    int nconn1 = 0;
    int nconn2 = 0;
    for (int j=0;j<test1.nbonds;j++)
    {
      if (test1.bonds[j][0] == i)
        clist1[6*i+nconn1++] = anumbers[test1.bonds[j][1]];
      else if (test1.bonds[j][1] == i)
        clist1[6*i+nconn1++] = anumbers[test1.bonds[j][0]];
      if (test2.bonds[j][0] == i)
        clist2[6*i+nconn2++] = anumbers[test2.bonds[j][1]];
      else if (test2.bonds[j][1] == i)
        clist2[6*i+nconn2++] = anumbers[test2.bonds[j][0]];
    }
  }
  if (cont)
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<5;j++)
    if (clist1[6*i+j]<clist1[6*i+j+1])
    {
      //printf("  if %i %i %i %i \n",i,j,clist1[6*i+j],clist1[6*i+j+1]);
      int tmp = clist1[6*i+j];
      clist1[6*i+j] = clist1[6*i+j+1];
      clist1[6*i+j+1] = tmp;
      j = -1;
    }
    for (int j=0;j<5;j++)
    if (clist2[6*i+j]<clist2[6*i+j+1])
    {
      //printf("  if %i %i %i %i \n",i,j,clist2[6*i+j],clist2[6*i+j+1]);
      int tmp = clist2[6*i+j];
      clist2[6*i+j] = clist2[6*i+j+1];
      clist2[6*i+j+1] = tmp;
      j = -1;
    }
  }
#if 0
  for (int i=0;i<natoms;i++)
  {
    printf(" atom %i has element connections: %i %i %i %i %i %i \n",i,clist1[6*i+0],clist1[6*i+1],clist1[6*i+2],clist1[6*i+3],clist1[6*i+4],clist1[6*i+5]);
    printf(" atom %i has element connections: %i %i %i %i %i %i \n",i,clist2[6*i+0],clist2[6*i+1],clist2[6*i+2],clist2[6*i+3],clist2[6*i+4],clist2[6*i+5]);
  }
#endif

  if (diff>0) cont = 0;
  int* fixed = new int[natoms];
  for (int i=0;i<natoms;i++)
    fixed[i] = test1.coordn[i];

  int* bondok = new int[test1.nbonds];
  for (int i=0;i<test1.nbonds;i++) bondok[i] = 0;

// look for same bonds
  int nsamecoordn = 0;
  int found = 0;
  if (cont)
  for (int i=0;i<test1.nbonds;i++)
  {
    if (test2.bond_exists(test1.bonds[i][0],test1.bonds[i][1]))
    {
//      printf(" found bond: %i %i \n",test1.bonds[i][0],test1.bonds[i][1]);
      fixed[test1.bonds[i][0]]--;
      fixed[test1.bonds[i][1]]--;
      bondok[i] = 1;
    }
  }

// heavy atoms not fully matched to bonds, check local attachment
  int nsame = 0;
  if (cont)
  for (int i=0;i<ntype;i++)
  if (atype[i]!=1)
  {
    for (int j=0;j<natoms;j++)
    if (anumbers[j]==atype[i] && fixed[j]>0)
    {
      int same1 = 1; 
      for (int l=0;l<6;l++)
      if (clist1[6*j+l]!=clist2[6*j+l])
        same1 = 0;
      if (same1)
      {
        printf(" atom %i same elemental coord sphere \n",j);
        fixed[j] = 0;
      }
      if (!same1)
      {
        int found = 0;
        for (int k=0;k<natoms;k++)
        if (j!=k && anumbers[k]==atype[i] && fixed[k]>0)
        {
          int same2 = 1;
          for (int l=0;l<6;l++)
          if (clist1[6*j+l]!=clist2[6*k+l])
            same2 = 0;
          if (same2)
          {
            found = 1;
            printf(" atoms %i %i could be swapped \n",j,k);
            fixed[k] = 0;
            break;
          }
        }
        if (!found) //CPMZ heavy atom elemental coord sphere criterion 
          diff++;
      } //if !same1
    } //loop over heavy atom j
  } //loop over elements i

  if (diff>0) cont = 0;
#if 0
  if (cont)
  for (int i=0;i<natoms;i++)
    printf(" test1.coordn[%i]: %i fixed[%i]: %i \n",i,test1.coordn[i],i,fixed[i]);
#endif

  if (cont)
  {
   // test1.print_bonds();
   // test2.print_bonds();
    test1.print_xyz();
    test2.print_xyz();
  }

#if 0
//not using?
  int* ne = new int[ntype];
  for (int i=0;i<ntype;i++)
    ne[i] = 0;
  int** elemlist = new int*[ntype];
  for (int i=0;i<ntype;i++)
    elemlist[i] = new int[natoms];
  for (int i=0;i<ntype;i++)
  for (int j=0;j<natoms;j++)
  if (anumbers[j] == atype[i])
    elemlist[i][ne[i]++] = j;
#endif


//rotate within same element then compare...
#if 0
//not fully implemented
  if (cont)
  for (int b=0;b<test1.nbonds;b++)
  if (!bondok[b])
  {
    int a1 = test1.bonds[b][0];
    int a2 = test1.bonds[b][1];
   // printf(" looking for bond %i(%s) %i(%s) \n",a1,anames[a1].c_str(),a2,anames[a2].c_str());
   // printf(" fixed[a1]: %i fixed[a2]: %i \n",fixed[a1],fixed[a2]);
    int a12 = -1;
    int a21 = -1;
    if (fixed[a1]>0 && anumbers[a1]==1)
    {
      a12 = a1;
      a21 = a2;
    }
    else if (fixed[a2]>0 && anumbers[a2]==1)
    {
      a12 = a2;
      a21 = a1;
    }
    if (a12!=-1)
    for (int i=0;i<natoms;i++)
    if (i!=a12 && fixed[i]>0 && anumbers[i]==anumbers[a12])
    {
      //printf(" swap candidate: %i %i \n",a12,i);
      if (test2.bond_exists(i,a21))
      {
        printf(" swap okay: %i %i becomes %i %i \n",a12,a21,i,a21);
//        fixed[i]--;
//        fixed[a21]--;
        bondok[b] = 1;
        break;
      }
    } //if candidate i for a12
  } //loop over not okay bonds

  if (cont)
  for (int i=0;i<test1.nbonds;i++)
  if (!bondok[i])
  {
    printf(" bond %i not found \n",i);
    diff++;
  }
  printf(" printing which atoms may be rotated \n");
  for (int i=0;i<natoms;i++)
    printf(" test1.coordn[%i]: %i fixed[%i]: %i \n",i,test1.coordn[i],i,fixed[i]);
#endif

  delete [] atype;
  delete [] clist1;
  delete [] clist2;
  delete [] fixed;
  delete [] bondok;
  delete [] coordn1;
  delete [] coordn2;
  delete [] xyz3;
  test1.freemem();
  test2.freemem();
  
  return diff;
}

void ZStruct::swap_atoms(double* xyz, int a1, int a2) {

  double* xyz3 = new double[3];

  xyz3[0] = xyz[3*a1+0];
  xyz3[1] = xyz[3*a1+1];
  xyz3[2] = xyz[3*a1+2];
  xyz[3*a1+0] = xyz[3*a2+0];
  xyz[3*a1+1] = xyz[3*a2+1];
  xyz[3*a1+2] = xyz[3*a2+2];
  xyz[3*a2+0] = xyz3[0];
  xyz[3*a2+1] = xyz3[1];
  xyz[3*a2+2] = xyz3[2];

  delete [] xyz3;
  return;
}

//uses coordination # only
int ZStruct::diff_structurec(int natoms, string* anames, int* anumbers, double* xyz1, double* xyz2) {
 
  int diff = 0;
  ICoord test1;
  ICoord test2;
  test1.init(natoms,anames,anumbers,xyz1);
  test2.init(natoms,anames,anumbers,xyz2);

  //printf(" testing diff_structurecq \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms) 
  {
    printf(" same coordination #s \n");
    diff = 0;
  }
  else
  {
    diff = natoms-nsamecoordn;
    printf(" different coordination# (by %i) \n",diff);
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}


//uses coordination # only
int ZStruct::diff_structurecq(int natoms, string* anames, int* anumbers, double* xyz1, double* xyz2) {
 
  int diff = 0;
  ICoord test1;
  ICoord test2;
  test1.init(natoms,anames,anumbers,xyz1);
  test2.init(natoms,anames,anumbers,xyz2);

  //printf(" testing diff_structurecq \n");

  int nsamecoordn = 0;
  for (int k=0;k<natoms;k++)
    if (test1.coordn[k]==test2.coordn[k]) nsamecoordn++;
//  printf(" nsamecoordn: %i \n",nsamecoordn);
  if (nsamecoordn==natoms) 
  {
    //printf(" same coordination #s \n");
    diff = 0;
  }
  else
  {
    diff = natoms-nsamecoordn;
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
}




int ZStruct::init(string xyzfile, string xyzlist, double emax, double emin){
  
  //printf(" setting oemax to: %1.4f \n",emaximum);
  oemax = emax;
  oemin = emin;
  init(xyzfile,xyzlist);

  return 0;
}

int ZStruct::init(string xyzfile, string xyzlist){

  printf(" initializing ZStruct ! \n");
  cout << " xyzfile: " << xyzfile << " xyzlist: " << xyzlist << endl;

  printf("\n\n -------------------------------------- \n");
  printf(" ---- 1. reading in database ---- \n");

  basis1.init(xyzlist);

  printf("\n done with database generation \n");

  printf("\n\n -------------------------------------- \n");
  printf(" ---- 2. reading starting geometry ---- \n");
  int done = geom0.init(xyzfile);

  geom0.print_ic();

  printf(" allocating pointers for isomers \n");
  geoms = new ICoord[NMAXSTRUCTS];
//  geoms[0] = geom0;
  geoms[0].init(geom0.natoms,geom0.anames,geom0.anumbers,geom0.coords);
  geomshadow.init(geom0.natoms,geom0.anames,geom0.anumbers,geom0.coords);

  Mopac mopt;
  mopt.alloc(geoms[0].natoms);
  mopt.reset(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);
  //geoms[0].seenergy = mopt.opt("scratch/xyzfile.xyzm_initial");
  //Assuming starting geometry is optimized
//  geoms[0].reset(mopt.xyz);

  geom1.alloc(geom0.natoms);
  geom1r.alloc(geom0.natoms);

  nrepeat = 0;
  repeat = new int[NMAXSTRUCTS];
  for (int i=0;i<NMAXSTRUCTS;i++)
    repeat[i] = -1;

  frozen = new int[geom0.natoms];
  for (int i=0;i<geom0.natoms;i++) frozen[i] = 0;

  nfail = 0;
  nmopacdiff = 0;
  nmopaclargediff = 0;

  for (int i=1;i<NMAXSTRUCTS;i++)
    geoms[i].alloc(geom0.natoms);
  for (int i=1;i<NMAXSTRUCTS;i++)
    geoms[i].reset(geom0.natoms,geoms[0].anames,geoms[0].anumbers,geoms[0].coords);

 return 0;
}

void ZStruct::write_initialxyz(int id, int nsave) {
  
  double* e = new double[3];
  e[0] = e[1] = e[2] = -1;
  for (int i=1;i<nsave;i++)
  {
    e[0] = geoms[id].seenergy; e[1] = geoms[i].seenergy;
    string nstr=StringTools::int2str(i,4,"0"); // 2,"0" is number of 0s in total
    string xyzfilename = "scratch/initial"+nstr+".xyz";
    print_double_xyz_save(xyzfilename,geoms[id].natoms,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,geoms[i].coords,e);
  }
  delete [] e;
  return;
}



void ZStruct::read_frozen(int natoms){ 
   
  string xyzfile = "frozen.xyz";
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    cout << "!!!!Error opening frozen xyz file!!!!" << endl;
    exit(-1);
  } 
  
 // cout <<"  -reading file..." << endl;
  
  int nfrozen = 0;

  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    nfrozen=atoi(line.c_str());
  }
  cout << " nfrozen: " << nfrozen << endl;

  int* flist = new int[nfrozen];  

  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<nfrozen;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    flist[i] = atoi(line.c_str());
    //vector<string> tok_line = StringTools::tokenize(line, " \t");
    //flist[i]=atoi(tok_line[0].c_string());
  }
  
  infile.close();
  for (int i=0;i<nfrozen;i++)
    frozen[flist[i]] = 1;

/*
  for (int i=0;i<natoms;i++)
    printf(" atom %i is frozen? %i \n",i,frozen[i]);  
  printf("\n");
*/

  delete [] flist;
}   


void ZStruct::alloc_elists() {

  elist = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) elist[j]=-1;
  elistts = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) elistts[j]=-1;
  dftelist = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) dftelist[j]=-1;
  dftelistlst = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) dftelistlst[j]=-1;
  dftelistgsm = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) dftelistgsm[j]=-1;
  dftelistts = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) dftelistts[j]=-1;

  return;
}

void ZStruct::free_elists() {

  delete [] elist;
  delete [] elistts;
  delete [] dftelist;
  delete [] dftelistlst;
  delete [] dftelistgsm;
  delete [] dftelistts;

  return;
}



//      string cmd = "./dircreate "+num+" "+"10000"+" "+num2;
//      ostringstream emax;
//      emax << oemax;
//      string emaxs = emax.str();
