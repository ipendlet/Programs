#include "zstruct.h"
#include "utils.h"
#include "dft.h"
#include "pgsm.h"
#include "constants.h"
#include "print.h"
using namespace std;


void ZStruct::go_mech_it(int nsteps) {

  if (nsteps>1)
  {
    printf("\n\n*****************************\n");
    printf("  Getting a mechanism! \n");
    printf("*****************************\n");

    printf(" oemax is: %1.4f \n",oemax);
    printf(" Steps remaining: %i \n",nsteps);
  }

  int* final_ints = new int[1000];
  int nkept = do_one_step(final_ints);

  printf("\n nkept: %i \n",nkept);

 // in the future, make a compendium of intermediates already examined
 //  to avoid repeats

  if (nkept > 30)
  {
    nkept=0;
    printf("\n too many structures! \n");
  }
 //preparing folders for next runs
  if (nsteps>1)
  for (int i=0;i<nkept;i++)
  {
    printf(" decided to run search on intermediate: %i with E: %1.1f \n",final_ints[i],dftelist[final_ints[i]]);

    string num =  StringTools::int2str(final_ints[i],4,"0");
    string num2 = StringTools::int2str(nsteps-1,1,"0");
    string cmd = "./dircreate "+num+" "+"10000"+" "+num2;
    system(cmd.c_str());

    cmd = "qsub "+num+"/zstruct.qsh";
    system(cmd.c_str());
    printf("  submitted %s \n",cmd.c_str());

  }

  system("rm final.xyzm_*");
  system("rm final.xyz_*");
  system("rm molecule*");
//  system("rm scratch/xyzfile*");

  if (nsteps<2)
    printf("\n Done moving forward! \n");
  else
    printf("\n Branch jobs submitted, zstruct exiting \n");


  delete [] final_ints;
  free_elists();

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

  int ngen;
#if !SKIP_ISOMERS
  ngen = generate_isomers(0,nbranch); 
#else
  ngen = niso = SKIP_ISOMERS;
  printf(" \n WARNING: isomers not generated: %i \n\n",niso);
  fflush(stdout);
#endif
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

#if 0
//  save_unique(1);

  printf(" STOPPING AFTER MOPAC! \n");
  return nkept;
#endif
#if 0
// Debug for branch submission
  printf("\n determing ints to continue \n");
  for (int i=1;i<nts;i++)
    if (savelist[i]) 
    {
      final_ints[nkept++]=i;
    }

  save_unique(0);

  printf(" STOPPING AFTER MOPAC! \n");
  return nkept;
#endif
 
  fflush(stdout);
#if 1
//  dft_sp_para();
  dft_opt_para();
#else
  printf(" skipping dft opt \n");
#endif
  fflush(stdout);


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
  int natoms = geoms[0].natoms0;
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
  max_check = emin + DFT_EMAX;
  if (max_check<0) max_check=10; // avoid exothermicity problems

  printf("\n max E for saving: %1.1f \n",max_check);

  for (int i=1;i<niso;i++) savelist[i] = 0;
  for (int i=1;i<nts;i++)
  {
//    printf(" checking: %i emin emints: %1.2f %1.2f elist elistts %1.2f %1.2f \n",i,emin,emints,elist[i],elistt$
    if (dftelist[i]<max_check)
      savelist[i] = 1;
  }

  printf("\n SKIPPING SAVE_UNIQUE \n\n");
//  save_unique(1);
  fflush(stdout);

  printf(" saved structure %i with E(DFT) (%1.4f) E(PM6) %1.1f \n",0,dftelist[0],elist[0]-elist[0]);
  for (int i=1;i<nts;i++) 
    if (savelist[i])
      printf(" saved structure %i with E(DFT) %1.2f E(PM6) %1.1f \n",i,dftelist[i],elist[i]-elist[0]);

  fflush(stdout);

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
  write_initialxyz(0,niso);


#if 1
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

#if 0
//skipping SE GSM
  nkeep = generate_paths_se(0,nts);
#endif


  printf(" \n recommending structures for DFT GSM \n");
  double dfteminlst = 25000;
  for (int i=1;i<nts;i++)
    if (dftelistlst[i]<dfteminlst && dftelistlst[i]>MINTSE) dfteminlst = dftelistlst[i];
  double max_check_ts = dfteminlst + LST_EMAX_TS;

  savelist[0] = 1;
#if 0
  for (int i=1;i<nts;i++)
  {
//    printf(" checking: %i emin emints: %1.2f %1.2f elist elistts %1.2f %1.2f \n",i,emin,emints,dftelist[i],dftelistlst[i]);
    if (dftelistlst[i]<max_check_ts && savelist[i])
    {
      savelist[i] = 1;
//CPMZ expand criteria
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
  printf("\n max lst ts energy criteria: %1.1f \n",max_check_ts);
#else
  printf(" skipping lst energy check \n");
#endif
  for (int i=0;i<nts;i++) 
    if (savelist[i])
      printf(" saved structure %i with E %1.1f E(LST/TS) %1.1f \n",i,dftelist[i],dftelistlst[i]);





  fflush(stdout);
#if 1
//  nkeep = generate_paths_dft();
  nkeep = gsm_dft_para();
#else
  printf(" skipping gsm dft \n");
#endif
  fflush(stdout);


  printf(" \n recommending structures for final DFT TS search \n");
  double dftemingsm = 25000;
  for (int i=1;i<nts;i++)
    if (dftelistgsm[i]<dftemingsm && dftelistgsm[i]>MINTSE) dftemingsm = dftelistgsm[i];
  max_check_ts = dftemingsm + GSM_EMAX_TS;

  for (int i=1;i<nts;i++)
  {
//    printf(" checking: %i emin emints: %1.2f %1.2f elist elistts %1.2f %1.2f \n",i,emin,emints,dftelist[i],dftelistgsm[i]);
    if (dftelistgsm[i]<max_check_ts && savelist[i])
    {
      savelist[i] = 1;
//CPMZ expand criteria
      if (close_value(dftelistgsm[i],0.0,EVSORIGTOL) || dftelistgsm[i]<-EVSORIGTOL)
      {
        printf(" WARNING: probably duplicate structure (or GSM failed) %i \n",i);
        savelist[i] = 0;
      }
    }
    else
      savelist[i] = 0;
  }
  nkeep = 0;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      nkeep++;
  printf("\n max ts energy criteria after GSM: %1.1f \n",max_check_ts);
  printf(" saved %i structures for TS search \n",nkeep);
  for (int i=0;i<nts;i++) 
    if (savelist[i])
      printf(" saved structure %i with E %1.1f ETS(LST) %1.1f ETS(GSM) %1.1f \n",i,dftelist[i],dftelistlst[i],dftelistgsm[i]);

  fflush(stdout);
#if 1
  nkeep = dft_ts_para();
#else
  printf(" skipping ts search \n");
#endif
  fflush(stdout);


  double dftemints = 25000;
  for (int i=1;i<nts;i++)
    if (savelist[i])
      if (dftelistts[i]<dftemints && dftelistts[i]>MINTSE) dftemints = dftelistts[i];

//CPMZ change criteria!
  max_check_ts = dftemints + 10;
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

  printf("\n continue list generated\n");
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
    start = max;
    string cmd = "./qmakel";
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
  dft1.init(geoms[0].natoms0,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);
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
  dft1.init(geoms[0].natoms0,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);
//  for (int i=0;i<maxprocs;i++)
//    dft1[i].init(geoms[0].natoms,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

//  string* rlist = new string[maxprocs];
  string* rclist = new string[nts]; 
  printf(" DFT preparing structure ");
#if QSTART
  printf("using qstart0 ");
#endif
  for (int i=0;i<nts;i++)
  {
    if (savelist[i])
    {
      printf(" %i",i);
      string nstr=StringTools::int2str(i,4,"0"); 
      string filename = "scratch/paradft"+nstr;
      //rlist[i%maxprocs] = filename;
      rclist[i] = filename;
      dft1.reset(geoms[i].natoms0,geoms[i].anumbers,geoms[i].anames,geoms[i].coords);
      dft1.sp_dnr(filename,nstr);
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
  dft1.init(geoms[0].natoms0,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

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
      dft1.reset(geoms[i].natoms0,geoms[i].anumbers,geoms[i].anames,geoms[i].coords);
      if (i!=0)
        dft1.opt_dnr(filename,nstr);
      else
        dft1.sp_dnr(filename,nstr);
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
      printf(" about to get_opt_energy %s \n",rclist[i].c_str());
      fflush(stdout);
      dftelist[i] = dft1.get_opt_energy(rclist[i]);
      dft1.get_structure(rclist[i],geoms[i].coords);
      printf(" xyz read in from DFT %i E is %1.6f \n",i,dftelist[i]);
      //cout << " xyz read in from DFT " << i << " E is " << dftelist[i] << endl;
      fflush(stdout);
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
  dft1.init(geoms[0].natoms0,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

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
      dft1.reset(geoms[i].natoms0,geoms[i].anumbers,geoms[i].anames,geoms[i].coordsts);
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
        printf(" final GSM energy for %i is: %1.2f \n",i,energyts);

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




// old code
void ZStruct::begin() {
  

  printf("\n ZStruct begin() \n");

  int max_steps = 2;
//  int nbranch = 100;
  int nunique;
  niso = 1;

  int** path = new int*[15]; //holds id of geoms
  for (int i=0;i<15;i++)
    path[i] = new int[max_steps];
  
//  path[0][0] = geom0.id = geoms[0].id = 0;

  int* olist = new int[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) olist[j]=-1;
  elist = new double[NMAXSTRUCTS];
  for (int j=0;j<NMAXSTRUCTS;j++) elist[j]=-1;

  int ngen;
  ngen = generate_isomers(0,5000); // id 0 refers to geoms[0]
//  nunique = remove_duplicates(ngen,olist);
  nunique = ngen;
  path[0][1] = 1;

  int xxt; // for random geom choice
  srand(5); 

  int steps;
  for (int i=1;i<=max_steps;i++)  
  {
    printf("\n\n beginning step %i \n",i);
    steps = i;
    for (int j=0;j<NMAXSTRUCTS;j++) olist[j]=-1;
    ngen = generate_isomers(path[0][i],1000);
//    nunique = remove_duplicates(ngen,olist);
    nunique = ngen;
    printf(" found %i unique structures \n",nunique);
    if (nunique<1)
    {
      printf(" failed to find new structure \n"); 
      break;
    }
    xxt = rand()%1000;
    printf(" xxt: %i \n",xxt);
//    path[0][i+1] = niso - xxt;
//    printf(" choosing geometry: %i \n",niso-nbranch);
    path[0][i+1] = niso-nunique; //this is done stupidly

  }

//CPMZ need to check if repeat[i] is correct
  for (int i=0;i<niso;i++)
    printf(" repeated %i: %i \n",i,repeat[i]);

  printf("\n saving all low-energy structures \n");
  system("rm final.xyz*");
  for (int i=0;i<niso;i++)
  {
    string nstr=StringTools::int2str(i,4,"0"); 
    string finfile = "final.xyz_"+nstr;
//    if (repeat[i]<0)
      geoms[i].print_xyz_save(finfile);
  }
  system("cat final.xyz_* > final.xyz");

  return;
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
    printf(" same coordination #s, checking connectivity \n");
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
    printf(" same coordination #s, checking connectivity \n");
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
  }

  test1.freemem();
  test2.freemem();
  
  return diff;
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




int ZStruct::init(string xyzfile, string xyzlist, double emaximum){
  
  //printf(" setting oemax to: %1.4f \n",emaximum);
  oemax = emaximum;
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

  printf(" after init(), printing ic \n");
  geom0.print_ic();

  printf(" getting solvent \n");
  frozen = new int[geom0.natoms+4];
  for (int i=0;i<geom0.natoms+4;i++) frozen[i] = 0;
  solvent = new int[geom0.natoms+4];
  for (int i=0;i<geom0.natoms+4;i++) solvent[i] = 0;
  int nsol = read_solvent(geom0.natoms);
  read_frozen(geom0.natoms);

#if 0
  printf(" opting initial structure \n");
//  int nsol0 = read_solvent(geom0.natoms0);
  geom0.opt(solvent);
  //geom0.opt();
  //geom0.print_xyz();
  double* e = new double[4];
  e[0] = 1; e[1] = 2;
  print_double_xyz_save("opt.xyz",geom0.natoms0,geom0.anames,geom0.anumbers,geom0.coords0,geom0.coords,e);

  printf(" ending early! \n");
  exit(1);
#endif
 
  printf(" allocating pointers for isomers \n");
  geoms = new ICoord[NMAXSTRUCTS];
//  geoms[0] = geom0;
  geoms[0].init(geom0.natoms0,geom0.anames,geom0.anumbers,geom0.coords);
  geomshadow.init(geom0.natoms0,geom0.anames,geom0.anumbers,geom0.coords);
  //printf("\n printing geoms[0] after init \n");
  //geoms[0].print_ic();

  Mopac mopt;
  mopt.alloc(geoms[0].natoms0);
  mopt.reset(geoms[0].natoms0,geoms[0].anumbers,geoms[0].anames,geoms[0].coords);

  addSolventMopac(geoms[0].natoms0,mopt);

//  printf(" before initial mopt \n");
  geoms[0].seenergy = mopt.opt("scratch/xyzfile.xyzm_initial0");
  //Assuming starting geometry is optimized
//  geoms[0].reset(mopt.xyz);

  geom1.alloc(geom0.natoms0);
  geom1r.alloc(geom0.natoms0);
//  geom1.skipTMadd = 1;
//  geom1r.skipTMadd = 1;

  nrepeat = 0;
  repeat = new int[NMAXSTRUCTS];
  for (int i=0;i<NMAXSTRUCTS;i++)
    repeat[i] = -1;

  nfail = 0;
  nmopacdiff = 0;
  nmopaclargediff = 0;

  for (int i=1;i<NMAXSTRUCTS;i++)
    geoms[i].alloc(geom0.natoms0);
//  for (int i=1;i<NMAXSTRUCTS;i++)
//    geoms[i].skipTMadd = 1;

 return 0;
}


void ZStruct::addSolventMopac(int natoms, Mopac mopt) {

#if NOMOOPT
  int nsolfrz = 0;
  int* solfrz = new int[natoms];
  for (int i=0;i<natoms;i++)
  if (frozen[i] && solvent[i])
  {
    solfrz[i] = 1;
    nsolfrz++;
  }
  else
    solfrz[i] = 0;
  for (int i=0;i<natoms;i++)
  for (int j=0;j<i;j++)
  {
    if (solfrz[i] && !solfrz[j] && geoms[0].bond_exists(i,j))
    {
      printf(" unfreezing: %i \n",i);
      solfrz[i] = 0;
      nsolfrz--;
    }
  }
  mopt.addSolvent(solfrz,nsolfrz);
  delete [] solfrz;
#endif

  return;
}

void ZStruct::write_initialxyz(int id, int nsave) {
  
  double* e = new double[3];
  e[0] = e[1] = e[2] = -1;
  for (int i=1;i<nsave;i++)
  {
    e[0] = geoms[id].seenergy; e[1] = geoms[i].seenergy;
    string nstr=StringTools::int2str(i,4,"0"); // 2,"0" is number of 0s in total
    string xyzfilename = "scratch/initial"+nstr+".xyz";
    print_double_xyz_save(xyzfilename,geoms[id].natoms0,geoms[id].anames,geoms[id].anumbers,geoms[id].coords,geoms[i].coords,e);
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



int ZStruct::read_solvent(int natoms){ 
   
  for (int i=0;i<natoms;i++)
    solvent[i] = 0;

  string xyzfile = "solvent.xyz";
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    cout << " No solvent xyz file found." << endl;
    return 0;
  } 
  
 // cout <<"  -reading file..." << endl;
  
  int nsol = 0;

  string line;
  bool success=true;
  success=getline(infile, line);
  if (success){
    int length=StringTools::cleanstring(line);
    nsol=atoi(line.c_str());
  }
  cout << " nsolvent: " << nsol << endl;

  int* slist = new int[nsol];  

  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<nsol;i++){
    success=getline(infile, line);
    int length=StringTools::cleanstring(line);
    slist[i] = atoi(line.c_str());
  }
  
  infile.close();
  for (int i=0;i<nsol;i++)
    solvent[slist[i]] = 1;

#if 0
  for (int i=0;i<natoms;i++)
    printf(" atom %i is solvent? %i \n",i,solvent[i]);  
  printf("\n");
#endif
#if 1
  printf(" solvent atoms:");
  for (int i=0;i<natoms;i++)
  if (solvent[i])
    printf(" %i",i);
  printf("\n");
#endif

  delete [] slist;

  return nsol;
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

