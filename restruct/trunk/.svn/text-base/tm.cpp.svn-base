#include "icoord.h"
#include "utils.h"


void ICoord::TM_setup() {

  coord_num();
  //printf(" in TM_setup \n");
  TM_planes_2();
  TM_angles(); // which calls remove TM_pairs
  //printf(" in TM_setup, ic_create_nobonds() \n");

#if 1
  // replace deleted triangles
  for (int i=0;i<ntri;i++)
  {
    //printf(" checking triangle: %i %i %i \n",tmtriangles[i][1],tmtriangles[i][2],tmtriangles[i][3]);
    //printf(" checking triangleb: %i %i %i %i \n",tmtriangles[i][5],tmtriangles[i][6],tmtriangles[i][7],tmtriangles[i][8]);
    //print_bonds();
    if (!bond_exists(tmtriangles[i][7],tmtriangles[i][8]))
    if (bond_exists(tmtriangles[i][1],tmtriangles[i][2]))
    if (bond_exists(tmtriangles[i][2],tmtriangles[i][3]))
    if (coordn[tmtriangles[i][8]]==tmtriangles[i][6])
    {
      printf(" reattaching triangle: %i %i \n",tmtriangles[i][7],tmtriangles[i][8]);
      bonds[nbonds][0]=tmtriangles[i][7];
      bonds[nbonds][1]=tmtriangles[i][8];
      nbonds++;
    }
  }
#endif
 
  //print_bonds();
  //printf(" about to ic_create_nobonds \n");
  ic_create_nobonds();
  make_angles_tm();

  //printf(" in TM_setup, after ic_create, printing ic's \n");
  //print_bonds();
  //print_ic();
  ///mm_init_tm();
  //printf(" after ic_create_nobonds\n");

  return;
} 




void ICoord::TM_angles() {
         
  ntmnoangles = 0;
  ntmlinangles = 0;
  nangles = 0;

  //printf("\n in TM_angles \n");
           
  double* d = new double[16];
  for (int i=0;i<16;i++) d[i]=-1.0;
            
//  for (int i=0;i<4;i++)
  int tm1 = -1;
  int tm2 = -1;
  for (int i=natoms0;i<natoms;i++)
    if (isTM(i)) 
     tm2 = i;

  int a1,a2,a3,a4; 
  a1=a2=a3=a4=-1;
  for (int i=0;i<natoms0;i++)
  {
    if (isTM(i) && nplanes)
    {
      printf(" circle at %i: %i %i %i %i \n",i,planes[0][0],planes[0][1],planes[0][2],planes[0][3]);
      tm1 = i;
      a1 = planes[0][0];
      a2 = planes[0][1];
      a3 = planes[0][2];
      a4 = planes[0][3];
    }     
  } //loop i over natoms0

  
 
// first find atoms in plane that are far from equilibrium

  if (nplanes)
  {
    d[0] = distance(tm1,a1);
    d[1] = distance(tm1,a2);
    d[2] = distance(tm1,a3);
    d[3] = distance(tm1,a4);
  }
  else 
    d[0] = d[1] = d[2] = d[3] = 0.0;

  int* m = new int[2]; //atoms moving onto TM center 
  m[0]=m[1]=-1;
  int* n = new int[4]; //atoms already on TM center 
  n[0]=n[1]=n[2]=n[3]=-1;

  int nfar = 0;
  int nnear = 0;
  if (nplanes)
  for (int i=0;i<4;i++)
  {
    if (d[i]>(getR(tm1)+getR(planes[0][i]))/2)
    {
      m[nfar++] = planes[0][i];
      //printf(" atom pair in plane is far: %i %i \n",tm1,m[nfar-1]);
    }
    else
    {
      n[nnear++] = planes[0][i];
      //printf(" atom pair in plane is near: %i %i \n",tm1,n[nnear-1]);
    }
  }

  //printf(" doing angles near square \n");
  for (int i=0;i<nnear;i++)
  {
    for (int j=0;j<i;j++)
    {
      if (angle_val(n[i],tm1,n[j])>160)
      {
        //printf(" found linear bond in already connected plane: %i %i %i \n",n[i],tm1,n[j]);
#if 0
        tmnoangles[ntmnoangles][0] = n[i];
        tmnoangles[ntmnoangles][1] = tm1;
        tmnoangles[ntmnoangles][2] = n[j];
        ntmnoangles++;
#else
        tmlinangles[ntmlinangles][0] = n[i];
        tmlinangles[ntmlinangles][1] = tm1;
        tmlinangles[ntmlinangles][2] = n[j];
        ntmlinangles++;
#endif
       //find complementary linear near-far pair
        for (int k1=0;k1<nnear;k1++)
        for (int k2=0;k2<nfar;k2++)
        if (n[k1]!=n[i] && n[k1]!=n[j])
        if (m[k2]!=n[i] && m[k2]!=n[j])
        {
          tmlinangles[ntmlinangles][0] = n[k1];
          tmlinangles[ntmlinangles][1] = tm1;
          tmlinangles[ntmlinangles][2] = m[k2];
          ntmlinangles++;
        }

#if 1
 //CPMZ maybe?
       //find complementary linear far-far pair
        for (int k1=0;k1<nfar;k1++)
        for (int k2=0;k2<nfar;k2++)
        if (m[k1]!=n[i] && m[k1]!=n[j])
        if (m[k2]!=n[i] && m[k2]!=n[j])
        {
          tmlinangles[ntmlinangles][0] = m[k1];
          tmlinangles[ntmlinangles][1] = tm1;
          tmlinangles[ntmlinangles][2] = m[k2];
          ntmlinangles++;
        }
#endif

       //find complementary linear pair
        for (int k1=0;k1<nnear;k1++)
        for (int k2=0;k2<k1;k2++)
        if (n[k1]!=n[i] && n[k1]!=n[j])
        if (n[k2]!=n[i] && n[k2]!=n[j])
        {
#if 0
          tmnoangles[ntmnoangles][0] = n[k1];
          tmnoangles[ntmnoangles][1] = tm1;
          tmnoangles[ntmnoangles][2] = n[k2];
          ntmnoangles++;
#else
          tmlinangles[ntmlinangles][0] = n[k1];
          tmlinangles[ntmlinangles][1] = tm1;
          tmlinangles[ntmlinangles][2] = n[k2];
          ntmlinangles++;
#endif
        }

      }
    }
    //if moved to plane
  }


  //printf(" doing angles above and below \n");
  if (nplanes)
  for (int i=0;i<natoms0;i++)
  {
    if (bond_exists(i,tm2))
    for (int j=0;j<i;j++)
    {
      if (bond_exists(j,tm2))
      {
        //printf(" found above and below pair : %i %i %i \n",i,tm1,j);
#if 0
        tmnoangles[ntmnoangles][0] = i;
        tmnoangles[ntmnoangles][1] = tm1;
        tmnoangles[ntmnoangles][2] = j;
        ntmnoangles++;
#else
        tmlinangles[ntmlinangles][0] = i;
        tmlinangles[ntmlinangles][1] = tm1;
        tmlinangles[ntmlinangles][2] = j;
        ntmlinangles++;
#endif
      }
    }
  }

#if 0
//CPMZ reconsider this
  //printf(" doing angles far from square \n");
  for (int i=0;i<nfar;i++)
  {
    d[0] = distance(m[i],a1); 
    d[1] = distance(m[i],a2); 
    d[2] = distance(m[i],a3); 
    d[3] = distance(m[i],a4); 

    printf(" working on %i with d: %1.1f %1.1f %1.1f %1.1f \n",m[i],d[0],d[1],d[2],d[3]);

    if (d[3]<0.1)
    {
      tmnoangles[ntmnoangles][0] = a4;
      tmnoangles[ntmnoangles][1] = tm1;
      if (d[0]>d[1] && d[0]>d[2])
        tmnoangles[ntmnoangles][2] = a1;
      else if (d[1]>d[0] && d[1]>d[2])
        tmnoangles[ntmnoangles][2] = a2;
      else if (d[2]>d[0] && d[2]>d[1])
        tmnoangles[ntmnoangles][2] = a3;
      ntmnoangles++;
    } // if adding to site a4 
    else if (d[2]<0.1)
    {
      tmnoangles[ntmnoangles][0] = a3;
      tmnoangles[ntmnoangles][1] = tm1;
      if (d[0]>d[1] && d[0]>d[3])
        tmnoangles[ntmnoangles][2] = a1;
      else if (d[1]>d[0] && d[1]>d[3])
        tmnoangles[ntmnoangles][2] = a2;
      else if (d[3]>d[0] && d[3]>d[1])
        tmnoangles[ntmnoangles][2] = a4;
      ntmnoangles++;
    } // if adding to site a3
    else if (d[1]<0.1)
    {
      tmnoangles[ntmnoangles][0] = a2;
      tmnoangles[ntmnoangles][1] = tm1;
      if (d[0]>d[2] && d[0]>d[3])
        tmnoangles[ntmnoangles][2] = a1;
      else if (d[2]>d[0] && d[2]>d[3])
        tmnoangles[ntmnoangles][2] = a3;
      else if (d[3]>d[0] && d[3]>d[2])
        tmnoangles[ntmnoangles][2] = a4;
      ntmnoangles++;
    } // if adding to site a2
    else if (d[0]<0.1)
    {
      tmnoangles[ntmnoangles][0] = a1;
      tmnoangles[ntmnoangles][1] = tm1;
      if (d[1]>d[2] && d[1]>d[3])
        tmnoangles[ntmnoangles][2] = a2;
      else if (d[2]>d[1] && d[2]>d[3])
        tmnoangles[ntmnoangles][2] = a3;
      else if (d[3]>d[1] && d[3]>d[2])
        tmnoangles[ntmnoangles][2] = a4;
      ntmnoangles++;
    } // if adding to site a1

  } //loop i over nfar
#endif

#if 0
  printf(" converting no angles to linear \n");
  for (int i=0;i<ntmnoangles;i++)
  {
    tmlinangles[ntmlinangles][0] = tmnoangles[i][0];
    tmlinangles[ntmlinangles][1] = tmnoangles[i][1];
    tmlinangles[ntmlinangles][2] = tmnoangles[i][2];
    ntmlinangles++;
  }
  ntmnoangles = 0;
#endif
#if 0
  for (int i=0;i<ntmnoangles;i++)
    printf(" found noangle: %i %i %i \n",tmnoangles[i][0],tmnoangles[i][1],tmnoangles[i][2]);
  for (int i=0;i<ntmlinangles;i++)
    printf(" found linangle: %i %i %i \n",tmlinangles[i][0],tmlinangles[i][1],tmlinangles[i][2]);
#endif

//other atoms connected to TM that are far away
  nfar=0;
  m[0]=m[1]=-1;
#if 0
  for (int i=0;i<nbonds;i++)
    if ( distance(bonds[i][0],bonds[i][1]) > (getR(bonds[i][0])+getR(bonds[i][1]))/2 )
    {
      //printf(" atom pair is far (%1.2f %1.2f): %i %i \n",bonds[i][0],bonds[i][1],distance(bonds[i][0],bonds[i][1]),(getR(bonds[i][0])+getR(bonds[i][1]))/2);
//      m[nfar++] = 
    }
    else
      printf(" atom pair is not far (%1.2f %1.2f): %i %i \n",bonds[i][0],bonds[i][1],distance(bonds[i][0],bonds[i][1]),(getR(bonds[i][0])+getR(bonds[i][1]))/2);
#endif



  //printf(" preparing for angle creation \n");
  remove_TM_pairs();
  //printf(" now making tm angles \n");
  //doing this after make_nobond in TM_setup
  //make_angles_tm();


 //for near atoms
#if 0
  for (int i=0;i<nangles;i++)
  {
    if (isTM(angles[i][1]))
      if (50<anglev[i] && anglev[i]<130)
        angled[i]=anglev[i];
  }  
#endif
#if 0
  printf(" angled array \n");
  for (int i=0;i<nangles;i++)
  {
    printf(" %1.1f",angled[i]);
  }
  printf("\n");
#endif

  delete [] d;
  delete [] m;
  delete [] n;
  
  return;
}


void ICoord::remove_TM_pairs() {

  //printf(" in remove_TM_pairs() \n");

  int tm1 = -1;
  int tmx = -1;
  for (int i=0;i<natoms;i++)
    if (oatomtm[i]>-1)
    {
      tm1 = oatomtm[i];
      tmx = i;
    }

  printf(" found TM pair: %i %i \n",tm1,tmx);

  if (tm1>-1)
  {
    for (int i=0;i<nbonds;i++)
    {
      if (bonds[i][0]==tmx)
        bonds[i][0]=tm1;
      else if (bonds[i][1]==tmx)
        bonds[i][1]=tm1;
    }

//CPMZ need to reassign all others as well? probably not
    natoms--;
    //printf(" natoms: %i \n",natoms);
    coord_num();  
    //print_ic();   

    //purge duplicate bonds
    for (int i=0;i<nbonds;i++)
    for (int j=0;j<i;j++)
    {
      int rem1 = -1;
      if (bonds[i][0]==bonds[j][0])
      {
        if (bonds[i][1]==bonds[j][1])
          rem1 = i;
      }
      else if (bonds[i][1]==bonds[j][1])
      {
        if (bonds[i][0]==bonds[j][0])
          rem1 = i;
      }
      else if (bonds[i][1]==bonds[j][0])
      {
        if (bonds[i][0]==bonds[j][1])
          rem1 = i;
      }
      else if (bonds[i][0]==bonds[j][1])
      {
        if (bonds[i][1]==bonds[j][0])
          rem1 = i;
      }
      if (rem1>-1)
      {
        for (int k=rem1;k<nbonds-1;k++) 
        {
          bonds[k][0] = bonds[k+1][0];
          bonds[k][1] = bonds[k+1][1];
        }
        printf(" removing duplicate bond: %i \n",rem1);
        nbonds--;
      }
    } //loop i,j over nbonds

  }
  else printf(" no TM to reassign \n");


  return;
}






//setup prior to bond changes

void ICoord::TM_planes() {

  if (natomstm==0) return;

  //printf("\n in TM_planes \n");
  nplanes = 0;
  for (int i=0;i<natoms;i++)
  {
    int s1,s2,s3;
    int b1,b2,b3;
    int found = 0;
    int a1,a2,a3,a4; 

    if (isTM(i))
    {
      //printf(" need to assign one plane, other atoms are type 2 \n");

      a1=a2=a3=a4=-1;

      s1=s2=s3=-1;
      b1=b2=b3=-1;

//CPMZ check this for "clockwise"
      for (int j1=0;j1<nangles;j1++)
      {
        if (found) break;
        if (angles[j1][1]==i)
        {
          a1 = angles[j1][0];
          a2 = angles[j1][2];
          s1 = j1;

          for (int j2=0;j2<nangles;j2++)
          {
            if (found) break;
            if (angles[j2][1]==i && (angles[j2][0]==a2 || angles[j2][2]==a2) && s1!=j2 && a3==-1)
            {
              if (angles[j2][0]==a2)
                a3 = angles[j2][2];
              else
                a3 = angles[j2][0];
              s1=j2;
              for (int j3=0;j3<nangles;j3++)
              {
                if (found) break;
                if (a1<0 || a2<0 || a3<0)
                  break;
                if (angles[j3][1]==i && ((angles[j3][0]==a3 && angles[j3][2]!=a2) || (angles[j3][2]==a3 && angles[j3][0]!=a2))
                   && s1!=j3 && a4==-1)
                {
                  if (angles[j3][0]==a3)
                    a4 = angles[j3][2];
                  else
                    a4 = angles[j3][0];
                  s1=j3;

                  //printf(" possible circle in plane: %i %i %i %i \n",a1,a2,a3,a4);
                  int t1,t2,t3,t4; 
                  t1=t2=t3=t4=0;
                  t1 = isImptor(i,a1,a2,a3);
                  t2 = isImptor(i,a1,a3,a4);
                  t3 = isImptor(i,a2,a3,a4);
                  //printf(" t1 t2 t3: %i %i %i \n",t1,t2,t3);
                  if (abs(t1)+abs(t2)+abs(t3)>2)
                  {
                    //printf(" found circle in plane: %i %i %i %i \n",a1,a2,a3,a4);
        
                    planes[nplanes][0] = a1;
                    planes[nplanes][1] = a2;
                    planes[nplanes][2] = a3;
                    planes[nplanes][3] = a4;
                    nplanes++;
                    found = 1;
                    break;
                  } // if plane found

                  a4 = -1;
                } // found angle 3
              } // loop over angle 3
              a3 = -1;
            } // found angle 2
          } // loop over angle 2
          a2 = -1;
        } // found angle 1
      } // loop over angle 1
    } // assign plane to TM



    if (isTM(i) && nplanes==0)
    {
      //printf(" haven't found any 4-planes, looking for 3-planes \n");
      int tm1 = i;
      int nc = 0;
      int* a = new int[8];
      a[0]=a[1]=a[2]=a[3]=a[4]=a[5]=a[6]=a[7]=-1;
      for (int j=0;j<nbonds;j++)
      {
        if (bonds[j][0]==tm1)
        {
          a[nc++] = bonds[j][1];
          //printf(" tm1 bonded to: %i \n",bonds[j][1]); 
        }
        else if (bonds[j][1]==tm1)
        {
          a[nc++] = bonds[j][0];
          //printf(" tm1 bonded to: %i \n",bonds[j][0]); 
        }

      }
      //printf(" found %i connections \n",nc);

      //printf(" connections: %i %i %i %i \n",a[0],a[1],a[2],a[3]);
      if (nc==4)
      {
        int* t = new int[8];
        t[0] = isImptor(i,a[0],a[1],a[2]);
        t[1] = isImptor(i,a[0],a[2],a[3]);
        t[2] = isImptor(i,a[1],a[2],a[3]);
        t[3] = isImptor(i,a[1],a[3],a[0]);
        t[4] = isImptor(i,a[2],a[3],a[0]);
        t[5] = isImptor(i,a[2],a[0],a[1]);
        t[6] = isImptor(i,a[3],a[0],a[1]);
        t[7] = isImptor(i,a[3],a[1],a[2]);
        //printf(" t1 t2 t3 t4 t5 t6 t7 t8: %i %i %i %i %i %i %i %i \n",t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7]);

        planes[nplanes][3] = -1;
        if (t[0]==1)       { planes[nplanes][0]=a[0]; planes[nplanes][1]=a[1]; planes[nplanes][2]=a[2]; }
        else if (t[0]==-1) { planes[nplanes][0]=a[0]; planes[nplanes][1]=a[2]; planes[nplanes][2]=a[1]; }
        else if (t[1]==1)  { planes[nplanes][0]=a[0]; planes[nplanes][1]=a[2]; planes[nplanes][2]=a[3]; }
        else if (t[1]==-1) { planes[nplanes][0]=a[0]; planes[nplanes][1]=a[3]; planes[nplanes][2]=a[2]; }
        else if (t[2]==1)  { planes[nplanes][0]=a[1]; planes[nplanes][1]=a[2]; planes[nplanes][2]=a[3]; }
        else if (t[2]==-1) { planes[nplanes][0]=a[1]; planes[nplanes][1]=a[3]; planes[nplanes][2]=a[2]; }
        else if (t[3]==1)  { planes[nplanes][0]=a[1]; planes[nplanes][1]=a[3]; planes[nplanes][2]=a[0]; }
        else if (t[3]==-1) { planes[nplanes][0]=a[1]; planes[nplanes][1]=a[0]; planes[nplanes][2]=a[3]; }
        else if (t[4]==1)  { planes[nplanes][0]=a[2]; planes[nplanes][1]=a[3]; planes[nplanes][2]=a[3]; }
        else if (t[4]==-1) { planes[nplanes][0]=a[2]; planes[nplanes][1]=a[0]; planes[nplanes][2]=a[0]; }
        else if (t[5]==1)  { planes[nplanes][0]=a[2]; planes[nplanes][1]=a[0]; planes[nplanes][2]=a[1]; }
        else if (t[5]==-1) { planes[nplanes][0]=a[2]; planes[nplanes][1]=a[1]; planes[nplanes][2]=a[0]; }
        else if (t[6]==1)  { planes[nplanes][0]=a[3]; planes[nplanes][1]=a[0]; planes[nplanes][2]=a[1]; }
        else if (t[6]==-1) { planes[nplanes][0]=a[3]; planes[nplanes][1]=a[1]; planes[nplanes][2]=a[0]; }
        else if (t[7]==1)  { planes[nplanes][0]=a[3]; planes[nplanes][1]=a[1]; planes[nplanes][2]=a[2]; }
        else if (t[7]==-1) { planes[nplanes][0]=a[3]; planes[nplanes][1]=a[2]; planes[nplanes][2]=a[1]; }
        for (int j=0;j<8;j++)
          if (t[j]!=0)
          {
            nplanes++; 
            printf(" 3-plane located (in 4 coord): %i %i %i %i \n",planes[0][0],planes[0][1],planes[0][2],planes[0][3]);
            break;
          }

        delete [] t;
      }
      else if (nc==3)
      {
        int t1 = 0;
        t1 = isImptor(i,a[0],a[1],a[2]);
        printf(" t1: %i for %i %i %i \n",t1,a[0],a[1],a[2]);

        planes[nplanes][3]=-1;
        if (t1==1)  { planes[nplanes][0]=a[0]; planes[nplanes][1]=a[1]; planes[nplanes][2]=a[2]; }
        else if (t1==-1) { planes[nplanes][0]=a[0]; planes[nplanes][1]=a[2]; planes[nplanes][2]=a[1]; }
        if (abs(t1))
        {
          nplanes++;
          printf(" 3-plane located (in 3 coord): %i %i %i %i \n",planes[0][0],planes[0][1],planes[0][2],planes[0][3]);
        }
      }

      delete [] a;
    } // if no planes (i.e. assign 3-plane)




    if (isTM(i) && coordn[i]>3)
    {
      //printf(" found TM %i with coordn (>4): %i \n",i,coordn[i]);

      int shift = -1;
      int* planenm = new int[nbonds];
      for (int j=0;j<nbonds;j++) planenm[j] = 0;

      for (int j=0;j<nbonds;j++)
        for (int k=0;k<4;k++)
          if ( (bonds[j][0]==planes[0][k] || bonds[j][1]==planes[0][k]) 
            && (isTM(bonds[j][0]) || isTM(bonds[j][1])) )
          {
            planenm[j] = 1;
           // printf(" bond %i in plane: %i %i \n",j,bonds[j][0],bonds[j][1]);
          }

      for (int j=0;j<natoms;j++)
        if (oatomtm[j]==i)
          shift = j;
     // printf(" moving out of plane bonds to %i \n",shift);
      for (int j=0;j<nbonds;j++)
        if (!planenm[j])
        {
         // printf(" bond not in plane: %i %i \n",bonds[j][0],bonds[j][1]);
          if (bonds[j][0]==i)
            bonds[j][0]=shift;
          else if (bonds[j][1]==i)
            bonds[j][1]=shift;
        }

      delete [] planenm;

    } // if coordn>4, moved extra atoms to copy TM


    if (isTM(i) && coordn[i]==4)
    {
//      printf(" found TM %i with coordn (=4): %i \n",i,coordn[i]);
      if (nplanes)
        printf(" already assigned one plane: %i %i %i %i \n",planes[0][0],planes[0][1],planes[0][2],planes[0][3]);

    } // coord # 4

    if (isTM(i) && coordn[i]==3)
    {
//      printf(" found TM %i with coord # 3  \n",i);
      if (nplanes)
        printf(" already assigned one plane: %i %i %i %i \n",planes[0][0],planes[0][1],planes[0][2],planes[0][3]);

    }



  } //loop i over natoms

  //printf(" recalculating coord_num  \n");
  coord_num();
  //print_ic();

//  for (int i=0;i<4;i++)
  //creating imptors for plane
  
  //printf(" setting up imptors for plane \n");
  //printf(" nimptors: %i \n",nimptor);
  //printf(" natoms0: %i \n",natoms0);
  for (int i=0;i<natoms0;i++)
  {
    int a1,a2,a3,a4;

    a1 = planes[0][0];
    a2 = planes[0][1];
    a3 = planes[0][2];
    a4 = planes[0][3];

    if (isTM(i) && nplanes)
    {  
      int t1,t2,t3,t4;
      t1=t2=t3=t4=0; 
      t1 = isImptor(i,a1,a2,a3);
      t2 = isImptor(i,a1,a4,a3);
      t3 = isImptor(i,a2,a3,a4);
      t4 = isImptor(i,a3,a4,a1);
      //printf(" t1 t2 t3 t4: %i %i %i %i \n",t1,t2,t3,t4);
 
      imptor[nimptor][2]=i;
      imptor[nimptor][0]=a1;
      if (t1>0)
      {
        imptor[nimptor][1]=a2;
        imptor[nimptor][3]=a3;
      }
      else
      {
        imptor[nimptor][1]=a3;
        imptor[nimptor][3]=a2;
      }
      nimptor++;

      if (a4!=-1)
      {
        imptor[nimptor][2]=i;
        imptor[nimptor][0]=a3;
        if (t4==1)
        {
          imptor[nimptor][1]=a4;
          imptor[nimptor][3]=a1;
        }
        else if (t4==-1)
        {
          imptor[nimptor][1]=a1;
          imptor[nimptor][3]=a4;
        }
        if (t4)
          nimptor++;

        imptor[nimptor][2]=i;
        imptor[nimptor][0]=a1;
        if (t2==1)
        {
          imptor[nimptor][1]=a4;
          imptor[nimptor][3]=a3;
        }
        else if (t2==-1)
        {
          imptor[nimptor][1]=a3;
          imptor[nimptor][3]=a4;
        }
        if (t2)
          nimptor++; 
      } // if a4!=1

    }
  } //loop i over natoms
  

  return;
}





void ICoord::TM_planes_2() {

  //printf("\n in TM_planes_2, nplanes: %i \n",nplanes);


//1. loop through atoms in plane, see if any removed
//2. loop through atoms in plane, see if any added
//3. add/replace etc
//4. set up imptors for plane

  int nintact = 0;
  int* intact = new int[4];
  int nremoved = 0;
  int* removed = new int[4];
  for (int i=0;i<natoms0;i++)
  {
    if (isTM(i))
    for (int j=0;j<4;j++)
    {
      if (bond_exists(i,planes[0][j]))
      {
        intact[nintact++] = planes[0][j];
        //printf(" found pre-existing bond %i %i \n",i,planes[0][j]);
      }
      else if (planes[0][j]!=-1)
      {
        removed[nremoved++] = planes[0][j];
        //printf(" found destroyed bond %i %i \n",i,planes[0][j]);
      } 

    } //loop j over plane atoms
  } //loop i over natoms

  //printf(" found %i intact bonds \n",nintact);
  //printf(" found %i destroyed bonds \n",nremoved);

  int nformed = 0;
  int* formed = new int[4];
  for (int i=0;i<natoms0;i++)
  {
    if (isTM(i))
    for (int j=0;j<natoms;j++)
    if (i!=j)
    {
      if (bond_exists(i,j))
      {
        int found = 0;
        for (int k=0;k<nintact;k++)
          if (intact[k]==j)
            found = 1;
        if (!found)
        {
          formed[nformed++] = j;
          printf(" found new bond %i %i \n",i,j);
        }
      }

    } //loop j over all atoms
  } //loop i over natoms

  for (int i=0;i<nformed;i++)
    if (tmtriangles[i][7]==formed[i] || tmtriangles[i][8]==formed[i])
    {
      for (int j=i;j<nformed-1;j++)
        formed[j] = formed[j+1];
      nformed--;
      printf(" skip adding %i in plane, was tmtri \n",formed[i]);
    }

  //printf(" found %i new bonds \n",nformed);
 


//  for (int i=0;i<nformed;i++)
//    if (oatomtm[formed[i]]>-1)
//      formed[i]

  int nplaced = 0;
  if (nformed>0)
  for (int i=0;i<natoms0;i++)
  {
    if (isTM(i) && nplaced<nremoved && nplaced<nformed)
    for (int j=0;j<4;j++)
    {
//CPMZ will need to place in nearest site
      if (planes[0][j]==removed[nplaced])
      {
        //printf(" replace %i in plane with %i \n",planes[0][j],formed[nplaced]);
        planes[0][j] = formed[nplaced++];
        break;
      }
    } //loop j over in plane

    if (isTM(i) && nplaced>=nremoved && nplaced<nformed)
    for (int j=0;j<4;j++)
    {
      if (planes[0][j]==-1)
      {
        //printf(" place atom %i in plane \n",formed[nplaced]);
       	planes[0][j] = formed[nplaced++];
        break;
      }
    } //loop j over in plane
  } //loop i over natoms

  //remaining plane spots should be empty
  for (int i=0;i<4;i++)
    for (int j=0;j<nremoved;j++)
      if (planes[0][i]==removed[j])
        planes[0][i] = -1;


  //printf(" setting up imptors for plane (in TM_plane_2) \n");
  nimptor = 0;
  for (int i=0;i<natoms0;i++)
  {
    int a1,a2,a3,a4;

    a1 = planes[0][0];
    a2 = planes[0][1];
    a3 = planes[0][2];
    a4 = planes[0][3];
  
#if 0
    if (oatomtm[a1]>-1) a1 = oatomtm[a1];
    if (oatomtm[a2]>-1) a2 = oatomtm[a2];
    if (oatomtm[a3]>-1) a3 = oatomtm[a3];
    if (oatomtm[a4]>-1) a4 = oatomtm[a4];
#endif

    if (isTM(i))
    {  
      if (a1!=-1 && a2!=-1 && a3!=-1)
      {
        imptor[nimptor][2]=i;
        imptor[nimptor][0]=a1;
        imptor[nimptor][1]=a2;
        imptor[nimptor][3]=a3;
        imptorv[nimptor] = torsion_val(a1,a2,i,a3);
        nimptor++;
      }

      if (a1!=-1 && a3!=-1 && a4!=-1)
      {
        imptor[nimptor][2]=i;
        imptor[nimptor][0]=a3;
        imptor[nimptor][1]=a4;
        imptor[nimptor][3]=a1;
        imptorv[nimptor] = torsion_val(a3,a4,i,a1);
        nimptor++;
      } // if a4!=1

      if (a2!=-1 && a3!=-1 && a4!=-1)
      {
        imptor[nimptor][2]=i;
        imptor[nimptor][0]=a2;
        imptor[nimptor][1]=a3;
        imptor[nimptor][3]=a4;
        imptorv[nimptor] = torsion_val(a2,a3,i,a4);
        nimptor++;
      } // if a4!=1

    }
  } //loop i over natoms

//  for (int i=0;i<nimptor;i++)
//    printf(" after tm add, imptor %i: %i %i %i %i: %1.1f \n",i+1,imptor[i][0],imptor[i][1],imptor[i][2],imptor[i][3],imptorv[i]);

  delete [] intact;
  delete [] removed;
  delete [] formed;

  return;
}

