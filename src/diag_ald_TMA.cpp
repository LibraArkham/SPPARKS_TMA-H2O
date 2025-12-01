/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_ald_TMA.h"
#include "app.h"
#include "app_ald_TMA.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;


enum{VACANCY,O,OH,Ala,OHAlaX3,OAlaX2,OAlaX2H2O,OAlaXOH,OAlaXOHH2O,OAlaX,OAlaXH2O,OAlaOH,OAlaOH2,AlaOH,AlaOH2,Alb,OHAlbX3,OAlbX2,OAlbX2H2O,OAlbXOH,OAlbXOHH2O,OAlbX,OAlbXH2O,OAlbOH,OAlbOH2,AlbOH,AlbOH2,H2O,QCM,EVENTS,ONE,TWO,THREE,FOUR,OHG,XG,TYPE1_TOTAL,TYPE2_TOTAL,TYPE3_TOTAL,TYPE4_TOTAL,ALL_EVENTS_TOTAL,OAla,OAlb,ALA_TOTAL,ALB_TOTAL};


/* ---------------------------------------------------------------------- */

DiagAldTMA::DiagAldTMA(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ald/TMA") != 0)
    error->all(FLERR,"Diag_style ald requires app_style ald");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style ald command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style ald command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagAldTMA::~DiagAldTMA()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::init()
{
  appaldTMA = (AppAldTMA *) app;
  
  int none = appaldTMA->none;
  int ntwo = appaldTMA->ntwo;
  int nthree = appaldTMA->nthree;
  int nfour = appaldTMA->nfour;

  for (int i = 0; i < nlist; i++) {
      if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"O") == 0) which[i] = O;
      else if (strcmp(list[i],"OH") == 0) which[i] = OH;
      else if (strcmp(list[i],"Ala") == 0) which[i] = Ala;
      else if (strcmp(list[i],"OHAlaX3") == 0) which[i] = OHAlaX3;
      else if (strcmp(list[i],"OAlaX2") == 0) which[i] = OAlaX2;
      else if (strcmp(list[i],"OAlaX2H2O") == 0) which[i] = OAlaX2H2O;
      else if (strcmp(list[i],"OAlaXOH") == 0) which[i] = OAlaXOH;
      else if (strcmp(list[i],"OAlaX") == 0) which[i] = OAlaX;
      else if (strcmp(list[i],"OAlaOH") == 0) which[i] = OAlaOH;
      else if (strcmp(list[i],"OAlaOH2") == 0) which[i] = OAlaOH2;
      else if (strcmp(list[i],"AlaOH") == 0) which[i] = AlaOH;
      else if (strcmp(list[i],"AlaOH2") == 0) which[i] = AlaOH2;
      else if (strcmp(list[i],"Alb") == 0) which[i] = Alb;
      else if (strcmp(list[i],"OHAlbX3") == 0) which[i] = OHAlbX3;
      else if (strcmp(list[i],"OAlbX2") == 0) which[i] = OAlbX2;
      else if (strcmp(list[i],"OAlbX2H2O") == 0) which[i] = OAlbX2H2O;
      else if (strcmp(list[i],"OAlbXOH") == 0) which[i] = OAlbXOH;
      else if (strcmp(list[i],"OAlbX") == 0) which[i] = OAlbX;
      else if (strcmp(list[i],"OAlbOH") == 0) which[i] = OAlbOH;
      else if (strcmp(list[i],"OAlbOH2") == 0) which[i] = OAlbOH2;
      else if (strcmp(list[i],"AlbOH") == 0) which[i] = AlbOH;
      else if (strcmp(list[i],"AlbOH2") == 0) which[i] = AlbOH2;
      else if (strcmp(list[i],"OAla") == 0) which[i] = OAla;
      else if (strcmp(list[i],"OAlb") == 0) which[i] = OAlb;
      else if (strcmp(list[i],"H2O") == 0) which[i] = H2O;
      else if (strcmp(list[i],"OAlaXOHH2O") == 0) which[i] = OAlaXOHH2O;
      else if (strcmp(list[i],"OAlaXH2O") == 0) which[i] = OAlaXH2O;
      else if (strcmp(list[i],"OAlbXOHH2O") == 0) which[i] = OAlbXOHH2O;
      else if (strcmp(list[i],"OAlbXH2O") == 0) which[i] = OAlbXH2O;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
      else if (strcmp(list[i],"QCM") == 0) which[i] = QCM;
      else if (strcmp(list[i],"OHG") == 0) which[i] = OHG;
      else if (strcmp(list[i],"XG") == 0) which[i] = XG;
      else if (strcmp(list[i],"type1_total") == 0) which[i] = TYPE1_TOTAL;
      else if (strcmp(list[i],"type2_total") == 0) which[i] = TYPE2_TOTAL;
      else if (strcmp(list[i],"type3_total") == 0) which[i] = TYPE3_TOTAL;
      else if (strcmp(list[i],"type4_total") == 0) which[i] = TYPE4_TOTAL;
      else if (strcmp(list[i],"all_events") == 0) which[i] = ALL_EVENTS_TOTAL;
      else if (strcmp(list[i],"ala_total") == 0) which[i] = ALA_TOTAL;
      else if (strcmp(list[i],"alb_total") == 0) which[i] = ALB_TOTAL;
      
      else if (list[i][0] == 'f') {
        which[i] = FOUR;
        int n = atoi(&list[i][1]);
        if (n < 1 || n > nfour)
          error->all(FLERR,"Invalid value setting in diag_style ald");
        index[i] = n - 1;
      }
      else if (list[i][0] == 's') {
        which[i] = ONE;
        int n = atoi(&list[i][1]);
        if (n < 1 || n > none) 
          error->all(FLERR,"Invalid value setting in diag_style ald");
        index[i] = n - 1;
      } else if (list[i][0] == 'd') {
        which[i] = TWO;
        int n = atoi(&list[i][1]);
        if (n < 1 || n > ntwo) 
          error->all(FLERR,"Invalid value setting in diag_style ald");
        index[i] = n - 1;
      } else if (list[i][0] == 'v') {
        which[i] = THREE;
        int n = atoi(&list[i][1]);
        if (n < 1 || n > nthree) 
          error->all(FLERR,"Invalid value setting in diag_style ald");
        index[i] = n - 1;
      } else error->all(FLERR,"Invalid value setting in diag_style ald");
  }

  siteflag = 1; 

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::compute()
{
  int sites[800], ivalue;
  if (siteflag) {
    // 初始化所有物种计数为0（补充缺失的）
    sites[VACANCY] = 0; sites[O] = 0; sites[OH] = 0; sites[Ala] = 0; sites[OHAlaX3] = 0;
    sites[OAlaX2] = 0; sites[OAlaX2H2O] = 0; sites[OAlaXOH] = 0; sites[OAlaXOHH2O] = 0; // 添加
    sites[OAlaX] = 0; sites[OAlaXH2O] = 0; // 添加
    sites[OAlaOH] = 0; sites[OAlaOH2] = 0; sites[AlaOH] = 0; sites[AlaOH2] = 0;
    sites[Alb] = 0; sites[OHAlbX3] = 0; sites[OAlbX2] = 0; sites[OAlbX2H2O] = 0;
    sites[OAlbXOH] = 0; sites[OAlbXOHH2O] = 0; // 添加
    sites[OAlbX] = 0; sites[OAlbXH2O] = 0; // 添加
    sites[OAlbOH] = 0; sites[OAlbOH2] = 0;
    sites[AlbOH] = 0; sites[AlbOH2] = 0; sites[OAla] = 0; sites[OAlb] = 0; sites[H2O] = 0;
    
    int *element = appaldTMA->element;
    int nlocal = appaldTMA->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == O) ivalue = sites[O];
    else if (which[i] == OH) ivalue = sites[OH];
    else if (which[i] == Ala) ivalue = sites[Ala];
    else if (which[i] == OHAlaX3) ivalue = sites[OHAlaX3];
    else if (which[i] == OAlaX2) ivalue = sites[OAlaX2];
    else if (which[i] == OAlaX2H2O) ivalue = sites[OAlaX2H2O];
    else if (which[i] == OAlaXOH) ivalue = sites[OAlaXOH];
    else if (which[i] == OAlaX) ivalue = sites[OAlaX];
    else if (which[i] == OAlaOH) ivalue = sites[OAlaOH];
    else if (which[i] == OAlaOH2) ivalue = sites[OAlaOH2];
    else if (which[i] == AlaOH) ivalue = sites[AlaOH];
    else if (which[i] == AlaOH2) ivalue = sites[AlaOH2];
    else if (which[i] == Alb) ivalue = sites[Alb];
    else if (which[i] == OHAlbX3) ivalue = sites[OHAlbX3];
    else if (which[i] == OAlbX2) ivalue = sites[OAlbX2];
    else if (which[i] == OAlbX2H2O) ivalue = sites[OAlbX2H2O];
    else if (which[i] == OAlbXOH) ivalue = sites[OAlbXOH];
    else if (which[i] == OAlbX) ivalue = sites[OAlbX];
    else if (which[i] == OAlbOH) ivalue = sites[OAlbOH];
    else if (which[i] == OAlbOH2) ivalue = sites[OAlbOH2];
    else if (which[i] == AlbOH) ivalue = sites[AlbOH];
    else if (which[i] == AlbOH2) ivalue = sites[AlbOH2];
    else if (which[i] == OAla) ivalue = sites[OAla];
    else if (which[i] == OAlb) ivalue = sites[OAlb];
    else if (which[i] == H2O) ivalue = sites[H2O];
    else if (which[i] == OAlaXOHH2O) ivalue = sites[OAlaXOHH2O];
    else if (which[i] == OAlaXH2O) ivalue = sites[OAlaXH2O];
    else if (which[i] == OAlbXOHH2O) ivalue = sites[OAlbXOHH2O];
    else if (which[i] == OAlbXH2O) ivalue = sites[OAlbXH2O];
    else if (which[i] == EVENTS) ivalue = appaldTMA->nevents;
    else if (which[i] == ONE) ivalue = appaldTMA->scount[index[i]];
    else if (which[i] == TWO) ivalue = appaldTMA->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appaldTMA->vcount[index[i]];
    else if (which[i] == FOUR) ivalue = appaldTMA->fcount[index[i]];
   else if (which[i] == QCM) ivalue = 0*sites[VACANCY] + 16*sites[O] + 17*sites[OH] + 27*sites[Ala] + 
                                  89*sites[OHAlaX3] + 73*sites[OAlaX2] + 91*sites[OAlaX2H2O] + 
                                  75*sites[OAlaXOH] + 93*sites[OAlaXOHH2O] + // 添加
                                  58*sites[OAlaX] + 76*sites[OAlaXH2O] + // 添加 
                                  60*sites[OAlaOH] + 77*sites[OAlaOH2] + 44*sites[AlaOH] + 61*sites[AlaOH2] + 
                                  27*sites[Alb] + 89*sites[OHAlbX3] + 73*sites[OAlbX2] + 91*sites[OAlbX2H2O] + 
                                  75*sites[OAlbXOH] + 93*sites[OAlbXOHH2O] + // 添加
                                  58*sites[OAlbX] + 76*sites[OAlbXH2O] + // 添加
                                  60*sites[OAlbOH] + 77*sites[OAlbOH2] + 44*sites[AlbOH] + 61*sites[AlbOH2] + 
                                  43*sites[OAla] + 43*sites[OAlb] + 18*sites[H2O];    else if (which[i] == OHG) ivalue = sites[OH]+sites[OHAlaX3]+sites[OAlaXOH]+sites[OAlaOH]+2*sites[OAlaOH2]+sites[AlaOH]+2*sites[AlaOH2]+sites[OHAlbX3]+sites[OAlbXOH]+sites[OAlbOH]+2*sites[OAlbOH2]+sites[AlbOH]+2*sites[AlbOH2];
    else if (which[i] == XG) ivalue = 3*sites[OHAlaX3]+2*sites[OAlaX2]+2*sites[OAlaX2H2O]+sites[OAlaXOH]+sites[OAlaX]+3*sites[OHAlbX3]+2*sites[OAlbX2]+2*sites[OAlbX2H2O]+sites[OAlbXOH]+sites[OAlbX]+sites[OAlaXOHH2O]+sites[OAlaXH2O]+sites[OAlbXOHH2O]+sites[OAlbXH2O];
    else if (which[i] == TYPE1_TOTAL) {
      // 计算所有Type 1事件的总数
      ivalue = 0;
      for (int j = 0; j < appaldTMA->none; j++) {
        ivalue += appaldTMA->scount[j];
      }
    }
    else if (which[i] == TYPE2_TOTAL) {
      // 计算所有Type 2事件的总数
      ivalue = 0;
      for (int j = 0; j < appaldTMA->ntwo; j++) {
        ivalue += appaldTMA->dcount[j];
      }
    }
    else if (which[i] == TYPE3_TOTAL) {
      // 计算所有Type 3事件的总数
      ivalue = 0;
      for (int j = 0; j < appaldTMA->nthree; j++) {
        ivalue += appaldTMA->vcount[j];
      }
    }
    else if (which[i] == TYPE4_TOTAL) {
      // 计算所有Type 4事件的总数
      ivalue = 0;
      for (int j = 0; j < appaldTMA->nfour; j++) {
        ivalue += appaldTMA->fcount[j];
      }
    }
    else if (which[i] == ALL_EVENTS_TOTAL) {
      // 计算所有类型事件的总数
      ivalue = 0;
      for (int j = 0; j < appaldTMA->none; j++) {
        ivalue += appaldTMA->scount[j];
      }
      for (int j = 0; j < appaldTMA->ntwo; j++) {
        ivalue += appaldTMA->dcount[j];
      }
      for (int j = 0; j < appaldTMA->nthree; j++) {
        ivalue += appaldTMA->vcount[j];
      }
      for (int j = 0; j < appaldTMA->nfour; j++) {
        ivalue += appaldTMA->fcount[j];
      }
    }
    else if (which[i] == ALA_TOTAL) {
  // 计算所有含Ala元素的物种总和
  ivalue = sites[Ala] + sites[OHAlaX3] + sites[OAlaX2] + sites[OAlaX2H2O] + 
           sites[OAlaXOH] + sites[OAlaXOHH2O] + // 添加
           sites[OAlaX] + sites[OAlaXH2O] + // 添加
           sites[OAlaOH] + sites[OAlaOH2] + 
           sites[AlaOH] + sites[AlaOH2] + sites[OAla];
}
else if (which[i] == ALB_TOTAL) {
  // 计算所有含Alb元素的物种总和
  ivalue = sites[Alb] + sites[OHAlbX3] + sites[OAlbX2] + sites[OAlbX2H2O] + 
           sites[OAlbXOH] + sites[OAlbXOHH2O] + // 添加
           sites[OAlbX] + sites[OAlbXH2O] + // 添加
           sites[OAlbOH] + sites[OAlbOH2] + 
           sites[AlbOH] + sites[AlbOH2] + sites[OAlb];
}

    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %6d ",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %6s ",list[i]);
    str += strlen(str);
  }
}