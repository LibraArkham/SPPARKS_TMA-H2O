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
/* ----------------------------------------------------------------------
   ALD application of HfO2 was developed by 
   mahdi shirazi: m.shirazi@tue.nl, TU/e department of applied physics,
   Simon D. Elliott: simon.elliott@schrodinger.com, Schrodinger Materials Science.
   This application is a part of SPPARKS and authors retian the above term.
   See the manual-app-ald and examples folders for more information.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(ald/TMA,AppAldTMA)

#else

#ifndef SPK_APP_ALD_TMA_H
#define SPK_APP_ALD_TMA_H

#include "app_lattice.h"



namespace SPPARKS_NS {

class AppAldTMA : public AppLattice {
  friend class DiagAldTMA;

 public:
  AppAldTMA(class SPPARKS *, int, char **);
  ~AppAldTMA();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);
  int species_to_enum(const char* species_name); // 添加这一行
 private:
  int engstyle;
  int *coord,*element;      // variables on each lattice site 
  int firsttime;
  int hello;
  double T1,T2,T3,T4;          // time period during ALD
  double cycle;
  int pressureOn;
  
  // 添加同z平面近邻相关变量
  int **same_z_neighbors;    // 每个位点的同z平面近邻列表
  int *num_same_z_neighbors; // 每个位点的同z平面近邻数量
  int max_same_z_neighbors;  // 每个位点最多保存的同z平面近邻数量
  double z_tolerance;        // z坐标比较的容差值

  int *esites;
  int *echeck;

  int none,ntwo,nthree,nfour;
  double *srate,*drate,*vrate,*frate;/* two type of reaction, therefore we need only two pointers here, I deleted trate,tcount,toutput */
  double *spropensity,*dpropensity,*vpropensity,*fpropensity;
  /* int *stype,**dtype,**ttype; we do not need any type, we have only one type of crystal that was red by read_sites*/
  int *sinput,**dinput,**vinput,**finput;
  int *soutput,**doutput,**voutput,**foutput;
  int comevent;
  double **comneigh;
  int *scount,*dcount,*vcount,*fcount;
  double *sA,*dA,*vA,*fA;
  int *sexpon,*dexpon,*vexpon,*fexpon;
  int *scoord,*dcoord,*vcoord,*fcoord;//coord options
  int *spresson,*dpresson,*vpresson,*fpresson; //pressure options

  struct Event {           // one event for an owned site
    int style;             // reaction style = SINGLE,DOUBLE,TRIPLE
    int which;             // which reaction of this type
    int jpartner,kpartner,gpartner; // neighbors of site I, it can be first or second 
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, int, double, int, int,int);
  void grow_reactions(int);
  void count_coord(int);
  void count_coordO(int);
  void remove_mask(int);
  void remove_mask_2(int);
  void put_mask(int);
  void put_mask_2(int);
  void update_coord(int,int,int,int,int);
  void output_event_propensities();  // 输出事件propensity的函数
  // 添加新的函数声明
  void precompute_same_z_neighbors();
  double horizontal_distance(int i, int j);
};

}

#endif
#endif
