#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald_TMA.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "comm_lattice.h" 

using namespace SPPARKS_NS;

enum{VACANCY,O,OH,Ala,OHAlaX3,OAlaX2,OAlaX2H2O,OAlaXOH,OAlaX,OAlaOH,OAlaOH2,AlaOH,AlaOH2,Alb,OHAlbX3,OAlbX2,OAlbX2H2O,OAlbXOH,OAlbX,OAlbOH,OAlbOH2,AlbOH,AlbOH2,OAla,OAlb,H2O};


#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppAldTMA::AppAldTMA(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  cycle = 0;
  pressureOn = 1;
  hello = 1;
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = nfour = 0;
  srate = drate = vrate = frate = NULL;
  spropensity = dpropensity = vpropensity = fpropensity =NULL;
  sinput = soutput = NULL;
  dinput = doutput = NULL;
  vinput = voutput = NULL;
  finput = foutput = NULL;
  comneigh = NULL;
  scount = dcount = vcount = fcount = NULL;
  sA = dA = vA = fA = NULL;
  scoord = dcoord = vcoord = fcoord = NULL;
  dcoord2 = NULL;
  vcoord2 = NULL; 
  fcoord2 = NULL;
  sexpon = dexpon = vexpon = fexpon = NULL;
  spresson = dpresson = vpresson = fpresson = NULL;
  same_z_neighbors = NULL;
  num_same_z_neighbors = NULL;
  z_tolerance = 1e-6;  // 设置默认容差值
  max_same_z_neighbors = 6;  // 每个位点最多保存6个同z平面近邻
}

/* ---------------------------------------------------------------------- */

AppAldTMA::~AppAldTMA()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(vrate);
  memory->sfree(frate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(vpropensity);
  memory->sfree(fpropensity);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->sfree(dinput);
  memory->sfree(doutput);
  memory->sfree(vinput);
  memory->sfree(voutput);
  memory->sfree(finput);
  memory->sfree(foutput);
  memory->sfree(comneigh);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(vcount);
  memory->sfree(fcount);
  memory->sfree(sA);
  memory->sfree(dA);
  memory->sfree(vA);
  memory->sfree(fA);
  memory->sfree(scoord);
  memory->sfree(dcoord);
  memory->sfree(vcoord);
  memory->sfree(fcoord);
  memory->sfree(sexpon);
  memory->sfree(dexpon);
  memory->sfree(vexpon);
  memory->sfree(fexpon);
  memory->sfree(spresson);
  memory->sfree(dpresson);
  memory->sfree(vpresson);
  memory->sfree(fpresson);
  memory->sfree(dcoord2);
  memory->sfree(vcoord2);
  memory->sfree(fcoord2);
  if (same_z_neighbors) {
    for (int i = 0; i < nlocal; i++)
      delete [] same_z_neighbors[i];
    delete [] same_z_neighbors;
  }
  delete [] num_same_z_neighbors;
}
/* ----------------------------------------------------------------------
   检查coord是否匹配，支持"all"参数
------------------------------------------------------------------------- */
bool AppAldTMA::coord_matches(int site_coord, int required_coord)
{
  // 如果required_coord是COORD_ALL，则匹配任何coord
  if (required_coord == COORD_ALL) return true;
  
  // 否则进行精确匹配
  return (site_coord == required_coord);
}
/* ----------------------------------------------------------------------
   解析coord字符串，支持数字和"all"关键字
   返回：数字值或COORD_ALL(-999)
------------------------------------------------------------------------- */
int AppAldTMA::parse_coord_value(const char* coord_str) {
  if (strcmp(coord_str, "all") == 0 || strcmp(coord_str, "ALL") == 0) {
    return COORD_ALL;
  } else {
    return atoi(coord_str);
  }
}
/* ---------------------------------------------------------------------- */
int AppAldTMA::species_to_enum(const char* species_name) {
  if (strcmp(species_name,"VAC") == 0) return VACANCY;
  else if (strcmp(species_name,"O") == 0) return O;
  else if (strcmp(species_name,"OH") == 0) return OH;
  else if (strcmp(species_name,"Ala") == 0) return Ala;
  else if (strcmp(species_name,"OHAlaX3") == 0) return OHAlaX3;
  else if (strcmp(species_name,"OAlaX2") == 0) return OAlaX2;
  else if (strcmp(species_name,"OAlaX2H2O") == 0) return OAlaX2H2O;
  else if (strcmp(species_name,"OAlaXOH") == 0) return OAlaXOH;
  else if (strcmp(species_name,"OAlaX") == 0) return OAlaX;
  else if (strcmp(species_name,"OAlaOH") == 0) return OAlaOH;
  else if (strcmp(species_name,"OAlaOH2") == 0) return OAlaOH2;
  else if (strcmp(species_name,"AlaOH") == 0) return AlaOH;
  else if (strcmp(species_name,"AlaOH2") == 0) return AlaOH2;
  else if (strcmp(species_name,"Alb") == 0) return Alb;
  else if (strcmp(species_name,"OHAlbX3") == 0) return OHAlbX3;
  else if (strcmp(species_name,"OAlbX2") == 0) return OAlbX2;
  else if (strcmp(species_name,"OAlbX2H2O") == 0) return OAlbX2H2O;
  else if (strcmp(species_name,"OAlbXOH") == 0) return OAlbXOH;
  else if (strcmp(species_name,"OAlbX") == 0) return OAlbX;
  else if (strcmp(species_name,"OAlbOH") == 0) return OAlbOH;
  else if (strcmp(species_name,"OAlbOH2") == 0) return OAlbOH2;
  else if (strcmp(species_name,"AlbOH") == 0) return AlbOH;
  else if (strcmp(species_name,"AlbOH2") == 0) return AlbOH2;
  else if (strcmp(species_name,"OAla") == 0) return OAla;
  else if (strcmp(species_name,"OAlb") == 0) return OAlb;
  else if (strcmp(species_name,"H2O") == 0) return H2O;

  return -1; // 返回-1表示未找到匹配的物种
}

void AppAldTMA::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      // Type I事件：支持coord="all"
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
      
      int input_species = species_to_enum(arg[1]);
      if (input_species < 0) error->all(FLERR,"Illegal event command");
      sinput[none] = input_species;

      int output_species = species_to_enum(arg[2]);
      if (output_species < 0) error->all(FLERR,"Illegal event command");
      soutput[none] = output_species;

      sA[none] = atof(arg[3]);
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = parse_coord_value(arg[6]);  // 使用新的解析函数
      spresson[none] = atoi(arg[7]);
      none++;
      
    } else if (rstyle == 2) {
      // Type II事件：支持两个coord参数，都可以是"all"
      if (narg != 12) error->all(FLERR,"Illegal event command - Type II needs 12 args");

      int input_species1 = species_to_enum(arg[1]);
      if (input_species1 < 0) error->all(FLERR,"Illegal event command");
      dinput[ntwo][0] = input_species1;

      int output_species1 = species_to_enum(arg[2]);
      if (output_species1 < 0) error->all(FLERR,"Illegal event command2");
      doutput[ntwo][0] = output_species1;

      int input_species2 = species_to_enum(arg[3]);
      if (input_species2 < 0) error->all(FLERR,"Illegal event command2");
      dinput[ntwo][1] = input_species2;

      int output_species2 = species_to_enum(arg[4]);
      if (output_species2 < 0) error->all(FLERR,"Illegal event command2");
      doutput[ntwo][1] = output_species2;

      dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = parse_coord_value(arg[8]);   // 使用新的解析函数
      dcoord2[ntwo] = parse_coord_value(arg[9]);  // 使用新的解析函数
      dpresson[ntwo] = atoi(arg[10]);
      ntwo++;
      
    } else if (rstyle == 3) {
      // Type III事件：支持两个coord参数，都可以是"all"
      if (narg != 12) error->all(FLERR,"Illegal event command31 - Type III needs 12 args");

      int input_species1 = species_to_enum(arg[1]);
      if (input_species1 < 0) error->all(FLERR,"Illegal event command32");
      vinput[nthree][0] = input_species1;

      int output_species1 = species_to_enum(arg[2]);
      if (output_species1 < 0) error->all(FLERR,"Illegal event command33");
      voutput[nthree][0] = output_species1;

      int input_species2 = species_to_enum(arg[3]);
      if (input_species2 < 0) error->all(FLERR,"Illegal event command34");
      vinput[nthree][1] = input_species2;

      int output_species2 = species_to_enum(arg[4]);
      if (output_species2 < 0) error->all(FLERR,"Illegal event command35");
      voutput[nthree][1] = output_species2;

      vA[nthree] = atof(arg[5]);
      vexpon[nthree] = atoi(arg[6]);
      vrate[nthree] = atof(arg[7]);
      vcoord[nthree] = parse_coord_value(arg[8]);   // 使用新的解析函数
      vcoord2[nthree] = parse_coord_value(arg[9]);  // 使用新的解析函数
      vpresson[nthree] = atoi(arg[10]);
      nthree++;
      
    } else if (rstyle == 4) {
      // Type IV事件：支持两个coord参数，都可以是"all"
      if (narg != 12) error->all(FLERR,"Illegal event command - Type IV needs 12 args");

      int input_species1 = species_to_enum(arg[1]);
      if (input_species1 < 0) error->all(FLERR,"Illegal event command");
      finput[nfour][0] = input_species1;

      int output_species1 = species_to_enum(arg[2]);
      if (output_species1 < 0) error->all(FLERR,"Illegal event command");
      foutput[nfour][0] = output_species1;

      int input_species2 = species_to_enum(arg[3]);
      if (input_species2 < 0) error->all(FLERR,"Illegal event command");
      finput[nfour][1] = input_species2;

      int output_species2 = species_to_enum(arg[4]);
      if (output_species2 < 0) error->all(FLERR,"Illegal event command");
      foutput[nfour][1] = output_species2;

      fA[nfour] = atof(arg[5]);
      fexpon[nfour] = atoi(arg[6]);
      frate[nfour] = atof(arg[7]);
      fcoord[nfour] = parse_coord_value(arg[8]);    // 使用新的解析函数
      fcoord2[nfour] = parse_coord_value(arg[9]);   // 使用新的解析函数
      fpresson[nfour] = atoi(arg[10]);
      nfour++;
    }
  }
  else if (strcmp(command,"pulse_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse time");
      T1 = atof(arg[0]);
      T3 = atof(arg[1]);
  }
  else if (strcmp(command,"purge_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal purge time");
      T2 = atof(arg[0]);
      T4 = atof(arg[1]);
  } else error->all(FLERR,"Unrecognized command38");
}
/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAldTMA::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */
void AppAldTMA::precompute_same_z_neighbors()
{
  // 分配内存
  num_same_z_neighbors = new int[nlocal];
  same_z_neighbors = new int*[nlocal];
  
  for (int i = 0; i < nlocal; i++) {
    // 初始化计数器
    num_same_z_neighbors[i] = 0;
    
    // 为每个位点分配内存来存储同z平面近邻
    same_z_neighbors[i] = new int[max_same_z_neighbors];
    
    // 获取位点i的z坐标
    double z_i = xyz[i][2];
    
    // 创建一个临时数组来存储所有与位点i在同一z平面上的位点及其水平距离
    int* temp_neighbors = new int[nlocal];
    double* temp_distances = new double[nlocal];
    int num_temp = 0;
    
    // 遍历所有位点，找出与位点i在同一z平面上的位点
    for (int j = 0; j < nlocal; j++) {
      if (i == j) continue;  // 跳过自身
      
      // 获取位点j的z坐标
      double z_j = xyz[j][2];
      
      // 检查z坐标是否相同（使用容差值）
      if (fabs(z_i - z_j) <= z_tolerance) {
        // 计算水平距离
        double dist = horizontal_distance(i, j);
        
        // 将位点j添加到临时数组中
        temp_neighbors[num_temp] = j;
        temp_distances[num_temp] = dist;
        num_temp++;
      }
    }
    
    // 对临时数组按照水平距离进行排序
    for (int j = 0; j < num_temp; j++) {
      for (int k = j + 1; k < num_temp; k++) {
        if (temp_distances[j] > temp_distances[k]) {
          // 交换距离
          double temp_dist = temp_distances[j];
          temp_distances[j] = temp_distances[k];
          temp_distances[k] = temp_dist;
          
          // 交换位点索引
          int temp_idx = temp_neighbors[j];
          temp_neighbors[j] = temp_neighbors[k];
          temp_neighbors[k] = temp_idx;
        }
      }
    }
    
    // 选择最近的max_same_z_neighbors个位点作为同z平面近邻
    for (int j = 0; j < num_temp && j < max_same_z_neighbors; j++) {
      same_z_neighbors[i][j] = temp_neighbors[j];
      num_same_z_neighbors[i]++;
    }
    
    // 释放临时数组
    delete [] temp_neighbors;
    delete [] temp_distances;
  }
  
  // 输出调试信息
  /*if (comm->me == 0) {
    if (screen) fprintf(screen, "预计算了%d个位点的同z平面近邻\n", nlocal);
    if (logfile) fprintf(logfile, "预计算了%d个位点的同z平面近邻\n", nlocal);
  }*/
}

double AppAldTMA::horizontal_distance(int i, int j)
{
  double dx = xyz[i][0] - xyz[j][0];
  double dy = xyz[i][1] - xyz[j][1];
  return sqrt(dx*dx + dy*dy);
}

void AppAldTMA::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity
    comneigh = memory->grow(comneigh,12*maxneigh,2,"app/ald:comneigh");
    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = (int *) memory->smalloc(12*maxneigh*sizeof(int),"app:esites");
    //esites = new int[12*maxneigh]; 
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAldTMA::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
  

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  precompute_same_z_neighbors();

  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app_ald");
  for (int m = 0; m < none; m++) {
    spropensity[m] = sA[m]*pow(temperature,sexpon[m])*exp(-srate[m]/(temperature*8.617333262e-5));
    scount[m] = 0;
  if (spropensity[m] == 0.0) error->warning(FLERR,"spropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = dA[m]*pow(temperature,dexpon[m])*exp(-drate[m]/(temperature*8.617333262e-5));
    dcount[m] = 0;
  if (dpropensity[m] == 0.0) error->warning(FLERR,"dpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nthree; m++) {
    vpropensity[m] = vA[m]*pow(temperature,vexpon[m])*exp(-vrate[m]/(temperature*8.617333262e-5));
    vcount[m] = 0;
  if (vpropensity[m] == 0.0) error->warning(FLERR,"vpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nfour; m++) {
    fpropensity[m] = fA[m]*pow(temperature,fexpon[m])*exp(-frate[m]/(temperature*8.617333262e-5));
    fcount[m] = 0;
  if (fpropensity[m] == 0.0) error->warning(FLERR,"fpropensity cannot be 0.0 for app_ald");
  }

    output_event_propensities();
  
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAldTMA::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAldTMA::site_propensity(int i)
{
  int j,k,m,l,o,p;

  clear_events(i);
  double proball = 0.0;

  // Type I events (单分子反应)
  for (m = 0; m < none; m++) {
    int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
    
    // 修改coord检查，支持"all"
    if (element[i] == sinput[m] && 
        coord_matches(coord[i], scoord[m]) &&  // 使用新的检查函数
        (spresson[m] == pressureOn || spresson[m] == 0) && 
        (coordi <= numneigh[i])) {
        
      add_event(i,1,m,spropensity[m],-1,-1,-1);
      proball += spropensity[m];
    }
  }

  // Type II events (双分子第二近邻反应)
  int nextneib = 1;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (int kk = 0; kk < numneigh[j]; kk++) {
      k = neighbor[j][kk];
      if (i == k) continue;
      for (m = 0; m < ntwo; m++) {
        int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
        int coordk = coord[k] < 0 ? (coord[k] % 10 + 10) : (coord[k] % 10);
        
        // 修改条件检查：同时检查两个物种的coord，支持"all"
        if ((element[i] == dinput[m][0] && element[k] == dinput[m][1]) && 
            (dpresson[m] == pressureOn || dpresson[m] == 0) && 
            coord_matches(coord[i], dcoord[m]) &&      // 第一个物种coord检查
            coord_matches(coord[k], dcoord2[m]) &&     // 第二个物种coord检查
            (coordi <= numneigh[i]) && (coordk < numneigh[k])) {
            
          comevent = 1;
          for (int ii = 0; ii < nextneib; ii++) {
            if (comneigh[ii][0] == k && comneigh[ii][1] == dpropensity[m]) comevent = 0;
          }
          if (comevent) {	
            add_event(i,2,m,dpropensity[m],-1,k,-1);
            proball += dpropensity[m];
            comneigh[nextneib][0] = k;
            comneigh[nextneib][1] = dpropensity[m];
            nextneib++;
          } 
        }
      }
    }
  }

  // Type III events (双分子第一近邻反应)
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < nthree; m++) {
      int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
      int coordj = coord[j] < 0 ? (coord[j] % 10 + 10) : (coord[j] % 10);
      
      // 修改条件检查：同时检查两个物种的coord，支持"all"
      if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && 
          coord_matches(coord[i], vcoord[m]) &&       // 第一个物种coord检查
          coord_matches(coord[j], vcoord2[m]) &&      // 第二个物种coord检查
          (vpresson[m] == pressureOn || vpresson[m] == 0) && 
          (coordi <= numneigh[i]) && (coordj < numneigh[j])) {
          
        add_event(i,3,m,vpropensity[m],j,-1,-1);
        proball += vpropensity[m];
      }
    }
  }

  // Type IV events (同z平面反应)
  for (int p = 0; p < num_same_z_neighbors[i]; p++) {
    int g = same_z_neighbors[i][p]; 
    for (m = 0; m < nfour; m++) {
      int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
      int coordg = coord[g] < 0 ? (coord[g] % 10 + 10) : (coord[g] % 10);
      
      // 修改条件检查：同时检查两个物种的coord，支持"all"
      if (element[i] == finput[m][0] && element[g] == finput[m][1] && 
          coord_matches(coord[i], fcoord[m]) &&       // 第一个物种coord检查
          coord_matches(coord[g], fcoord2[m]) &&      // 第二个物种coord检查
          (fpresson[m] == pressureOn || fpresson[m] == 0) && 
          (coordi <= numneigh[i]) && (coordg < numneigh[g])) {
          
        add_event(i,4,m,fpropensity[m],-1,-1,g);
        proball += fpropensity[m];
      }
    }
  }

  // Type V event (空事件)
  add_event(i,5,0,0.1,-1,-1,-1);
  proball += 0.1;

  return proball;
}

/* ----------------------------------------------------------------------
   输出所有事件的propensity值到终端和日志文件
------------------------------------------------------------------------- */
void AppAldTMA::output_event_propensities()
{
  // 输出到屏幕
  if (screen) {
    fprintf(screen, "\n=== Event Propensities at Temperature %.2f K ===\n", temperature);
    
    // 输出Type I事件 (s1, s2, ...)
    if (none > 0) {
      fprintf(screen, "Type I Events: ");
      for (int m = 0; m < none; m++) {
        fprintf(screen, "s%d:%.3e ", m+1, spropensity[m]);
        if ((m+1) % 5 == 0 && m != none-1) fprintf(screen, "\n               ");
      }
      fprintf(screen, "\n");
    }
    
    // 输出Type II事件 (d1, d2, ...)
    if (ntwo > 0) {
      fprintf(screen, "Type II Events: ");
      for (int m = 0; m < ntwo; m++) {
        fprintf(screen, "d%d:%.3e ", m+1, dpropensity[m]);
        if ((m+1) % 5 == 0 && m != ntwo-1) fprintf(screen, "\n                ");
      }
      fprintf(screen, "\n");
    }
    
    // 输出Type III事件 (v1, v2, ...)
    if (nthree > 0) {
      fprintf(screen, "Type III Events: ");
      for (int m = 0; m < nthree; m++) {
        fprintf(screen, "v%d:%.3e ", m+1, vpropensity[m]);
        if ((m+1) % 5 == 0 && m != nthree-1) fprintf(screen, "\n                 ");
      }
      fprintf(screen, "\n");
    }
    
    // 输出Type IV事件 (f1, f2, ...)
    if (nfour > 0) {
      fprintf(screen, "Type IV Events: ");
      for (int m = 0; m < nfour; m++) {
        fprintf(screen, "f%d:%.3e ", m+1, fpropensity[m]);
        if ((m+1) % 5 == 0 && m != nfour-1) fprintf(screen, "\n                ");
      }
      fprintf(screen, "\n");
    }
    
    fprintf(screen, "===============================================\n\n");
  }
  
  // 输出到日志文件
  if (logfile) {
    fprintf(logfile, "\n=== Event Propensities at Temperature %.2f K ===\n", temperature);
    
    // 输出Type I事件
    if (none > 0) {
      fprintf(logfile, "Type I Events: ");
      for (int m = 0; m < none; m++) {
        fprintf(logfile, "s%d:%.3e ", m+1, spropensity[m]);
        if ((m+1) % 5 == 0 && m != none-1) fprintf(logfile, "\n               ");
      }
      fprintf(logfile, "\n");
    }
    
    // 输出Type II事件
    if (ntwo > 0) {
      fprintf(logfile, "Type II Events: ");
      for (int m = 0; m < ntwo; m++) {
        fprintf(logfile, "d%d:%.3e ", m+1, dpropensity[m]);
        if ((m+1) % 5 == 0 && m != ntwo-1) fprintf(logfile, "\n                ");
      }
      fprintf(logfile, "\n");
    }
    
    // 输出Type III事件
    if (nthree > 0) {
      fprintf(logfile, "Type III Events: ");
      for (int m = 0; m < nthree; m++) {
        fprintf(logfile, "v%d:%.3e ", m+1, vpropensity[m]);
        if ((m+1) % 5 == 0 && m != nthree-1) fprintf(logfile, "\n                 ");
      }
      fprintf(logfile, "\n");
    }
    
    // 输出Type IV事件
    if (nfour > 0) {
      fprintf(logfile, "Type IV Events: ");
      for (int m = 0; m < nfour; m++) {
        fprintf(logfile, "f%d:%.3e ", m+1, fpropensity[m]);
        if ((m+1) % 5 == 0 && m != nfour-1) fprintf(logfile, "\n                ");
      }
      fprintf(logfile, "\n");
    }
    
    fprintf(logfile, "===============================================\n\n");
    fflush(logfile);  // 确保立即写入文件
  }
}
/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppAldTMA::site_event(int i, class RandomPark *random)
{
  int j,k,g,m,n,mm,jj;
  //int elcoord = element[i];
  
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;
  g = events[ievent].gpartner;
  int elcoord_i = element[i];
  int elcoord_j = element[j];
  int elcoord_k = element[k];
  int elcoord_g = element[g];



  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
    } 
  else if (rstyle == 2 && j == -1) {
    element[i] = doutput[which][0];
    element[k] = doutput[which][1];
    dcount[which]++;
    }
  else if (rstyle == 3 && k == -1) {
    element[i] = voutput[which][0];
    element[j] = voutput[which][1];
    vcount[which]++;
    }
  else if (rstyle == 4) {
    element[i] = foutput[which][0];
    element[g] = foutput[which][1];
    fcount[which]++;
  }
  
  else if (rstyle == 5) {}
  else { error->all(FLERR,"Illegal execution event"); }

  //update_coord(elcoord,i,j,k,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 oxygen pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;}
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }

  

  int nsites = 0;
  int isite = i2site[i];


// mask 
   if (rstyle == 1) {
    if ((elcoord_i == OHAlaX3 || elcoord_i == OHAlbX3) && element[i] == OH) {
      remove_mask(i);
      remove_mask_2(i);
    }
    else if ((elcoord_i == OAlaX2H2O || OAlbX2H2O) && (element[i] == OAlaXOH || OAlbXOH)) {
      remove_mask(i);
      remove_mask_2(i);
    }
  }
  
  else if (rstyle == 3) {
    if (elcoord_i == OH && (element[i] == OHAlaX3 || element[i] == OHAlbX3 )) {
      put_mask(i);
      put_mask_2(i);
     }
     
  }

  else if (rstyle == 4) {
    if ((elcoord_g == OAlaX2 || OAlbX2) && (element[g] == OAlaX || OAlbX)) {
      remove_mask(g);
      remove_mask_2(g);
    }
    }
  
  count_coord(i);
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;
for (int n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
    count_coord(m);
    propensity[isite] = site_propensity(m);
    esites[nsites++] = isite;
    // std::cout<<"nsites:"<<nsites<<std::endl;
    echeck[isite] = 1;
    }
    for (int jj = 0; jj< numneigh[m];jj++) {
    mm = neighbor[m][jj];
    isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
            count_coord(mm);
            propensity[isite] = site_propensity(mm);
            esites[nsites++] = isite;
            // std::cout<<"nsites:"<<nsites<<std::endl;
            echeck[isite] = 1;
        }
        for (int ss = 0; ss< numneigh[mm]; ss++) {
            int s = neighbor[mm][ss];
            isite = i2site[s];
            if (isite >= 0 && echeck[isite] == 0) {
              count_coord(s);
              propensity[isite] = site_propensity(s);
              esites[nsites++] = isite;
              // std::cout<<"nsites:"<<nsites<<std::endl;
              echeck[isite] = 1;
            }
            for (int gg = 0; gg< numneigh[s]; gg++) {
              int g = neighbor[s][gg];
              isite = i2site[g];
              if (isite >= 0 && echeck[isite] == 0) {
                count_coord(g);
                propensity[isite] = site_propensity(g);
                esites[nsites++] = isite;
                // std::cout<<"nsites:"<<nsites<<std::endl;
                echeck[isite] = 1;
              }
            }
        }
    }
  }



  solve->update(nsites,esites,propensity);
   // clear echeck array

  for (m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAldTMA::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAldTMA::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner, int gpartner)
{
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].gpartner = gpartner;
  events[freeevent].propensity = propensity;

  if ( propensity == 0 ) error->all(FLERR,"propensity in add_event wrong app ald");
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single and double
------------------------------------------------------------------------- */

void AppAldTMA::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/ald:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/ald:spropensity");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/ald:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/ald:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/ald:scount");
    sA = (double *) 
      memory->srealloc(sA,n*sizeof(double),"app/ald:sA");
    sexpon = (int *) 
      memory->srealloc(sexpon,n*sizeof(int),"app/ald:sexpon");
    scoord = (int *) 
      memory->srealloc(scoord,n*sizeof(int),"app/ald:scoord");
    spresson = (int *) 
      memory->srealloc(spresson,n*sizeof(int),"app/ald:spresson");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/ald:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/ald:dpropensity");
    dinput = memory->grow(dinput,n,2,"app/ald:dinput");
    doutput = memory->grow(doutput,n,2,"app/ald:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/ald:dcount");
    dA = (double *) 
      memory->srealloc(dA,n*sizeof(double),"app/ald:dA");
    dexpon = (int *) 
      memory->srealloc(dexpon,n*sizeof(int),"app/ald:dexpon");
    dcoord = (int *) 
      memory->srealloc(dcoord,n*sizeof(int),"app/ald:dcoord");
    dpresson = (int *) 
      memory->srealloc(dpresson,n*sizeof(int),"app/ald:dpresson");
    dcoord2 = (int *) 
      memory->srealloc(dcoord2,n*sizeof(int),"app/ald:dcoord2");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    vrate = (double *)
      memory->srealloc(vrate,n*sizeof(double),"app/ald:vrate");
    vpropensity = (double *)
      memory->srealloc(vpropensity,n*sizeof(double),"app/ald:vpropensity");
    vinput = memory->grow(vinput,n,2,"app/ald:vinput");
    voutput = memory->grow(voutput,n,2,"app/ald:voutput");
    vcount = (int *)
      memory->srealloc(vcount,n*sizeof(int),"app/ald:vcount");
    vA = (double *)
      memory->srealloc(vA,n*sizeof(double),"app/ald:vA");
    vexpon = (int *)
      memory->srealloc(vexpon,n*sizeof(int),"app/ald:vexpon");
    vcoord = (int *)
      memory->srealloc(vcoord,n*sizeof(int),"app/ald:vcoord");
    vpresson = (int *)
      memory->srealloc(vpresson,n*sizeof(int),"app/ald:vpresson");
    vcoord2 = (int *)
      memory->srealloc(vcoord2,n*sizeof(int),"app/ald:vcoord2");
     }
  else if (rstyle == 4) {
    int n = nfour + 1;
    frate = (double *)
      memory->srealloc(frate,n*sizeof(double),"app/ald:frate");
    fpropensity = (double *)
      memory->srealloc(fpropensity,n*sizeof(double),"app/ald:fpropensity");
    finput = memory->grow(finput,n,2,"app/ald:finput");
    foutput = memory->grow(foutput,n,2,"app/ald:foutput");
    fcount = (int *)
      memory->srealloc(fcount,n*sizeof(int),"app/ald:fcount");
    fA = (double *)
      memory->srealloc(fA,n*sizeof(double),"app/ald:fA");
    fexpon = (int *)
      memory->srealloc(fexpon,n*sizeof(int),"app/ald:fexpon");
    fcoord = (int *)
      memory->srealloc(fcoord,n*sizeof(int),"app/ald:fcoord");
    fpresson = (int *)
      memory->srealloc(fpresson,n*sizeof(int),"app/ald:fpresson");
    fcoord2 = (int *)
      memory->srealloc(fcoord2,n*sizeof(int),"app/ald:fcoord2");
  }
}

/* ----------------------------------------------------------------------
   count c.n after densification
------------------------------------------------------------------------- */
void AppAldTMA::count_coord(int i) // i: Oxygen species(does not necessarily hold)
{
 if (0 <= coord[i]) coord[i] = 0;
  else if (-10 < coord[i] && coord[i] < 0) coord[i] = -10;
  else if (-20 < coord[i] && coord[i] < -10) coord[i] = -20;
  else if (-30 < coord[i] && coord[i] < -20) coord[i] = -30;
  else if (-40 < coord[i] && coord[i] < -30) coord[i] = -40;
  else if (-50 < coord[i] && coord[i] < -40) coord[i] = -50;

  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    if (element[nn] != VACANCY) coord[i]++;
}
}

/* ---------------------------------------------------- ------------------
   count c.n of oxygen before adsorption
------------------------------------------------------------------------- */
void AppAldTMA::count_coordO(int i)
{
    int fullO = 0;
    int emptyO = 0;
    int totalS = 0;

    int isite = i2site[i];
    int nsites = 0;

	for (int m = 0; m < numneigh[i]; m++) {
		int mm = neighbor[i][m];
		for (int s = 0; s < numneigh[mm]; s++) {
			int ss = neighbor[mm][s];
			isite = i2site[ss];
			if (i==ss)  continue;
			if (isite >= 0 && echeck[isite] == 0) {
			  if ( element[ss] >= O && element[ss] <= AlbOH2 ) {fullO++;}
			  else if (element[ss] == VACANCY) {emptyO++;}
		          esites[nsites++] = isite;
		          echeck[isite] = 1;
		        }
    
		}
	}
   totalS = fullO+emptyO;
   if ( float(fullO) > 4*totalS/5 and coord[i] > -20) {coord[i] += -20;} // decrease the coord of the oxygen site to render it inactive for adsorption
   for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   put mask for affected sites
------------------------------------------------------------------------- */
void AppAldTMA::put_mask(int i) {
  int nsites = 0;
  int isite = i2site[i];
  echeck[isite] = 1;
  esites[nsites++] = isite;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    int isite = i2site[nn];
    // std::cout<<"isite:"<<isite<<std::endl;
    // std::cout<<"echeck:"<<echeck[isite]<<std::endl;
    if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nn] -= 0;
        // std::cout << "put mask " << i << std::endl;
    }
    for (int ss = 0; ss < numneigh[nn]; ss++){
      int nnn = neighbor[nn][ss];
      int isite = i2site[nnn];
      if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nnn] -= 50;
          }
      
  }
  }
  // for (int s = 0; s < numneigh[i]; s++) {
  //   int nn = neighbor[i][s];
  //   coord[nn] -= 10;
  // }
  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  // std::cout << "put mask " << i << std::endl;
}

/* ----------------------------------------------------------------------
   将 mask 放置到同 z 平面的最近 6 个近邻上
   使用预先计算的 same_z_neighbors 数组
------------------------------------------------------------------------- */
void AppAldTMA::put_mask_2(int i) {
  int nsites = 0;
  int isite = i2site[i];
  
  // 标记中心位点
  echeck[isite] = 1;
  esites[nsites++] = isite;
  
  // 遍历预先计算的同 z 平面近邻
  // 这些近邻已经按 xy 平面距离排序
  for (int j = 0; j < num_same_z_neighbors[i]; j++) {
    int nn = same_z_neighbors[i][j];
    isite = i2site[nn];
    
    if (isite >= 0 && echeck[isite] == 0) {
      // 标记位点为已处理
      echeck[isite] = 1;
      esites[nsites++] = isite;
      
      // 对位点的 coord 值减去 10
      coord[nn] -= 10;
    }
  }
  
  // 重置 echeck 数组，为下次计算做准备
  for (int m = 0; m < nsites; m++) {
    echeck[esites[m]] = 0;
    esites[m] = 0;
  }
}

void AppAldTMA::remove_mask(int i) {
  int nsites = 0;
  int isite = i2site[i];
  echeck[isite] = 1;
  esites[nsites++] = isite;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    int isite = i2site[nn];
    // std::cout<<"isite:"<<isite<<std::endl;
    // std::cout<<"echeck:"<<echeck[isite]<<std::endl;
    if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nn] += 0;
        // std::cout << "put mask " << i << std::endl;
    }
    for (int ss = 0; ss < numneigh[nn]; ss++){
      int nnn = neighbor[nn][ss];
      int isite = i2site[nnn];
      if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nnn] += 50;
          }
      
  }
  }
  // for (int s = 0; s < numneigh[i]; s++) {
  //   int nn = neighbor[i][s];
  //   coord[nn] += 10;
  // }
  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  // std::cout << "remove mask " << i << std::endl;
}

void AppAldTMA::remove_mask_2(int i) {
  int nsites = 0;
  int isite = i2site[i];
  
  // 标记中心位点
  echeck[isite] = 1;
  esites[nsites++] = isite;
  
  // 遍历预先计算的同 z 平面近邻
  // 这些近邻已经按 xy 平面距离排序
  for (int j = 0; j < num_same_z_neighbors[i]; j++) {
    int nn = same_z_neighbors[i][j];
    isite = i2site[nn];
    
    if (isite >= 0 && echeck[isite] == 0) {
      // 标记位点为已处理
      echeck[isite] = 1;
      esites[nsites++] = isite;
      
      // 对位点的 coord 值减去 10
      coord[nn] += 10;
    }
  }
  
  // 重置 echeck 数组，为下次计算做准备
  for (int m = 0; m < nsites; m++) {
    echeck[esites[m]] = 0;
    esites[m] = 0;
  }
}