#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> /* required for randomize() and random() */
#include <time.h>
#include <unistd.h>
#include "nrutil.h"
#include "nr.h"    /* numerical recipes: for GJ elim */
#include "neusglob.h"
#include "aryeh.h"
#include "myCoordServer.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif


#define DEBUG

#ifdef DEBUG
#define COMPILE_STAMP fprintf(stdout, "%s was compiled on %s at %s\n", __FILE__ , __DATE__ , __TIME__ )
#elif
#define COMPILE_STAMP
#endif

#ifdef DEBUG
#define OUTPUT( string ) fprintf(stdout, string );
#elif
#define OUTPUT( string )
#endif

#ifdef DEBUG
#define LINE_STAMP fprintf(stdout, "%d:  line %d of %s  \n" , me,__LINE__ , __FILE__ )
#elif
#define LINE_STAMP
#endif


/* super-global variables */

ServedCoords Coords_weight;
ServedCoords Coords_wadj;
ServedCoords Coords_wavg;
ServedCoords Coords_elaps;
ServedCoords Coords_tau;
ServedCoords Coords_stravg;
ServedCoords Coords_z;
ServedCoords Coords_count;
ServedCoords Coords_narr;
ServedCoords Coords_on;
ServedCoords Coords_navg;
ServedCoords Coords_bdedit;
ServedCoords Coords_ntau;
ServedCoords Coords_prob;
ServedCoords Coords_nprob;
ServedCoords Coords_thist;
ServedCoords Coords_full;
ServedCoords Coords_nlist;
ServedCoords Coords_twlist;
ServedCoords Coords_pts;
ServedCoords Coords_wlist;
ServedCoords Coords_from;
ServedCoords Coords_clock;
ServedCoords Coords_mytim;
ServedCoords Coords_nwstack;

double ***weight, ***wavg, ***tau;
opoint_t **stravg;
opoint_t ***z;
int ***count, ***on, **navg, ***ntau;
int ****thist;
int allind[2];
int stredit=0;  /* stredit = 0 (hasn't been done yet) 
	         stredit = 1 (it's done already)
		 stredit = -1 (it's in progress) */

/* global variables */

double *twavg1, *twavg2, *tweight1, *tweight2, *ttau1, *ttau2;
double *tstravg1, *tstravg2, *tz1, *tz2, *tpt, *lfrac;
double globtim=0.;
int nglobtim=0;
int *tint1, *tint2, *ton1, *ton2, *tprob1, *tprob2, *tnprob1, *tnprob2, *tthist1, *tthist2;
int *tnwstack1, *tnwstack2;
int verb=2;
int me, nproc, status, beadsp2, beadsp3, lockind, lo[2], hi[2], ld[1];
double tstepsq, tempw[1];

char crname[30]="cutRNAin.dat";      // crystal structure file
char secname[30]="cutRNAsec.con";    // secondary contacts
char tertname[30]="cutRNAtert.con";  // tertiary contacts
char emname[30]="emerg.dat";         // emergency points
char sname[30]="isolv";         // emergency points
char filename[30];         // emergency points
double **r0, **r0d6, **r0d12, **nextstep, **r0dists, **rdists, **vsolv, **rsolv, **f1;
int **pos, **pos_poly;
double ***v_cmax, ***v_cmay, ***v_cmaz;
double ***v_cmx, ***v_cmy, ***v_cmz; 
xpoint_t **rvecs;
int ***Delta;
int ***box, ***box_poly;
int ****lab, ***nlab;
double polymass=300.;
double solvmass=18.;
double calpha[2], salpha[2], var, var2, solvstep, phan[NDIM_c], n_av;
double phanr0 = 1.;
double phanr06;
int N=503200, maxd, strmfrq;
int Lx=400, Ly=170, Lz=185;
double blen=5.;              /* length of box in each dir */
double cutoff=3;
int nx,ny,nz;
double blinv;
double kf, R0, ooR02, eph1_12, eph2_12, ep1_6, sigma1, sigma16, sigma2, sigma26;
void wrxyz(int step);
int globn;


const int ndim=NDIM_c, mxlist=MXLIST_c, nop=NOP_c, npart=NPART_c;
int state, bas, mytim, polyerr=0;
opoint_t ocoor;
point_t coor;
int myclock;
point_t ptconv(int lat, int dir, int bead);
opoint_t ptavg(opoint_t x, opoint_t y);
opoint_t projx(point_t x);
int bascheck();
int which(int lat, int dir);
int whichx(int lat, int dir, point_t x);
void gaexit(int sig);
void gtst(int lat, int dir);
void addpoint (int lat, int dir, int bead, point_t x, double wt, int from1, int from2);
void fixwt(int lat, int dir1, int dir2, int reg1, int reg2);
double calcsum(int lat, int dir);
double ogtdist(opoint_t px, opoint_t py);
void repar(int lat, int dir);
void nobad(int lat, int dir, int bead);
void delpoint(int lat, int dir, int bead, int pt);
void wrstring(), getstring();
int getprob(int lat, int dir, int bead, int i, int j, int n1, int n2);
int getnprob(int lat, int dir, int bead, int i, int j);
void chksec(int lat, int dir, int bead);
void chkprim(int lat, int dir, int bead);
int chkdir(int lat, int dir, int bead);
int rstart[2], tempi[1], tempi2[1], endtraj;

void **alloc2d(int varsize, int n, int p) ; 
void ***alloc3d(int varsize, int n, int p, int q) ;
void ****alloc4d(int varsize, int n, int p, int q, int r) ;
void *****alloc5d(int varsize, int n, int p, int q, int r, int s) ;
void *******alloc7d(int varsize, int n, int p, int q, int r, int s, int t, int u) ;
double gssran(long *idum) ;
double mypow2(double x);
double mypow7(double x);
double mypow8(double x);
double mypow14(double x);
double ran2(long *idum) ;
FILE * mainout, *xyzout;
char fname[30];
char outname[30];
char ratname[30];
char wtname[30];
FILE * outPath;

/*-------------------------------------------------------------------------*/

int main(int argc,char** argv){

  int i, j, k, l, m, dim, size, temp, done, worker, mtag, n1, n2, part, lim;
  int op, dir1, dir2, from, to, rem, limit, pt, bead1, bead2;
  int lat, dir, bead, tcount;
  int rc, maxsize, ind, flag, mcount, it1, it2;
  int cycle, endcycle, gotmin;
  unsigned int iseed = (unsigned int)time(NULL);
  void back(int lat, int dir, int bead);
  void createservers(), initvar();
  int strinit();
  void move();
  void string(), stack(int lat, int dir, int bead), wrpath(int ind);
  int getmin(int *lat, int *dir, int *bead);
  void strread(), wtread(), fxread(), probread();
  void wtstack(int lat, int dir, int bead);
  void updtwts(int ind);
  void wrextend(int index);
  void rddelta(), rdcryst(), rdsolv(), rdemerg();
  double rate[2][2], sum1, sum2, t1, t2, telaps;
  int wtupdt_s;
  FILE * outFile;
  FILE * ratFile, *filein;

  ld[0] = 0;

  if (seed ==0) {
      seed = iseed;
  }

  tstepsq=tstep*tstep;
  beadsp2 = beads*beads;
  beadsp3 = beadsp2*beads;
  wtupdt_s = wtupdt/100;

  //COMPILE_STAMP;

#ifdef _OPENMP
  int desired = MPI_THREAD_MULTIPLE;
  int provided;
  MPI_Init_thread(&argc, &argv, desired, &provided);
#else
  MPI_Init(&argc, &argv);
#endif

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
//  printf("proc %d: MPI_Init done\n",me); fflush(stdout);
  
  GA_Initialize();
//  printf("proc %d: GA_Initialize done\n",me); fflush(stdout);
  
  MA_init(MT_DBL, 128*1024*1024, 16*1024*1024);
//  printf("proc %d: MA_init done\n",me); fflush(stdout);

  createservers();

  GA_Sync();
  if (me == 0) {
    printf("Starting..\n"); fflush(stdout);
  }
  nx=Lx/(int)blen;
  ny=Ly/(int)blen+1;
  nz=Lz/(int)blen;
  blinv = 1./blen;

  var=1./(beta*solvmass);
  var2=1./(beta*polymass);
  strmfrq = (int)(1./tstep);       /* stream every 4 ps (1 time unit) */
  solvstep = tstep*strmfrq;
  phan[0] = 20.;
  phan[1] = 0.;
  phan[2] = 92.5;
  n_av=N/(float)(nx*(ny-1)*nz);	

  kf = 20;
  R0 = 2;
  ooR02 = 1./(R0*R0);
  eph1_12 = 0.7*12;
  eph2_12 = 0.7*12;
  ep1_6 = 1*6;
  sigma1 = 7;
  sigma2 = 3.5;
  sigma16 = mypow8(sigma1)/mypow2(sigma1);
  sigma26 = mypow8(sigma2)/mypow2(sigma2);
  phanr06 = mypow8(phanr0)/mypow2(phanr0);
  calpha[0]=0.;
  salpha[0]=1.;
  calpha[1]=0.;
  salpha[1]=1.;

  if ((verb >= 3) || ((verb >=2) && (me == 0))) {
      sprintf(outname,"%s%d%s","out",me,".dat");
      mainout = fopen(outname,"w");
      fprintf(mainout,"rank = %d\n", me);
      fclose(mainout);
  }

  if (me == 0) {
      sprintf(wtname,"%s","weight.dat");
      ratFile = fopen(wtname,"w");
      fclose(ratFile);
      
      sprintf(ratname,"%s","rate.dat");
      ratFile = fopen(ratname,"w");
      fclose(ratFile);
  }

  seed *= me+1;
  for (i=0;i<100;i++){
    sum1 = ran2(&seed);
  }
  
  lim = 2*beads;
  twavg1 = (double *) malloc(2*beads*sizeof(double));
  twavg2 = (double *) malloc(2*beads*sizeof(double));
  tweight1 = (double *) malloc(2*beads*sizeof(double));
  tweight2 = (double *) malloc(2*beads*sizeof(double));
  ttau1 = (double *) malloc(2*beads*sizeof(double));
  ttau2 = (double *) malloc(2*beads*sizeof(double));
  tstravg1 = (double *) malloc(beads*nop*sizeof(double));
  tstravg2 = (double *) malloc(beads*nop*sizeof(double));
  tz1 = (double *) malloc(2*beads*nop*sizeof(double));
  tz2 = (double *) malloc(2*beads*nop*sizeof(double));
  tpt = (double *) malloc(2*ndim*npart*sizeof(double));
  lfrac = (double *) malloc(2*2*beads*2*beads*sizeof(double));

  tint1 = (int *) malloc(2*beads*sizeof(int));
  tint2 = (int *) malloc(2*beads*sizeof(int));
  ton1 = (int *) malloc(2*beads*sizeof(int));
  ton2 = (int *) malloc(2*beads*sizeof(int));
  tnwstack1 = (int *) malloc(2*beads*sizeof(int));
  tnwstack2 = (int *) malloc(2*beads*sizeof(int));
  lim = 2*beads*2*beads*2*beads;
  tprob1 = (int *) malloc(lim*sizeof(int));
  tprob2 = (int *) malloc(lim*sizeof(int));
  lim = 2*beads*2*beads;
  tnprob1 = (int *) malloc(lim*sizeof(int));
  tnprob2 = (int *) malloc(lim*sizeof(int));
  tthist1 = (int *) malloc(2*beads*tres*sizeof(int));
  tthist2 = (int *) malloc(2*beads*tres*sizeof(int));

  z = (opoint_t ***) alloc3d(sizeof(opoint_t), 2, 2, beads);
  on = (int ***) alloc3d(sizeof(int), 2, 2, beads);
  weight = (double ***) alloc3d(sizeof(double), 2, 2, beads);
  wavg = (double ***) alloc3d(sizeof(double), 2, 2, beads);
  stravg= (opoint_t **) alloc2d(sizeof(opoint_t), 2, beads);
  navg = (int **) alloc2d(sizeof(int), 2, beads);
  count= (int ***) alloc3d(sizeof(int), 2, 2, beads);
  thist = (int ****) alloc4d(sizeof(int), 2, 2, beads, tres);
  if (wtalg == 2) {
    tau = (double ***) alloc3d(sizeof(double), 2, 2, beads);
    ntau = (int ***) alloc3d(sizeof(int), 2, 2, beads);
  }

  box = (int ***) alloc3d(sizeof(int), nx, ny, nz);
  box_poly = (int ***) alloc3d(sizeof(int), nx, ny, nz);
  
  maxd = 5;
  lab = (int ****) alloc4d(sizeof(int),nx,ny,nz,maxd);
  nlab = (int ***) alloc3d(sizeof(int),nx,ny,nz);

  r0 = (double **) alloc2d(sizeof(double), npart, 3);
  rvecs = (xpoint_t **) alloc2d(sizeof(xpoint_t), npart, npart);
  r0dists = (double **) alloc2d(sizeof(double), npart, npart);
  r0d6 = (double **) alloc2d(sizeof(double), npart, npart);
  r0d12 = (double **) alloc2d(sizeof(double), npart, npart);
  rdists = (double **) alloc2d(sizeof(double), npart, npart);
  nextstep = (double **) alloc2d(sizeof(double), npart, 3);
  Delta = (int ***) alloc3d(sizeof(int), npart, npart,4);

  rsolv = (double **) alloc2d(sizeof(double), N, 3);
  vsolv = (double **) alloc2d(sizeof(double), N, 3);

  v_cmax = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmay = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmaz = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmx = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmy = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmz = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  pos = (int **) alloc2d(sizeof(int), N, 3);
  pos_poly = (int **) alloc2d(sizeof(int), npart, 3);

  f1 = (double **) alloc2d(sizeof(double), npart, 3);

  rdcryst();
  rddelta();
  rdsolv();
  rdemerg();
 
  initvar();
  
  if (me == 0) {
    printf("proc %d is initializing the string\n",me); fflush(stdout);
    
    /* initialize string */

    if (strcmp(flname,"NOREAD") == 0) {
	if (verb >= 2) {
	    printf("Automatically initializing string\n");
	}
	if (strinit())
	    gaexit(2) ;  
    } else {
	if (verb >= 2) {
	    printf("Using string configuration from %s\n", flname);
	}
	strread();
    }
    wrpath(0);
    
    /* send string to server */

    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (i=0;i<beads;i++) {
	for (j=0;j<nop;j++) {
	  tz1[mcount] = z[0][dir][i].x[j];
	  tz2[mcount] = z[1][dir][i].x[j];
	  mcount++;
	}
      }
    }
    //LINE_STAMP;
    PutToServer_dbl(&Coords_z,0,-1,tz1);
    //LINE_STAMP;
    PutToServer_dbl(&Coords_z,1,-1,tz2);
    
    /* read in fprob */

    if (strcmp(pname,"NOREAD") != 0) {
      if (verb >= 2) {
	printf("Reading prob data from %s\n", pname);
      }
      probread();
    }

    /* initialize weights */
    
    if (strcmp(wname,"NOREAD") == 0) {
      if (verb >= 2) {
	printf("Automatically initializing weights\n");
      }
      mcount = 0;
      for (j=0;j<2;j++) {
	for (k=0;k<beads;k++) {
	  tweight1[mcount] = 1./(float)(2*beads);
	  tweight2[mcount] = 1./(float)(2*beads-2);
	  mcount++;
	}
      }
      //LINE_STAMP;
      PutToServer_dbl(&Coords_weight,0,-1,tweight1);
      //LINE_STAMP;
      PutToServer_dbl(&Coords_weight,1,-1,tweight2);
      //LINE_STAMP;
      PutToServer_dbl(&Coords_wadj,0,-1,tweight1);
      //LINE_STAMP;
      PutToServer_dbl(&Coords_wadj,1,-1,tweight2);
    } else {
      if (verb >= 2) {
	printf("Using weight configuration from %s\n", wname);
      }
      wtread();
    }

    if (strcmp(fxname,"NOREAD") != 0) {  /* read in fluxes */
      if (verb >= 2) {
	printf("Using flux configuration from %s\n", fxname);
      }
      fxread();
    }
    printf("proc %d is done with the string\n",me); fflush(stdout);
  }

  /* check point */
  GA_Sync();

  if (me == 0) {
      t1 = MPI_Wtime();
  }
  getstring();

  /* ---------------start of dynamics loop-------------------------- */
  
  for (cycle=1;cycle<=T;cycle++) { 
    
    if (me == 0) {
      printf("Starting cycle %d of %d\n",cycle,T);
    
      mcount = 0;
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<=beads-1;bead++) {
	  tint1[mcount] = 1;
	  mcount++;
	}
      }
      //LINE_STAMP;
      PutToServer_int(&Coords_count,0,-1,tint1);  // set all counts to 1
      //LINE_STAMP;
      PutToServer_int(&Coords_count,1,-1,tint1);
      PutToServer_int(&Coords_on,0,-1,tint1);    // set all regions to 'on'
      PutToServer_int(&Coords_on,1,-1,tint1);
    }
    GA_Sync();

    endcycle = 0;
    while (!endcycle) {
      
      /* start a new trajectory */
      tcount = getmin(&lat,&dir,&bead);   /* returns (lat,dir,bead) for the region with the lowest counter */
      
      if (lat == -1) {
	endcycle = 1;
      } else {
	
	if ((verb >= 3) || ((verb >= 2) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"Running on %d %d %d, start count = %d\n",lat,dir,bead,tcount);
	  fclose(mainout);
	}

	if (verb >= 3) {
	  sprintf(filename,"flow%d.xyz",me);
	  xyzout = fopen(filename,"w");
	  fprintf(xyzout,"%i \n polymer movie\n", npart);
	}
	
	if (((tcount != 1) || (strcmp(fxname,"NOREAD") != 0))|| (cycle!=1)) {
	  back(lat,dir,bead);  /* initialize trajectory */
	} else {
	  coor = ptconv(lat,dir,bead);
	  ocoor = projx(coor);
	  gtst(lat,dir);
	  rstart[0] = -1;
	  rstart[1] = -1;
	  mytim = 0;
	  myclock = 0;
	}

	endtraj = 0;
	while (!endtraj) {
	    
	  /* move */
	    
	  move();

	  if (!polyerr) {

	    mytim += 1;
	    myclock += 1;
		
	    /* stack */
		
	    if ((mytim-1)%xyzfrq ==0) {
	      if (verb >= 3) {
		mainout = fopen(outname,"a");
		fprintf(mainout,"proc %d: mytim = %d writing to xyz \n",me, mytim);
		fclose(mainout);
		wrxyz(mytim-1);
	      }
	    }
		
	    if (tcount%stkfrq == 0) {
	      stack(lat,dir,bead);
	    }
		
	    if (wtalg == 1) { /* local wt alg */
	      if ((tcount + cycle*every)%wtupdt_s == 0) {
		wtstack(lat,dir,bead);
	      }
	    }
		
	    if (mytim%chkfrq == 0) {
		    
	      /* check */
		    
	      if (!chkdir(lat,dir,bead)) {
		chksec(lat,dir,bead);
		chkprim(lat,dir,bead);
	      }
	    }
	    tcount++;
	    if (tcount%globfrq == 0) {
	      ind = dir*beads + bead;
	      allind[0] = lat;
	      allind[1] = ind;
	      tempi[0] = globfrq + NGA_Read_inc(Coords_elaps.ga,allind,globfrq);
	      tempi[0] = globfrq + NGA_Read_inc(Coords_count.ga,allind,globfrq);
	      if (tempi[0] > every) {
		if (verb >= 3) {
		  mainout = fopen(outname,"a");
		  fprintf(mainout,"Ended traj: count\n");
		  fclose(mainout);
		}
		/* write coordinates to (someone's) fluxlist */
		ind = dir*beads + bead;
		//LINE_STAMP;                                    
		GetFromServer_dbl(&Coords_weight,lat,ind,tempw);
		
		if (!chkdir(lat,dir,bead)) {
		  chkprim(lat,dir,bead);
		}
		if (!endtraj) {
		  if (rstart[0] != -1) {
		    addpoint(lat,dir,bead,coor,tempw[0],rstart[0],rstart[1]);
		  }
		  endtraj = 1;
		}
	      }
	    }
	  } else {   /* if polyerr */
	    endtraj = 1;
	    polyerr = 0;
	  }
	}  /* end of trajectory loop */
	ind = dir*beads + bead;
	allind[0] = lat;
	allind[1] = ind;
	tempi[0] = tcount%globfrq + NGA_Read_inc(Coords_elaps.ga,allind,tcount%globfrq);
	tempi[0] = tcount%globfrq + NGA_Read_inc(Coords_count.ga,allind,tcount%globfrq);
	if ((verb >= 3) || ((verb >=2) && (me == 0))) {
	    mainout = fopen(outname,"a");
	    fprintf(mainout,"Count on %d %d %d, is now = %d\n",lat,dir,bead,tempi[0]);
	    fclose(mainout);
	}
      }
    }       /* end of cycle */

      
    if (phase == 1) {  /* move the string */
      
      if (me == 0) {
	
	if (verb >= 2) {
	  printf("Moving string..\n",cycle);
	}
	string();
	wrpath(cycle);
      }
      /* check point */
      GA_Sync();
      
      getstring();
      
      for (lat=0;lat<=1;lat++) {
	for (dir=0;dir<=1;dir++) {
	  for (bead=0;bead<=beads-1;bead++) {
	    ind = dir*beads + bead;
	    if ((ind%nproc) == me) {
	      if (verb >=3) {
		printf("proc %d is doing ind = %d of %d\n",me,ind,beads+beads-1);
	      }
	      nobad(lat,dir,bead);
	    }
	  }
	}
      }
      
      /* checkpoint */
      GA_Sync();

    }
    
    /* compute the rate / update weights / write thist */

    GA_Sync();
    if (me == 0) {
      printf("proc %d: TCOB\n",me); fflush(stdout);

      /* compute rate */

      //LINE_STAMP;
      GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
      //LINE_STAMP;
      GetFromServer_dbl(&Coords_weight,1,-1,tweight2);
      
      mcount = 0;
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<=beads-1;bead++) {
	  weight[0][dir][bead] = tweight1[mcount];
	  weight[1][dir][bead] = tweight2[mcount];
	  mcount++;
	}
      }
      
      for (dir=0;dir<=1;dir++) {
	for (lat=0;lat<=1;lat++) {
	  rate[lat][dir] = 0.;
	  for (bead=0;bead<=beads-1;bead++) {
	    ind = dir*beads + bead;
	    //LINE_STAMP;
	    GetFromServer_int(&Coords_elaps,lat,ind,tempi);
	    telaps = (float)tempi[0];
	    //LINE_STAMP;
	    GetFromServer_int(&Coords_narr,lat,ind,tempi);

	    if (telaps != 0.) {
	      rate[lat][dir] += tempi[0]/telaps;
	    }
	  }
	  rate[lat][dir] /= calcsum(lat,dir);
	}
      }
      ratFile = fopen(ratname,"a");
      fprintf(ratFile,"%e %e\n",0.5*(rate[0][0]+rate[1][0]),0.5*(rate[0][1]+rate[1][1]));
      fclose(ratFile);

      /* update weights */

      if ((cycle*every)%wtupdt == 0) {
	updtwts((cycle*every)/wtupdt);
      }

      /* write out thist */

      if ((cycle*every)%wrfrq == 0) {
	wrextend((cycle*every)/wrfrq);
      }
      printf("proc %d: done TCOB\n",me); fflush(stdout);
    }

    /* checkpoint */
    GA_Sync();

    /* do other end-of-cycle calcs: projections, etc. */
  }
  
  /* ---------------end of dynamics loop---------------------------- */
  
  /* set the restart point: write fluxes, weights, etc. */
  
  if (me == 0) {
    t2 = MPI_Wtime();
    printf("total time of dynamics steps: %f\n",t2-t1); fflush(stdout);
    printf("total time of locking steps: %f\n",globtim); fflush(stdout);
    printf("total number of locking steps: %d\n",nglobtim); fflush(stdout);
    printf("time per locking step: %f\n",globtim/(float)nglobtim); fflush(stdout);
    GetFromServer_int(&Coords_count,0,-1,tint1);  
    GetFromServer_int(&Coords_count,1,-1,tint2);
    
    mcount = 0;
    tcount = 0;
    for (i=0;i<2;i++) {
      for (j=0;j<beads;j++) {
	tcount += tint1[mcount];
	tcount += tint2[mcount];
	mcount++;
      }
    }
    printf("total number of dynamics steps: %d\n",tcount); fflush(stdout);
    printf("time per dynamics step: %f\n",(t2-t1)/(float)tcount); fflush(stdout);
    

  /* write weights */

    sprintf(outname,"%s","fwts.dat");
    outFile = fopen(outname,"w");
    for (i=0;i<=1;i++) {
      for (j=0; j<=1; j++) {
	for (k=0; k<=beads-1; k++) {
	  fprintf(outFile,"%e\n",weight[i][j][k]);
	}
      }
    }
      fclose(outFile);

    /* write fluxes */
  
    sprintf(outname,"%s","fflux.dat");
    outFile = fopen(outname,"w");
    for (i=0;i<=1;i++) {
      for (j=0; j<=1; j++) {
	for (k=0; k<=beads-1; k++) {
	  if ((i == 0) || (k!=beads-1)) {
	    fprintf(outFile,"%d %d %d\n",i,j,k);
	    ind = j*beads + k;
	    //LINE_STAMP;
	    GetFromServer_int(&Coords_full,i,ind,tempi);
	    if (tempi[0] == 1) {
	      limit = mxlist;
	    } else {
	      //LINE_STAMP;
	      GetFromServer_int(&Coords_nlist,i,ind,tempi);
	      limit = tempi[0];
	    }

	    //LINE_STAMP;
	    GetFromServer_dbl(&Coords_twlist,i,ind,tempw);
	    fprintf(outFile,"%d %e\n", limit, tempw[0]);

	    for (pt=0; pt<limit; pt++) {
	      ind = 2*mxlist*beads*j + 2*mxlist*k + 2*pt;
	      //LINE_STAMP;
	      GetFromServer_int(&Coords_from,i,ind,tempi);
	      //LINE_STAMP;
	      GetFromServer_int(&Coords_from,i,ind+1,tempi2);
	      fprintf(outFile,"%d %d ", tempi[0], tempi2[0]);

	      ind = mxlist*beads*j + mxlist*k + pt;
	      //LINE_STAMP;
	      GetFromServer_dbl(&Coords_wlist,i,ind,tempw);
	      //LINE_STAMP;
	      GetFromServer_int(&Coords_clock,i,ind,tempi);
	      fprintf(outFile,"%e %d ",tempw[0],tempi[0]);
	      GetFromServer_int(&Coords_mytim,i,ind,tempi);
	      fprintf(outFile,"%d ",tempi[0]);

	      lo[0] = i;
	      hi[0] = i;
	      lo[1] = 2*npart*ndim*mxlist*beads*j + 2*npart*ndim*mxlist*k + 2*npart*ndim*pt;
	      hi[1] = 2*npart*ndim*mxlist*beads*j + 2*npart*ndim*mxlist*k + 2*npart*ndim*(pt+1) -1;
	      NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);

	      mcount = 0;
	      for (part=0; part<npart;part++) {
		for (dim=0; dim<=ndim-1; dim++) {
		  fprintf(outFile,"%e ", tpt[mcount]);
		  mcount++;
		  fprintf(outFile,"%e ", tpt[mcount]);
		  mcount++;
		}
	      }
	      fprintf(outFile,"\n");
	    }
	  }
	}
      }
    }
    fclose(outFile);

    /* write probs */
  
    GetFromServer_int(&Coords_nprob,0,-1,tnprob1);
    //LINE_STAMP;
    GetFromServer_int(&Coords_nprob,1,-1,tnprob2);
    //LINE_STAMP;
    GetFromServer_int(&Coords_prob,0,-1,tprob1);
    //LINE_STAMP;
    GetFromServer_int(&Coords_prob,1,-1,tprob2);

    sprintf(outname,"%s","fprob.dat");
    outFile = fopen(outname,"w");
    for (lat=0;lat<=1;lat++) {
      for (dir=0; dir<=1; dir++) {
	for (bead=0; bead<=beads-1; bead++) {
	  for (n1=0; n1<2; n1++) {
	    for (n2=0; n2<=beads-1; n2++) {
	      it1 = getnprob(lat,dir,bead,n1,n2);
	      if (it1 != 0) {
		fprintf(outFile,"%d %d %d %d %d %d\n",lat,dir,bead,n1,n2,it1);
		for (i=0; i<2; i++) {
		  for (j=0; j<=beads-1; j++) {
		    it2 = getprob(lat,dir,bead,n1,n2,i,j);
		    if (it2 != 0) {
		      fprintf(outFile,"%d %d %d\n",i,j,it2);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    fclose(outFile);
		    
    /* write results */
  }
  
  /* check point */
//  printf("proc %d: final check\n",me); fflush(stdout);
  GA_Sync();

  status = DestroyCoordServer(&Coords_wavg);
  status = DestroyCoordServer(&Coords_wadj);
  status = DestroyCoordServer(&Coords_weight);
  status = DestroyCoordServer(&Coords_elaps);
  status = DestroyCoordServer(&Coords_tau);

  status = DestroyCoordServer(&Coords_count);
  status = DestroyCoordServer(&Coords_narr);
  status = DestroyCoordServer(&Coords_on);
  status = DestroyCoordServer(&Coords_bdedit);
  status = DestroyCoordServer(&Coords_ntau);
  status = DestroyCoordServer(&Coords_navg);

  status = DestroyCoordServer(&Coords_stravg);
  status = DestroyCoordServer(&Coords_z);

  status = DestroyCoordServer(&Coords_prob);
  status = DestroyCoordServer(&Coords_nprob);
  status = DestroyCoordServer(&Coords_thist);
  status = DestroyCoordServer(&Coords_full);
  status = DestroyCoordServer(&Coords_nlist);
  status = DestroyCoordServer(&Coords_twlist);
  status = DestroyCoordServer(&Coords_pts);
  status = DestroyCoordServer(&Coords_wlist);
  status = DestroyCoordServer(&Coords_from);
  status = DestroyCoordServer(&Coords_clock);
  status = DestroyCoordServer(&Coords_mytim);
  status = DestroyCoordServer(&Coords_nwstack);

  if (me == 0) GA_Print_stats();
  
  GA_Terminate();
  MPI_Finalize();
  
  return(0);
}




/*=========================================================
  END OF MAIN
  ============================================================*/

/*-------------------------------------------------------------------*/
void back (int lat, int dir, int bead) {

  /* This subroutine resets the position of a walker
     to a saved point, or if none exist, to the position
     of the bead (guestimating for undefined coordinates) */
  
  double sum, a, ttwlist, oran, t1, t2;
  int i, j, newreg, ind, limit, k, mcount;
  int tnlist, tfull, part, dim, good;
  point_t point;
  FILE *filein;


  lockind = getlockind(0,lat,dir,bead);
  //LINE_STAMP;
  if (me==0) {
    t1 = MPI_Wtime();
  }
  GA_Lock(lockind);
  if (me==0) {
    t2 = MPI_Wtime();
    globtim += t2-t1;
    nglobtim++;
  }

  good = 0;
  while (!good) {
    good = 1;
    
    ind = dir*beads + bead;
    //LINE_STAMP;
    GetFromServer_int(&Coords_nlist,lat,ind,tempi);
    tnlist = tempi[0];
    
    //LINE_STAMP;
    GetFromServer_int(&Coords_full,lat,ind,tempi);
    tfull = tempi[0];
    
    if (tnlist != 0) {
      sum = 0.;
      oran = ran2(&seed);
      a = oran;
      //LINE_STAMP;
      GetFromServer_dbl(&Coords_twlist,lat,ind,tempw);
      ttwlist = tempw[0];
      a *= ttwlist;
      
      if (ttwlist <= 0) {
	printf("proc %d:  twlist < 0!! %e, lat,dir,bead = %d %d %d\n", me, ttwlist, lat, dir, bead);
	printf("proc %d:  nlist = %d; full = %d\n", me, tnlist, tfull);
	if (tfull) {
	  limit = mxlist;
	} else {
	  limit = tnlist;
	}
	for (i=0;i<=limit-1;i++) {
	  ind = mxlist*beads*dir + mxlist*bead + i;
	  //LINE_STAMP;
	  GetFromServer_dbl(&Coords_wlist,lat,ind,tempw);
	  printf("i = %d; wt = %e\n",i,tempw[0]);
	}
	gaexit(9);
      }

      j = 1;

      while (sum < a) {
	if ((!tfull) && (j > tnlist)) {  /* nlist is the actual number of points in the list */
	  sum = 0.;                                                  /* so subtract 1 to use as an array index */
	  printf("proc %d:  Error!  Weights not deleted!\n",me);
	  printf("Random Number: %e\n", a/ttwlist);
	  printf("dir = %d\n", dir);
	  printf("nlist(%d,%d)=%d \n", lat,bead,tnlist);
	  for (i=0;i<tnlist;i++) {
	    ind = mxlist*beads*dir + mxlist*bead + i;
	    //LINE_STAMP;
	    GetFromServer_dbl(&Coords_wlist,lat,ind,tempw);
	    sum = sum + tempw[0];
	  }
	  printf("j = %d\n", j);
	  printf("Manual sum: %e\n", sum);
	  printf("twlist = %e\n", ttwlist);
	  printf("Resetting twlist..\n");
	  j = tnlist+1;
	  printf("Use j = %d\n", j);
	  tempw[0] = sum;
	  ind = dir*beads + bead;
	  //LINE_STAMP;
	  GA_Init_fence();
	  PutToServer_dbl(&Coords_twlist,lat,ind,tempw);
	  GA_Fence();
	  break;
	} else if ((tfull) && (j > mxlist)) {
	  sum = 0.;
	  printf("proc %d:  Error!  Weights not deleted!\n", me);
	  printf("Random Number: %e\n", a/ttwlist);
	  printf("dir = %d\n", dir);
	  printf("List full\n");
	  for (i=0;i<mxlist;i++) {
	    ind = mxlist*beads*dir + mxlist*bead + i;
	    //LINE_STAMP;
	    GetFromServer_dbl(&Coords_wlist,lat,ind,tempw);
	    sum = sum + tempw[0];
	  }
	  printf("Manual sum: %e\n", sum);
	  printf("twlist = %e\n", ttwlist);
	  printf("Resetting twlist..\n");
	  j = mxlist+1;
	  printf("Use j = %d\n", j);
	  tempw[0] = sum;
	  ind = dir*beads + bead;
	  //LINE_STAMP;
	  GA_Init_fence();
	  PutToServer_dbl(&Coords_twlist,lat,ind,tempw);
	  GA_Fence();
	  break;
	}
	ind = mxlist*beads*dir + mxlist*bead + j-1;
	//LINE_STAMP;
	GetFromServer_dbl(&Coords_wlist,lat,ind,tempw);
	sum += tempw[0];
	j++;
      }

      lo[0] = lat;
      hi[0] = lat;
      lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(j-2);
      hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(j-1) -1;
      NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);

      mcount = 0;
      for (part=0; part<npart;part++) {
	for (dim=0; dim<=ndim-1; dim++) {
	  coor.x[part][dim] = tpt[mcount];
	  mcount++;
	  coor.v[part][dim] = tpt[mcount];
	  mcount++;
	}
      }

      ind = 2*mxlist*beads*dir + 2*mxlist*bead + 2*(j-2);
      //LINE_STAMP;
      GetFromServer_int(&Coords_from,lat,ind,tempi);
      rstart[0] = tempi[0];
      //LINE_STAMP;
      GetFromServer_int(&Coords_from,lat,ind+1,tempi);
      rstart[1] = tempi[0];
      
      ind = mxlist*beads*dir + mxlist*bead + (j-2);
      //LINE_STAMP;
      GetFromServer_int(&Coords_clock,lat,ind,tempi);
      myclock = tempi[0];
      GetFromServer_int(&Coords_mytim,lat,ind,tempi);
      mytim = tempi[0];
    } else {
      j = 0;
      ind = dir*beads + bead;
      //LINE_STAMP;
      GetFromServer_int(&Coords_count,lat,ind,tempi);
      if (tempi[0] != 0) {
	  if (verb >= 2) {
	      printf("proc %d:  %d  Reset  %d, %d, %d, without a flux point!\n", me, tempi[0],lat+1, dir+1, bead+1); 
	  }
      }
      rstart[0] = -1;
      rstart[1] = -1;
      coor = ptconv(lat,dir,bead);
      mytim = 0;
    }
    
    
    /* get ocoor, state, bas, check the validity of the point */
    
    ocoor = projx(coor);
    
    gtst(lat,dir);
    
    bas = bascheck();
    if (dir + bas == 1) {
      printf("proc %d:  Error! Started in opposite basin!\n",me);
      if (j == 0) {
	printf("Turning %d %d %d off.\n", lat, dir, bead);
	
	ind = dir*beads + bead;
	tempi[0] = 0;
	//LINE_STAMP;
	PutToServer_int(&Coords_on,lat,ind,tempi);
      } else {
	printf("Bad flux point j = %d\n", j);
	good = 0; 
	//	gaexit(97);
      }
    }
    newreg = which(lat,dir);
    if (newreg != bead) {
      if (verb >= 2) {
	printf("proc %d:  Error:  Replica reinitialized in wrong region!\n",me);
	printf("which(%d,%d) = %d\n", lat, dir, newreg);
	printf("opoint: %e %e %e\n",ocoor.x[0],z[lat][dir][bead],z[lat][dir][newreg]);
	printf("rank is %d; j is %d; lat,dir,bead are (%d,%d,%d) \n", me, j, lat,dir,bead);
	printf("tnlist: %d, oran: %e\n",tnlist,oran);
	printf("from: %d %d\n",rstart[0],rstart[1]);
	printf("deleting point %d\n",j-2);
	if (j >= 2) {
	  if (tfull) {
	    globn = mxlist;
	  } else {
	    globn = tnlist;
	  }
	  delpoint(lat,dir,bead,j-2);
	  ind = dir*beads + bead;
	  tempi[0] = globn;
	  PutToServer_int(&Coords_nlist,lat,ind,tempi);
	  tempi[0] = 0;
	  PutToServer_int(&Coords_full,lat,ind,tempi);
	}
	good = 0;
      }
	//  gaexit(96);
    }
  }
  GA_Unlock(lockind);
}
/*-------------------------------------------------------------------------*/
 void fixwt(int lat, int dir1, int dir2, int reg1, int reg2) {
   
   double temp, w[1], tempe[1], telaps, telaps0;
   int i, i2;

   i = dir1*beads + reg1;
   i2 = dir2*beads + reg2;
   //LINE_STAMP;
   GetFromServer_dbl(&Coords_wadj,lat,i,w);

   //LINE_STAMP;
   GetFromServer_int(&Coords_elaps,lat,i,tempi);
   telaps = (float)tempi[0];

   //LINE_STAMP;
   GetFromServer_int(&Coords_elaps,0,0,tempi);
   telaps0 = (float)tempi[0];

   if (telaps != 0.) {
     temp = s*w[0]*telaps0/telaps;
   } else {
     temp = 0.;
   }

   w[0] -= temp;
   //LINE_STAMP;
   PutToServer_dbl(&Coords_wadj,lat,i,w);

   lockind = getlockind(2,lat,dir2,reg2);
   GA_Lock(lockind);
   GA_Init_fence();
   GetFromServer_dbl(&Coords_wadj,lat,i2,w);
   w[0] += temp;
   PutToServer_dbl(&Coords_wadj,lat,i2,w);
   GA_Fence();
   GA_Unlock(lockind);

 }
/*-------------------------------------------------------------------------*/
void wtstack(int lat, int dir, int bead) {
   
  int i;
  double tempwadj[1], tempwavg[1];
  
  i = dir*beads + bead;

  if (wtalg != 1) {
    printf("Error! wtalg = %d\n",wtalg);
    gaexit(10);
  }

  //LINE_STAMP;
  GetFromServer_dbl(&Coords_wadj,lat,i,tempwadj);

  lockind = getlockind(2,lat,dir,bead);
  GA_Lock(lockind);
  GA_Init_fence();
  GetFromServer_dbl(&Coords_wavg,lat,i,tempwavg);
  tempwavg[0] += tempwadj[0];
  PutToServer_dbl(&Coords_wavg,lat,i,tempwavg);
  GA_Fence();
  GA_Unlock(lockind);

  allind[0] = lat;
  allind[1] = i;
  tempi[0] = 1 + NGA_Read_inc(Coords_nwstack.ga,allind,1);
}
/*-------------------------------------------------------------------------*/
void updtwts(int ind) {
   
  void getfrac(int lat, int dir, int bead);
  double lim2;
  int tgt,lat,dir,bead,good,lim,i, mcount,j,n1,n2,mcount2;
  FILE * wtout, *tmpfile;
  char tmpname[30];
  int wtmat(int lat, int sdir, int sbead);

  lim = 2*beads;
  
  if (wtalg == 1) {

    /* set weight to wavg, reset wavg */
    
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_wavg,0,-1,twavg1);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_wavg,1,-1,twavg2);
    //LINE_STAMP;
    GetFromServer_int(&Coords_nwstack,0,-1,tnwstack1);
    //LINE_STAMP;
    GetFromServer_int(&Coords_nwstack,1,-1,tnwstack2);
    
    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	if (tnwstack1[mcount]) {
	  tweight1[mcount] = twavg1[mcount]/(float)tnwstack1[mcount];
	}
	if (tnwstack2[mcount]) {
	  tweight2[mcount] = twavg2[mcount]/(float)tnwstack2[mcount];
	}
	tnwstack1[mcount] = 0;
	tnwstack2[mcount] = 0;
	twavg1[mcount] = 0.;
	twavg2[mcount] = 0.;
	
	weight[0][dir][bead] = tweight1[mcount]; 
	weight[1][dir][bead] = tweight2[mcount]; 
	mcount++;
      }
    }

    //LINE_STAMP;
    PutToServer_dbl(&Coords_wavg,0,-1,twavg1);
    //LINE_STAMP;
    PutToServer_dbl(&Coords_wavg,1,-1,twavg2);
    //LINE_STAMP;
    PutToServer_dbl(&Coords_weight,0,-1,tweight1);
    //LINE_STAMP;
    PutToServer_dbl(&Coords_weight,1,-1,tweight2);
    //LINE_STAMP;
    PutToServer_int(&Coords_nwstack,0,-1,tnwstack1);
    //LINE_STAMP;
    PutToServer_int(&Coords_nwstack,1,-1,tnwstack2);

  } else {

    /* perform matrix calculations, reset counters */
    
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_weight,1,-1,tweight2);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_tau,0,-1,ttau1);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_tau,1,-1,ttau2);
    //LINE_STAMP;
    GetFromServer_int(&Coords_ntau,0,-1,tint1);
    //LINE_STAMP;
    GetFromServer_int(&Coords_ntau,1,-1,tint2);
    //LINE_STAMP;
    GetFromServer_int(&Coords_nprob,0,-1,tnprob1);
    //LINE_STAMP;
    GetFromServer_int(&Coords_nprob,1,-1,tnprob2);
    //LINE_STAMP;
    GetFromServer_int(&Coords_prob,0,-1,tprob1);
    //LINE_STAMP;
    GetFromServer_int(&Coords_prob,1,-1,tprob2);

    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	weight[0][dir][bead] = tweight1[mcount];
	weight[1][dir][bead] = tweight2[mcount];
	tau[0][dir][bead] = ttau1[mcount];
	tau[1][dir][bead] = ttau2[mcount];
	ntau[0][dir][bead] = tint1[mcount];
	ntau[1][dir][bead] = tint2[mcount];
	mcount++;
      }
    }

    for (lat=0;lat<=1;lat++) {
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<=beads-1;bead++) {
	  if ((lat == 0) || (bead != beads-1)) {
	    wavg[lat][dir][bead] = 0.;
	    getfrac(lat,dir,bead);
	  }
	}
      }
    }

    good = 1;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	if (good) {
	  if (!wtmat(1,dir,bead)) {
	    printf("Weight update not completed!\n");
	    good = 0;
	  }
	}
      }
      for (bead=0;bead<=beads-2;bead++) {
	if (good) {
	  if (!wtmat(2,dir,bead)) {
	    printf("Weight update not completed!\n");
	    good = 0;
	  }
	}
      }
    }
    if (good) {
      for (lat=0;lat<=1;lat++) {
	if (lat == 0) {
	  lim = 2*beads;
	} else {
	  lim = 2*(beads-1);
	}
	lim2 = 1.0/(float)lim;
	for (dir=0;dir<=1;dir++) {
	  for (bead=0;bead<=beads-1;bead++) {
	    if ((lat == 0) || (bead != beads-1)) {
	      weight[lat][dir][bead] = wavg[lat][dir][bead]*lim2;
	    }
	  }
	}
	if (verb >= 3) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"Weight update completed\n");
	  fclose(mainout);
	}
      }
      
      mcount = 0;
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<=beads-1;bead++) {
	  tweight1[mcount] = weight[0][dir][bead]; 
	  tweight2[mcount] = weight[1][dir][bead]; 
	  mcount++;
	}
      }
      //LINE_STAMP;
      PutToServer_dbl(&Coords_weight,0,-1,tweight1);
      //LINE_STAMP;
      PutToServer_dbl(&Coords_weight,1,-1,tweight2);
    }
  }
    
  /* write weights */

  if (me == 0) {
      wtout = fopen(wtname,"a");
      fprintf(wtout,"%d ",ind*wtupdt);
      for (lat=0;lat<=1;lat++) {
	  for (dir=0;dir<=1;dir++) {
	      for (bead=0;bead<=beads-1;bead++) {
		  if ((lat == 0) || (bead != beads-1)) {
		      fprintf(wtout,"%e ",weight[lat][dir][bead]);
		  }
	      }
	  }
      }
      fprintf(wtout,"\n");
      fclose(wtout);
  }
}
/*-------------------------------------------------------------------------*/
int wtmat(int lat, int stdir, int stbead) {

  /* this function calculates the weights */

  double **fmat, **b, *vec, sum, avtau, tw, wt, tn, tprob;
  double temp1;
  int lim, i, j, dir, bead, regi, dirind[2*beads], regind[2*beads];
  int ndir, ni, odir, good, sdir, sbead, ind, mcount;
  FILE * matout;
  good = 1;

  sdir = stdir + 1;
  sbead = stbead + 1;  /* to be compatible with Fortran indexing */

  if (lat == 1) {
    lim = 2*beads;
  } else {
    lim = 2*(beads-1);
  }

  fmat = (double **) alloc2d(sizeof(double), lim, lim);
  b = (double **) alloc2d(sizeof(double), lim, 1);
  if ((vec = (double *) malloc(lim * sizeof(double))) == NULL) {
    printf ("no memory left (vec)!\n") ;
    gaexit(1) ;
  }

  for(i=0;i<=lim-1;i++) {
    for(j=0;j<=lim-1;j++) {
      fmat[i][j] = 0.;
    }
  }

  for (i=1;i<=lim-1;i++) {
    sum = 0.;
    if (i+sbead <= lim/2) {
      dir = sdir;
      regi = i+sbead;
    } else if (i+sbead <= lim) {
      if (sdir == 1) {
	dir = 2;
      } else {
	dir = 1;
      }
      regi = i+sbead-lim/2;
    } else {
      dir = sdir;
      regi = i+sbead-lim;
    }
    if (dir == 2) {
      regi = lim/2-regi+1;
    }
    dirind[i-1] = dir;
    regind[i-1] = regi;

    avtau = tau[lat-1][dir-1][regi-1]/(float)ntau[lat-1][dir-1][regi-1];
    ind = (dir-1)*beads + (regi-1);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_twlist,lat-1,ind,tempw);
    tw = tempw[0];
    temp1 = 1./(avtau*tw);
    //tw = flist[lat-1][dir-1][regi-1].twlist;
    if (dir == 1) {
      if (regi != lim/2) {
	ndir = 1;
	ni = regi+1;
      } else {
	ndir = 2;
	ni = lim/2;
      }
    } else {
      if (regi != 1) {
	ndir = 2;
	ni = regi-1;
      } else {
	ndir = 1;
	ni = 1;
      }
    }
    if (i == lim-1) {
      dirind[i] = ndir;
      regind[i] = ni;
    }
    for (odir=1;odir<=2;odir++) {
      for (j=1;j<=lim/2;j++) {

	if (getnprob(lat-1,dir-1,regi-1,odir-1,j-1) != 0) {
	  tprob = getprob(lat-1,dir-1,regi-1,odir-1,j-1,ndir-1,ni-1)/(float)getnprob(lat-1,dir-1,regi-1,odir-1,j-1);
	  if (tprob != 0.) {
	    mcount = 4*beadsp2*lat + 2*beadsp2*dir + 2*beads*regi + beads*odir + j;
	    tn = lfrac[mcount];
	    sum += tprob*tn*temp1;
	  }
	}
      }
    }

    if (sum == 0.) {
      printf("1: Incomplete statistics!  Need longer wtupdt!\n");
      printf("%d %d %d %d %d %d %d\n", sdir, sbead, i, regind[i-1], dirind[i-1], ndir, ni);
      good = 0;
      free(fmat);
      free(b);
      free(vec);
      return good;
    }

    fmat[i-1][i-1] = sum;

    if (abs(sum) > 1) {
      printf("sum error! %d %e %e %e\n",i, avtau, tw, sum);
    }
    avtau = tau[lat-1][ndir-1][ni-1]/(float)ntau[lat-1][ndir-1][ni-1];
    mcount = 4*beadsp2*lat + 2*beadsp2*ndir + 2*beads*ni + beads*dir + regi;
    tn = lfrac[mcount];

    ind = (ndir-1)*beads + (ni-1);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_twlist,lat-1,ind,tempw);
    sum = tn/(avtau*tempw[0]);
    //sum = tn/(avtau*flist[lat-1][ndir-1][ni-1].twlist);

    if (sum == 0.) {
      printf("2: Incomplete statistics!  Need longer wtupdt!\n");
      printf("%d %d %d %d %d %e\n", lat, ndir, ni, dir, regi, tn);
      good = 0;
      free(fmat);
      free(b);
      free(vec);
      return good;
    }
    fmat[i][i-1] = -sum;
  }

  for (i=1;i<=lim;i++) { /* write 1s for normalization */
    fmat[i-1][lim-1] = 1.;
  }

  if (verb >= 2) {
    matout = fopen("fmat.dat","w");
    for (i=0;i<=lim-1;i++) {
      for (j=0;j<=lim-1;j++) {
	fprintf(matout,"%e ",fmat[j][i]);
      }
      fprintf(matout,"\n");
    }
    fclose(matout);
  }

  /* invert matrix */
  
  for (j=0;j<=lim-1;j++) {
    b[j][0] = j;
  }

  gaussj(fmat,lim,b,1);

  for (i=0;i<=lim-1;i++) {
    wt = weight[lat-1][dirind[i]-1][regind[i]-1];
    wavg[lat-1][dirind[i]-1][regind[i]-1] += wt + wfrac*(fmat[lim-1][i] - wt);

    if (fmat[lim-1][i] < 0.0) {
      printf("negative weights! %d %d %d\n", lat-1, dirind[i], regind[i]);
      good = 0;
      free(fmat);
      free(b);
      free(vec);
      return good;
    }
  }

  free(fmat);
  free(b);
  free(vec);
  return good;
}
/*-------------------------------------------------------------------------*/
void wrextend(int index) {

  /* this function writes out a theta histogram */

  int lat, dir, bead, ind, mcount, mcount2, tsum;
  double x[tres], sum, ext, ootsum;
  FILE * tout;
  char tname[30];
  
  for (ind=0;ind<=tres-1;ind++) {
    x[ind] = 0.;
  }

  //LINE_STAMP;
  GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
  //LINE_STAMP;
  GetFromServer_dbl(&Coords_weight,1,-1,tweight2);
  //LINE_STAMP;
  GetFromServer_int(&Coords_thist,0,-1,tthist1);
  //LINE_STAMP;
  GetFromServer_int(&Coords_thist,1,-1,tthist2);
 
  mcount = 0;
  mcount2 = 0;
  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      weight[0][dir][bead] = tweight1[mcount];
      weight[1][dir][bead] = tweight2[mcount];
      mcount++;
      for (ind=0;ind<tres;ind++) {
	thist[0][dir][bead][ind] = tthist1[mcount2];
	thist[1][dir][bead][ind] = tthist2[mcount2];
	mcount2++;
      }
    }
  }

  for (lat=0;lat<=1;lat++) {
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	tsum = 0;
	for (ind=0;ind<=tres-1;ind++) {
	  tsum += thist[lat][dir][bead][ind];
	}
	if (tsum) {
	    ootsum = 1./(float)tsum;
	    for (ind=0;ind<=tres-1;ind++) {
		if (thist[lat][dir][bead][ind]) {
		    x[ind] += thist[lat][dir][bead][ind]*weight[lat][dir][bead]*ootsum;
		}
	    }
	}
      }
    }
  }

  mcount2 = 0;
  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      for (ind=0;ind<tres;ind++) {
	tthist1[mcount2] = 0;
	tthist2[mcount2] = 0;
	mcount2++;
      }
    }
  }

  //LINE_STAMP;
  PutToServer_int(&Coords_thist,0,-1,tthist1);
  //LINE_STAMP;
  PutToServer_int(&Coords_thist,1,-1,tthist2);

  sum = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += x[ind];
  }

  sprintf(tname,"%s%d%s","thist",index,".dat");
  tout = fopen(tname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin + (tmax-tmin)*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x[ind]/sum);
  }
  fclose(tout);
}
/*-------------------------------------------------------------------------*/
void gtst(int lat, int dir) {

  /* this function updates the state variable for 
     the walker in the (lat,dir,bead) region */
  
  int nlat, i, j, off;
  double dist, lodist;

  if (lat == 0) {
    off = 1;
    nlat = 1;
  } else {
    off = 0;
    nlat = 0;
  }

  for (i=0;i<=beads-off-1;i++) {
    dist = 0.;
    dist = ogtdist(z[nlat][dir][i],ocoor);
    if (i == 0) {
      lodist = dist;
      state = i;
    }else if (dist < lodist) {
      lodist = dist;
      state = i;
    }
  }
}

/*-------------------------------------------------------------------------*/
int which(int lat, int dir) {

  /* this function returns the region that the walker
   in (lat,dir,bead) is in, presently */
       
   int i, j, off, temp;
  double dist, lodist;

  if (lat == 0) {
    off = 0;
  } else {
    off = 1;
  }

  for (i=0;i<=beads-off-1;i++) {
    dist = ogtdist(z[lat][dir][i],ocoor);
    if (i == 0) {
      lodist = dist;
      temp = i;
    } else if (dist < lodist) {
      lodist = dist;
      temp = i;
    }
  }
  return temp;
}

/*-------------------------------------------------------------------------*/
int whichx(int lat, int dir, point_t x) {
   
   /* this function returns the region that the point 'x'
      is in, presently */
   
   int i, j, off, temp;
   double dist, lodist;
   opoint_t ox;
   
   ox = projx(x);
   
   if (lat == 0) {
     off = 0;
   } else {
     off = 1;
   }
   
   for (i=0;i<=beads-off-1;i++) {
     dist = ogtdist(z[lat][dir][i],ox);
     if (i == 0) {
       lodist = dist;
       temp = i;
     } else if (dist < lodist) {
       lodist = dist;
       temp = i;
     }
  }
   return temp;
 }
/*-------------------------------------------------------------------------*/
double ogtdist(opoint_t px, opoint_t py) {
   
  /* this function returns the distance between two 
     points in order parameter space */
   
   int i;
   double temp;
   
   temp = 0.;
   
   for (i=0;i<nop;i++) {
     temp = temp + (px.x[i] - py.x[i])*(px.x[i] - py.x[i]);
   }
   temp = sqrt(temp);
   
   return temp;
 }
/*-------------------------------------------------------------------------*/
void addpoint (int lat, int dir, int bead, point_t x, double wt, int from1, int from2) {

  /* this function adds the point 'x' to the (lat,dir,bead) fluxlist
     along with the weight 'wt' */
  
  int tgt, n, i, ind, tfull, wind, part, dim, ptind, mcount;
  double tw[1], ttw[1], t1,t2;
  
  // printf("proc %d: adding point to [%d,%d,%d] w = %e\n",me,lat,dir,bead,wt); fflush(stdout);

  lockind = getlockind(0,lat,dir,bead);
  //LINE_STAMP;
  if (me==0) {
	t1 = MPI_Wtime();
    }
  GA_Lock(lockind);
  if (me==0) {
      t2 = MPI_Wtime();
      globtim += t2-t1;
      nglobtim++;
  }
  GA_Init_fence();
  ind = dir*beads + bead;
  
  allind[0] = lat;
  allind[1] = ind;
  tempi[0] = 1 + NGA_Read_inc(Coords_nlist.ga,allind,1);

  //LINE_STAMP;
  GetFromServer_int(&Coords_full,lat,ind,tempi);  
  tfull = tempi[0];
  if (tfull == 0) {   /* not full */
    //LINE_STAMP;
    GetFromServer_int(&Coords_nlist,lat,ind,tempi);  
    if (tempi[0] == mxlist) {
      tempi2[0] = 1;
      //LINE_STAMP;
      PutToServer_int(&Coords_full,lat,ind,tempi2);
    }
  } else {
    //LINE_STAMP;
    GetFromServer_int(&Coords_nlist,lat,ind,tempi);  
    if (tempi[0] == mxlist+1) {
      tempi[0] = 1;
      //LINE_STAMP;
      PutToServer_int(&Coords_nlist,lat,ind,tempi);
    }
    n = tempi[0] - 1;

    wind = mxlist*beads*dir + mxlist*bead + n;
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_wlist,lat,wind,tw);
    //LINE_STAMP;

    GetFromServer_dbl(&Coords_twlist,lat,ind,ttw);
    ttw[0] -= tw[0];
    PutToServer_dbl(&Coords_twlist,lat,ind,ttw);
    
    //    flist[lat][dir][bead].twlist -= flist[lat][dir][bead].wlist[n];
  }
  
  //LINE_STAMP;
  GetFromServer_int(&Coords_nlist,lat,ind,tempi);
  n = tempi[0] - 1;
  wind = mxlist*beads*dir + mxlist*bead + n;
  tw[0] = wt;
  //LINE_STAMP;
  PutToServer_dbl(&Coords_wlist,lat,wind,tw);
  //flist[lat][dir][bead].wlist[n] = wt;

  GetFromServer_dbl(&Coords_twlist,lat,ind,ttw);
  ttw[0] += wt;
  PutToServer_dbl(&Coords_twlist,lat,ind,ttw);

  //flist[lat][dir][bead].twlist += wt;

  mcount = 0;
  for (part=0; part<npart;part++) {
    for (dim=0; dim<=ndim-1; dim++) {
      tpt[mcount] = x.x[part][dim];  // position
      mcount++;
      tpt[mcount] = x.v[part][dim];  // velocity
      mcount++;
    }
  }	    
  
  lo[0] = lat;
  hi[0] = lat;
  lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*n;
  hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(n+1) -1;
  NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);
  //  flist[lat][dir][bead].pts[n] = x; 

  tempi[0] = from1;
  ind = 2*mxlist*beads*dir + 2*mxlist*bead + 2*n;
  //LINE_STAMP;
  PutToServer_int(&Coords_from,lat,ind,tempi);
  // flist[lat][dir][bead].from[n][0] = from1;

  tempi[0] = from2;
  //LINE_STAMP;
  PutToServer_int(&Coords_from,lat,ind+1,tempi);
  // flist[lat][dir][bead].from[n][1] = from2;

  ind = mxlist*beads*dir + mxlist*bead + n;
  tempi[0] = myclock;
  //LINE_STAMP;
  PutToServer_int(&Coords_clock,lat,ind,tempi);
  tempi[0] = mytim;
  PutToServer_int(&Coords_mytim,lat,ind,tempi);
  // flist[lat][dir][bead].clock[n] = myclock;
  
  ind = dir*beads + bead;
  //LINE_STAMP;
  GetFromServer_int(&Coords_on,lat,ind,tempi);
  
  if (!tempi[0]) {
    printf("proc %d:  Turning %d %d %d back on.\n", me, lat, dir, bead);
    tempi[0] = 1;
    //LINE_STAMP;
    PutToServer_int(&Coords_on,lat,ind,tempi);
  }

  GA_Fence();
  GA_Unlock(lockind);

 }
/*-------------------------------------------------------------------------*/
opoint_t ptavg(opoint_t x, opoint_t y) {

  /* this function returns the average of two points in order parameter space */

  opoint_t temp;
  int i;

  for (i=0;i<=nop-1;i++) {
    temp.x[i] = 0.5*(x.x[i] + y.x[i]);
  }

  return temp;
}
/*-------------------------------------------------------------------------*/
 double calcsum(int lat, int dir) { 

   double temp;
   int lim, i;

   temp = 0.;
   if (lat == 1) {
     lim = beads;
   } else {
     lim = beads-1;
   }
   for (i=0;i<=lim-1;i++) {
     temp = temp + weight[lat][dir][i];
   }
   
   return temp;
 }
/*-------------------------------------------------------------------------*/
 void stack(int lat, int dir, int bead) { 
   
   /* this function adds the current state of the replicas to stravg */
   
   int i, j, index, ind;
   double temps[1];
   point_t rt;

   if (phase==1) {
     if (lat == 0) {
   
       allind[0] = dir;
       allind[1] = bead;
       tempi[0] = 1 + NGA_Read_inc(Coords_navg.ga,allind,1);
       
       lockind = getlockind(3,lat,dir,bead);
       GA_Lock(lockind);
       GA_Init_fence();
       
       for (j=0;j<=nop-1;j++) {
	 index = bead*nop + j;
	 
	 GetFromServer_dbl(&Coords_stravg,lat,index,temps);
	 temps[0] += ocoor.x[j];
	 PutToServer_dbl(&Coords_stravg,lat,index,temps);
       }
       GA_Fence();
       GA_Unlock(lockind);
     }
   }

   index = tres*(ocoor.x[0]-tmin)/(tmax-tmin);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index! %f\n",ocoor.x[0]);
   } else {
     allind[0] = lat;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_thist.ga,allind,1);
   }
 }
/*-------------------------------------------------------------------------*/
void string() {

  /* this function uses the average statistics of the beads to move the 
     string images */
  
  /* move zs (1st lattice, both directions, not the endpoints) */
  
  int lat, dir, i, j, n, bead, last, next, index, mcount, index2;
  double temp;

  /* get stravg from server */

  //LINE_STAMP;
  GetFromServer_dbl(&Coords_stravg,0,-1,tstravg1);
  //LINE_STAMP;
  GetFromServer_dbl(&Coords_stravg,1,-1,tstravg2);
  //LINE_STAMP;
  GetFromServer_int(&Coords_navg,0,-1,tint1);
  //LINE_STAMP;
  GetFromServer_int(&Coords_navg,1,-1,tint2);

  index = 0;
  index2 = 0;
  for (i=0;i<beads;i++) {
    navg[0][i] = tint1[index2];
    navg[1][i] = tint2[index2];
    index2++;
    for (j=0;j<nop;j++) {
      stravg[0][i].x[j] = tstravg1[index];
      stravg[1][i].x[j] = tstravg2[index];
      index++;
    }
  } 

  lat = 0;
  for (dir=0;dir<=1;dir++) {
    for (i=1;i<=beads-2;i++) {
      for (j=0;j<=nop-1;j++) {
	z[lat][dir][i].x[j] += (stravg[dir][i].x[j]/(float)navg[dir][i] - z[lat][dir][i].x[j])*frac;
      }
    }
  }
  
  /* move 1st lattice endpoints */
  
  n = beads-1;
  for (j=0;j<=nop-1;j++) {
    z[lat][0][0].x[j] += (0.5*(stravg[0][0].x[j]/(float)navg[0][0]+stravg[1][0].x[j]/(float)navg[1][0])-z[lat][0][0].x[j])*frac;
    z[lat][0][n].x[j] += (0.5*(stravg[0][n].x[j]/(float)navg[0][n]+stravg[1][n].x[j]/(float)navg[1][n])-z[lat][0][n].x[j])*frac;
  }
  z[lat][1][0] = z[lat][0][0];
  z[lat][1][n] = z[lat][0][n];

  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      for (j=0;j<=nop-1;j++){
	stravg[dir][bead].x[j] = 0.;
      }
      navg[dir][bead] = 0;
    }
  }

  GA_Zero(Coords_stravg.ga);
  GA_Zero(Coords_navg.ga);

  /* smooth */
  
  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"presmooth:\n");
    wrstring();
    fclose(mainout);
  }

  for (dir=0;dir<=1;dir++) {
    for (bead=1;bead<=beads-2;bead++) {
      last = bead-1;
      next = bead+1;
      for (j=0;j<=nop-1;j++) {
	z[0][dir][bead].x[j] += kappa*(z[0][dir][next].x[j]+z[0][dir][last].x[j]-2.*z[0][dir][bead].x[j]);
      }
    }
  }
  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"prerepar:\n");
    wrstring();
    fclose(mainout);
  }

  repar(0,0);
  repar(0,1);

  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"postrepar:\n");
    wrstring();
    fclose(mainout);
  }

  for (dir=0;dir<=1;dir++) {
    for (i=1;i<=beads-1;i++) {
      z[1][dir][i-1] = ptavg(z[0][dir][i],z[0][dir][i-1]);      
    }
  }  

  /* send string to server */

  mcount = 0;
  for (dir=0;dir<=1;dir++) {
    for (i=0;i<beads;i++) {
      for (j=0;j<nop;j++) {
	tz1[mcount] = z[0][dir][i].x[j];
	tz2[mcount] = z[1][dir][i].x[j];
	mcount++;
      }
    }
  }
  //LINE_STAMP;
  PutToServer_dbl(&Coords_z,0,-1,tz1);
  //LINE_STAMP;
  PutToServer_dbl(&Coords_z,1,-1,tz2);
}
/*-------------------------------------------------------------------------*/
void wrstring() {
  
  int j, k;

  for (j=0;j<=1;j++) {
    for (k=0;k<=beads-1;k++) {
      fprintf(mainout,"%e %e\n",z[0][j][k].x[0],z[0][j][k].x[1]);
    }
  }

}
/*-------------------------------------------------------------------------*/
void wrpath(int ind) {

  int j, k, op;

  sprintf(fname,"%s%d%s","path", ind, "_1.dat");
  outPath = fopen(fname,"w");
  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      for (op=0;op<=nop-1;op++) {
	fprintf(outPath,"%e ",z[0][j][k].x[op]);
      }
      fprintf(outPath,"\n");
    }
  }
  fclose(outPath);

  if (verb >=2) {
    sprintf(fname,"%s%d%s","path", ind, "_2.dat");
    outPath = fopen(fname,"w");
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-2; k++) {
	for (op=0;op<=nop-1;op++) {
	  fprintf(outPath,"%e ",z[1][0][k].x[op]);
	}
	fprintf(outPath,"\n");
      }
    }
    fclose(outPath); 
  }
}
/*-------------------------------------------------------------------------*/
void repar(int lat, int dir) {

  opoint_t oldz[beads];
  int i, last, next, j, test;
  double sum, dist[beads-1], even, targt;

  for (i=0;i<=beads-1;i++) {
    oldz[i] = z[lat][dir][i];
  }

  sum = 0.;
  for(i=0;i<=beads-2;i++) { 
    next = i+1;
    dist[i] = ogtdist(oldz[next],oldz[i]);
    sum += dist[i];
  }
  even = sum/(float)(beads-1);

  for (i=1;i<=beads-2;i++) {
    targt = i*even;

    test = 0;
    sum = 0.;

    while (sum <= targt) {
      test++;

      if (test == beads) {
	printf("proc %d: Error in repar!\n", me);
	gaexit(9);
      }
      sum += dist[test-1];
    }
    next = test;  /* image belongs between test and test-1 */
    test--;
    for (j=0;j<=nop-1;j++) {
      z[lat][dir][i].x[j] = oldz[next].x[j] + (oldz[test].x[j]-oldz[next].x[j])*(sum-targt)/dist[test];
    }
  }
}
/*-------------------------------------------------------------------------*/
void nobad(int lat, int dir, int bead) {

  /* this function looks through the flux list to see if any flux entry
     points are no longer in the Voronoi regions */

  int i, j, k, l, limit, lb, ub, chkall, tfull,ptind, part, dim, ind, mcount;
  point_t x;

  ind = dir*beads + bead;
  //LINE_STAMP;
  GetFromServer_int(&Coords_nlist,lat,ind,tempi);
  globn = tempi[0];
  //LINE_STAMP;
  GetFromServer_int(&Coords_full,lat,ind,tempi);
  tfull = tempi[0];

  if ((globn !=0) || (tfull)) {
    if (tfull) {
      limit = mxlist;
      tempi[0] = mxlist;
      //LINE_STAMP;
      PutToServer_int(&Coords_nlist,lat,ind,tempi);
      globn = mxlist;
      // flist[lat][dir][bead].nlist = mxlist;

    }else{
      limit = globn;
    }
    
    k = 0;
    chkall = 0;
    while (k < globn) {
      k++;
      lo[0] = lat;
      hi[0] = lat;
      lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(k-1);
      hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*k -1;
      NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);

      mcount = 0;
      for (part=0; part<npart;part++) {
	for (dim=0; dim<=ndim-1; dim++) {
	  x.x[part][dim] = tpt[mcount];
	  mcount++;
	  x.v[part][dim] = tpt[mcount];
	  mcount++;
	}
      } // load up x

      if (whichx(lat,dir,x) != bead) {
	if (tfull) {
	  tfull = 0;
	  tempi[0] = 0;
	  //LINE_STAMP;
	  PutToServer_int(&Coords_full,lat,ind,tempi);
	}
	delpoint(lat,dir,bead,k-1);
	k--;       /* recheck kth element */
      }
    }

    // publish globn
    tempi[0] = globn;
    //LINE_STAMP;
    PutToServer_int(&Coords_nlist,lat,ind,tempi);
  }
}
/*-------------------------------------------------------------------------*/
void delpoint(int lat, int dir, int bead, int pt) {
  
  /* this function deletes the 'pt' point from the (lat,dir,bead) fluxlist */
  
  int n, i, part,ind,wind;
  int ptind, ptind2, dim, mcount;
  double tw[1], ttw[1];
  
  n = globn - 1;
  ind = beads*dir + bead;
  if (pt != n) {
    wind = mxlist*beads*dir + mxlist*bead + pt;
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_wlist,lat,wind,tw);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_twlist,lat,ind,ttw);
    ttw[0] -= tw[0];
    //LINE_STAMP;
    PutToServer_dbl(&Coords_twlist,lat,ind,ttw);
    //flist[lat][dir][bead].twlist -= flist[lat][dir][bead].wlist[pt];

    wind = mxlist*beads*dir + mxlist*bead + n;
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_wlist,lat,wind,tw);
    wind = mxlist*beads*dir + mxlist*bead + pt;
    //LINE_STAMP;
    PutToServer_dbl(&Coords_wlist,lat,wind,tw);
    //flist[lat][dir][bead].wlist[pt] = flist[lat][dir][bead].wlist[n]; /* replace kth element with last element */

    tw[0] = 0.;
    wind = mxlist*beads*dir + mxlist*bead + n;
    //LINE_STAMP;
    PutToServer_dbl(&Coords_wlist,lat,wind,tw);
    //flist[lat][dir][bead].wlist[n] = 0.;

    lo[0] = lat;
    hi[0] = lat;
    lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*n;
    hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(n+1) -1;
    NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);
    lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*pt;
    hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(pt+1) -1;
    NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);
    
    mcount = 0;
    for (part=0; part<npart;part++) {
      for (dim=0; dim<=ndim-1; dim++) {
	tpt[mcount] = 0.;
	mcount++;
	tpt[mcount] = 0.;
	mcount++;
      }
    }
    lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*n;
    hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(n+1) -1;
    NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);

    //flist[lat][dir][bead].pts[pt] = flist[lat][dir][bead].pts[n];  then set .pts[n] to zero

    wind = 2*mxlist*beads*dir + 2*mxlist*bead + 2*n;
    GetFromServer_int(&Coords_from,lat,wind,tempi);
    GetFromServer_int(&Coords_from,lat,wind+1,tempi2);
    wind = 2*mxlist*beads*dir + 2*mxlist*bead + 2*pt;
    PutToServer_int(&Coords_from,lat,wind,tempi);
    PutToServer_int(&Coords_from,lat,wind+1,tempi2);
    tempi[0] = 0;
    wind = 2*mxlist*beads*dir + 2*mxlist*bead + 2*n;
    PutToServer_int(&Coords_from,lat,wind,tempi);
    PutToServer_int(&Coords_from,lat,wind+1,tempi);
    // flist[lat][dir][bead].from[pt][0] = flist[lat][dir][bead].from[n][0];  then set [n] to zero

    globn--;
    // flist[lat][dir][bead].nlist--;

  } else {

    wind = mxlist*beads*dir + mxlist*bead + n;
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_wlist,lat,wind,tw);
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_twlist,lat,ind,ttw);
    ttw[0] -= tw[0];
    PutToServer_dbl(&Coords_twlist,lat,ind,ttw);

    // flist[lat][dir][bead].twlist -= flist[lat][dir][bead].wlist[n];

    tw[0] = 0.;
    //LINE_STAMP;
    PutToServer_dbl(&Coords_wlist,lat,wind,tw);
    // flist[lat][dir][bead].wlist[n] = 0.;
    
    mcount = 0;
    for (part=0; part<npart;part++) {
      for (dim=0; dim<=ndim-1; dim++) {
	tpt[mcount] = 0.;
	mcount++;
	tpt[mcount] = 0.;
	mcount++;
      }
    }
    lo[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*n;
    hi[1] = 2*npart*ndim*mxlist*beads*dir + 2*npart*ndim*mxlist*bead + 2*npart*ndim*(n+1) -1;
    NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);
    //set point to zero
    
    globn--;
    // flist[lat][dir][bead].nlist--;
  }
}
/*-------------------------------------------------------------------------*/
void getfrac(int lat, int dir, int bead) {

  /* this function saves the amount of weight on the (lat,dir,bead) 
     flux list that came from region (olat,dir2,bead2) [in the lfrac array],
     for all dir2 and bead2 */

  double temp, *tmpwt;
  int i, j, ind, tfull, wind, *tmpfrom, mcount, mcount2, limit;

  ind = beads*dir + bead;
  //LINE_STAMP;
  GetFromServer_int(&Coords_full,lat,ind,tempi);
  tfull = tempi[0];
  if (tfull) {
    limit = mxlist;
  } else {
    //LINE_STAMP;
    GetFromServer_int(&Coords_nlist,lat,ind,tempi);
    limit = tempi[0];
  }
  mcount = 4*beadsp2*lat + 2*beadsp2*dir + 2*beads*bead;
  for (i=0;i<2;i++) {
    for (j=0;j<beads;j++) {
      lfrac[mcount] = 0.;
      mcount++;
    }
  }

  tmpfrom = (int *) malloc(2*limit*sizeof(int));
  tmpwt = (double *) malloc(limit*sizeof(double));

  lo[0] = lat;
  hi[0] = lat;
  lo[1] = 2*mxlist*beads*dir + 2*mxlist*bead;
  hi[1] = lo[1] + 2*limit - 1;
  NGA_Get(Coords_from.ga,lo,hi,tmpfrom,ld);

  lo[1] = mxlist*beads*dir + mxlist*bead;
  hi[1] = lo[1] + limit - 1;
  NGA_Get(Coords_wlist.ga,lo,hi,tmpwt,ld);

  mcount = 0;
  for (i=0;i<=limit-1;i++) {
    tempi[0] = tmpfrom[mcount];
    mcount++;
    tempi2[0] = tmpfrom[mcount];
    mcount++;

    mcount2 = 4*beadsp2*lat + 2*beadsp2*dir + 2*beads*bead + beads*tempi[0] + tempi2[0];
    lfrac[mcount2] += tmpwt[i];
  }

  free(tmpfrom);
  free(tmpwt);

}
/*-------------------------------------------------------------------------*/
int getmin(int *lat, int *dir, int *bead) {
  
  /* this function returns the (lat,dir,bead) for the region with the 
     lowest counter */
  
  int low, i, j, k, mcount, ind, start, limit, mind;
  double rnd;

  //LINE_STAMP;
  GetFromServer_int(&Coords_count,0,-1,tint1);
  //LINE_STAMP;
  GetFromServer_int(&Coords_count,1,-1,tint2);
  //LINE_STAMP;
  GetFromServer_int(&Coords_on,0,-1,ton1);
  //LINE_STAMP;
  GetFromServer_int(&Coords_on,1,-1,ton2);


  low = every;
  *lat = -1;
  *dir = 0;
  *bead = 0;

  limit = 4*beads-2;
  start = limit*(me+1)/nproc;
  
  for(ind=0;ind<limit;ind++) {
    mcount = ind + start;
    if (mcount >= limit) {
      mcount -= limit;
    }
    i = mcount/(2*beads);
    if (i==0) {
      j = mcount/beads;
    } else {
      j = (mcount-(2*beads))/(beads-1);
    }
    if (i==0) {
      k = mcount - j*beads;
    } else {
      k = mcount-(2*beads) - j*(beads-1);
    }
    mind = j*beads + k;

    if (i == 0) {
      count[i][j][k] = tint1[mind];
      on[i][j][k] = ton1[mind];
    } else {
      count[i][j][k] = tint2[mind];
      on[i][j][k] = ton2[mind];
    }
    if (on[i][j][k]) {
      if (count[i][j][k] < low) {
	low = count[i][j][k];
	*lat = i;
	*dir = j;
	*bead = k;
      }
    }
  }
  
  if (low >= every) {
    *lat = -1;           /* all regions have run for 'every' steps */
  }

  return low;
}
/*-------------------------------------------------------------------------*/
void wtread() {

  /* this function reads the initial values of the weights
     from the input file 'wname' */

  int i, j, k, op, dir, mcount;
  FILE * wtin;
  float temp;
  
  wtin = fopen(wname,"r");

  for (i=0;i<=1;i++) {
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-1; k++) {
	fscanf(wtin,"%e",&temp);
	weight[i][j][k] = temp;
      }
    }
  }
  fclose(wtin);

  mcount = 0;
  for (j=0;j<2;j++) {
    for (k=0;k<beads;k++) {
      tweight1[mcount] = weight[0][j][k];
      tweight2[mcount] = weight[1][j][k];
      twavg1[mcount] = 0.;
      mcount++;
    }
  }
  //LINE_STAMP;
  PutToServer_dbl(&Coords_weight,0,-1,tweight1);
  //LINE_STAMP;
  PutToServer_dbl(&Coords_weight,1,-1,tweight2);
  //LINE_STAMP;
  PutToServer_dbl(&Coords_wadj,0,-1,tweight1);
  //LINE_STAMP;
  PutToServer_dbl(&Coords_wadj,1,-1,tweight2);
  //LINE_STAMP;
  PutToServer_dbl(&Coords_wavg,0,-1,twavg1);
  //LINE_STAMP;
  PutToServer_dbl(&Coords_wavg,1,-1,twavg1);

}
/*-------------------------------------------------------------------------
--------------------------------------------
the following functions are system-dependent
--------------------------------------------
-------------------------------------------------------------------------*/

int strinit() {

  /* this function initializes the string */

  int lat, dir, bead, op, i;

  /* string 1: forward & backward */

  for (i=0;i<=beads-1;i++) {
    z[0][0][i].x[0] = bmin + bwidth + (bmax-bmin-2*bwidth)*(i)/(float)(beads-1);
    z[0][1][i].x[0] = z[0][0][i].x[0];
  }

  for (dir=0;dir<=1;dir++) {
    for (i=1;i<=beads-1;i++) {
      z[1][dir][i-1] = ptavg(z[0][dir][i],z[0][dir][i-1]);
    }
  }

  return 0;
}
/*-------------------------------------------------------------------------*/
void strread() {

  /* this function reads the initial string from the input
     file 'flname' */
  printf("Error! strread not written yet!\n");
  gaexit(2);
}
/*-------------------------------------------------------------------------*/
void fxread() {

  /* this function reads the initial values of the fluxes
     from the input file 'fxname' */

  int i, j, k, op, dir, ti, tj, tk, limit, pt, dim, part, ind, mcount;
  FILE * fxin;
  float temp;
  float tw;
  
  fxin = fopen(fxname,"r");

  for (i=0;i<=1;i++) {
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-1; k++) {
	if ((i == 0) || (k!=beads-1)) {
	  fscanf(fxin,"%d %d %d",&ti,&tj,&tk);
	  if (((ti != i) || (tj != j)) || (tk != k)) {
	    printf("Error!\n");
	    gaexit(2);
	  }
	  
	  fscanf(fxin,"%d %f", &tempi[0], &tw);
	  tempw[0] = tw;
	  ind = j*beads + k;
	  //LINE_STAMP;
	  PutToServer_int(&Coords_nlist,i,ind,tempi);
	  //LINE_STAMP;
	  PutToServer_dbl(&Coords_twlist,i,ind,tempw);
	  
	  if (tempi[0] == mxlist) {
	    tempi2[0] = 1;
	  } else {
	    tempi2[0] = 0;
	  }
	  //LINE_STAMP;
	  PutToServer_int(&Coords_full,i,ind,tempi2);
	  limit = tempi[0];
	  
	  for (pt=0; pt<=limit-1; pt++) {
	    fscanf(fxin,"%d %d ", &tempi[0], &tempi2[0]);
	    ind = 2*mxlist*beads*j + 2*mxlist*k + 2*pt;
	    //LINE_STAMP;
	    PutToServer_int(&Coords_from,i,ind,tempi);
	    //LINE_STAMP;
	    PutToServer_int(&Coords_from,i,ind+1,tempi2);
	    
	    fscanf(fxin,"%f %d ",&tw,&tempi[0]);
	    tempw[0] = tw;
	    ind = mxlist*beads*j + mxlist*k + pt;
	    //LINE_STAMP;
	    PutToServer_dbl(&Coords_wlist,i,ind,tempw);
	    //LINE_STAMP;
	    PutToServer_int(&Coords_clock,i,ind,tempi);
	    fscanf(fxin,"%d ",&tempi[0]);
	    PutToServer_int(&Coords_mytim,i,ind,tempi);
	    
	    mcount = 0;
	    for (part=0; part<npart;part++) {
	      for (dim=0; dim<=ndim-1; dim++) {
		fscanf(fxin,"%f ", &tw);
		tpt[mcount] = tw;
		mcount++;
		fscanf(fxin,"%f ", &tw);
		tpt[mcount] = tw;
		mcount++;
	      }
	    }	    

	    lo[0] = i;
	    hi[0] = i;
	    lo[1] = 2*npart*ndim*mxlist*beads*j + 2*npart*ndim*mxlist*k + 2*npart*ndim*pt;
	    hi[1] = 2*npart*ndim*mxlist*beads*j + 2*npart*ndim*mxlist*k + 2*npart*ndim*(pt+1) -1;
	    NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);
	  }
	}
      }
    }
  }
  fclose(fxin);
}
/*-------------------------------------------------------------------------*/
void probread() {

  /* this function reads the initial values of prob and nprob
     from the input file 'pname' */

  int i, j, lat, dir, bead, n1, n2, tn, tp, isum;
  FILE * fxin;
  float temp;
  float tw;
  
  fxin = fopen(pname,"r");

  while(fscanf(fxin,"%d %d %d %d %d %d", &(lat), &(dir), &(bead), &(n1), &(n2), &(tn)) != EOF ){
    allind[0] = lat;
    allind[1] = 2*beadsp2*dir + 2*beads*bead + beads*n1 + n2;
    tempi[0] = tn + NGA_Read_inc(Coords_nprob.ga,allind,tn);
    isum = 0;
    while(isum < tn) {
      fscanf(fxin,"%d %d %d", &i, &j, &tp);
      isum += tp;
      allind[0] = lat;
      allind[1] = 4*beadsp3*dir + 4*beadsp2*bead + 2*beadsp2*n1 + 2*beads*n2 + beads*i + j;
      tempi[0] = tp + NGA_Read_inc(Coords_prob.ga,allind,tp);
    }
  }
  fclose(fxin);
}
/*--------------------------------------------------------------------*/
point_t ptconv(int lat, int dir, int bead) {

  /* this function returns a point in the full space that 
     corresponds to z(lat,dir,bead) in o.p. space
     It is used only for restarting the replica when
     no flux points are available */

  int i, j, part, good, odir;
  point_t temp;
  
  good = 0;

  if (dir == 1) {
      odir = 0;
  } else {
      odir = 1;
  }

  for (i=0;i<nemerg;i++) {
    ocoor.x[0] = emerg[i].op;
    if ((which(lat,dir) == bead) && (bascheck() != odir)) {
      for (part=0;part<npart;part++) {
	for (j=0;j<ndim;j++) {
	  temp.x[part][j] = emerg[i].x[part][j];
	  temp.v[part][j] = emerg[i].v[part][j];
	}
      }
      good = 1;
      break;
    }
  }

  //printf("emerg[%d].op = %f\n",i,emerg[i].op);
  if (!good) {
    printf("ptconv couldn't find a point for (%d,%d,%d) = %e\n",lat,dir,bead,z[lat][dir][bead].x[0]);
    gaexit(2);
  }

  return temp;
}
/*-------------------------------------------------------------------------*/

int bascheck() {
  
  /* this function checks if the walker is in a basin */

  int temp, i;
  double dist;

  temp = -1;

  for (i=0;i<=1;i++) {
    dist = ogtdist(ocoor,basin[i]);
    if (dist < bwidth) {
      temp = i;
    }
  }  

  return temp;
}
/*-------------------------------------------------------------------------*/
opoint_t projx(point_t x) {
  
  /* this function returns the projection of coor onto the order parameter space */
  
  int i;
  opoint_t temp;
  double mn, mx;
  
  for (i=0;i<npart;i++) {
    if (i==0) {
      mn = x.x[i][0];
      mx = x.x[i][0];
    } else {
      if (x.x[i][0] < mn) {
	mn = x.x[i][0];
      }
      if (x.x[i][0] > mx) {
	mx = x.x[i][0];
      }
    }
  }

  temp.x[0] = mx-mn;
  
  return temp;
}
/*-------------------------------------------------------------------------*/
void move() {
  
  /* moves coor, updates ocoor */

  int i;
  point_t oldx, rt;
  double psi[2],phi[2],fx,fy,temp1,temp2,dvdx,dvdy;
  void dostep(), stream(), dorotate();
  
  dostep();

  if (myclock%strmfrq == 0) {
    stream();
    dorotate();
  }
  ocoor = projx(coor);  /* update o.p. projection */
}
/*-------------------------------------------------------------------------*/
void dostep(){
  int i, j, k, x;
  double a, f2[npart][ndim],temp;
  point_t f, force();
    
#pragma omp parallel for default(shared) private(i,x,a) schedule(static)
  for(i=0; i< npart; i++){
    for(x=0; x< 3; x++){	
      /* calculate a(t) = F/m */        
      // don't divide by polymass though, since the internal polymer forces are scaled assuming m = 1
      coor.x[i][x] += 0.5*tstepsq*f1[i][x] + tstep*coor.v[i][x];
    }
  }
  
#pragma omp parallel for default(shared) private(i,j,k,temp) schedule(static)
  for(i=0; i< npart; i++){
    for(j=i; j< npart; j++){
      temp = 0.;
      for (k=0;k<3;k++) {
	temp += (coor.x[i][k] - coor.x[j][k])*(coor.x[i][k] - coor.x[j][k]);
      }
      rdists[i][j] = sqrt(temp);
      rdists[j][i] = rdists[i][j];
      
      for (k=0;k<ndim;k++) {
	rvecs[i][j].x[k] = coor.x[i][k]-coor.x[j][k];
	rvecs[j][i].x[k] = -rvecs[i][j].x[k];
      }
    }
  }
  
  f = force();

#pragma omp parallel for default(shared) private(i,x) schedule(static)
  for(i=0; i< npart; i++){
    for(x=0; x< 3; x++){
      f2[i][x] = f.x[i][x];
      coor.v[i][x] += 0.5*(f1[i][x]+f2[i][x])*tstep;
      f1[i][x] = f2[i][x];
    }
  }
}
/*-------------------------------------------------------------------------*/
point_t force(){
  int i, j, p1, p2;
  xpoint_t f_1, vcmult(double a, xpoint_t x), vplus(xpoint_t x, xpoint_t y), trvec;
  double dvdr, tempr, pre;
  point_t f;

#pragma omp parallel default(shared) private(i,j,p1,p2,pre,f_1,tempr,trvec)
  {
#pragma omp for schedule(static)
      for (i=0;i<npart;i++) {
	  for (j=0;j<3;j++) {
	      f.x[i][j] = 0.;
	  }
      }

#pragma omp for schedule(dynamic)
      for (p1=0;p1<npart-1;p1++) {
	  for (p2=p1+1;p2<npart;p2++) {
	      
	      for (j=0;j<3;j++) {
		  f_1.x[j] = 0.;
	      }
	      
	      /* compute force from p2 acting on p1 */
	      
	      if (p2-p1 >= 3) {
		  
		  /* sec and tert LJ forces */
		  
		  if (Delta[p1][p2][0] == 1 ){  // seclist
		      pre = eph1_12*( r0d12[p1][p2]/mypow14(rdists[p1][p2]) - r0d6[p1][p2]/mypow8(rdists[p1][p2]));
		      f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
		  }
		  if (Delta[p1][p2][1] == 1 ){  // terlist
		      pre = eph2_12*( r0d12[p1][p2]/mypow14(rdists[p1][p2]) - r0d6[p1][p2]/mypow8(rdists[p1][p2]));
		      f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
		  }
		  
		  /*Repulsive lj forces- 3 or more away (non-native)*/
		  
		  if( cutoff > 0 ){  // !seclist and !terlist
		      if( rdists[p1][p2] < cutoff*sigma1 && Delta[p1][p2][0] == 0 && Delta[p1][p2][1] == 0){
			  pre = ep1_6*sigma16/mypow8(rdists[p1][p2]);
			  f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
		      }
		  } else {
		      if( Delta[p1][p2][0] == 0 && Delta[p1][p2][1] == 0){
			  pre = ep1_6*sigma16/mypow8(rdists[p1][p2]);
			  f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
		      }
		  }
	      }

	      /* repulsive LJ forces - 2 away */

	      if (p2-p1 == 2) {
		  pre = ep1_6*sigma26/mypow8(rdists[p1][p2]);
		  f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
	      }
	      
	      /* connection force */
	      
	      if (p2-p1 == 1) {
		  if(abs(rdists[p1][p2]-r0dists[p1][p2]) < R0) {  
		      pre = -kf*(rdists[p1][p2] - r0dists[p1][p2])/(rdists[p1][p2]*(1.-mypow2(rdists[p1][p2]-r0dists[p1][p2])*ooR02));
		      f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
		  } else {
		      pre = -kf*(rdists[p1][p2]-r0dists[p1][p2])/rdists[p1][p2];
		      f_1 = vplus(f_1,vcmult(pre,rvecs[p1][p2]));
		  }
	      }
	      
	      /* add forces to f matrix */
	      
	      for (j=0;j<3;j++) {
		  f.x[p1][j] += f_1.x[j];
		  f.x[p2][j] -= f_1.x[j];
	      }
	  }
	  
	  if (coor.x[p1][1] < cutoff) { // add LJ repulsion if bead is close to bottom wall
	    f.x[p1][1] += ep1_6*phanr06/mypow7(coor.x[p1][1]);
	  }

	  if (p1 == 0) {  // add connection force to phantom bead
	      tempr = 0.;
	      for (j=0;j<3;j++) {
		  f_1.x[j] = 0.;
		  trvec.x[j] = coor.x[p1][j] - phan[j];
		  tempr += mypow2(trvec.x[j]);
	      }
	      tempr = sqrt(tempr);
	      if(abs(tempr-phanr0) < R0) {  
		  pre = -kf*(tempr - phanr0)/(tempr*(1.-mypow2(tempr-phanr0)*ooR02));
		  f_1 = vplus(f_1,vcmult(pre,trvec));
	      } else {
		  pre = -kf*(tempr-phanr0)/tempr;
		  f_1 = vplus(f_1,vcmult(pre,trvec));
	      }
	      
	      for (j=0;j<3;j++) {
		  f.x[p1][j] += f_1.x[j];
	      }
	  }
      }
  } /* -- end of OpenMP parallel region */
  return f;
}
/*-----------------------------------------------------------------------------*/
xpoint_t vcmult(double a, xpoint_t x) {
  int i;
  xpoint_t temp;

  for (i=0;i<3;i++) {
    temp.x[i] = x.x[i]*a;
  }
  return temp;

}
/*-----------------------------------------------------------------------------*/
xpoint_t vplus(xpoint_t x, xpoint_t y) {
  int i;
  xpoint_t temp;

  for (i=0;i<3;i++) {
    temp.x[i] = x.x[i]+y.x[i];
  }
  return temp;
}
/*-----------------------------------------------------------------------------*/
void rdcryst() {

  FILE *cfile;
  int i, j;
  double val,dx,dy,dz;

  cfile= fopen(crname,"rt");
  i=0;
  while(fscanf(cfile, "%lf", &(val) ) != EOF ){
    if( i >= npart*3 ){ printf("cryst structure too big?\n"); gaexit(1);}
    r0[i/3][i%3] = val;
    i++;
  }
  
  for(i=0; i< npart; i++){
    for(j=0; j< npart; j++){
	dx = r0[i][0]-r0[j][0];
	dy = r0[i][1]-r0[j][1];
	dz = r0[i][2]-r0[j][2];
	r0dists[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
	r0d6[i][j] = mypow2(mypow2(r0dists[i][j]))*mypow2(r0dists[i][j]);
	r0d12[i][j] = mypow2(r0d6[i][j]);
    }
  }
}
/*-----------------------------------------------------------------------------*/
void rddelta() {

  /*---------Read in Delta values, or write them out----------------------*/
  /*secondary contacts; if blank, none read*/ 

  FILE *dfile;
  int i, j, val2;
  
  dfile = fopen(secname,"rt");
  if(dfile == NULL){
    printf("No secondary contacts read, setting them to zero\n");
    for(i=0; i< npart; i++){
      for(j=0; j< npart; j++){
	Delta[i][j][0] = 0;
      }
    }
  }else{
    i=0;
    while(fscanf(dfile, "%i", &(val2) ) != EOF ){
      i++;
    }
    if( i != npart*npart ){
      printf("Firt contact list is different size than polymer: %i %i\n", i, npart*npart);
      gaexit(0);
    }
    fclose(dfile);
    dfile=fopen(secname,"rt");
    i=0;
    while(fscanf(dfile, "%i", &(val2) ) != EOF){
      Delta[i/npart][i%npart][0] = (int) val2;
      i++;
    }
  }
  
  /*	tertiary contacts, if blank, none read*/
  dfile = fopen(tertname,"rt");
  if(dfile == NULL){
    printf("No tertiary concacts read, setting them to zero\n");	
    for(i=0; i< npart; i++){
      for(j=0; j< npart; j++){
	Delta[i][j][1] = 0; 
      }
    }
  }else{
    i=0;	
    while(fscanf(dfile, "%i", &(val2) ) != EOF ){		
      i++;								
    }		
    if( i != npart*npart ){			
      printf("Second contact list is different size than polymer: %i %i\n", i, npart*npart);
      gaexit(0);
    }
    fclose(dfile);
    dfile = fopen(tertname,"rt");			
    i=0;							
    while(fscanf(dfile, "%i", &(val2) ) != EOF){
      Delta[i/npart][i%npart][1] = (int) val2;
      i++;			
    }
  }
}
/*-----------------------------------------------------------------------------*/
void rdsolv(){

  FILE *filein;
  int j, k, good, ind;
  double temp, nudge=1.e-4, rx;
  char stemp[30];

  filein = fopen(sname,"r");

  for (j=0; j<N; j++) {
    for (k=0; k<3; k++) {
      fscanf(filein,"%s",&stemp);
      rsolv[j][k] = atof(stemp);
    }
    for (k=0; k<3; k++) {
      fscanf(filein,"%s",&stemp);
      vsolv[j][k] = atof(stemp);
    }
    good = 0;
    while (!good) {
      good = 1;
      if(rsolv[j][0]<0){
	rsolv[j][0]+=Lx;
	good = 0;
      } else if(rsolv[j][0]>Lx){
	rsolv[j][0]-=Lx;
	good = 0;
      }
      if(rsolv[j][1]<0){
	rsolv[j][1]+=Ly;
	good = 0;
      } else if (rsolv[j][1]>Ly){
	rsolv[j][1]-=Ly;
	good = 0;
      }
      if (rsolv[j][1] == 0) {
	rsolv[j][1] += nudge;
      } else if(rsolv[j][1] == Ly) {
	rsolv[j][1] -= nudge;
      }
      if(rsolv[j][2]<0){
	rsolv[j][2]+=Lz;
	good = 0;
      } else if(rsolv[j][2]>Lz){
	rsolv[j][2]-=Lz;
	good = 0;
      }
    }
  }
  fclose(filein);
}
/*-----------------------------------------------------------------------------*/
void rdemerg(){

  FILE *filein;
  int j, part, dim;
  char stemp[30];

  filein = fopen(emname,"r");

  for (j=0;j<nemerg;j++) {
    fscanf(filein,"%s",&stemp);
    emerg[j].op = atof(stemp);
    for (part=0;part<npart;part++) {
      for (dim=0;dim<ndim;dim++) {
	fscanf(filein,"%s",&stemp);
	emerg[j].x[part][dim] = atof(stemp);
      }
      for (dim=0;dim<ndim;dim++) {
	emerg[j].v[part][dim] = sqrt(-var2*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
      }
    }
  }
  fclose(filein);
}
/*-----------------------------------------------------------------------------*/
void stream() {
  
  /* This function streams the solvent particles forward */
  
  int i, j, good;
  double temp3;
  
  for(i=0;i<N;i++){
    temp3=rsolv[i][1]+vsolv[i][1]*solvstep;         /*first move along y coordinate as it has fixed wall*/
    if(temp3<0.){
      for(j=0;j<3;j++){
	vsolv[i][j]=-vsolv[i][j];
      }
    } else if (temp3 > Ly) {
      vsolv[i][1]=-vsolv[i][1];
    } else {
      for(j=0;j<3;j++){
	rsolv[i][j]+=vsolv[i][j]*solvstep;
      }
    }
    good = 0;
    while (!good) {
      good = 1;
      if(rsolv[i][0]<0){
	rsolv[i][0]+=Lx;
	good = 0;
      } else if(rsolv[i][0]>Lx){
	rsolv[i][0]-=Lx;
	vsolv[i][0] += grav*solvstep;
	good = 0;
      }
      if((rsolv[i][1]<0) || (rsolv[i][1]>Ly)){
	printf("Error (Ly): i=%d!\n",i);
	gaexit(2);
      }
      if(rsolv[i][2]<0){
	rsolv[i][2]+=Lz;
	good = 0;
      } else if(rsolv[i][2]>Lz){
	rsolv[i][2]-=Lz;
	good = 0;
      }
    }
  }
}
/*-----------------------------------------------------------------------------*/
void dorotate() {

    double deltax[ndim], rho, psi, temp1, rndvec[ndim], odisp[ndim], rnd[3];
    int i, j, k, l, m, n, good, rr, ind, ii, err, errj, errk, errl, erri;
    double temp4, lim, lim2;
    xpoint_t distchg(double odisp[3], double rndvec[3], int rr), ndisp;
    FILE *filein;

#pragma omp parallel for default(shared) private(i) schedule(static)
    for (i=0;i<3;i++) {
	deltax[i] = ran2(&seed);
    }
      
    /*----------------------------What box is each particle in?-------------------------------*/

#pragma omp parallel for default(shared) private(j,k,l) schedule(static)
    for (j=0;j<nx;j++) {
	for (k=0;k<ny;k++) {
	    for (l=0;l<nz;l++) {
		nlab[j][k][l] = 0;
	    }
	}
    }
 
    good = 0;
    while (!good) {
	good = 1;
#pragma omp parallel for default(shared) private(i,j,k,l) schedule(static)
	for(i=0;i<N;i++){
	    if (good) {
		if(rsolv[i][0]>deltax[0]){
		    pos[i][0]=(int)(blinv*rsolv[i][0]-deltax[0])+1;
		    if (pos[i][0] == nx) {
			pos[i][0] = 0;
		    }
		} else {
		    pos[i][0]=0;
		}
		j = pos[i][0];
		if(rsolv[i][1]>deltax[1]){
		    pos[i][1]=(int)(blinv*rsolv[i][1]-deltax[1])+1;
		} else {
		    pos[i][1]=0;
		}
		k = pos[i][1];
		if(rsolv[i][2]>deltax[2]){
		    pos[i][2]=(int)(blinv*rsolv[i][2]-deltax[2])+1;
		    if (pos[i][2] == nz) {
			pos[i][2] = 0;
		    }
		} else {
		    pos[i][2]=0;
		}
		l = pos[i][2];
		if (nlab[j][k][l] < maxd) {
		    lab[j][k][l][nlab[j][k][l]] = i;
		    nlab[j][k][l]++;
		} else {
		    good = 0;
		}
	    }
	} /* end of parallel region */
	if (!good) {
	    maxd += 5;
	    free(lab);
	    //printf("Increasing maxd to %d\n",maxd);
	    lab = (int ****) alloc4d(sizeof(int),nx,ny,nz,maxd);
	    for (j=0;j<nx;j++) {
		for (k=0;k<ny;k++) {
		    for (l=0;l<nz;l++) {
			nlab[j][k][l] = 0;
		    }
		}
	    }
	}
    }
    
    /*----------------------------What box is each polymer bead in?-------------------------*/
    
    err = 0;
#pragma omp parallel for default(shared) private(i) schedule(static)
    for(i=0;i<npart;i++){
	if(coor.x[i][0]>deltax[0]){
	    pos_poly[i][0]=(int)(blinv*coor.x[i][0]-deltax[0])+1;
	} else {
	    pos_poly[i][0]=0;
	    if ((pos_poly[i][0] < 0) || (pos_poly[i][0] >= nx)) {
		err = 1;
		erri = i;
	    }
	}
	if(coor.x[i][1]>deltax[1]){
	    pos_poly[i][1]=(int)(blinv*coor.x[i][1]-deltax[1])+1;
	} else {    
	    pos_poly[i][1]=0;
	    if ((pos_poly[i][1] < 0) || (pos_poly[i][1] >=ny)) {
		err = 1;
		erri = i;
	    }
	}
	if(coor.x[i][2]>deltax[2]){
	    pos_poly[i][2]=(int)(blinv*coor.x[i][2]-deltax[2])+1;
	} else {
	    pos_poly[i][2]=0;
	    if ((pos_poly[i][2] < 0) || (pos_poly[i][2] >=nz)) {
		err = 1;
		erri = i;
	    }
	}
    } /* end of parallel region */
    
    if (err) {
	printf("Polymer outside of box! x[%d] = (%e,%e,%e)\n",erri,pos_poly[erri][0],pos_poly[erri][1],pos_poly[erri][2]);	      
	polyerr = 1;
    }

    if (!polyerr) {

    /*-------------------------How many particles are in each box-------------------------------*/
 
	err = 0;
#pragma omp parallel default(shared) private(i,j,k,l)
	{
#pragma omp for schedule(static) nowait
	    for(j=0;j<nx;j++){
		for(k=0;k<ny;k++){
		    for(l=0;l<nz;l++){
			v_cmax[j][k][l]=0.;
			v_cmay[j][k][l]=0.;
			v_cmaz[j][k][l]=0.;
			box[j][k][l]=0;	
			box_poly[j][k][l]=0;
		    }
		}
	    }
	
#pragma omp for schedule(static)
	    for(i=0;i<N;i++){
		j = pos[i][0];
		k = pos[i][1];
		l = pos[i][2];
		if ((((j<0) || (j>=nx)) || ((k<0) || (k>=ny))) || ((l<0) || (l>=nz))) {
		    err = 1;
		}
		if (!err) {
		    box[j][k][l]++;
		}
	    }
	} /* end of parallel region */
	
	if (err) {
	    printf("Solvent outside of box! %d %d %d %d\n",i,j,k,l);
	    gaexit(2);
	}
	
	/*Fill up boxes at y=0 and with fake velocities*/
	
#pragma omp parallel for default(shared) private(i,k,temp4) schedule(static)      
	for(i=0;i<nx;i++){
	    for(k=0;k<nz;k++){
		while ((int)n_av>box[i][0][k]) {
		    box[i][0][k]++;
		    
		    temp4=sqrt(-var*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
		    v_cmax[i][0][k]+=solvmass*temp4;
		    
		    temp4=sqrt(-var*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
		    v_cmay[i][0][k]+=solvmass*temp4;
		    
		    temp4=sqrt(-var*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
		    v_cmaz[i][0][k]+=solvmass*temp4;
		}
	    }
	}
	
	/*---------------------How many polymer beads are in each box------------------------------*/
	
#pragma omp parallel for default(shared) private(i,j,k,l) schedule(static)      
	for(i=0;i<npart;i++){
	    j=pos_poly[i][0];
	    k=pos_poly[i][1];
	    l=pos_poly[i][2];
	    box_poly[j][k][l]++;
	    
	}  /* end of parallel region */
	
    
/*-------------------------Calculate Velocity of Center of Mass--------------------------------*/

#pragma omp parallel for default(shared) private(i,j,k,m) schedule(static)      
	for(m=0;m<N;m++){
	    i=pos[m][0];
	    j=pos[m][1];
	    k=pos[m][2];
	    v_cmax[i][j][k]+=solvmass*vsolv[m][0];
	    v_cmay[i][j][k]+=solvmass*vsolv[m][1];
	    v_cmaz[i][j][k]+=solvmass*vsolv[m][2];
	}
	
#pragma omp parallel for default(shared) private(i,j,k,n) schedule(static)
	for(n=0;n<npart;n++){
	    i= pos_poly[n][0];
	    j= pos_poly[n][1];
	    k= pos_poly[n][2];
	    
	    v_cmax[i][j][k]+= polymass*coor.v[n][0];
	    v_cmay[i][j][k]+= polymass*coor.v[n][1];
	    v_cmaz[i][j][k]+= polymass*coor.v[n][2];
	}
	
#pragma omp parallel for default(shared) private(i,j,k,lim,lim2) schedule(static)      
	for(i=0;i<nx;i++){
	    for(j=0;j<ny;j++){
		for(k=0;k<nz;k++){
		    lim = solvmass*box[i][j][k]+polymass*box_poly[i][j][k];
		    if (lim != 0.) {
			lim2 = 1./lim;
			v_cmx[i][j][k]=v_cmax[i][j][k]*lim2;
			v_cmy[i][j][k]=v_cmay[i][j][k]*lim2;
			v_cmz[i][j][k]=v_cmaz[i][j][k]*lim2;
		    }
		}
	    }
	}
	
/*------------------------------Calculate new velocities------------------------------------------*/	
	
	for(i=0;i<nx;i++){
	    for(j=0;j<ny;j++){
		for(k=0;k<nz;k++){
		    if (box[i][j][k]+box_poly[i][j][k] > 0) {
			
#pragma omp parallel default(shared) private(ii,ind,odisp,ndisp)
			{
#pragma omp for schedule(static)
			    for (ii=0;ii<3;ii++) {
				rnd[ii]=ran2(&seed);
			    }
			    rho = 2.*rnd[0] - 1.;
			    psi=twoPi*rnd[1];               /*rotation angle*/
			    temp1 = sqrt(1.-mypow2(rho));
			    rndvec[0] = cos(psi)*temp1;
			    rndvec[1] = sin(psi)*temp1;
			    rndvec[2] = rho;
			    if (rnd[2] > 0.5) {
				rr = 0;
			    } else {
				rr = 1;
			    }
			    
			    // solvents
			    
#pragma omp for schedule(static)	      
			    for (ii=0;ii<nlab[i][j][k];ii++) {
				ind = lab[i][j][k][ii];
				odisp[0] = vsolv[ind][0] - v_cmx[i][j][k];
				odisp[1] = vsolv[ind][1] - v_cmy[i][j][k];
				odisp[2] = vsolv[ind][2] - v_cmz[i][j][k];
				ndisp = distchg(odisp,rndvec,rr);
				vsolv[ind][0] = v_cmx[i][j][k] + ndisp.x[0];
				vsolv[ind][1] = v_cmy[i][j][k] + ndisp.x[1];
				vsolv[ind][2] = v_cmz[i][j][k] + ndisp.x[2];
			    }
			    
			    // polymer beads
			    
#pragma omp for schedule(static)	      
			    for (ind=0;ind<npart;ind++) {
				if (((pos_poly[ind][0]==i)&&(pos_poly[ind][1]==j))&&(pos_poly[ind][2]==k)) {
				    odisp[0] = coor.v[ind][0] - v_cmx[i][j][k];
				    odisp[1] = coor.v[ind][1] - v_cmy[i][j][k];
				    odisp[2] = coor.v[ind][2] - v_cmz[i][j][k];
				    ndisp = distchg(odisp,rndvec,rr);
				    coor.v[ind][0] = v_cmx[i][j][k] + ndisp.x[0];
				    coor.v[ind][1] = v_cmy[i][j][k] + ndisp.x[1];
				    coor.v[ind][2] = v_cmz[i][j][k] + ndisp.x[2];
				}
			    }
			} /* end of parallel region */
		    }
		}
	    }
	}
    }
}
/*--------------------------------------*/
xpoint_t distchg(double odisp[3], double rndvec[3], int rr) {
  xpoint_t temp, cross(double vperp[3], double rndvec[3]), tempv;
  double t1, dot(double odisp[3], double rndvec[3]), vpar[3], vperp[3];
  double r1, r2, r3 ,r4;
  int i;

  t1 = dot(odisp,rndvec);
  for (i=0;i<3;i++) {
    vpar[i] = rndvec[i]*t1;
    vperp[i] = odisp[i]-vpar[i];
  }
  tempv = cross(vperp,rndvec);
  for(i=0;i<3;i++) {
    temp.x[i] = calpha[rr]*vperp[i] + salpha[rr]*tempv.x[i] + vpar[i];
  }
  return temp;
}
/*--------------------------------------*/
xpoint_t cross(double a[3], double b[3]) {
  xpoint_t temp;
  
  temp.x[0] = a[1]*b[2] - a[2]*b[1];
  temp.x[1] = a[2]*b[0] - a[0]*b[2];
  temp.x[2] = a[0]*b[1] - a[1]*b[0];
  return temp;
}
/*--------------------------------------*/
double dot(double a[3], double b[3]) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
/*---------------------------------------*/
void gaexit(int sig) {
  printf("Error!! Exiting early with signal (%d)\n",sig); fflush(stdout);

  GA_Print_stats();
  GA_Terminate();
  MPI_Finalize();
  exit(sig);
}
/*----------------------Random generator-----------------------------------------------------------*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum){
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if(idum<=0){
    if(-(*idum)<1) *idum=1;
    else *idum=-(*idum);
    idum2=(*idum);
    for(j=NTAB+7;j>=0;j--){
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if(*idum<0) *idum+=IM1;
      if(j<NTAB) iv[j]=*idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if(*idum<0) *idum+=IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2<0) idum2+=IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j]=*idum;
  if(iy<1) iy+=IMM1;
  if((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
}
/*------------------------------------------------------*/
void createservers() {

  int MaxNum = 2*beads;

//  printf("proc %d is creating coord servers...\n",me); fflush(stdout);
  status = CreateCoordServer_dbl(&Coords_weight, &MaxNum, "test");
  status = CreateCoordServer_dbl(&Coords_wadj, &MaxNum, "test");
  status = CreateCoordServer_dbl(&Coords_wavg, &MaxNum, "test");
  status = CreateCoordServer_dbl(&Coords_tau, &MaxNum, "test");
  status = CreateCoordServer_dbl(&Coords_twlist, &MaxNum, "test");
  
  status = CreateCoordServer_int(&Coords_elaps, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_count, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_narr, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_on, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_bdedit, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_ntau, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_full, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nlist, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nwstack, &MaxNum, "test");
  
  MaxNum = beads;
  status = CreateCoordServer_int(&Coords_navg, &MaxNum, "test");
  
  MaxNum = beads*nop;
  status = CreateCoordServer_dbl(&Coords_stravg, &MaxNum, "test");
  
  MaxNum = 2*beads*nop;
  status = CreateCoordServer_dbl(&Coords_z, &MaxNum, "test");
  
  MaxNum = 2*beads*2*beads*2*beads;
  status = CreateCoordServer_int(&Coords_prob, &MaxNum, "test");
  MaxNum = 2*beads*2*beads;
  status = CreateCoordServer_int(&Coords_nprob, &MaxNum, "test");
  MaxNum = 2*beads*tres;
  status = CreateCoordServer_int(&Coords_thist, &MaxNum, "test");
  
  MaxNum = 2*beads*mxlist*2*npart*ndim;
  status = CreateCoordServer_dbl(&Coords_pts, &MaxNum, "test");
  MaxNum = 2*beads*mxlist;
  status = CreateCoordServer_dbl(&Coords_wlist, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_clock, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_mytim, &MaxNum, "test");
  MaxNum = 2*beads*mxlist*2;
  status = CreateCoordServer_int(&Coords_from, &MaxNum, "test");

  MaxNum = 4*2*2*beads;
  status = GA_Create_mutexes(MaxNum);

//  printf("all coord servers created by proc %d\n",me); fflush(stdout);
}

void initvar() {

  int i, j, k, ind, pt, op, mcount, dir1, bead1, dir2, bead2, *tmpfrom;

  /* initialize arrays */

  basin[0].x[0] = bmin;
  basin[1].x[0] = bmax;

  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"Initializing variables...\n");
    fclose(mainout);
  }

  for(i=0;i<npart;i++){
    for(j=0;j<3;j++){
      f1[i][j]=0.;
    }
  }

  if (me ==0) {
    mcount=0;
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-1; k++) {
	tint1[mcount] = 1;
	mcount++;
      }
    }
    //LINE_STAMP;
    PutToServer_int(&Coords_on,0,-1,tint1);
    //LINE_STAMP;
    PutToServer_int(&Coords_on,1,-1,tint1);
      
    tmpfrom = (int *) malloc(2*beads*mxlist*sizeof(int));

    mcount=0;
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-1; k++) {
	for (ind=0; ind<mxlist; ind++) {
	  for (i=0; i<2; i++) {
	    tmpfrom[mcount] = -1;
	    mcount++;
	  }
	}
      }
    }
    PutToServer_int(&Coords_from,0,-1,tmpfrom);
    PutToServer_int(&Coords_from,1,-1,tmpfrom);
    free(tmpfrom);
  }
      
  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      navg[j][k] = 0;
      for (op=0; op<=nop-1;op++) {
	stravg[j][k].x[op] = 0.;     /* super-globals */
      }
      for (i=0; i<=1; i++) {
	for (pt=0; pt<=tres-1; pt++) {
	  thist[i][j][k][pt] = 0;
	}

	for (op=0; op<=nop-1; op++) {
	  z[i][j][k].x[op] = 0.;
	}
	if (wtalg == 2) {
	  tau[i][j][k] = 0.;
	  ntau[i][j][k] = 0;
	}

      }
    }
  }
  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"Done initialization\n");
    fclose(mainout);
  }
}

void getstring() {

  int dir, i, j, mcount;

/* get string from server */

  //LINE_STAMP;
  GetFromServer_dbl(&Coords_z,0,-1,tz1);
  //LINE_STAMP;
  GetFromServer_dbl(&Coords_z,1,-1,tz2);

  mcount = 0;
  for (dir=0;dir<=1;dir++) {
    for (i=0;i<beads;i++) {
      for (j=0;j<nop;j++) {
	z[0][dir][i].x[j] = tz1[mcount];
	z[1][dir][i].x[j] = tz2[mcount];
	mcount++;
	if (verb >= 3) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"%d %d %f %f\n",dir, i,tz1[mcount-1],tz2[mcount-1]);
	  fclose(mainout);
	}
      }
    }
  }
}
int getprob(int lat, int dir, int bead, int i, int j, int n1, int n2) {
  
  int ind;
  int temp;

  if ((((lat > 1) || (dir > 1)) || ((i > 1) || (n1 > 1))) || (((bead >= beads) || (j >= beads)) || (n2 >= beads))) {
    printf("Error in getprob!  Array indices out of bounds\n"); fflush(stdout);
    printf("%d %d %d %d %d %d %d\n",lat,dir,bead,i,j,n1,n2); fflush(stdout);
    gaexit(98);
  }
  ind = 4*beadsp3*dir + 4*beadsp2*bead + 2*beadsp2*i + 2*beads*j + beads*n1 + n2;
  if (lat == 0) {
    temp = tprob1[ind];
  } else {
    temp = tprob2[ind];
  }
  return(temp);
}
int getnprob(int lat, int dir, int bead, int i, int j) {

  int ind;
  int temp;

  if ((((lat > 1) || (dir > 1)) || (i > 1)) || ((bead >= beads) || (j >= beads)))  {
    printf("Error in getnprob!  Array indices out of bounds\n"); fflush(stdout);
    printf("%d %d %d %d %d\n",lat,dir,bead,i,j); fflush(stdout);
    gaexit(97);
  }
  ind = 2*beadsp2*dir + 2*beads*bead + beads*i + j; 

  if (lat == 0) {
    temp = tnprob1[ind];
  } else {
    temp = tnprob2[ind];
  }
  return(temp);
}
void wrxyz(int step) {
  
  /* This function writes a point to the xyz function open in 'xyzout' */

  int j, k;
  
  if (step/xyzfrq != 0) {
    fprintf(xyzout,"\n %d \n",npart);
  }
  for(j=0; j< npart; j++){
    if(j<9){
      fprintf(xyzout," atom00%i", j+1);
      for(k=0; k<3; k++){
	if(coor.x[j][k] >= 0 ){
	  fprintf(xyzout,"   %e", coor.x[j][k]);
	}
	if(coor.x[j][k] < 0 ){
	  fprintf(xyzout,"  %e", coor.x[j][k]);
	}
      }
      fprintf(xyzout,"\n");
    }
    if(j<99 && j >= 9){
      fprintf(xyzout," atom0%i", j+1);
      for(k=0; k<3; k++){
	if(coor.x[j][k] >= 0 ){
	  fprintf(xyzout,"   %e", coor.x[j][k]);
	}
	if(coor.x[j][k] < 0 ){
	  fprintf(xyzout,"  %e", coor.x[j][k]);
	}
      }
      fprintf(xyzout,"\n");
    }
    if(j>=99){
      fprintf(xyzout," atom%i", j+1);
      for(k=0; k<3; k++){
	if(coor.x[j][k] >= 0 ){
	  fprintf(xyzout,"   %e", coor.x[j][k]);
	}
	if(coor.x[j][k] < 0 ){
	  fprintf(xyzout,"  %e", coor.x[j][k]);
	}
      }
      fprintf(xyzout,"\n");
    }
  }
}
double mypow2(double x) {
    return(x*x);
}
double mypow7(double x) {
  double x2, x4;

    x2 = x*x;
    x4 = x2*x2;
    return(x4*x2*x);
}
double mypow8(double x) {
    double temp1, temp2;

    temp1 = x*x;
    temp2 = temp1*temp1;
    return(temp2*temp2);
}
double mypow14(double x) {
    double xp2, xp4;

    xp2 = x*x;
    xp4 = xp2*xp2;
    return(xp4*xp4*xp4*xp2);
}
int getlockind(int var, int lat, int dir, int bead) {

  return var*2*2*beads + lat*2*beads + dir*beads + bead;
}

int chkdir(int lat, int dir, int bead) {
  
  /* returns 1 if there is a basin crossing, 0 otherwise */
  int temp, tbas, odir, n1, n2, olat, newreg, ind;

  temp = 0;
  tbas = bascheck();
		    
  if (tbas + dir == 1) {  /* either tbas = 0, dir = 1 or tbas = 1, dir = 0 */
    ind = dir*beads + bead;
    temp = 1;
			
    /* basin crossing */
			
    allind[0] = lat;
    allind[1] = dir*beads + bead;
    tempi[0] = 1 + NGA_Read_inc(Coords_narr.ga,allind,1);
			
    if (dir == 1) {
      odir = 0;
    }else{
      odir = 1;
    }
    if (verb >= 3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Ended traj: dir cross\n");
      fclose(mainout);
    }
    newreg = whichx(lat,odir,coor);
    if (wtalg == 1) {
      fixwt(lat,dir,odir,bead,newreg); /* transfer weight */
    } else if (rstart[0] != -1) {
      allind[0] = lat;
      allind[1] = dir*beads + bead;
      /*start lock*/
      
      lockind = getlockind(1,lat,dir,bead);
      GA_Lock(lockind);
      GA_Init_fence();
      GetFromServer_dbl(&Coords_tau,lat,ind,tempw);
      tempw[0] += mytim;
      PutToServer_dbl(&Coords_tau,lat,ind,tempw);
      GA_Fence();
      GA_Unlock(lockind);
      
      /*end lock*/
      
      tempi[0] = 1 + NGA_Read_inc(Coords_ntau.ga,allind,1);
      
      n1 = rstart[0];
      n2 = rstart[1];
      
      allind[0] = lat;
      allind[1] = 4*beadsp3*dir + 4*beadsp2*bead + 2*beadsp2*n1 + 2*beads*n2 + beads*odir + newreg;
      tempi[0] = 1 + NGA_Read_inc(Coords_prob.ga,allind,1);			  
      allind[0] = lat;
      allind[1] = 2*beadsp2*dir + 2*beads*bead + beads*n1 + n2;
      tempi[0] = 1 + NGA_Read_inc(Coords_nprob.ga,allind,1);
      
      //nprob[lat][dir][bead][n1][n2]++;
      //prob[lat][dir][bead][n1][n2][odir][newreg]++;
    }
    if (lat == 1) {
      olat =0;
    } else {
      olat =1;
    }
    newreg = whichx(olat,odir,coor);
    ind = dir*beads + bead;
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_weight,lat,ind,tempw);
    mytim = 0;
    addpoint(olat,odir,newreg,coor,tempw[0],dir,state);  /* send point to other flux list */
    endtraj = 1;
  }
  return temp;
}

void chksec(int lat, int dir, int bead) {

  int oldst, olat, ind;

  oldst = state;
  gtst(lat,dir);
  if (state != oldst) {
    if (lat == 1) {
      olat = 0;
    }else{
      olat = 1;
    }
    if (verb >= 3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Adding flux point to %d %d %d\n",olat,dir,state);
      fclose(mainout);
    }
    ind = dir*beads + bead;
    //LINE_STAMP;
    GetFromServer_dbl(&Coords_weight,lat,ind,tempw);
    mytim = 0;
    addpoint(olat,dir,state,coor,tempw[0],dir,oldst);
  }
}

void chkprim(int lat, int dir, int bead) {

  int newreg, ind, n1, n2;
  
  /* prim crossing */
  newreg = which(lat,dir);
  if (newreg != bead) {
    ind = dir*beads + bead;

    if (wtalg == 1) {
      fixwt(lat,dir,dir,bead,newreg); /* fix weight */
    } else if (rstart[0] != -1) {
      allind[0] = lat;
      allind[1] = dir*beads + bead;
      
      lockind = getlockind(1,lat,dir,bead);
      GA_Lock(lockind);
      GA_Init_fence();
      GetFromServer_dbl(&Coords_tau,lat,ind,tempw);
      tempw[0] += mytim;
      PutToServer_dbl(&Coords_tau,lat,ind,tempw);
      GA_Fence();
      GA_Unlock(lockind);

      tempi[0] = 1 + NGA_Read_inc(Coords_ntau.ga,allind,1); 
      
      n1 = rstart[0];
      n2 = rstart[1];
		    
      allind[1] = 4*beadsp3*dir + 4*beadsp2*bead + 2*beadsp2*n1 + 2*beads*n2 + beads*dir + newreg;
      tempi[0] = 1 + NGA_Read_inc(Coords_prob.ga,allind,1); 
      
      allind[1] = 2*beadsp2*dir + 2*beads*bead + beads*n1 + n2;
      tempi[0] = 1 + NGA_Read_inc(Coords_nprob.ga,allind,1); 
		    
      //nprob[lat][dir][bead][n1][n2]++;
      //prob[lat][dir][bead][n1][n2][dir][newreg]++;
    }
    endtraj = 1;
    if (verb >= 3) {
      mainout = fopen(outname,"a");
      ind = dir*beads + bead;
      fprintf(mainout,"Ended traj: prim cross\n");
      fclose(mainout);
    }
  }
}
