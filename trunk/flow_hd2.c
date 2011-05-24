#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> /* required for randomize() and random() */
#include <time.h>
#include <unistd.h>
#include "neusglob.h"
#include "../aryeh.h"
#include "../myCoordServer.h"

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
ServedCoords Coords_count;
ServedCoords Coords_narr;
ServedCoords Coords_on;
ServedCoords Coords_bdedit;
ServedCoords Coords_nflux;
ServedCoords Coords_thist;
ServedCoords Coords_thistsig;
ServedCoords Coords_thistfr;
ServedCoords Coords_full;
ServedCoords Coords_nlist;
ServedCoords Coords_twlist;
ServedCoords Coords_pts;
ServedCoords Coords_wlist;
ServedCoords Coords_from;
ServedCoords Coords_clock;
ServedCoords Coords_mytim;
ServedCoords Coords_nwstack;
ServedCoords Coords_regind;
ServedCoords Coords_nreg;

double myweight=0.;
double **weight, **wavg;
opoint_t cryst;
int **count, **on, **navg, **fmat;
int ***thist, ***thist2, ***thist3;
int allind[2];
int stredit=0;  /* stredit = 0 (hasn't been done yet) 
	         stredit = 1 (it's done already)
		 stredit = -1 (it's in progress) */

/* global variables */

int nfxfile = 60;
int nlockreg = 5;
int nreg[2];
int nsec, ntert;
double *twavg1, *tweight1, *ttau1;
double *tpt, *lfrac;
double globtim=0.;
int nglobtim=0;
int *tint1, *tint2, *ton1, *ton2, *tnflux1, *tnflux2, *tthist1, *tthist2, *tthist3, *tthist5, ****nin, **nout;
int *tnwstack1, *tnwstack2, *tmpfrom, **allcnt, **allon, **secind, **tertind, *tregind, ***myreg;
int verb=2;
int me, nproc, status, beadsp2, beadsp3, lockind, lo[2], hi[2], ld[1];
double tstepsq, tempw[1], wcasig, wcacut;
double sigfacsq=25.;
double rsigfacsq=4.;

char crname[70]="../cutRNAin.dat";      // crystal structure file
char secname[70]="../cutRNAsec.con";    // secondary contacts
char tertname[70]="../cutRNAtert.con";  // tertiary contacts
char emname[30]="ptlist.dat";         // emergency points
char sname[30]="isolv";         // emergency points
char filename[30];         // emergency points
double **r0, **r02, **r0d6, **r0d12, **nextstep, **r0dists, **rdists, **vsolv, **rsolv, **f1;
int **pos, **pos_poly;
double ***v_cmax, ***v_cmay, ***v_cmaz;
double ***v_cmx, ***v_cmy, ***v_cmz; 
xpoint_t **rvecs;
int ***Delta;
int ***box, ***box_poly;
int ****lab, ***nlab;
double polymass=300.;
double solvmass=36.;
double calpha[2], salpha[2], var, var2, solvstep, phan[NDIM_c], n_av;
double phanr0 = 1.;
double phanr06;
int N=503200, maxd, strmfrq;
int Lx=384, Ly=192, Lz=384;
double blen=8.;              /* length of box in each dir */
double kcalconv = 6695.5;
double cutoff=3.;
int nx,ny,nz;
double blinv;
double kf, R0, ooR02, eph1_12, eph2_12, ep1_6, sigma1, sigma16, sigma2, sigma26;
void wrxyz();
int globn;


const int ndim=NDIM_c, mxlist=MXLIST_c, nop=NOP_c, npart=NPART_c;
int onpart=262;
int state, bas, mytim=0, sofar, polyerr=0;
opoint_t ocoor;
point_t coor;
int myclock=0, tmpswitch;
void scalflux();
void wtmat3();
opoint_t projx(point_t x);
double projxsig(point_t x);
void writereg();
void wrtflux(int find, int cycle);
void printiox(point_t x, int dir);
void getop(int i, int j, int *op);
int bascheck();
int which(int dir,int bead);
int whichx(int dir, point_t x);
void gaexit(int sig);
void addpoint (int dir, int bead, point_t x, double wt, int from1, int from2);
void fixwt(int dir1, int dir2, int reg1, int reg2);
double calcsum(int dir);
double ogtdist(opoint_t px, opoint_t py);
void repar(int dir);
void nobad(int dir, int bead);
void getregdef();
void chkprim(int dir, int bead);
int chkdir(int dir, int bead);
int findreg(iopoint_t iox, int dir);
int getlockind(int var, int dir, int bead);
int getlistind(int dir, int bead, int from1, int from2);
int rstart[2], tempi[1], tempi2[1], endtraj;

void **alloc2d(int varsize, int n, int p) ; 
void ***alloc3d(int varsize, int n, int p, int q) ;
void ****alloc4d(int varsize, int n, int p, int q, int r) ;
void *****alloc5d(int varsize, int n, int p, int q, int r, int s) ;
void *******alloc7d(int varsize, int n, int p, int q, int r, int s, int t, int u) ;
double gssran(long *idum) ;
double mypow2(double x);
double mypow4(double x);
double mypow6(double x);
double mypow7(double x);
double mypow8(double x);
double mypow14(double x);
double ran2(long *idum) ;
FILE * mainout, *xyzout;
char fname[30];
char outname[30];
char ratname[30];
char onname[30], twname[30];
char wtname[30];
FILE * outPath, *outFile;
double pmin=1.e-40;

/*-------------------------------------------------------------------------*/

int main(int argc,char** argv){

  int i, j, k, l, m, dim, size, temp, done, worker, mtag, n1, n2, part, lim;
  int op, dir1, dir2, from, to, rem, limit, pt, bead1, bead2, tnlist;
  int dir, bead, tcount, listind;
  int rc, maxsize, ind, flag, mcount, it1, it2;
  int cycle, endcycle, gotmin;
  unsigned int iseed = (unsigned int)time(NULL);
  int back(int dir, int bead);
  void createservers(), initvar();
  void move();
  void string(), stack(int dir, int bead), wrpath(int ind);
  int getmin(int *dir, int *bead);
  void strread(), wtread(), fxread(int find), nfluxread();
  void wtstack(int dir, int bead);
  void updtwts(int ind);
  void wrextend(int index);
  void rddelta(), rdcryst(), rdsolv(), rdemerg();
  double rate[2], sum1, sum2, t1, t2, telaps, alpha;
  int wtupdt_s;
  FILE * ratFile, *filein, *onfile, *twfile;

  nreg[0] = 0;
  nreg[1] = 0;

  ld[0] = 0;
  tmpswitch = 0;

  if (seed ==0) {
      seed = iseed;
  }

  tstepsq=tstep*tstep;
  beadsp2 = beads*beads;
  beadsp3 = beadsp2*beads;
  wtupdt_s = wtupdt*every/100;

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
  
  GA_Initialize_ltd(800*1024*1024*sizeof(double));

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
  strmfrq = (int)(1.5/tstep);       /* stream every 6 ps (1.5 time units) */
  solvstep = tstep*strmfrq;
  phan[0] = 120.;
  phan[1] = 25.;
  phan[2] = 192;
  n_av=N/(float)(nx*(ny-1)*nz);	

  kf = 20*kcalconv;
  R0 = 2;
  ooR02 = 1./(R0*R0);
  eph1_12 = 0.7*12*bfrac*kcalconv;
  eph2_12 = 0.7*12*bfrac*kcalconv;
  ep1_6 = 1*6*bfrac*kcalconv;
  sigma1 = 7;
  sigma2 = 3.5;
  sigma16 = mypow8(sigma1)/mypow2(sigma1);
  sigma26 = mypow8(sigma2)/mypow2(sigma2);
  phanr06 = mypow8(phanr0)/mypow2(phanr0);
  alpha = 0.243*Pi;
  calpha[0]=cos(alpha);
  salpha[0]=sin(alpha);
  calpha[1]=cos(-alpha);
  salpha[1]=sin(-alpha);
  wcasig = 2.0;
  wcacut = wcasig*pow(2.,1./6.);

  if ((verb >= 3) || ((verb >=2) && (me == 0))) {
      sprintf(outname,"%s%d%s","out",me,".dat");
      mainout = fopen(outname,"w");
      fprintf(mainout,"rank = %d\n", me);
      fclose(mainout);
  }
  //  sprintf(outname,"%s","out0.dat");
  if (me == 0) {
      sprintf(wtname,"%s","weight.dat");
      ratFile = fopen(wtname,"w");
      fclose(ratFile);
      
      sprintf(ratname,"%s","rate.dat");
      ratFile = fopen(ratname,"w");
      fclose(ratFile);

      sprintf(onname,"%s","on.dat");
      onfile = fopen(onname,"w");
      fclose(onfile);

      sprintf(twname,"%s","tw.dat");
      twfile = fopen(twname,"w");
      fclose(twfile);
  }

  seed *= me+1;
  for (i=0;i<100;i++){
    sum1 = ran2(&seed);
  }
  
  lim = 2*beads;
  twavg1 = (double *) malloc(2*beads*sizeof(double));
  tweight1 = (double *) malloc(2*beads*sizeof(double));
  tpt = (double *) malloc(3*ndim*npart*sizeof(double));
  lfrac = (double *) malloc(2*2*beads*2*beads*sizeof(double));
  tmpfrom = (int *) malloc(2*beads*mxlist*2*sizeof(int));
  tregind = (int *) malloc(beads*nop*sizeof(int));

  tint1 = (int *) malloc(2*beads*sizeof(int));
  ton1 = (int *) malloc(2*beads*sizeof(int));
  tnwstack1 = (int *) malloc(2*beads*sizeof(int));
  tthist1 = (int *) malloc(2*beads*tres*sizeof(int));
  tthist3 = (int *) malloc(2*beads*tres*sizeof(int));
  tthist5 = (int *) malloc(2*beads*tres*sizeof(int));
  tnflux1 = (int *) malloc(2*beads*2*beads*sizeof(int));
  tnflux2 = (int *) malloc(2*beads*2*beads*sizeof(int));

  nout = (int **) alloc2d(sizeof(int), 2, beads);
  nin = (int ****) alloc4d(sizeof(int), 2, beads, 2, beads);
  myreg = (int ***) alloc3d(sizeof(int), 2, beads, nop);
  on = (int **) alloc2d(sizeof(int), 2, beads);
  weight = (double **) alloc2d(sizeof(double), 2, beads);
  wavg = (double **) alloc2d(sizeof(double), 2, beads);
  navg = (int **) alloc2d(sizeof(int), 2, beads);
  lim = 2*beads;
  allcnt = (int **) alloc2d(sizeof(int), 2, lim);
  allon = (int **) alloc2d(sizeof(int), 2, lim);
  fmat = (int **) alloc2d(sizeof(int), lim, lim);
  count= (int **) alloc2d(sizeof(int), 2, beads);
  thist = (int ***) alloc3d(sizeof(int), 2, beads, tres);
  thist2 = (int ***) alloc3d(sizeof(int), 2, beads, tres);
  thist3 = (int ***) alloc3d(sizeof(int), 2, beads, tres);

  box = (int ***) alloc3d(sizeof(int), nx, ny, nz);
  box_poly = (int ***) alloc3d(sizeof(int), nx, ny, nz);
  
  maxd = 5;
  lab = (int ****) alloc4d(sizeof(int),nx,ny,nz,maxd);
  nlab = (int ***) alloc3d(sizeof(int),nx,ny,nz);

  r0 = (double **) alloc2d(sizeof(double), npart, 3);
  r02 = (double **) alloc2d(sizeof(double), npart, npart);
  rvecs = (xpoint_t **) alloc2d(sizeof(xpoint_t), npart, npart);
  r0dists = (double **) alloc2d(sizeof(double), npart, npart);
  r0d6 = (double **) alloc2d(sizeof(double), npart, npart);
  r0d12 = (double **) alloc2d(sizeof(double), npart, npart);
  rdists = (double **) alloc2d(sizeof(double), npart, npart);
  nextstep = (double **) alloc2d(sizeof(double), npart, 3);
  Delta = (int ***) alloc3d(sizeof(int), onpart, onpart,4);

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

  initvar();

  if (splitread) {
    if (me < nfxfile) {
      GA_Init_fence();
      fxread(me);
      GA_Fence();
    }
  } else {
    if (me == 0) {
      fxread(-1);
    }
  }

  GA_Sync();

  if (me == 0) {
    GA_Init_fence();
    
    if (strcmp(fxname,"NOREAD") != 0) {  /* read in fluxes */
      writereg();
    } else {
      rdemerg();
    }

    /* read in nflux data */

    if (strcmp(nfname,"NOREAD") != 0) {
      if (verb >= 2) {
	printf("Reading nflux data from %s\n", nfname); fflush(stdout);
      }
      nfluxread();
      if (wtalg ==3) {
	mcount=0;
	for (j=0;j<2;j++) {
	  for (k=0;k<beads;k++) {
	    tweight1[mcount] = 1./(float)(2*beads);
	    mcount++;
	  }
	}
	GA_Init_fence();
	PutToServer_dbl(&Coords_weight,0,-1,tweight1);
	GA_Fence();
	updtwts(0);
	scalflux();
      }
    }

    /* initialize weights */
    
    if (strcmp(nfname,"NOREAD") == 0) {
      if (strcmp(wname,"NOREAD") == 0) {
	if (verb >= 2) {
	  printf("Automatically initializing weights\n"); fflush(stdout);
	}
	mcount = 0;
	for (j=0;j<2;j++) {
	  for (k=0;k<beads;k++) {
	    if (wtalg == 0) {
	      tweight1[mcount] = 0.;
	    } else {
	      tweight1[mcount] = 1./(float)(2*beads);
	    }
	    mcount++;
	  }
	}
	if (wtalg ==0) {
	  tweight1[0] = 0.5;
	  tweight1[2*beads-1] = 0.5;
	  wtalg = 1;
	}
	PutToServer_dbl(&Coords_weight,0,-1,tweight1);
	if (wtalg == 4) {
	  mcount = 0;
	  for (j=0;j<2;j++) {
	    for (k=0;k<beads;k++) {
	      tweight1[mcount] = 0.;
	      mcount++;
	    }
	  }
	}
	PutToServer_dbl(&Coords_wadj,0,-1,tweight1);
      } else {
	if (verb >= 2) {
	  printf("Using weight configuration from %s\n", wname); fflush(stdout);
	}
	wtread();
      }
    }

    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	ind = dir*beads + bead;
	GetFromServer_dbl(&Coords_twlist,0,ind,tempw);
	t1 = tempw[0];
	GetFromServer_dbl(&Coords_twlist,1,ind,tempw);
	if (t1 + tempw[0] > 0.) {
	  tint1[mcount] = 1;
	} else {
	  tint1[mcount] = 0;
	}
	mcount++;
      }
    }
    PutToServer_int(&Coords_on,0,-1,tint1);    // set regions to 'on' or 'off'  

    printf("proc %d is done with the string\n",me); fflush(stdout);
    GA_Fence();
  }

  /* check point */
  if (verb >=3) {
    t2 = MPI_Wtime();
    mainout = fopen(outname,"a");
    fprintf(mainout,"%f: Bef sync2\n",t2); 
    fclose(mainout);
  }
  GA_Sync();
  if (verb >=3) {
    t2 = MPI_Wtime();
    mainout = fopen(outname,"a");
    fprintf(mainout,"%f: Sync2\n",t2);
    fclose(mainout);
  }

  sprintf(filename,"flow_ufold.xyz");
    
  if (me==0) {
    xyzout = fopen(filename,"w");
    fprintf(xyzout,"%i \n polymer movie\n", npart);
    fclose(xyzout);
  }

  if (me == 0) {
    t1 = MPI_Wtime();
  }
  getregdef();

  /* ---------------start of dynamics loop-------------------------- */
  
  for (cycle=1;cycle<=T;cycle++) { 
    
    if (me == 0) {
      printf("Starting cycle %d of %d\n",cycle,T);
    
      mcount = 0;
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<=beads-1;bead++) {
	  tint1[mcount] = 0;
	  mcount++;
	}
      }
      GA_Init_fence();
      PutToServer_int(&Coords_count,0,-1,tint1);  // set all counts to 0
      GA_Fence();
     
    }
    if (verb >=3) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: Bef beginning of cycle\n",t2);
      fclose(mainout);
    }
    GA_Sync();
    if (verb >=3) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: Beginning of cycle\n",t2);
      fclose(mainout);
    }

    endcycle = 0;
    while (!endcycle) {
      
      /* start a new trajectory */
      tcount = getmin(&dir,&bead);   /* returns (dir,bead) for the region with the lowest counter */
      
      if (dir == -1) {
	endcycle = 1;
      } else {
	
	if ((verb >= 3) || ((verb >= 2) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"Running on %d %d, start count = %d\n",dir,bead,tcount);
	  fclose(mainout);
	}

	endtraj = 0;
	if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  t2 = MPI_Wtime();
	  fprintf(mainout,"before back (%f)\n",t2);
	  fclose(mainout);
	}
	if (!back(dir,bead)) {  /* initialize trajectory */
	  endtraj = 1;
	} else {
	  ind = dir*beads + bead;
	  GetFromServer_int(&Coords_count,0,ind,tempi);
	  if (tempi[0] > every) {
	    endtraj = 1;
	  }
	}
	if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  t2 = MPI_Wtime();
	  fprintf(mainout,"after back (%f)\n",t2);
	  fclose(mainout);
	}
	while (!endtraj) { 
	  if (verb >= 4) {
	    mainout = fopen(outname,"a");
	    fprintf(mainout,"updating count and elaps...\n");
	    fclose(mainout);
	  }
	  if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
	    mainout = fopen(outname,"a");
	    t2 = MPI_Wtime();
	    fprintf(mainout,"starting dynamics (%f)\n",t2);
	    fclose(mainout);
	  }

	  allind[0] = 0;
	  allind[1] = dir*beads + bead;
	  tcount = chkfrq + NGA_Read_inc(Coords_count.ga,allind,chkfrq);
	  tempi[0] = chkfrq + NGA_Read_inc(Coords_elaps.ga,allind,chkfrq);

	  /* move */
	  for (i=0;i<chkfrq;i++) {  
	    if (!polyerr) {
	      if (verb >= 4) {
		mainout = fopen(outname,"a");
		fprintf(mainout,"move %d of %d...\n",i+1,chkfrq);
		fclose(mainout);
	      }
	      move();
	      mytim++;
	      myclock++;
	    }
	  }

	  if (!polyerr) {

	    /* check */	
	    if (verb >= 4) {
	      mainout = fopen(outname,"a");
	      fprintf(mainout,"chkdir...\n");
	      fclose(mainout);
	    }	    
	    if (!chkdir(dir,bead)) {
	      if (verb >= 4) {
		mainout = fopen(outname,"a");
		fprintf(mainout,"chkprim...\n");
		fclose(mainout);
	      }
	      chkprim(dir,bead);
	    }
	
	    /* stack */
		
	    if (tcount%stkfrq == 0) {
	      if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
		mainout = fopen(outname,"a");
		t2 = MPI_Wtime();
		fprintf(mainout,"in stack (%f)\n",t2);
		fclose(mainout);
	      }
	      if (verb >= 4) {
		mainout = fopen(outname,"a");
		fprintf(mainout,"stacking...\n");
		fclose(mainout);
	      }
	      stack(dir,bead);
	      if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
		mainout = fopen(outname,"a");
		t2 = MPI_Wtime();
		fprintf(mainout,"out stack (%f)\n",t2);
		fclose(mainout);
	      }
	    }
		
	    if (wtalg == 1) { /* local wt alg */
	      if ((tcount + cycle*every)%wtupdt_s == 0) {
		if (verb >= 4) {
		  mainout = fopen(outname,"a");
		  fprintf(mainout,"weight stacking...\n");
		  fclose(mainout);
		}
		wtstack(dir,bead);
		if (verb >= 4) {
		  mainout = fopen(outname,"a");
		  fprintf(mainout,"done weight stacking...\n");
		  fclose(mainout);
		}
	      }
	    }
		
	    if ((tcount >= every) && (!endtraj)) {
	      if (verb >= 3) {
		mainout = fopen(outname,"a");
		fprintf(mainout,"Ended traj: count\n");
		fclose(mainout);
	      }
	      /* write coordinates to your fluxlist */

	      addpoint(dir,bead,coor,myweight,rstart[0],rstart[1]);
	      endtraj = 1;
	    }
	  } else {   /* if polyerr */
	    if (verb >= 3) {
	      mainout = fopen(outname,"a");
	      fprintf(mainout,"polyerr!\n");
	      fclose(mainout);
	    }
	    endtraj = 1;
	    polyerr = 0;
	    printf("polyerr: (%d,%d) from (%d,%d)\n",dir,bead,rstart[0],rstart[1]);
	    //	    gaexit(2);
	  }
	}  /* end of trajectory loop */
	if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  t2 = MPI_Wtime();
	  fprintf(mainout,"stopping dynamics (%f)\n",t2);
	  fclose(mainout);
	}
      }
    }       /* end of cycle */
    
      
    /* compute the rate / update weights / write thist */
    if (verb >=3) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: Bef TCOB\n",t2);
      fclose(mainout);
    }
    GA_Sync();
    if (verb >=3) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: TCOB\n",t2);
      fclose(mainout);
    }
    if (me == 0) {
      printf("proc %d: TCOB\n",me); fflush(stdout);

      /* update weights */

      if (cycle%wtupdt == 0) {
	updtwts(cycle/wtupdt);
      }
    }
    GA_Sync();
    
    if (me == 0) {
    
      /* write 'on' data */
      GetFromServer_int(&Coords_on,0,-1,ton1);
      
      onfile = fopen(onname,"a");
      mcount = 0;
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<beads;bead++) {
	  fprintf(onfile,"%d",ton1[mcount]);
	  mcount++;
	}
	fprintf(onfile,"\n");
      }
      fprintf(onfile,"\n");
      fclose(onfile);

      /*write twlist data */
      GetFromServer_dbl(&Coords_twlist,0,-1,tweight1);
      GetFromServer_dbl(&Coords_twlist,1,-1,twavg1);
      
      twfile = fopen(twname,"a");
      mcount = 0;
      for (dir=0;dir<=1;dir++) {
	for (bead=0;bead<beads;bead++) {
	  tweight1[mcount] += twavg1[mcount];
	  if ((tweight1[mcount] != 0.) && (ton1[mcount])) {
	    fprintf(twfile,"%5d ",(int)(log10(tweight1[mcount])));
	  } else {
	    fprintf(twfile,"      ");
	  }
	  mcount++;
	}
	fprintf(twfile,"\n");
      }
      fprintf(twfile,"\n");
      fclose(twfile);

      /* compute rate */

      if (cycle%wrfrq == 0) {
	GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
      
	mcount = 0;
	for (dir=0;dir<=1;dir++) {
	  for (bead=0;bead<=beads-1;bead++) {
	    weight[dir][bead] = tweight1[mcount];
	    mcount++;
	  }
	}
      
	for (dir=0;dir<=1;dir++) {
	  rate[dir] = 0.;
	  for (bead=0;bead<=beads-1;bead++) {
	    ind = dir*beads + bead;
	    GetFromServer_int(&Coords_elaps,0,ind,tempi);
	    telaps = (float)tempi[0];
	    GetFromServer_int(&Coords_narr,0,ind,tempi);
	    
	    if (telaps != 0.) {
	      rate[dir] += tempi[0]*weight[dir][bead]/telaps;
	    }
	  }
	  rate[dir] /= calcsum(dir);
	}
	ratFile = fopen(ratname,"a");
	fprintf(ratFile,"%e %e\n",rate[0],rate[1]);
	fclose(ratFile);

	/* write out thist */

	wrextend(cycle/wrfrq);
            
	/* write nflux */

	if (wtalg != 1) {
	  /*	  GetFromServer_int(&Coords_nflux,0,-1,tnflux1);
	  GetFromServer_int(&Coords_nflux,1,-1,tnflux2);
	
	  sprintf(fname,"%s","tnflux1.dat");
	  outFile = fopen(fname,"w");
	  lim = 2*beads;
	  k = 0;
	  for (i=0;i<lim;i++) {
	    for (j=0;j<lim;j++) {
	      fprintf(outFile,"%d ",tnflux1[k]);
	      k++;
	    }
	    fprintf(outFile,"\n");
	  }
	  fclose(outFile);

	  sprintf(fname,"%s","tnflux2.dat");
	  outFile = fopen(fname,"w");
	  lim = 2*beads;
	  k = 0;
	  for (i=0;i<lim;i++) {
	    for (j=0;j<lim;j++) {
	      fprintf(outFile,"%d ",tnflux2[k]);
	      k++;
	    }
	    fprintf(outFile,"\n");
	  }
	  fclose(outFile);*/
	
	  printf("proc %d: done TCOB\n",me); fflush(stdout);
	}
      }
    } else if (me <= nfxfile) {

      /* write fluxes, scale */
      if (cycle%wrfrq == 0) {
	if (tmpswitch==0) {
	  tmpswitch=1;
	} else {
	  tmpswitch=0;
	}
	wrtflux(me,cycle);
      }
    } else if (me == nfxfile+1) {
      writereg();
    }

    /* checkpoint */
    if (verb >=3) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: Bef end of cycle\n",t2);
      fclose(mainout);
    }
    GA_Sync();
    if (verb >=3) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: End of cycle\n",t2);
      fclose(mainout);
    }
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
    
    mcount = 0;
    tcount = 0;
    for (i=0;i<2;i++) {
      for (j=0;j<beads;j++) {
	tcount += tint1[mcount];
	mcount++;
      }
    }
    printf("total number of dynamics steps: %d\n",tcount); fflush(stdout);
    printf("time per dynamics step: %f\n",(t2-t1)/(float)tcount); fflush(stdout);
    

  /* write weights */

    sprintf(fname,"%s","fwts.dat");
    outFile = fopen(fname,"w");
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-1; k++) {
	fprintf(outFile,"%e\n",weight[j][k]);
      }
    }
    fclose(outFile);

    /* write fluxes */
  
    sprintf(fname,"%s","fflux.dat");
    outFile = fopen(fname,"w");
    for (j=0; j<=1; j++) {
      for (k=0; k<=beads-1; k++) {
	fprintf(outFile,"%d %d\n",j,k);
	for (listind=0;listind<2;listind++) {
	  ind = j*beads + k;
	  //LINE_STAMP;
	  GetFromServer_int(&Coords_full,listind,ind,tempi);
	  if (tempi[0] == 1) {
	    limit = mxlist;
	  } else {
	    GetFromServer_int(&Coords_nlist,listind,ind,tempi);
	    limit = tempi[0];
	  }
	  
	  GetFromServer_dbl(&Coords_twlist,listind,ind,tempw);
	  fprintf(outFile,"%d %e\n", limit, tempw[0]);
	  
	  for (pt=0; pt<limit; pt++) {
	    ind = 2*mxlist*beads*j + 2*mxlist*k + 2*pt;
	    //LINE_STAMP;
	    GetFromServer_int(&Coords_from,listind,ind,tempi);
	    //LINE_STAMP;
	    GetFromServer_int(&Coords_from,listind,ind+1,tempi2);
	    fprintf(outFile,"%d %d ", tempi[0], tempi2[0]);
	    
	    ind = mxlist*beads*j + mxlist*k + pt;
	    //LINE_STAMP;
	    GetFromServer_dbl(&Coords_wlist,listind,ind,tempw);
	    GetFromServer_int(&Coords_clock,listind,ind,tempi);
	    fprintf(outFile,"%e %d ",tempw[0],tempi[0]);
	    GetFromServer_int(&Coords_mytim,listind,ind,tempi);
	    fprintf(outFile,"%d ",tempi[0]);
	    
	    lo[0] = listind;
	    hi[0] = listind;
	    lo[1] = 3*npart*ndim*mxlist*beads*j + 3*npart*ndim*mxlist*k + 3*npart*ndim*pt;
	    hi[1] = 3*npart*ndim*mxlist*beads*j + 3*npart*ndim*mxlist*k + 3*npart*ndim*(pt+1) -1;
	    NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);
	    
	    mcount = 0;
	    for (part=0; part<npart;part++) {
	      for (dim=0; dim<=ndim-1; dim++) {
		fprintf(outFile,"%e %e %e ", tpt[mcount], tpt[mcount+1], tpt[mcount+2]);
		mcount += 3;
	      }
	    }
	    fprintf(outFile,"\n");
	  }
	}
      }
    }
    fclose(outFile);

    if (wtalg != 1) {
    /* write nflux */
  
      GetFromServer_int(&Coords_nflux,0,-1,tnflux1);
      GetFromServer_int(&Coords_nflux,1,-1,tnflux2);
      
      sprintf(fname,"%s","fnflux1.dat");
      outFile = fopen(fname,"w");
      lim = 2*beads;
      k = 0;
      for (i=0;i<lim;i++) {
	for (j=0;j<lim;j++) {
	  fprintf(outFile,"%d ",tnflux1[k]);
	  k++;
	}
	fprintf(outFile,"\n");
      }
      fclose(outFile);
      
      sprintf(fname,"%s","fnflux2.dat");
      outFile = fopen(fname,"w");
      lim = 2*beads;
      k = 0;
      for (i=0;i<lim;i++) {
	for (j=0;j<lim;j++) {
	  fprintf(outFile,"%d ",tnflux2[k]);
	  k++;
	}
	fprintf(outFile,"\n");
      }
      fclose(outFile);
      
    }
		    
    /* write results */
  }
  
  /* check point */
//  printf("proc %d: final check\n",me); fflush(stdout);
  GA_Sync();

  status = DestroyCoordServer(&Coords_wavg);
  status = DestroyCoordServer(&Coords_wadj);
  status = DestroyCoordServer(&Coords_weight);
  status = DestroyCoordServer(&Coords_elaps);

  status = DestroyCoordServer(&Coords_count);
  status = DestroyCoordServer(&Coords_narr);
  status = DestroyCoordServer(&Coords_on);
  status = DestroyCoordServer(&Coords_bdedit);
  status = DestroyCoordServer(&Coords_nflux);


  status = DestroyCoordServer(&Coords_thist);
  status = DestroyCoordServer(&Coords_thistsig);
  status = DestroyCoordServer(&Coords_thistfr);
  status = DestroyCoordServer(&Coords_full);
  status = DestroyCoordServer(&Coords_nlist);
  status = DestroyCoordServer(&Coords_twlist);
  status = DestroyCoordServer(&Coords_pts);
  status = DestroyCoordServer(&Coords_wlist);
  status = DestroyCoordServer(&Coords_from);
  status = DestroyCoordServer(&Coords_clock);
  status = DestroyCoordServer(&Coords_mytim);
  status = DestroyCoordServer(&Coords_nwstack);
  status = DestroyCoordServer(&Coords_regind);
  status = DestroyCoordServer(&Coords_nreg);

  if (me == 0) GA_Print_stats();
  
  GA_Terminate();
  MPI_Finalize();
  
  return(0);
}




/*=========================================================
  END OF MAIN
  ============================================================*/

/*-------------------------------------------------------------------*/
int back (int dir, int bead) {

  /* This subroutine resets the position of a walker
     to a saved point, or if none exist, to the position
     of the bead (guestimating for undefined coordinates) */
  
  double sum, a, ttwlist, oran, tlist1, tlist2, t2, tmpwlist[mxlist];
  int i, j, newreg, ind, limit, k, mcount, tfrom2[2], gotit, listind;
  int tnlist, tfull, part, dim, good, wind;
  point_t point;
  FILE *filein;

  good = 0;
  gotit = 1;
  while (!good) {
    good = 1;

    ind = dir*beads + bead;    
    GetFromServer_dbl(&Coords_twlist,0,ind,tempw);
    tlist1 = tempw[0];
    GetFromServer_dbl(&Coords_twlist,1,ind,tempw);
    tlist2 = tempw[0];
    oran = (tlist1+tlist2)*ran2(&seed);
    if (oran < tlist1) {
      listind = 0;
    } else {
      listind = 1;
    }

    lockind = getlockind(0,dir,bead);
    if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout,"(%f) getting nf/f",t2);
      fclose(mainout);
    }
    GA_Lock(lockind);
    if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout,"(%f)",t2);
      fclose(mainout);
    }
    GA_Init_fence();
    if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout,"(%f)",t2);
      fclose(mainout);
    }
    GetFromServer_int(&Coords_nlist,listind,ind,tempi);
    tnlist = tempi[0];
    GetFromServer_int(&Coords_full,listind,ind,tempi);
    tfull = tempi[0];
    if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout,"(%f)",t2);
      fclose(mainout);
    }
    GA_Fence();
    if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout,"(%f)",t2);
      fclose(mainout);
    }
    GA_Unlock(lockind);
    if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout," (%f)\n",t2);
      fclose(mainout);
    }
    
    ttwlist = 0.;
    if (tnlist != 0) {
      sum = 0.;
      oran = ran2(&seed);
      a = oran;

      // get complete wlist
      lo[0] = listind;
      hi[0] = listind;
      lo[1] = mxlist*beads*dir + mxlist*bead;
      hi[1] = mxlist*beads*dir + mxlist*(bead + 1) - 1;

      if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
	mainout = fopen(outname,"a");
	t2 = MPI_Wtime();
	fprintf(mainout,"(%f) weight ops",t2);
	fclose(mainout);
      }

      lockind = getlockind(1+listind,dir,bead);
      GA_Lock(lockind);
      if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
	mainout = fopen(outname,"a");
	t2 = MPI_Wtime();
	fprintf(mainout," (%f)",t2);
	fclose(mainout);
      }
      GA_Init_fence();
      if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
	mainout = fopen(outname,"a");
	t2 = MPI_Wtime();
	fprintf(mainout," (%f)",t2);
	fclose(mainout);
      }
      NGA_Get(Coords_wlist.ga,lo,hi,tmpwlist,ld);

      // get twlist

      GetFromServer_dbl(&Coords_twlist,listind,ind,tempw);
      ttwlist = tempw[0];

      if (tfull) {
	limit = mxlist;
      } else {
	limit = tnlist;
      }

      if (ttwlist < 0) {
	printf("proc %d:  twlist < 0!! %e, dir,bead = %d %d\n", me, ttwlist, dir, bead);
	printf("proc %d:  nlist = %d; full = %d\n", me, tnlist, tfull);
	sum = 0.;
	for (i=0;i<=limit-1;i++) {
	  printf("i = %d; wt = %e\n",i,tmpwlist[i]);
	  sum += tmpwlist[i];
	}
	tempw[0] = sum;
	ind = dir*beads + bead;
	PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	ttwlist = tempw[0];
      }

      //      printf("In back: dir,bead = %d, %d; tnlist=%d; ttwlist=%f\n",dir,bead,tnlist,ttwlist);
      if (ttwlist !=0.) {

	a *= ttwlist;
	j = 1;
	while (sum < a) {
	  if ((!tfull) && (j > tnlist)) {  /* nlist is the actual number of points in the list */
	    sum = 0.;                                                  /* so subtract 1 to use as an array index */
	    printf("proc %d:  Error!  Weights not deleted!\n",me);
	    printf("Random Number: %e\n", a/ttwlist);
	    printf("dir = %d\n", dir);
	    printf("nlist(%d)=%d \n", bead,tnlist);
	    for (i=0;i<tnlist;i++) {
	      sum = sum + tmpwlist[i];
	    }
	    printf("j = %d\n", j);
	    printf("Manual sum: %e\n", sum);
	    printf("twlist = %e\n", ttwlist);
	    printf("Resetting twlist..\n");
	    j = tnlist+1;
	    printf("Use j = %d\n", j);
	    tempw[0] = sum;
	    ind = dir*beads + bead;
	    
	    PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	    break;
	  } else if ((tfull) && (j > mxlist)) {
	    sum = 0.;
	    printf("proc %d:  Error!  Weights not deleted!\n", me);
	    printf("Random Number: %e\n", a/ttwlist);
	    printf("dir = %d\n", dir);
	    printf("List full\n");
	    for (i=0;i<mxlist;i++) {
	      sum = sum + tmpwlist[i];
	    }
	    printf("Manual sum: %e\n", sum);
	    printf("twlist = %e\n", ttwlist);
	    printf("Resetting twlist..\n");
	    j = mxlist+1;
	    printf("Use j = %d\n", j);
	    tempw[0] = sum;
	    ind = dir*beads + bead;
	    
	    PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	    break;
	  }
	  sum += tmpwlist[j-1];
	  j++;
	}

	wind = mxlist*beads*dir + mxlist*bead + j-2;      
	myweight = spltfac*tmpwlist[j-2];
	tempw[0] = (1.-spltfac)*tmpwlist[j-2];
	if (tempw[0] < pmin) {
	  tempw[0] = 0.;
	  myweight = tmpwlist[j-2];
	}
	tmpwlist[j-2] = tempw[0];
	if (myweight < 0) {
	  printf("proc %d: Error! myweight = %f (dir,bead,listind = %d,%d,%d)\n",me,myweight,dir,bead,listind);
	  gaexit(20);
	}
	
	PutToServer_dbl(&Coords_wlist,listind,wind,tempw);
	GetFromServer_dbl(&Coords_twlist,listind,ind,tempw);
	tempw[0] -= myweight;
	PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	if (verb >=4) {
	  t2 = MPI_Wtime();
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"%f: Subtracting weight from (%d,%d,%d):  %e - %e = %e\n",t2,dir,bead,listind,tempw[0]+myweight,myweight,tempw[0]); 
	  fclose(mainout);
	}
	sum = 0.;
	for (i=0;i<limit;i++) {
	  sum += tmpwlist[i];
	}
	if (fabs(sum - tempw[0]) > 0.001*sum) {
	  tempw[0] = sum;
	  PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	}
	if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  t2 = MPI_Wtime();
	  fprintf(mainout," (%f)",t2);
	  fclose(mainout);
	}
	GA_Fence();
	if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  t2 = MPI_Wtime();
	  fprintf(mainout," (%f)",t2);
	  fclose(mainout);
	}
	GA_Unlock(lockind);
	if ((verb >= 4) || ((verb >= 2) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  t2 = MPI_Wtime();
	  fprintf(mainout," (%f)\n",t2);
	  fclose(mainout);
	}
	
	lo[0] = listind;
	hi[0] = listind;
	lo[1] = 3*npart*ndim*mxlist*beads*dir + 3*npart*ndim*mxlist*bead + 3*npart*ndim*(j-2);
	hi[1] = 3*npart*ndim*mxlist*beads*dir + 3*npart*ndim*mxlist*bead + 3*npart*ndim*(j-1) -1;
	NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);
	
	mcount = 0;
	for (part=0; part<npart;part++) {
	  for (dim=0; dim<=ndim-1; dim++) {
	    coor.x[part][dim] = tpt[mcount];
	    mcount++;
	    coor.v[part][dim] = tpt[mcount];
	    mcount++;
	    f1[part][dim] = tpt[mcount];
	    mcount++;
	  }
	}
	
	lo[1] = 2*mxlist*beads*dir + 2*mxlist*bead + 2*(j-2);
	hi[1] = lo[1] + 1;
	NGA_Get(Coords_from.ga,lo,hi,tfrom2,ld);
	rstart[0] = tfrom2[0];
	rstart[1] = tfrom2[1];
	
	ind = mxlist*beads*dir + mxlist*bead + (j-2);
	//LINE_STAMP;
	GetFromServer_int(&Coords_clock,listind,ind,tempi);
	myclock = tempi[0];
	GetFromServer_int(&Coords_mytim,listind,ind,tempi);
	mytim = tempi[0];
	
	/* get ocoor, state, bas, check the validity of the point */
	
	ocoor = projx(coor);
	
	bas = bascheck();
	if (dir + bas == 1) {
	  printf("proc %d:  Error! Started in opposite basin!\n",me);
	  printf("Bad flux point j = %d\n", j);
	  good = 0; 
	}
	newreg = which(dir,bead);
	if (newreg != bead) {
	  if (verb >= 2) {
	    printf("proc %d:  Error:  Replica reinitialized in wrong region!\n",me);
	    printf("which(%d) = %d\n", dir, newreg);
	    printiox(coor, dir);
	    printf("rank is %d; j is %d; dir,bead are (%d,%d) \n", me, j, dir,bead);
	    printf("tnlist: %d, oran: %e\n",tnlist,oran);
	    printf("from: %d %d\n",rstart[0],rstart[1]);
	    good = 0;
	  }
	  //  gaexit(96);
	}
	if (!good) {
	  if (j >= 2) {
	    printf("deleting point %d\n",j-2);
	    wind = mxlist*beads*dir + mxlist*bead + j-2;
	    tempw[0] = 0.;
	    PutToServer_dbl(&Coords_wlist,listind,wind,tempw);
	    
	    ind = dir*beads + bead;
	    lockind = getlockind(1+listind,dir,bead);
	    GA_Lock(lockind);
	    GA_Init_fence();
	    GetFromServer_dbl(&Coords_twlist,listind,ind,tempw);
	    tempw[0] -= tmpwlist[j-2];
	    PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	    GA_Fence();
	    GA_Unlock(lockind);
	  }
	}
      } else {
	tempi[0] = 0;
	ind = dir*beads + bead;
	printf("proc %d: twlist = 0, turning (%d,%d) off!\n",me,dir,bead);
	if ((verb >= 3) || ((verb >= 2) && (me ==0))) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"twlist = 0, turning (%d,%d) off!\n",dir,bead);
	  fclose(mainout);
	}
	PutToServer_int(&Coords_on,0,ind,tempi);
	GA_Fence();
	GA_Unlock(lockind);
	gotit = 0;
      }
    } else {
      if ((verb >= 3) || ((verb >= 2) && (me ==0))) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"in region (%d,%d), no points in list %d (%f,%f,%f)!\n",dir,bead,listind,tlist1,tlist2,oran);
	fclose(mainout);
      }
      gotit = 0;
      tempi[0] = 0;
      ind = dir*beads + bead;
      printf("proc %d: twlist = 0, turning (%d,%d) off!\n",me,dir,bead);
      if ((verb >= 3) || ((verb >= 2) && (me ==0))) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"twlist = 0, turning (%d,%d) off!\n",dir,bead);
	fclose(mainout);
      }
      PutToServer_int(&Coords_on,0,ind,tempi);
    }
  }
  return gotit;
}
/*-------------------------------------------------------------------------*/
 void fixwt(int dir1, int dir2, int reg1, int reg2) {
   
   double temp, w[1], tempe[1], telaps, telaps0, iw;
   int i, i2;

   i = dir1*beads + reg1;
   i2 = dir2*beads + reg2;
   //LINE_STAMP;
   GetFromServer_dbl(&Coords_wadj,0,i,w);

   //LINE_STAMP;
   //GetFromServer_int(&Coords_elaps,0,i,tempi);
   //telaps = (float)tempi[0];

   //LINE_STAMP;
   //GetFromServer_int(&Coords_elaps,0,0,tempi);
   //telaps0 = (float)tempi[0];

   //if (telaps != 0.) {
     //temp = s*w[0]*telaps0/telaps;
   temp = s*w[0];
     /*} else {
     temp = 0.;
     }*/

   w[0] -= temp;
   //LINE_STAMP;
   PutToServer_dbl(&Coords_wadj,0,i,w);
   iw = w[0];

   lockind = getlockind(2,dir2,reg2);
   GA_Lock(lockind);
   GA_Init_fence();
   GetFromServer_dbl(&Coords_wadj,0,i2,w);
   w[0] += temp;
   PutToServer_dbl(&Coords_wadj,0,i2,w);
   GA_Fence();
   GA_Unlock(lockind);

 }
/*-------------------------------------------------------------------------*/
void wtstack(int dir, int bead) {
   
  int i;
  double tempwadj[1], tempwavg[1];
  
  i = dir*beads + bead;

  if (wtalg != 1) {
    printf("Error! wtalg = %d\n",wtalg);
    gaexit(10);
  }

  //LINE_STAMP;
  GetFromServer_dbl(&Coords_wadj,0,i,tempwadj);

  lockind = getlockind(2,dir,bead);
  GA_Lock(lockind);
  GA_Init_fence();
  GetFromServer_dbl(&Coords_wavg,0,i,tempwavg);
  tempwavg[0] += tempwadj[0];
  PutToServer_dbl(&Coords_wavg,0,i,tempwavg);
  GA_Fence();
  GA_Unlock(lockind);

  allind[0] = 0;
  allind[1] = i;
  tempi[0] = 1 + NGA_Read_inc(Coords_nwstack.ga,allind,1);
}
/*-------------------------------------------------------------------------*/
void updtwts(int ind) {
   
  void getfrac(int dir, int bead);
  double lim2, sum;
  int tgt,dir,bead,good,lim,i, mcount,j,n1,n2,mcount2, nav;
  FILE * wtout, *tmpfile;
  char tmpname[30];

  GA_Init_fence();
  lim = 2*beads;
  
  if (wtalg == 1) {

    /* set weight to wavg, reset wavg */
    
    GetFromServer_dbl(&Coords_wavg,0,-1,twavg1);
    GetFromServer_int(&Coords_nwstack,0,-1,tnwstack1);
    
    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	if (tnwstack1[mcount]) {
	  tweight1[mcount] = twavg1[mcount]/(float)tnwstack1[mcount];
	}
	tnwstack1[mcount] = 0;
	twavg1[mcount] = 0.;
	
	weight[dir][bead] = tweight1[mcount]; 
	mcount++;
      }
    }

    PutToServer_dbl(&Coords_wavg,0,-1,twavg1);
    PutToServer_dbl(&Coords_weight,0,-1,tweight1);
    PutToServer_int(&Coords_nwstack,0,-1,tnwstack1);

  } else if (wtalg == 2) {
    printf("Error!  wtalg = 2!\n");
    gaexit(2);
  } else if (wtalg == 3) {

    /* perform matrix calculations, reset counters */
    
    GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
    GetFromServer_int(&Coords_nflux,0,-1,tnflux1);
    GetFromServer_int(&Coords_nflux,1,-1,tnflux2);
    
    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	weight[dir][bead] = tweight1[mcount];
	mcount++;
      }
    }
    
    mcount = 0;
    for (i=0;i<lim;i++) {
      for (j=0;j<lim;j++) {
	if (tnflux1[mcount] != tnflux2[mcount]) {
	  printf("Error!  tnflux discrepancy! i,j = %d,%d: %d, %d\n",i,j,tnflux1[mcount],tnflux2[mcount]);
	  if ((i>0)&&(j>0)) {
	    n1 = abs(tnflux1[mcount]-fmat[i-1][j-1]);
	    n2 = abs(tnflux2[mcount]-fmat[i-1][j-1]);
	  } else {
	    n1 = abs(tnflux1[mcount]-fmat[i+1][j+1]);
	    n2 = abs(tnflux2[mcount]-fmat[i+1][j+1]);
	  }
	  if (n1 < n2) { // overwrite tnflux2 with tnflux1
	    tnflux2[mcount] = tnflux1[mcount];
	    PutToServer_int(&Coords_nflux,1,-1,tnflux2);
	  } else { // overwrite tnflux1 with tnflux2
	    tnflux1[mcount] = tnflux2[mcount];
	    PutToServer_int(&Coords_nflux,0,-1,tnflux1);
	  }
	  printf("tnfluxes are now %d,%d \n", tnflux1[mcount],tnflux2[mcount]);
	}
	fmat[i][j] = tnflux1[mcount];
	mcount++;
      }
    }

    for (i=0;i<lim;i++) {
      fmat[i][i] = 0;
      for (j=0;j<lim;j++) {
	if (i != j) {
	  fmat[i][i] -= fmat[j][i];
	}
      }
    }

    wtmat3();
    if (verb >= 3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Weight update completed\n");
      fclose(mainout);
    }
    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	tweight1[mcount] = weight[dir][bead]; 
	mcount++;
      }
    }
    PutToServer_dbl(&Coords_weight,0,-1,tweight1);
  } else {
    GetFromServer_dbl(&Coords_wadj,0,-1,twavg1);
    GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
    GetFromServer_int(&Coords_on,0,-1,ton1);
    mcount = 0;
    sum = 0.;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	sum += twavg1[mcount];
	mcount++;
      }
    }
    mcount = 0;
    for (dir=0;dir<=1;dir++) {
      for (bead=0;bead<=beads-1;bead++) {
	weight[dir][bead] = 0.;
	tweight1[mcount] = 0.;
	if (ton1[mcount]) {
	  tweight1[mcount] += twavg1[mcount]/sum;
	  weight[dir][bead] = tweight1[mcount];
	}
	twavg1[mcount] = 0.;
	mcount++;
      }
    }
    PutToServer_dbl(&Coords_wadj,0,-1,twavg1);
    PutToServer_dbl(&Coords_weight,0,-1,tweight1);
  }
    
  GA_Fence();

  /* write weights */

  sum = 0.;
  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      sum += weight[dir][bead];
    }
  }

  wtout = fopen(wtname,"a");
  fprintf(wtout,"%d ",ind*every*wtupdt);
  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      fprintf(wtout,"%e ",weight[dir][bead]/sum);
    }
  }
  fprintf(wtout,"\n");
  fclose(wtout);
}
/*-------------------------------------------------------------------------*/
void wtmat3() {

  int i, j, k, lim, ind1, n1, n2, ind2, ind3, offset, nreps=2000;
  double sum;
  FILE *matout, *fout2;

  lim = 2*beads;

  if (verb >= 2) {
    fout2 = fopen("imat.dat","w");
    for (i=0;i<=lim-1;i++) {
      for (j=0;j<=lim-1;j++) {
	if (fmat[i][j] > 0) {
	  fprintf(fout2,"%s","+");
	} else if (fmat[i][j] < 0) {
	fprintf(fout2,"%s","-");
	} else {
	  fprintf(fout2,"%s"," ");
	}
      }
      fprintf(fout2,"\n");
    }
    fclose(fout2);
  }

  if (verb >= 2) {
    matout = fopen("wtconv.dat","w");
  }
  
  for (k=0;k<nreps;k++) {
    if (verb >= 2) {
      fprintf(matout,"%d ", k);
    }
    //    offset = (int)lim*ran2(&seed);    
    offset = 0;
    for (ind3=0;ind3<lim;ind3++) {
      ind1 = ind3 + offset;
      if (ind1 >= lim) {
	ind1 -= lim;
      }
      i = ind1*2/lim;
      if (i) {
	j = ind1-lim/2;
      } else {
	j = ind1;
      }
      weight[i][j] = 0.;
      for (ind2=0;ind2<lim;ind2++) {
	if ((ind2 != ind1) && (fmat[ind1][ind2])) {
	  n1 = ind2*2/lim;
	  if (n1) {
	    n2 = ind2-lim/2;
	  } else {
	    n2 = ind2;
	  }
	  weight[i][j] -= (float)fmat[ind1][ind2]*weight[n1][n2]/(float)fmat[ind1][ind1];
	}
      }
    }
    if (verb >=2) {
      for (i=0;i<2;i++) {
	for (j=0;j<beads;j++) {
	  fprintf(matout,"%e ",weight[i][j]);
	}
      }
      fprintf(matout,"\n");
    }
  }
  if (verb >=2) {
    fclose(matout);
  }
  sum = 0.;
  for (i=0;i<2;i++) {
    for (j=0;j<beads;j++) {
      sum += weight[i][j];
    }
  }
  for (i=0;i<2;i++) {
    for (j=0;j<beads;j++) {
      weight[i][j] = weight[i][j]/sum;
    }
  }
}
void wrextend(int index) {
    
  /* this function writes out the conditional and total theta histograms 
    for the normal order parameter and the sig-modified o.p.*/

  int dir, bead, ind, mcount, mcount2, tsum;
  double x[tres], x2[tres], x3[tres], sum, sum2, sum3, ext, ootsum;
  FILE * tout;
  char tmpname[30];
  
  GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
  GetFromServer_int(&Coords_thist,0,-1,tthist1);
  GetFromServer_int(&Coords_thistsig,0,-1,tthist3);
  GetFromServer_int(&Coords_thistfr,0,-1,tthist5);
  
  mcount = 0;
  mcount2 = 0;
  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      weight[dir][bead] = tweight1[mcount];
      mcount++;
      for (ind=0;ind<tres;ind++) {
	thist[dir][bead][ind] = tthist1[mcount2];
	thist2[dir][bead][ind] = tthist3[mcount2];
	thist3[dir][bead][ind] = tthist5[mcount2];
	mcount2++;
      }
    }
  }

  /* write 1st dir histogram */
  
  for (ind=0;ind<=tres-1;ind++) {
    x[ind] = 0.;
    x2[ind] = 0.;
    x3[ind] = 0.;
  }

  dir = 0;
  for (bead=0;bead<=beads-1;bead++) {
      
    /* thist */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist[dir][bead][ind]) {
	  x[ind] += thist[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }

    /* thistsig */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist2[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist2[dir][bead][ind]) {
	  x2[ind] += thist2[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }

    /* thistfr */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist3[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist3[dir][bead][ind]) {
	  x3[ind] += thist3[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }
  }

  sum = 0.;
  sum2 = 0.;
  sum3 = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += x[ind];
    sum2 += x2[ind];
    sum3 += x3[ind];
  }

  sprintf(tmpname,"%s%d%s","thist1_",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[0] + (tmax[0]-tmin[0])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x[ind]/sum);
  }
  fclose(tout);

  sprintf(tmpname,"%s%d%s","t2hist1_",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[1] + (tmax[1]-tmin[1])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x2[ind]/sum2);
  }
  fclose(tout);

  sprintf(tmpname,"%s%d%s","t3hist1_",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[2] + (tmax[2]-tmin[2])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x3[ind]/sum3);
  }
  fclose(tout);

  /* write 2nd dir histograms */
  
  for (ind=0;ind<=tres-1;ind++) {
    x[ind] = 0.;
    x2[ind] = 0.;
    x3[ind] = 0.;
  }

  dir = 1;
  for (bead=0;bead<=beads-1;bead++) {
      
    /* thist */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist[dir][bead][ind]) {
	  x[ind] += thist[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }

    /* thistsig */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist2[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist2[dir][bead][ind]) {
	  x2[ind] += thist2[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }

    /* thistfr */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist3[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist3[dir][bead][ind]) {
	  x3[ind] += thist3[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }
  }

  sum = 0.;
  sum2 = 0.;
  sum3 = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += x[ind];
    sum2 += x2[ind];
    sum3 += x3[ind];
  }

  sprintf(tmpname,"%s%d%s","thist2_",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[0] + (tmax[0]-tmin[0])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x[ind]/sum);
  }
  fclose(tout);

  sprintf(tmpname,"%s%d%s","t2hist2_",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[1] + (tmax[1]-tmin[1])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x2[ind]/sum2);
  }
  fclose(tout);

  sprintf(tmpname,"%s%d%s","t3hist2_",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[2] + (tmax[2]-tmin[2])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x3[ind]/sum3);
  }
  fclose(tout);

  /* write total histograms */
  
  for (ind=0;ind<=tres-1;ind++) {
    x[ind] = 0.;
    x2[ind] = 0.;
    x3[ind] = 0.;
  }

  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      
      /* thist */
      tsum = 0;
      for (ind=0;ind<=tres-1;ind++) {
	tsum += thist[dir][bead][ind];
      }
      if (tsum) {
	ootsum = 1./(float)tsum;
	for (ind=0;ind<=tres-1;ind++) {
	  if (thist[dir][bead][ind]) {
	    x[ind] += thist[dir][bead][ind]*weight[dir][bead]*ootsum;
	  }
	}
      }
      
      /* thistsig */
      tsum = 0;
      for (ind=0;ind<=tres-1;ind++) {
	tsum += thist2[dir][bead][ind];
      }
      if (tsum) {
	ootsum = 1./(float)tsum;
	for (ind=0;ind<=tres-1;ind++) {
	  if (thist2[dir][bead][ind]) {
	    x2[ind] += thist2[dir][bead][ind]*weight[dir][bead]*ootsum;
	  }
	}
      }

      /* thistfr */
      tsum = 0;
      for (ind=0;ind<=tres-1;ind++) {
	tsum += thist3[dir][bead][ind];
      }
      if (tsum) {
	ootsum = 1./(float)tsum;
	for (ind=0;ind<=tres-1;ind++) {
	  if (thist3[dir][bead][ind]) {
	    x3[ind] += thist3[dir][bead][ind]*weight[dir][bead]*ootsum;
	  }
	}
      }
    }
  }
    
  sum = 0.;
  sum2 = 0.;
  sum3 = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += x[ind];
    sum2 += x2[ind];
    sum3 += x3[ind];
  }

  sprintf(tmpname,"%s%d%s","thist",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[0] + (tmax[0]-tmin[0])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x[ind]/sum);
  }
  fclose(tout);

  sprintf(tmpname,"%s%d%s","t2hist",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[1] + (tmax[1]-tmin[1])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x2[ind]/sum2);
  }
  fclose(tout);

  sprintf(tmpname,"%s%d%s","t3hist",index,".dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin[2] + (tmax[2]-tmin[2])*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,x3[ind]/sum3);
  }
  fclose(tout);

  mcount2 = 0;





  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      for (ind=0;ind<tres;ind++) {
	tthist1[mcount2] = 0;
	mcount2++;
      }
    }
  }

  PutToServer_int(&Coords_thist,0,-1,tthist1);
  PutToServer_int(&Coords_thistsig,0,-1,tthist1);
  PutToServer_int(&Coords_thistfr,0,-1,tthist1);

}
/*-------------------------------------------------------------------------*/
int which(int dir,int bead) {

  /* this function returns the region that the walker
   in (dir,bead) is in, presently */
       
  int i, j, temp, good;
  double dist, lodist;
  iopoint_t iox;

  good = 1;
   for (i=0;i<nop;i++) {
     iox.x[i] = (int)((cryst.x[i]-ocoor.x[i])/opstep[dir]);
     if (iox.x[i] != myreg[dir][bead][i]) {
       good = 0;
     }
   }

   if (!good) {
     temp = findreg(iox,dir);
   } else {
     temp = bead;
   }

  return temp;
}

/*-------------------------------------------------------------------------*/
int whichx(int dir, point_t x) {
   
   /* this function returns the region that the point 'x'
      is in, presently */
   
   int i, j, temp;
   double dist, lodist;
   opoint_t ox;
   iopoint_t iox;
   
   ox = projx(x);

   for (i=0;i<nop;i++) {
     iox.x[i] = (int)((cryst.x[i]-ox.x[i])/opstep[dir]);
   }

   temp = findreg(iox,dir);

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
void addpoint (int dir, int bead, point_t x, double wt, int from1, int from2) {

  /* this function adds the point 'x' to the (dir,bead) fluxlist
     along with the weight 'wt' */
  
  int tgt, n, i, ind, tfull, wind, part, dim, ptind, mcount, tn, tfrom2[2], listind;
  double tw[1], ttw[1], t1,t2, extra, newextra, tmpwlist[mxlist], sum, oosum;
  
  ind = dir*beads + bead;

  if (((from1 <0)||(from1>=beads))||((from2<0)||(from2>=beads))) {
    listind = 0;
  } else {
    listind = getlistind(dir,bead,from1,from2);
  }

  allind[0] = listind;
  allind[1] = ind;
  
  lockind = getlockind(0,dir,bead);
  GA_Lock(lockind);
  GA_Init_fence();
  GetFromServer_int(&Coords_full,listind,ind,tempi);  
  tfull = tempi[0];

  tn = 1 + NGA_Read_inc(Coords_nlist.ga,allind,1);
  if (tn == mxlist+1) {
    tempi[0] = 1;
    PutToServer_int(&Coords_nlist,listind,ind,tempi);
    tn = 1;
    if (!tfull) {   /* not full */
      tempi2[0] = 1;
      PutToServer_int(&Coords_full,listind,ind,tempi2);
      tfull = 1;
    }
  }
  GA_Fence();
  GA_Unlock(lockind);     // you now have guaranteed unique values of tn
  // and correct values of tfull
    
  wind = mxlist*beads*dir + mxlist*bead + (tn - 1);

  if (tfull) {

    lockind = getlockind(1+listind,dir,bead);
    GA_Lock(lockind);
    GA_Init_fence();
    
    lo[0] = listind;
    hi[0] = listind;
    lo[1] = mxlist*beads*dir + mxlist*bead;
    hi[1] = mxlist*beads*dir + mxlist*(bead + 1) - 1;
    
    NGA_Get(Coords_wlist.ga,lo,hi,tmpwlist,ld);

    extra = tmpwlist[tn-1];
    tmpwlist[tn-1] = wt;
    sum = 0.;
    for (i=0;i<mxlist;i++) {
      sum += tmpwlist[i];
    }
    oosum = 1./sum;
    sum = sum + extra;
    oosum = sum * oosum;
    for (i=0;i<mxlist;i++) {
      tmpwlist[i] = tmpwlist[i]*oosum;
    }
    NGA_Put(Coords_wlist.ga,lo,hi,tmpwlist,ld);
    ttw[0] = sum;
    PutToServer_dbl(&Coords_twlist,listind,ind,ttw);
    if (verb >=4) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: Adding weight to (%d,%d,%d):  %e + %e = %e\n",t2,dir,bead,listind,ttw[0]-wt,wt,ttw[0]); 
      fclose(mainout);
    }

    // flist[dir][bead].twlist += wt - flist[dir][bead].wlist[n];
    GA_Fence();
    GA_Unlock(lockind);

    //flist[dir][bead].wlist[n] = wt;      

  } else {
    
    lockind = getlockind(1+listind,dir,bead);
    GA_Lock(lockind);
    GA_Init_fence();

    tw[0] = wt;
    PutToServer_dbl(&Coords_wlist,listind,wind,tw);
    //flist[dir][bead].wlist[n] = wt;

    GetFromServer_dbl(&Coords_twlist,listind,ind,ttw);
    ttw[0] += wt;
    PutToServer_dbl(&Coords_twlist,listind,ind,ttw);
    if (verb >=4) {
      t2 = MPI_Wtime();
      mainout = fopen(outname,"a");
      fprintf(mainout,"%f: Adding weight to (%d,%d,%d):  %e + %e = %e\n",t2,dir,bead,listind,ttw[0]-wt,wt,ttw[0]); 
      fclose(mainout);
    }
    GA_Fence();
    GA_Unlock(lockind);

    //flist[dir][bead].twlist += wt;
  }

  mcount = 0;
  for (part=0; part<npart;part++) {
    for (dim=0; dim<=ndim-1; dim++) {
      tpt[mcount] = x.x[part][dim];  // position
      mcount++;
      tpt[mcount] = x.v[part][dim];  // velocity
      mcount++;
      tpt[mcount] = f1[part][dim];  // old force
      mcount++;
    }
  }	    
    
  lo[0] = listind;
  hi[0] = listind;
  lo[1] = 3*npart*ndim*mxlist*beads*dir + 3*npart*ndim*mxlist*bead + 3*npart*ndim*(tn-1);
  hi[1] = 3*npart*ndim*mxlist*beads*dir + 3*npart*ndim*mxlist*bead + 3*npart*ndim*tn -1;
  NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);
  //  flist[dir][bead].pts[n] = x; 

  lo[1] = 2*mxlist*beads*dir + 2*mxlist*bead + 2*(tn-1);
  hi[1] = lo[1] + 1;
  tfrom2[0] = from1;
  tfrom2[1] = from2;
  NGA_Put(Coords_from.ga,lo,hi,tfrom2,ld);
  // flist[dir][bead].from[n][0] = from1;
  // flist[dir][bead].from[n][1] = from2;

  ind = mxlist*beads*dir + mxlist*bead + (tn-1);
  tempi[0] = myclock;
  PutToServer_int(&Coords_clock,listind,ind,tempi);
  tempi[0] = mytim;
  PutToServer_int(&Coords_mytim,listind,ind,tempi);
  // flist[dir][bead].clock[n] = myclock;
  
  ind = dir*beads + bead;
  GetFromServer_int(&Coords_on,0,ind,tempi);
  
  if (!tempi[0]) {
    printf("proc %d:  Turning %d %d back on.\n", me, dir, bead);
    tempi[0] = 1;
    PutToServer_int(&Coords_on,0,ind,tempi);
  }
}
/*-------------------------------------------------------------------------*/
 double calcsum(int dir) { 

   double temp;
   int lim, i;

   temp = 0.;
   for (i=0;i<beads;i++) {
     temp = temp + weight[dir][i];
   }
   
   return temp;
 }
/*-------------------------------------------------------------------------*/
 void stack(int dir, int bead) { 
   
   /* this function stacks weights for wtalg=4, and stacks histogram values for thist */
   
   int i, j, index, ind;
   double temps[1], mx, mn, tw[1];
   point_t rt;
   double sum;

   lockind = getlockind(3,dir,bead);
   ind = dir*beads + bead;
   if (wtalg == 4) {
     GA_Lock(lockind);
     GA_Init_fence();
     GetFromServer_dbl(&Coords_wadj,0,ind,tw);
     tw[0] += myweight;
     PutToServer_dbl(&Coords_wadj,0,ind,tw);
     GA_Fence();
     GA_Unlock(lockind);     
   }

   sum = 0.;
   for (i=0;i<nop;i++) {
     sum += ocoor.x[i];
   }
   index = tres*(sum-tmin[0])/(tmax[0]-tmin[0]);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index (norm)! %f\n",sum);
   } else {
     allind[0] = 0;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_thist.ga,allind,1);
   }

   sum = projxsig(coor);
   index = tres*(sum-tmin[1])/(tmax[1]-tmin[1]);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index (sig)! %f\n",sum);
   } else {
     allind[0] = 0;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_thistsig.ga,allind,1);
   }

   temps[0] = 0.;
   for (i=0;i<ndim;i++) {
     temps[0] += (coor.x[0][i]-coor.x[npart-1][i])*(coor.x[0][i]-coor.x[npart-1][i]);
   }
   temps[0] = sqrt(temps[0]);
   index = tres*(temps[0]-tmin[2])/(tmax[2]-tmin[2]);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index (fr)! %f\n",temps[0]);
   } else {
     allind[0] = 0;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_thistfr.ga,allind,1);
   }
 }
/*-------------------------------------------------------------------------*/
void getfrac(int dir, int bead) {

  /* this function saves the amount of weight on the (dir,bead) 
     flux list that came from region (dir2,bead2) [in the lfrac array],
     for all dir2 and bead2 */

  double temp, *tmpwt;
  int i, j, ind, tfull, wind, *tmpfrom2, mcount, mcount2, limit;

  ind = beads*dir + bead;
  //LINE_STAMP;
  GetFromServer_int(&Coords_full,0,ind,tempi);
  tfull = tempi[0];
  if (tfull) {
    limit = mxlist;
  } else {
    //LINE_STAMP;
    GetFromServer_int(&Coords_nlist,0,ind,tempi);
    limit = tempi[0];
  }
  mcount = 2*beadsp2*dir + 2*beads*bead;
  for (i=0;i<2;i++) {
    for (j=0;j<beads;j++) {
      lfrac[mcount] = 0.;
      mcount++;
    }
  }
  if (limit != 0) {
    tmpfrom2 = (int *) malloc(2*limit*sizeof(int));
    tmpwt = (double *) malloc(limit*sizeof(double));

    lo[0] = 0;
    hi[0] = 0;
    lo[1] = 2*mxlist*beads*dir + 2*mxlist*bead;
    hi[1] = lo[1] + 2*limit - 1;
    if (verb >= 3) {
      printf("proc %d: NGA from limits: (%d,%d) -> (%d,%d); (dir,bead) = (%d,%d)\n",me,lo[0],lo[1],hi[0],hi[1],dir,bead); fflush(stdout);
    }
    NGA_Get(Coords_from.ga,lo,hi,tmpfrom2,ld);
    
    lo[1] = mxlist*beads*dir + mxlist*bead;
    hi[1] = lo[1] + limit - 1;
    if (verb >= 3) {
      printf("proc %d: NGA wt limits: (%d,%d) -> (%d,%d); (dir,bead) = (%d,%d)\n",me,lo[0],lo[1],hi[0],hi[1],dir,bead); fflush(stdout);
    }
    NGA_Get(Coords_wlist.ga,lo,hi,tmpwt,ld);
    
    mcount = 0;
    for (i=0;i<=limit-1;i++) {
      tempi[0] = tmpfrom2[mcount];
      mcount++;
      tempi2[0] = tmpfrom2[mcount];
      mcount++;
      
      mcount2 = 2*beadsp2*dir + 2*beads*bead + beads*tempi[0] + tempi2[0];
      lfrac[mcount2] += tmpwt[i];
    }

    free(tmpfrom2);
    free(tmpwt);
  }
}
/*-------------------------------------------------------------------------*/
int getmin(int *dir, int *bead) {
  
  /* this function returns the (dir,bead) for the region with the 
     lowest counter */
  
  int low, i, j, k, mcount, ind, start, limit, mind, sum, gotit, ti, tj, tk;
  double rnd;

  GetFromServer_int(&Coords_count,0,-1,tint1);
  GetFromServer_int(&Coords_on,0,-1,ton1);

  low = every;
  *dir = -1;
  *bead = 0;

  mcount = 0;
  sum = 0;
  for (j=0;j<2;j++) {
    for (k=0;k<beads;k++) {
      count[j][k] = tint1[mcount];
      on[j][k] = ton1[mcount];
      if (on[j][k] && (count[j][k] < every)) {
	sum += every-count[j][k];
      }
      mcount++;
    }
  }

  if (sum > 0) { // sum is the total number of steps outstanding
    rnd = (float)sum*ran2(&seed);
    
    sum = 0;
    gotit = 0;
    
    for (j=0;j<2;j++) {
      for (k=0;k<beads;k++) {
	if (!gotit) {
	  if (on[j][k] && (count[j][k] < every)) {
	    sum += every-count[j][k];
	    if ((float)sum > rnd) { 
	      gotit = 1;
	      tj = j;
	      tk = k;
	    }
	  }
	}
      }
    }
    if (!gotit) {
      printf("proc %d: sum error (sum,rnd) = (%d,%f) \n",me,sum,rnd); fflush(stdout);
      tj = 1;
      tk = beads-2;
    }
    *dir = tj;
    *bead = tk;
    low = count[tj][tk];
  } else {
    *dir = -1;           /* all regions have run for 'every' steps */
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

  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      fscanf(wtin,"%e",&temp);
      weight[j][k] = temp;
    }
  }
  fclose(wtin);

  mcount = 0;
  for (j=0;j<2;j++) {
    for (k=0;k<beads;k++) {
      tweight1[mcount] = weight[j][k];
      twavg1[mcount] = 0.;
      mcount++;
    }
  }

  PutToServer_dbl(&Coords_weight,0,-1,tweight1);
  PutToServer_dbl(&Coords_wadj,0,-1,tweight1);
  PutToServer_dbl(&Coords_wavg,0,-1,twavg1);

}
/*-------------------------------------------------------------------------
--------------------------------------------
the following functions are system-dependent
--------------------------------------------
-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
void fxread(int find) {

  /* this function reads the initial values of the fluxes
     from the input file 'fxname' */

  int j, k, op, dir, ti, tj, tk, limit, pt, dim, part, ind, mcount, listind, from1, from2;
  int scale;
  double fac, wt;
  FILE * fxin;
  float temp;
  float tw;
  point_t tx;
  char fxname1[30];
  
  if (find == -1) {
    sprintf(fxname1,"%s",fxname);
  } else {    
    sprintf(fxname1,"%s%d%s",fxname,find,".dat");
  }
  fxin = fopen(fxname1,"r");

  while(fscanf(fxin,"%d %d",&tj,&tk) != EOF ){
      
    if (verb >= 4) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"reading %d %d\n",tj,tk);
      fclose(mainout);
    }

    for (listind=0;listind<2;listind++) {
	  
      fscanf(fxin,"%d %le", &tempi[0], &wt);
      limit = tempi[0];

      if (verb >= 4) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"reading %d points from listind %d \n",tempi[0],listind);
	fclose(mainout);
      }
	
      for (pt=0; pt<=limit-1; pt++) {
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"point %d..",pt);
	  fclose(mainout);
	}
	fscanf(fxin,"%d %d ", &tempi[0], &tempi2[0]);
	from1 = tempi[0];
	from2 = tempi2[0];
	fscanf(fxin,"%le %d ",&wt,&tempi[0]);
	fscanf(fxin,"%d ",&tempi[0]);
	
	mcount = 0;
	for (part=0; part<npart;part++) {
	  for (dim=0; dim<=ndim-1; dim++) {
	    fscanf(fxin,"%f ", &tw);
	    tpt[mcount] = tw;
	    tx.x[part][dim] = tw;
	    mcount++;
	    
	    fscanf(fxin,"%f ", &tw);
	    tpt[mcount] = tw;
	    tx.v[part][dim] = tw;
	    mcount++;
	    
	    fscanf(fxin,"%f ", &tw);
	    tpt[mcount] = tw;
	    mcount++;
	    f1[part][dim] = tw;
	  }
	}	    
	  
	dir = 1;
	ocoor = projx(tx);
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"ocoor..");
	  fclose(mainout);
	}
	bas = bascheck();
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"bas..");
	  fclose(mainout);
	}
	if (bas != -1) {
	  dir = bas;
	}
	k = whichx(dir,tx);
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"dir=%d, bead = %d..\n",k,dir);
	  fclose(mainout);
	}
	if (wt > pmin) {
	  //	  printf("%d: Adding point from rdflux to %d %d with wt = %e\n",me,dir,k,wt);
	  addpoint(dir,k,tx,wt,from1,from2);
	}
      }
    }
  }

  fclose(fxin);
}
/*-------------------------------------------------------------------------*/
void nfluxread() {

  /* this function reads the initial values of nflux
     from the input file 'nfname' */

  int i, j, k, dir, bead, n1, n2, tn, tp, isum, lim;
  FILE * fxin;
  float temp;
  float tw;
  
  fxin = fopen(nfname,"r");
  lim = 2*beads;
  k = 0;
  for (i=0;i<lim;i++) {
    for (j=0;j<lim;j++) {
      fscanf(fxin,"%d ",&tn);
      tnflux1[k] = tn;
      k++;
    }
  }
  PutToServer_int(&Coords_nflux,0,-1,tnflux1);
  PutToServer_int(&Coords_nflux,1,-1,tnflux1);

  fclose(fxin);
}
/*-------------------------------------------------------------------------*/

int bascheck() {
  
  /* this function checks if the walker is in a basin */

  int temp, i, j, sum;
  double dist;

  temp = -1;

  sum = 0;

  for (i=0;i<nop;i++) {
    sum += ocoor.x[i];
  }

  if (sum < basin[0]) {
    temp = 0;
  } else if (sum > basin[1]) {
    temp = 1;
  }

  return temp;
}
/*-------------------------------------------------------------------------*/
double projxsig(point_t x) {
  
  /* this function returns the projection of coor onto the order parameter space, using sigfacsq */
  
  int i, j, k, ind1, ind2, segi, segj;
  opoint_t temp;
  double mn, mx, r2, sum, tdub;
  
  /* loop over contacts */

  sum = 0.;
  for (i=0;i<nsec;i++) {
    ind1 = secind[i][0];
    ind2 = secind[i][1];
    if ((ind1 < npart) && (ind2 < npart)) {
      r2 = 0.;
      for (j=0;j<ndim;j++) {
	r2 += (x.x[ind1][j]-x.x[ind2][j])*(x.x[ind1][j]-x.x[ind2][j]);
      }
      if (r2 < sigfacsq*r02[ind1][ind2]) {
	tdub = 1.;
      } else if (r2 < 4.*sigfacsq*r02[ind1][ind2]) {
	tdub = mypow4(sigfacsq*r02[ind1][ind2]/r2);
      } else {
	tdub = 0.;
      }
      sum += tdub;
    }
  }
  for (i=0;i<ntert;i++) {
    ind1 = tertind[i][0];
    ind2 = tertind[i][1];
    if ((ind1 < npart) && (ind2 < npart)) {
      r2 = 0.;
      for (j=0;j<ndim;j++) {
	r2 += (x.x[ind1][j]-x.x[ind2][j])*(x.x[ind1][j]-x.x[ind2][j]);
      }
      if (r2 < sigfacsq*r02[ind1][ind2]) {
	tdub = 1.;
      } else if (r2 < 4.*sigfacsq*r02[ind1][ind2]) {
	tdub = mypow4(sigfacsq*r02[ind1][ind2]/r2);
      } else {
	tdub = 0.;
      }
      sum += tdub;
    }
  }

  return sum;
}
/*-------------------------------------------------------------------------*/
opoint_t projx(point_t x) {
  
  /* this function returns the projection of coor onto the order parameter space, using rsigfacsq */
  
  int i, j, k, ind1, ind2, op;
  opoint_t temp;
  double mn, mx, r2, sum, tdub;
  
  for (i=0;i<nop;i++) {
    temp.x[i] = 0.;
  }

  /* loop over contacts */

  for (i=0;i<nsec;i++) {
    ind1 = secind[i][0];
    ind2 = secind[i][1];
    if ((ind1 < npart) && (ind2 < npart)) {
      getop(ind1,ind2,&op);
      r2 = 0.;
      for (j=0;j<ndim;j++) {
	r2 += (x.x[ind1][j]-x.x[ind2][j])*(x.x[ind1][j]-x.x[ind2][j]);
      }
      if (r2 < rsigfacsq*r02[ind1][ind2]) {
	tdub = 1.;
      } else if (r2 < 4.*rsigfacsq*r02[ind1][ind2]) {
	tdub = mypow4(rsigfacsq*r02[ind1][ind2]/r2);
      } else {
	tdub = 0.;
      }
      temp.x[op] += tdub;
    }
  }
  for (i=0;i<ntert;i++) {
    ind1 = tertind[i][0];
    ind2 = tertind[i][1];
    if ((ind1 < npart) && (ind2 < npart)) {
      getop(ind1,ind2,&op);
      r2 = 0.;
      for (j=0;j<ndim;j++) {
	r2 += (x.x[ind1][j]-x.x[ind2][j])*(x.x[ind1][j]-x.x[ind2][j]);
      }
      if (r2 < rsigfacsq*r02[ind1][ind2]) {
	tdub = 1.;
      } else if (r2 < 4.*rsigfacsq*r02[ind1][ind2]) {
	tdub = mypow4(rsigfacsq*r02[ind1][ind2]/r2);
      } else {
	tdub = 0.;
      }
      temp.x[op] += tdub;
    }
  }

  return temp;
}
/*-------------------------------------------------------------------------*/
void move() {
  
  /* moves coor, updates ocoor */

  int i;
  point_t oldx, rt;
  double psi[2],phi[2],fx,fy,temp1,temp2,dvdx,dvdy,t2;
  void dostep(), stream(), dorotate();
  
  if (verb >= 4) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"dostep...");
    fclose(mainout);
  }
  dostep();
  if (!polyerr) {
    if (myclock%strmfrq == 0) {
      if (verb >= 4) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"stream...");
	fclose(mainout);
      }
      stream();
      if (verb >= 4) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"dorotate...");
	fclose(mainout);
      }
      dorotate();
    }
    if (verb >= 4) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"project...\n");
      fclose(mainout);
    }
    ocoor = projx(coor);  /* update o.p. projection */
  }
}
/*-------------------------------------------------------------------------*/
void dostep(){
  int i, j, k, x;
  double a, f2[npart][ndim],temp,x1,x2,x3;
  point_t f, force();
    
#pragma omp parallel for default(shared) private(i,x,a) schedule(static)
  for(i=0; i< npart; i++){
    for(x=0; x< 3; x++){	
      /* calculate a(t) = F/m */        
      // forces are pre-divided by polymass to save calcs
      coor.x[i][x] += 0.5*tstepsq*f1[i][x]/polymass + tstep*coor.v[i][x];
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

#pragma omp parallel for default(shared) private(i,x,x1,x2,x3) schedule(static)
  for(i=0; i< npart; i++){
    for(x=0; x< 3; x++){
      f2[i][x] = f.x[i][x];
      coor.v[i][x] += 0.5*(f1[i][x]+f2[i][x])*tstep/polymass;
      f1[i][x] = f2[i][x];
    }
    x1=coor.x[i][0];
    x2=coor.x[i][1];
    x3=coor.x[i][2];
    if(((x1<0.)||(x1>Lx)) || (((x2<0.)||(x2>Ly)) || ((x3<0.)||(x3>Lz)))) {
      printf("polymer bead %d out of bounds! (%f,%f,%f)\n",i,x1,x2,x3);
      polyerr = 1;
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
		  
		  if( cutoff > 0. ){  // !seclist and !terlist
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

  for (p1=0;p1<npart;p1++) {
    if (coor.x[p1][1] < wcacut) { // add WCA repulsion if bead is close to bottom wall
      pre = mypow6(wcasig/coor.x[p1][1]);
      f.x[p1][1] += 4.*ep1_6/coor.x[p1][1]*(2.*pre*pre - pre);
    }
  }

  /* pre-scale the forces */
  /*
  for (p1=0;p1<npart;p1++) {
    for (p2=0;p2<ndim;p2++) {
      f.x[p1][p2] = f.x[p1][p2]/polymass;
    }
    }*/

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
  for(i=0;i<npart*3;i++) {
    fscanf(cfile, "%lf", &(val));
    r0[i/3][i%3] = val;
  }
  
  for(i=0; i< npart; i++){
    for(j=0; j< npart; j++){
	dx = r0[i][0]-r0[j][0];
	dy = r0[i][1]-r0[j][1];
	dz = r0[i][2]-r0[j][2];
	r0dists[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
	r02[i][j] = mypow2(r0dists[i][j]);
	r0d6[i][j] = mypow2(r02[i][j])*mypow2(r0dists[i][j]);
	r0d12[i][j] = mypow2(r0d6[i][j]);
    }
  }
}
/*-----------------------------------------------------------------------------*/
void rddelta() {

  /*---------Read in Delta values, or write them out----------------------*/
  /*secondary contacts; if blank, none read*/ 

  FILE *dfile;
  int i, j, val2, ti, ind1, ind2, op;
  
  nsec = 0;
  ntert = 0;

  dfile = fopen(secname,"rt");
  if(dfile == NULL){
    printf("No secondary contacts read, setting them to zero\n");
    for(i=0; i< onpart; i++){
      for(j=0; j< onpart; j++){
	Delta[i][j][0] = 0;
      }
    }
  }else{
    i=0;
    while(fscanf(dfile, "%i", &(val2) ) != EOF ){
      if (val2 == 1) {
	nsec++;
      }
      i++;
    }
    nsec /= 2;
    
    secind = (int **) alloc2d(sizeof(int), nsec, 2);
    
    if( i != onpart*onpart ){
      printf("First contact list is different size than polymer: %i %i\n", i, onpart*onpart);
      gaexit(0);
    }
    fclose(dfile);
    dfile=fopen(secname,"rt");
    i=0;
    ti=0;
    while(fscanf(dfile, "%i", &(val2) ) != EOF){
      ind1 = i/onpart;
      ind2 = i%onpart;
      Delta[ind1][ind2][0] = (int) val2;
      if ((val2 == 1) && (ind1 > ind2)) {
	secind[ti][0] = ind1;
	secind[ti][1] = ind2;
	ti++;
      }
      i++;
    }
  }
  
  /*	tertiary contacts, if blank, none read*/
  dfile = fopen(tertname,"rt");
  if(dfile == NULL){
    printf("No tertiary concacts read, setting them to zero\n");	
    for(i=0; i< onpart; i++){
      for(j=0; j< onpart; j++){
	Delta[i][j][1] = 0; 
      }
    }
  }else{
    i=0;	
    while(fscanf(dfile, "%i", &(val2) ) != EOF ){		
      if (val2 == 1) {
	ntert++;
      }
      i++;								
    }		
    ntert /= 2;
    
    tertind = (int **) alloc2d(sizeof(int), ntert, 2);

    if( i != onpart*onpart ){			
      printf("Second contact list is different size than polymer: %i %i\n", i, onpart*onpart);
      gaexit(0);
    }
    fclose(dfile);
    dfile = fopen(tertname,"rt");			
    i=0;
    ti=0;							
    while(fscanf(dfile, "%i", &(val2) ) != EOF){
      ind1 = i/onpart;
      ind2 = i%onpart;
      Delta[ind1][ind2][1] = (int) val2;
      if ((val2 == 1) && (ind1 > ind2)) {
	tertind[ti][0] = ind1;
	tertind[ti][1] = ind2;
	ti++;
      }
      i++;			
    }
  }

  for (i=0;i<nop;i++) {
    cryst.x[i] = 0.;
  }

  for (i=0;i<npart;i++) {
    for (j=0;j<npart;j++) {
      getop(i,j,&op);
      cryst.x[op] += (float)(Delta[i][j][0] + Delta[i][j][1]);
    }
  }

  if (me==0) {
    printf("crystal contacts:\n");
    for (i=0;i<nop;i++) {
      printf("%d %d %f \n",i,j,cryst.x[i]);
    }
  }
}
/*-----------------------------------------------------------------------------*/
void rdsolv(){

  FILE *filein;
  int j, k, good, ind;
  double temp, nudge=1.e-4, rx;
  float tw;
  filein = fopen(sname,"r");

  for (j=0; j<N; j++) {
    for (k=0; k<3; k++) {
      fscanf(filein,"%f",&tw);
      rsolv[j][k] = tw;
    }
    for (k=0; k<3; k++) {
      fscanf(filein,"%f",&tw);
      vsolv[j][k] = tw;
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
  int j, part, dim, notyet[2][beads], i, k, dir, bas;
  float tw;
  opoint_t tx;
  point_t px;

  filein = fopen(emname,"r");

  for (j=0;j<nemerg;j++) {
    for (part=0;part<npart;part++) {
      for (dim=0;dim<ndim;dim++) {
	fscanf(filein,"%e ",&tw);
	px.x[part][dim] = tw;
	fscanf(filein,"%e ",&tw);
	px.v[part][dim] = tw;
	fscanf(filein,"%e ",&tw);
	f1[part][dim] = tw;
      }
    }

    dir = 1;
    ocoor = projx(px);
    bas = bascheck();
    if (bas != -1) {
      dir = bas;
    }
    k = whichx(dir,px);
    printf("Adding point from rdemerg to %d %d\n",dir,k);
    addpoint(dir,k,px,1./(float)nemerg,-1,-1);
  }
  fclose(filein);
}
/*-----------------------------------------------------------------------------*/
void stream() {
  
  /* This function streams the solvent particles forward */
  
  int i, j, good;
  double temp3, tstar, sum, fac;
  
  sum =0.;
  for(i=0;i<N;i++){
    
    temp3=rsolv[i][1]+vsolv[i][1]*solvstep;         /*first move along y coordinate as it has fixed wall*/

    if (temp3 < 0.) {
      // what time did it cross?; evolve
      tstar = -rsolv[i][1]/vsolv[i][1];

      for(j=0;j<3;j++){
	rsolv[i][j]+=vsolv[i][j]*tstar;
      }

      // pick random velocities
      for (j=0;j<3;j++) {
	//	vsolv[i][j]=sqrt(-var*2.*log(ran2(&seed)))*cos(2*pi*ran2(&seed));
	vsolv[i][j] = -vsolv[i][j];
      }
      if (vsolv[i][1] < 0.) {
	vsolv[i][1] *= -1.;
      }
      for(j=0;j<3;j++){
	rsolv[i][j]+=vsolv[i][j]*(solvstep-tstar);
      }
    } else if (temp3 > Ly) {
      // what time did it cross?; evolve
      tstar = (Ly-rsolv[i][1])/vsolv[i][1];

      for(j=0;j<3;j++){
	rsolv[i][j]+=vsolv[i][j]*tstar;
      }

      // pick random velocities
      for (j=0;j<3;j++) {
	//	vsolv[i][j]=sqrt(-var*2.*log(ran2(&seed)))*cos(2*pi*ran2(&seed));
	vsolv[i][j] = -vsolv[i][j];
      }
      if (vsolv[i][1] > 0.) {
	vsolv[i][1] *= -1.;
      }
      for(j=0;j<3;j++){
	rsolv[i][j]+=vsolv[i][j]*(solvstep-tstar);
      }
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
	good = 0;
      }
      if((rsolv[i][1]<0) || (rsolv[i][1]>Ly)){
	printf("Error (Ly): i=%d, x,y,z=(%f,%f,%f)!\n",i,rsolv[i][0],rsolv[i][1],rsolv[i][2]);
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
    for (j=0;j<3;j++) {
      sum += vsolv[i][j]*vsolv[i][j];
    }
  }
  sum = (sum/(float)N)*solvmass;
  fac = sqrt(3./(beta*sum));
  for (i=0;i<N;i++) {
    for (j=0;j<3;j++) {
      vsolv[i][j] *= fac;
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
		if(blinv*rsolv[i][0]>deltax[0]){
		    pos[i][0]=(int)(blinv*rsolv[i][0]-deltax[0])+1;
		    if (pos[i][0] == nx) {
			pos[i][0] = 0;
		    }
		} else {
		    pos[i][0]=0;
		}
		j = pos[i][0];
		if(blinv*rsolv[i][1]>deltax[1]){
		    pos[i][1]=(int)(blinv*rsolv[i][1]-deltax[1])+1;
		} else {
		    pos[i][1]=0;
		}
		k = pos[i][1];
		if(blinv*rsolv[i][2]>deltax[2]){
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
	    //	    printf("proc %d: Increasing maxd to %d\n",me,maxd);
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
	if(blinv*coor.x[i][0]>deltax[0]){
	    pos_poly[i][0]=(int)(blinv*coor.x[i][0]-deltax[0])+1;
	} else {
	    pos_poly[i][0]=0;
	}
	if(blinv*coor.x[i][1]>deltax[1]){
	  pos_poly[i][1]=(int)(blinv*coor.x[i][1]-deltax[1])+1;
	} else {    
	  pos_poly[i][1]=0;
	}
	if(blinv*coor.x[i][2]>deltax[2]){
	    pos_poly[i][2]=(int)(blinv*coor.x[i][2]-deltax[2])+1;
	} else {
	    pos_poly[i][2]=0;
	}

	if ((pos_poly[i][0] < 0) || (pos_poly[i][0] >= nx)) {
	  err = 1;
	  erri = i;
	}
	if ((pos_poly[i][1] < 0) || (pos_poly[i][1] >=ny)) {
	  err = 1;
	  erri = i;
	}
	if ((pos_poly[i][2] < 0) || (pos_poly[i][2] >=nz)) {
	  err = 1;
	  erri = i;
	}

    } /* end of parallel region */


    if (err) {
	printf("Polymer outside of box! x[%d] = (%d,%d,%d)\n",erri,pos_poly[erri][0],pos_poly[erri][1],pos_poly[erri][2]);	      
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
	
	/*Fill up boxes at y=0 and y=Ly with fake velocities*/
	
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
	    while ((int)n_av>box[i][ny-1][k]) {
	      box[i][ny-1][k]++;
	      
	      temp4=sqrt(-var*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
	      v_cmax[i][ny-1][k]+=solvmass*temp4;
	      
	      temp4=sqrt(-var*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
	      v_cmay[i][ny-1][k]+=solvmass*temp4;
	      
	      temp4=sqrt(-var*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
	      v_cmaz[i][ny-1][k]+=solvmass*temp4;
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
				if ((j!=0)&&(j!=ny-1)) {
				  vsolv[ind][0] += solvstep*grav;  
				}
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
  printf("Terminating GA..\n"); fflush(stdout);
  status = DestroyCoordServer(&Coords_wavg);
  status = DestroyCoordServer(&Coords_wadj);
  status = DestroyCoordServer(&Coords_weight);
  status = DestroyCoordServer(&Coords_elaps);

  status = DestroyCoordServer(&Coords_count);
  status = DestroyCoordServer(&Coords_narr);
  status = DestroyCoordServer(&Coords_on);
  status = DestroyCoordServer(&Coords_bdedit);
  status = DestroyCoordServer(&Coords_nflux);


  status = DestroyCoordServer(&Coords_thist);
  status = DestroyCoordServer(&Coords_thistsig);
  status = DestroyCoordServer(&Coords_thistfr);
  status = DestroyCoordServer(&Coords_full);
  status = DestroyCoordServer(&Coords_nlist);
  status = DestroyCoordServer(&Coords_twlist);
  status = DestroyCoordServer(&Coords_pts);
  status = DestroyCoordServer(&Coords_wlist);
  status = DestroyCoordServer(&Coords_from);
  status = DestroyCoordServer(&Coords_clock);
  status = DestroyCoordServer(&Coords_mytim);
  status = DestroyCoordServer(&Coords_nwstack);
  status = DestroyCoordServer(&Coords_regind);
  status = DestroyCoordServer(&Coords_nreg);

  GA_Terminate();
  printf("Finalizing MPI..\n"); fflush(stdout);
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
  status = CreateCoordServer_dbl(&Coords_twlist, &MaxNum, "test");
  
  status = CreateCoordServer_int(&Coords_elaps, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_count, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_narr, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_on, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_bdedit, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_full, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nlist, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nwstack, &MaxNum, "test");
  
  MaxNum = 2*beads*2*beads;
  status = CreateCoordServer_int(&Coords_nflux, &MaxNum, "test");
  
  MaxNum = 2*beads*tres;
  status = CreateCoordServer_int(&Coords_thist, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_thistsig, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_thistfr, &MaxNum, "test");
  
  MaxNum = 2*beads*mxlist*3*npart*ndim;
  status = CreateCoordServer_dbl(&Coords_pts, &MaxNum, "test");
  MaxNum = 2*beads*mxlist;
  status = CreateCoordServer_dbl(&Coords_wlist, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_clock, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_mytim, &MaxNum, "test");
  MaxNum = 2*beads*mxlist*2;
  status = CreateCoordServer_int(&Coords_from, &MaxNum, "test");

  MaxNum = beads*nop;
  status = CreateCoordServer_int(&Coords_regind, &MaxNum, "test");

  MaxNum = 1;
  status = CreateCoordServer_int(&Coords_nreg, &MaxNum, "test");

  MaxNum = 5*2*nlockreg;
  status = GA_Create_mutexes(MaxNum);

  /*  lock indices:   one for each region, and each of the 5 variables ->

0: nlist/full
1/2: twlist/wlist
3: stravg / wrxyz
4: regind   */

//  printf("all coord servers created by proc %d\n",me); fflush(stdout);
}

void initvar() {

  int i, j, k, ind, pt, op, mcount, dir1, bead1, dir2, bead2, dir, bead;

  /* initialize arrays */

  basin[0] = bmin;
  basin[1] = bmax;

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
  }
      
  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      navg[j][k] = 0;
      for (pt=0; pt<=tres-1; pt++) {
	thist[j][k][pt] = 0;
      }
    }
  }
  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"Done initialization\n");
    fclose(mainout);
  }
}

void getregdef() {

  int dir, i, j, mcount, reg;

  /* get region definitions from server */

  for (dir=0;dir<2;dir++) {
    GetFromServer_int(&Coords_nreg,dir,0,tempi);
    GetFromServer_int(&Coords_regind,dir,-1,tregind);

    nreg[dir] = tempi[0];

    mcount = 0; 

    for (reg=0;reg<tempi[0];reg++) {
      for (i=0;i<nop;i++) {
	myreg[dir][reg][i] = tregind[mcount];
	mcount++;
      }
    }
  }
}
void wrxyz() {
  
  /* This function writes a point to the xyz function open in 'xyzout' */

  int j, k;
  lockind = getlockind(3,0,0);
  GA_Lock(lockind);
  GA_Init_fence();

  xyzout = fopen(filename,"a");

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
  fprintf(xyzout,"\n %d \n",npart);
  fclose(xyzout);
  GA_Fence();
  GA_Unlock(lockind);
}
double mypow2(double x) {
    return(x*x);
}
double mypow4(double x) {
  double x2;
  x2 = x*x;
    return(x2*x2);
}
double mypow6(double x) {
  double x2, x4;

    x2 = x*x;
    x4 = x2*x2;
    return(x4*x2);
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
int getlockind(int var, int dir, int bead) {
  
  int tbead;
  tbead = bead*nlockreg/beads;
  if (tbead > nlockreg) {
    tbead = nlockreg;
  } else if (tbead < 0) {
    tbead = 0;
  }
  
  return var*2*nlockreg + dir*nlockreg + tbead;
}

int chkdir(int dir, int bead) {
  
  /* returns 1 if there is a basin crossing, 0 otherwise */
  int temp, tbas, odir, n1, n2, newreg, ind, ind1, ind2;
  double t2;

  temp = 0;
  tbas = bascheck();
		    
  if (tbas + dir == 1) {  /* either tbas = 0, dir = 1 or tbas = 1, dir = 0 */
    if (dir == 1) {
      wrxyz();
    }

    ind = dir*beads + bead;
    temp = 1;
			
    /* basin crossing */
			
    allind[0] = 0;
    allind[1] = dir*beads + bead;
    tempi[0] = 1 + NGA_Read_inc(Coords_narr.ga,allind,1);
			
    if (dir == 1) {
      odir = 0;
    }else{
      odir = 1;
    }
   
    newreg = whichx(odir,coor);
    if (verb >= 3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Ended traj: dir cross into %d\n",newreg);
      fclose(mainout);
    }
    if (wtalg == 1) {
      fixwt(dir,odir,bead,newreg); /* transfer weight */
    }

    ind1 = dir*beads + bead;
    ind2 = odir*beads + newreg;
    allind[0] = 0;
    allind[1] = ind2*2*beads + ind1;
    
    tempi[0] = 1 + NGA_Read_inc(Coords_nflux.ga,allind,1);  
    allind[0] = 1;
    tempi[0] = 1 + NGA_Read_inc(Coords_nflux.ga,allind,1);  
 
    mytim = 0;
    addpoint(odir,newreg,coor,myweight,dir,bead);  /* send point to other flux list */
    endtraj = 1;
  }

  return temp;
}

void chkprim(int dir, int bead) {

  int newreg, ind, n1, n2, ind1, ind2;
  double t2;
  
  if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
    mainout = fopen(outname,"a");
    t2 = MPI_Wtime();
    fprintf(mainout,"in chkprim (%f)\n",t2);
    fclose(mainout);
  }

  /* prim crossing */
  if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
    mainout = fopen(outname,"a");
    t2 = MPI_Wtime();
    fprintf(mainout,"(%f) which",t2);
    fclose(mainout);
  }
  newreg = which(dir,bead);
  if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
    mainout = fopen(outname,"a");
    t2 = MPI_Wtime();
    fprintf(mainout," (%f)\n",t2);
    fclose(mainout);
  }
  if (newreg != bead) {
    ind = dir*beads + bead;
    
    if (verb >=3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Ended traj: prim cross into %d\n",newreg);
      fclose(mainout);
    }

    if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout,"(%f) addpoint",t2);
      fclose(mainout);
    }
    addpoint(dir,newreg,coor,myweight,dir,bead);  /* send point to other flux list */
    if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
      mainout = fopen(outname,"a");
      t2 = MPI_Wtime();
      fprintf(mainout," (%f)\n",t2);
      fclose(mainout);
    }

    if (wtalg == 1) {
      fixwt(dir,dir,bead,newreg); /* fix weight */
    }

    ind1 = dir*beads + bead;
    ind2 = dir*beads + newreg;
    allind[0] = 0;
    allind[1] = ind2*2*beads + ind1;
    
    tempi[0] = 1 + NGA_Read_inc(Coords_nflux.ga,allind,1);  
    allind[0] = 1;
    tempi[0] = 1 + NGA_Read_inc(Coords_nflux.ga,allind,1);  

    endtraj = 1;
    if (verb >= 3) {
      mainout = fopen(outname,"a");
      ind = dir*beads + bead;
      fprintf(mainout,"Ended traj: prim cross\n");
      fclose(mainout);
    }
  }

  if ((verb >= 4) || ((verb >= 3) && (me ==0))) {
    mainout = fopen(outname,"a");
    t2 = MPI_Wtime();
    fprintf(mainout,"out chkprim (%f)\n",t2);
    fclose(mainout);
  }
}
void wrtflux(int find, int cycle) {

  int j,k,pt,part,dim,mcount,ind,limit,scale,listind;
  double extra, fac,tsum[2];

  /* write fluxes, divide flux weights */
  
  scale = 0;
  if ((cycle%wtupdt == 0) && (wtalg == 3)) {
    scale = 1;
    GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
    mcount = 0;
    for (j=0;j<2;j++) {
      for (k=0;k<beads;k++) {
	weight[j][k] = tweight1[mcount];
	mcount++;
      }
    }
  }

  sprintf(fname,"%s%d%s%d%s","tflux", tmpswitch, "_", find-1, ".dat");
  outFile = fopen(fname,"w");
  for (j=0;j<2;j++) {
    for (k=find-1; k<beads; k=k+nfxfile) {
      
      fprintf(outFile,"%d %d\n",j,k);
      ind = j*beads + k;
      
      GetFromServer_dbl(&Coords_twlist,0,ind,tempw);
      tsum[0] = tempw[0];
      GetFromServer_dbl(&Coords_twlist,1,ind,tempw);
      tsum[1] = tempw[0];

      for (listind=0;listind<2;listind++) {

	ind = j*beads + k;
	GetFromServer_int(&Coords_full,listind,ind,tempi);
	if (tempi[0] == 1) {
	  limit = mxlist;
	} else {
	  //LINE_STAMP;
	  GetFromServer_int(&Coords_nlist,listind,ind,tempi);
	  limit = tempi[0];
	}
	
	//LINE_STAMP;
	if (scale) {
	  printf("(%d,%d,%d) weight(%e) tsum(%e,%e)\n",j,k,listind,weight[j][k],tsum[0],tsum[1]);
	  fac = weight[j][k]/(tsum[0]+tsum[1]);
	  tempw[0] = tsum[listind]*fac;
	  PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	  fprintf(outFile,"%d %e\n", limit, tempw[0]);
	} else {
	  fprintf(outFile,"%d %e\n", limit, tsum[listind]);
	}
      
	for (pt=0; pt<limit; pt++) {
	  ind = 2*mxlist*beads*j + 2*mxlist*k + 2*pt;
	  GetFromServer_int(&Coords_from,listind,ind,tempi);
	  GetFromServer_int(&Coords_from,listind,ind+1,tempi2);
	  
	  fprintf(outFile,"%d %d ", tempi[0], tempi2[0]);
	  
	  ind = mxlist*beads*j + mxlist*k + pt;
	  
	  GetFromServer_dbl(&Coords_wlist,listind,ind,tempw);
	  if (scale) {
	    tempw[0] = tempw[0]*fac;
	    PutToServer_dbl(&Coords_wlist,listind,ind,tempw);
	  }
	  
	  GetFromServer_int(&Coords_clock,listind,ind,tempi);
	  fprintf(outFile,"%e %d ",tempw[0],tempi[0]);
	  GetFromServer_int(&Coords_mytim,listind,ind,tempi);
	  fprintf(outFile,"%d ",tempi[0]);
	  
	  lo[0] = listind;
	  hi[0] = listind;
	  lo[1] = 3*npart*ndim*mxlist*beads*j + 3*npart*ndim*mxlist*k + 3*npart*ndim*pt;
	  hi[1] = 3*npart*ndim*mxlist*beads*j + 3*npart*ndim*mxlist*k + 3*npart*ndim*(pt+1) -1;
	  NGA_Get(Coords_pts.ga,lo,hi,tpt,ld);
	  
	  mcount = 0;
	  for (part=0; part<npart;part++) {
	    for (dim=0; dim<=ndim-1; dim++) {
	      fprintf(outFile,"%e %e %e ", tpt[mcount],tpt[mcount+1],tpt[mcount+2]);
	      mcount += 3;
	    }
	  }
	  fprintf(outFile,"\n");
	}
      }
    }
  }
  fclose(outFile);
}
void getop(int i, int j, int *op) {

  int ti, tj;

  if (i < j) {  // switch
    tj = i;
    ti = j;
  } else {
    ti = i;
    tj = j;
  } // now tj is always less than ti

  if (tj <= 20) {
    if (ti <= 30) {
      *op = 0;
    } else if (ti <= 217) {
      *op = 1;
    } else {
      *op = 2;
    }
  } else if (tj <= 57) {
    if (ti <= 137) {
      *op = 3;
    } else {
      *op = 4;
    }
  } else if (tj <= 126) {
    if (ti <= 182) {
      *op = 5;
    } else {
      *op = 6;
    }
  } else if (tj <= 186) {
    if (ti <= 212) {
      *op = 7;
    } else {
      *op = 8;
    }
  } else {
    *op = 9;
  }

}
int findreg(iopoint_t iox, int dir) {

  int reg, i, j, good, gotit, mcount, overwt, target, ind;

  if (verb>=4) {
    printf("Finding reg for point:\n");
    for (i=0;i<nop;i++) {
      printf("%d ",iox.x[i]);
    }    
    printf("\n");
  }

  gotit = 0;
  if (nreg[dir] > 0) {
    for (reg=0;reg<nreg[dir];reg++) {
      good = 1;
      for (i=0;i<nop;i++) {
	if (iox.x[i] != myreg[dir][reg][i]) {
	  good = 0;
	}
      }
      if (verb>=4) {
	printf("checking point %d:\n",reg);
	for (i=0;i<nop;i++) {
	  printf("%d ",myreg[dir][reg][i]);
	}    
	printf("\n");
      }
	
      if (good) {
	gotit = reg+1;
	break;
      }
    }
  }

  if (!gotit) {
    GetFromServer_int(&Coords_nreg,dir,0,tempi);
    GetFromServer_int(&Coords_regind,dir,-1,tregind);
    
    nreg[dir] = tempi[0];
    
    mcount = 0;  //reload myreg, check again
    gotit = 0;

    if (tempi[0] > 0) {
      for (reg=0;reg<tempi[0];reg++) {
	good = 1;
	for (i=0;i<nop;i++) {
	  myreg[dir][reg][i] = tregind[mcount];
	  if (iox.x[i] != myreg[dir][reg][i]) {
	    good = 0;
	  }
	  mcount++;
	}
	if (good) {
	  gotit = reg+1;
	}
      }
    }
  }

  if (!gotit) {
    // lock it up
    lockind = getlockind(4,0,dir);
    GA_Lock(lockind);
    GA_Init_fence();
    
    GetFromServer_int(&Coords_nreg,dir,0,tempi);
    GetFromServer_int(&Coords_regind,dir,-1,tregind);
    
    nreg[dir] = tempi[0];
    
    mcount = 0;  //reload myreg, check again
    gotit = 0;

    if (tempi[0] > 0) {
      for (reg=0;reg<tempi[0];reg++) {
	good = 1;
	for (i=0;i<nop;i++) {
	  myreg[dir][reg][i] = tregind[mcount];
	  if (iox.x[i] != myreg[dir][reg][i]) {
	    good = 0;
	  }
	  mcount++;
	}
	if (good) {
	  gotit = reg+1;
	}
      }
    }
    
    if (!gotit) {
      // add the region
      tempi[0]++;
      if (tempi[0] > beads) {
	printf("Exceeded region quota!\n");
	writereg();
	gaexit(77);
      } else {
	printf("%d: Adding region %d %d: ",me,dir,tempi[0]-1);
	nreg[dir] = tempi[0];
	gotit = tempi[0];
	for (i=0;i<nop;i++) {
	  tregind[mcount] = iox.x[i];
	  myreg[dir][gotit-1][i] = iox.x[i];
	  printf("%d ",iox.x[i]);
	  mcount++;
	}
	printf("\n");
	PutToServer_int(&Coords_regind,dir,-1,tregind);	
	PutToServer_int(&Coords_nreg,dir,0,tempi);
	ind = dir*beads + gotit-1;
	tempi[0] = 1;
	PutToServer_int(&Coords_on,0,ind,tempi);
      }
    }
    GA_Fence();
    GA_Unlock(lockind);
  }

  gotit--;

  if (verb>=4) {
    printf("returning %d\n",gotit);
  }

  return gotit;
}
void printiox(point_t x, int dir) {

  opoint_t ox;
  iopoint_t iox;
  int i, j;

  ox = projx(x);

  for (i=0;i<nop;i++) {
    iox.x[i] = (int)((cryst.x[i]-ox.x[i])/opstep[dir]);
    printf("%d ",iox.x[i]);
  }
  printf("\n");

}
void scalflux() {

  int j,k,pt,part,dim,mcount,ind,limit,scale,listind;
  double extra, fac,tsum[2];

  /* write fluxes, divide flux weights */
  
  GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
  mcount = 0;
  for (j=0;j<2;j++) {
    for (k=0;k<beads;k++) {
      weight[j][k] = tweight1[mcount];
      mcount++;
    }
  }

  for (j=0; j<2; j++) {
    for (k=0; k<=beads-1; k++) {
      ind = j*beads + k;

      GetFromServer_dbl(&Coords_twlist,0,ind,tempw);
      tsum[0] = tempw[0];
      GetFromServer_dbl(&Coords_twlist,1,ind,tempw);
      tsum[1] = tempw[0];

      for (listind=0;listind<2;listind++) {

	ind = j*beads + k;
	GetFromServer_int(&Coords_full,listind,ind,tempi);
	if (tempi[0] == 1) {
	  limit = mxlist;
	} else {
	  //LINE_STAMP;
	  GetFromServer_int(&Coords_nlist,listind,ind,tempi);
	  limit = tempi[0];
	}
      
	fac = weight[j][k]/(tsum[0]+tsum[1]);
	tempw[0] = tsum[listind]*fac;
	PutToServer_dbl(&Coords_twlist,listind,ind,tempw);

	for (pt=0; pt<limit; pt++) {
	  ind = mxlist*beads*j + mxlist*k + pt;
	  GetFromServer_dbl(&Coords_wlist,listind,ind,tempw);
	  tempw[0] = tempw[0]*fac;
	  PutToServer_dbl(&Coords_wlist,listind,ind,tempw);
	}
	
      }
    }
  }
}
void writereg() {
  
  int dir, reg, i, j, sum, mxcon, mncon, con, mnind, mcount;
  FILE *fout;

  GetFromServer_dbl(&Coords_twlist,0,-1,tweight1);
  GetFromServer_dbl(&Coords_twlist,1,-1,twavg1);
  GetFromServer_int(&Coords_on,0,-1,ton1);

  mcount = 0;
  for (dir=0;dir<=1;dir++) {
    for (reg=0;reg<beads;reg++) {
      tweight1[mcount] += twavg1[mcount];
      mcount++;
    }
  }
  
  getregdef();
  fout = fopen("tregs.dat","w");

  for (dir=0;dir<2;dir++) {
    fprintf(fout,"dir = %d\n",dir);
    if (nreg[dir] > 0) {

      // find mxcon, mncon

      for (reg=0;reg<nreg[dir];reg++) {
	sum = 0;
	for (i=0;i<nop;i++) {
	  sum += myreg[dir][reg][i];
	}
	if (reg == 0) {
	  mxcon = sum;
	  mncon = sum;
	  mnind = reg;
	} 
	if (sum > mxcon) {
	  mxcon = sum;
	}
	if (sum < mncon) {
	  mncon = sum;
	  mnind = reg;
	}
      }

      for (con=mxcon;con>=mncon;con--) {

	mcount = dir*beads;
	for (reg=0;reg<nreg[dir];reg++) {
	  sum = 0;
	  for (i=0;i<nop;i++) {
	    sum += myreg[dir][reg][i];
	  }
	
	  if (sum == con) {
	    fprintf(fout,"%3d %3d  ",reg,con);
	    for (i=0;i<nop;i++) {
	      fprintf(fout,"%d ",myreg[dir][reg][i]-myreg[dir][mnind][i]);
	    }
	    if (ton1[mcount]) {
	      fprintf(fout,"%5d ",(int)(log10(tweight1[mcount])));
	    } else {
	      fprintf(fout,"%s ","XXXXX");
	    }
	    fprintf(fout,"\n");
	  }
	  mcount++;
	}
      }
    }
  }
  fclose(fout);
}
int getlistind(int dir, int bead, int from1, int from2) {

  int sum1, sum2, i, temp;

  sum1 = 0;
  for (i=0;i<nop;i++) {
    sum1 += myreg[dir][bead][i];
  }
  sum2 = 0;
  for (i=0;i<nop;i++) {
    sum2 += myreg[from1][from2][i];
  }
  if (sum1 <= sum2) {
    temp = 0;
  } else {
    temp = 1;
  }

  return temp;

}
