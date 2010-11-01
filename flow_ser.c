#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> /* required for randomize() and random() */
#include <time.h>
#include <unistd.h>
#include "../nrutil.h"
#include "../nr.h"    /* numerical recipes: for GJ elim */
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

ServedCoords Coords_prob;
ServedCoords Coords_nprob;
ServedCoords Coords_stravg;
ServedCoords Coords_z;
ServedCoords Coords_narr;
ServedCoords Coords_navg;
ServedCoords Coords_thist;
ServedCoords Coords_nlist;
ServedCoords Coords_pts;
ServedCoords Coords_wlist;
ServedCoords Coords_rcount;
ServedCoords Coords_reg;
ServedCoords Coords_extrawt;
ServedCoords Coords_nin;
ServedCoords Coords_avwalk;
ServedCoords Coords_nout;
ServedCoords Coords_weight;

opoint_t **stravg;
opoint_t **z, *ocoor;
point_t *coor;
int **navg;
int allind[2];
int stredit=0;  /* stredit = 0 (hasn't been done yet) 
	         stredit = 1 (it's done already)
		 stredit = -1 (it's in progress) */

/* global variables */

int nsec, ntert, nav, tnlist;
double *tthist1, *tthist2, alf=1., bet=1.;
double *tstravg1, *tstravg2, *tz1, *tz2, *tpt, *tw1, *tw2;
double **wavg;
double globtim=0.;
int nglobtim=0;
int *ti1, *ti2, *tr1, *tr2, *treg1, *treg2, *tin1, *tin2, *tprob1, *tprob2, *tnprob1, *tnprob2;
int **secind, **tertind, ****nin, **nout, **avw;
int verb=2;
int me, nproc, status, beadsp2, beadsp3, lockind, lo[2], hi[2], ld[1], *myrcount;
double tstepsq, tempw[1], wcasig, wcacut, *myweight;
double sigfacsq=4.;

char crname[70]="../cutRNAin.dat";      // crystal structure file
char secname[70]="../cutRNAsec.con";    // secondary contacts
char tertname[70]="../cutRNAtert.con";  // tertiary contacts
char emname[30]="../emerg05.dat";         // emergency points
char ptname[30]="allpt.dat";         // points
char iwtname[30]="iweight";         // weights
char sname[30]="isolv";         // emergency points
char filename[30];         // emergency points
double **r0, **r02, **r0d6, **r0d12, **nextstep, **r0dists, **rdists, ***vsolv, ***rsolv, ***f1;
int **pos, **pos_poly;
int *done, alldone;
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
int N=503200, maxd, strmfrq, cycle;
int Lx=384, Ly=192, Lz=384;
double blen=8.;              /* length of box in each dir */
double kcalconv = 6695.5;
double cutoff=3.;
int nx,ny,nz;
double blinv, tw[1];
double kf, R0, ooR02, eph1_12, eph2_12, ep1_6, sigma1, sigma16, sigma2, sigma26;
void wrxyz(int step, int jj);
int globn;


const int ndim=NDIM_c, nop=NOP_c, npart=NPART_c;
int onpart=262;
int state, bas, mytim=0, sofar, polyerr=0;
int per=1, *dir, *bead;
int myclock=0, tmpswitch;

int getprob(int dir, int bead, int i, int j, int n1, int n2);
int getnprob(int dir, int bead, int i, int j);
opoint_t projx(point_t x);
int bascheck(int jj);
int which(int jj);
int whichx(int xdir, point_t x);
void getweight (int jj);
int updtwt(int ind);
void gaexit(int sig);
double ogtdist(opoint_t px, opoint_t py);
void getstring();
void chkprim(int jj);
int chkdir(int jj);
void addpoint (int jj);
int **rstart, tempi[1], tempi2[1], endtraj, mxlist;

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
char wtname[30];
char ioname[30];
char ptoname[30];
FILE * outPath, *outFile, *rfile;

/*-------------------------------------------------------------------------*/

int main(int argc,char** argv){

  int i, j, k, l, m, dim, size, temp, worker, mtag, n1, n2, part, lim, nreg;
  int op, dir1, dir2, to, rem, limit, pt, bead1, bead2, odir, obead, wtsucc=0;
  int tcount, ti[1], tlat, tbead, idir, ibead;
  int rc, maxsize, ind, flag, mcount, it1, it2;
  int endcycle, gotmin;
  float tf1, tf2, tf3;
  unsigned int iseed = (unsigned int)time(NULL);
  void createservers(), initvar();
  int strinit();
  void move();
  void stack(), wrpath(int ind);
  void strread();
  void wrextend(int index);
  void rddelta(), rdcryst(), rdsolv(), rdemerg();
  double rate[2], sum1, sum2, t1, t2, telaps, alpha, *tgt, *rrnd;
  FILE * ratFile, *filein, *ffile, *ptout;
    
  ld[0] = 0;
  tmpswitch = 0;

  if (seed ==0) {
      seed = iseed;
  }

  tstepsq=tstep*tstep;
  beadsp2 = beads*beads;
  beadsp3 = beadsp2*beads;

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
  //printf("proc %d: MPI_Init done\n",me); fflush(stdout);
  
  GA_Initialize();
  //printf("proc %d: GA_Initialize done\n",me); fflush(stdout);
  
  MA_init(MT_DBL, 128*1024*1024, 16*1024*1024);
  //printf("proc %d: MA_init done\n",me); fflush(stdout);

  mxlist = per*nproc/2;

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

  sprintf(ptoname,"%s%d%s","pt",me,".dat");
  ptout = fopen(ptoname,"w");

  if ((verb >= 3) || ((verb >=2) && (me == 0))) {
      sprintf(outname,"%s%d%s","out",me,".dat");
      mainout = fopen(outname,"w");
      fprintf(mainout,"rank = %d\n", me);
      fclose(mainout);
  }
  //  sprintf(outname,"%s","out0.dat");
  sprintf(wtname,"%s","mdweight.dat");
  sprintf(ratname,"%s","mdrate.dat");

  if (me == 0) {
      ratFile = fopen(wtname,"w");
      fclose(ratFile);
      
      ratFile = fopen(ratname,"w");
      fclose(ratFile);

  }

  seed *= me+1;
  for (i=0;i<100;i++){
    sum1 = ran2(&seed);
  }
  
  lim = 2*beads;
  tstravg1 = (double *) malloc(beads*nop*sizeof(double));
  tstravg2 = (double *) malloc(beads*nop*sizeof(double));
  tz1 = (double *) malloc(beads*nop*sizeof(double));
  tz2 = (double *) malloc(beads*nop*sizeof(double));
  tpt = (double *) malloc(3*ndim*npart*sizeof(double));
  tthist1 = (double *) malloc(tres*sizeof(double));
  tthist2 = (double *) malloc(tres*sizeof(double));

  treg1 = (int *) malloc(per*nproc*sizeof(int));
  treg2 = (int *) malloc(per*nproc*sizeof(int));
  tr1 = (int *) malloc(beads*sizeof(int));
  tr2 = (int *) malloc(beads*sizeof(int));
  ti1 = (int *) malloc(beads*sizeof(int));
  ti2 = (int *) malloc(beads*sizeof(int));
  tw1 = (double *) malloc(beads*sizeof(double));
  tw2 = (double *) malloc(beads*sizeof(double));
  tin1 = (int *) malloc(2*beads*beads*sizeof(int));
  tin2 = (int *) malloc(2*beads*beads*sizeof(int));
  wavg = (double **) alloc2d(sizeof(double), 2, beads);
  nin = (int ****) alloc4d(sizeof(int), 2, beads, 2, beads);
  rstart = (int **) alloc2d(sizeof(int), per, 2);
  nout = (int **) alloc2d(sizeof(int), 2, beads);
  avw = (int **) alloc2d(sizeof(int), 2, beads);

  myweight = (double *) malloc(per*sizeof(double));
  done = (int *) malloc(per*sizeof(int));
  rrnd = (double *) malloc(per*sizeof(double));
  tgt = (double *) malloc(per*sizeof(double));
  myrcount = (int *) malloc(per*sizeof(int));
  dir = (int *) malloc(per*sizeof(int));
  bead = (int *) malloc(per*sizeof(int));
  coor = (point_t *) malloc(per*sizeof(point_t));
  ocoor = (opoint_t *) malloc(per*sizeof(opoint_t));
  z = (opoint_t **) alloc2d(sizeof(opoint_t), 2, beads);
  stravg= (opoint_t **) alloc2d(sizeof(opoint_t), 2, beads);
  navg = (int **) alloc2d(sizeof(int), 2, beads);
  lim = 2*beads*2*beads*2*beads;
  tprob1 = (int *) malloc(lim*sizeof(int));
  lim = 2*beads*2*beads;
  tnprob1 = (int *) malloc(lim*sizeof(int));
  lim = 2*beads;
  
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

  rsolv = (double ***) alloc3d(sizeof(double), per, N, 3);
  vsolv = (double ***) alloc3d(sizeof(double), per, N, 3);

  v_cmax = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmay = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmaz = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmx = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmy = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  v_cmz = (double ***) alloc3d(sizeof(double), nx, ny, nz);
  pos = (int **) alloc2d(sizeof(int), N, 3);
  pos_poly = (int **) alloc2d(sizeof(int), npart, 3);

  f1 = (double ***) alloc3d(sizeof(double), per, npart, 3);

  rdcryst();
  rddelta();
  rdsolv();

  initvar();
  if (verb >= 4) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"after initvar\n");
    fclose(mainout);
  }
  if (strcmp(ptname,"NOREAD") == 0) {
    printf("Damn!\n");
    gaexit(2);
  } else {

    mcount = 0;
    mainout = fopen(ptname,"r");
    lim = 16;
    if (lim < nproc*per) {
      printf("not enough points!\n");
      gaexit(2);
    } else {
      for (k=1;k<=lim;k++) {		
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"reading k = %d",k);
	  fclose(mainout);
	}
	if ((mcount >= me*per) && (mcount < (me+1)*per)) {
	  ind = mcount - me*per;
	  if ((ind < 0) || (ind >= per)) {
	    printf("Error!  me = %d, ind = %d!",me, ind);
	    gaexit(2);
	  }
	  dir[ind] = 1;
	  bead[ind] = 19;
	  rstart[ind][0] = 0;
	  rstart[ind][1] = 19;
	  /*	  fscanf(mainout,"%d %d ", &tempi[0], &tempi2[0]);
	  fscanf(mainout,"%f %d ",&tw,&tempi[0]);
	  fscanf(mainout,"%d ",&tempi[0]);*/

	  for (part=0; part<npart;part++) {
	    for (dim=0; dim<=ndim-1; dim++) {
	      fscanf(mainout,"%e %e %e ",&tf1,&tf2,&tf3);
	      coor[ind].x[part][dim] = tf1;
	      coor[ind].v[part][dim] = tf2;
	      f1[ind][part][dim] = tf3;
	    }
	  }
	} else {
	  /*	  fscanf(mainout,"%d %d ", &tempi[0], &tempi2[0]);
	  fscanf(mainout,"%f %d ",&tw,&tempi[0]);
	  fscanf(mainout,"%d ",&tempi[0]);*/

	  for (part=0; part<npart;part++) {
	    for (dim=0; dim<=ndim-1; dim++) {
	      fscanf(mainout,"%e %e %e ",&tf1,&tf2,&tf3);
	    }
	  }
	}
	mcount++;
      }
    /*    for (ind=0;ind<per;ind++) {
      dir[ind] = 1;
      bead[ind] = 19;
      rstart[ind][0] = 0;
      rstart[ind][1] = 19;
    }
    for (i=0;i<npart;i++) {
    //    fscanf(mainout,"%s",&temps);
      for (j=0;j<3;j++) {
	fscanf(mainout,"%e ",&tf1);
	for (ind=0;ind<per;ind++) {
	  coor[ind].x[i][j] = tf1;
	  }*/

	/*	fscanf(mainout,"%e ",&tf1);
	for (ind=0;ind<per;ind++) {
	  coor[ind].v[i][j] = tf1;
	}

	fscanf(mainout,"%e ",&tf1);
	for (ind=0;ind<per;ind++) {
	  f1[ind][i][j] = tf1;
	  }*/
      fclose(mainout);
      printf("proc %d done reading\n",me);
    }
  }


  if (verb >=4) {
    mainout = fopen(outname,"a");
    for (ind=0;ind<per;ind++) {
      fprintf(mainout,"walker %d:\n",ind);
      fprintf(mainout,"positions:\n");
      for (part=0; part<npart;part++) {
	for (dim=0; dim<=ndim-1; dim++) {
	  fprintf(mainout,"%e ",coor[ind].x[part][dim]);
	}
      }
      fprintf(mainout,"\n");
      fprintf(mainout,"velocities:\n");
      for (part=0; part<npart;part++) {
	for (dim=0; dim<=ndim-1; dim++) {
	  fprintf(mainout,"%e ",coor[ind].v[part][dim]);
	}
      }
      fprintf(mainout,"\n");
      fprintf(mainout,"forces:\n");
      for (part=0; part<npart;part++) {
	for (dim=0; dim<=ndim-1; dim++) {
	  fprintf(mainout,"%e ",f1[ind][part][dim]);
	}
      }
      fprintf(mainout,"\n");
    }
    fclose(mainout);
  }
  
  if (me == 0) {
    GA_Init_fence();
    printf("proc %d is initializing the string\n",me); fflush(stdout);
    
    /* initialize string */

    if (strcmp(flname,"NOREAD") == 0) {
	if (verb >= 2) {
	    printf("Automatically initializing string\n"); fflush(stdout);
	}
	if (strinit())
	    gaexit(2) ;  
    } else {
	if (verb >= 2) {
	    printf("Using string configuration from %s\n", flname); fflush(stdout);
	}
	strread();
    }
    wrpath(0);
    
    /* send string to server */

    mcount = 0;
    for (i=0;i<beads;i++) {
      for (j=0;j<nop;j++) {
	tz1[mcount] = z[0][i].x[j];
	tz2[mcount] = z[1][i].x[j];
	mcount++;
      }
    }
    PutToServer_dbl(&Coords_z,0,-1,tz1);
    PutToServer_dbl(&Coords_z,1,-1,tz2);
    
    printf("proc %d is done with the string\n",me); fflush(stdout);
    GA_Fence();
  }

  /* check point */

  GA_Sync();

  for (j=0;j<per;j++) {
    sprintf(filename,"%s%d%s%d%s","flow",me,"_",j,".dat");
  
    xyzout = fopen(filename,"w");
    fprintf(xyzout,"%i \n polymer movie\n", npart);
    fclose(xyzout);
  }

  if (me == 0) {
    t1 = MPI_Wtime();
  }
  getstring();
  for (i=0;i<per;i++) {
    done[i] = 0;
  }
  alldone = 0;

  /* ---------------start of dynamics loop-------------------------- */
  cycle = 0;
  while (!alldone) {
    cycle++;
  
    /*    if (me == 0) {
      printf("Starting cycle %d of %d\n",cycle,T);
    
      }*/

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

      /* stack */
      for (j=0;j<per;j++) {
	if (!done[j]) {
	  if (cycle%xyzfrq ==0) {
	    if (verb >= 3) {
	      mainout = fopen(outname,"a");
	      fprintf(mainout,"proc %d: tcount = %d writing to xyz \n",me, tcount);
	      fclose(mainout);
	    }
	    sprintf(filename,"%s%d%s%d%s","flow",me,"_",j,".dat");
	    wrxyz(cycle,j);
	    /*	    ptout = fopen(ptoname,"a");
	    for (part=0;part<npart;part++) {
	      for (dim=0;dim<ndim;dim++) {
		fprintf(ptout,"%e ",coor[j].x[part][dim]);
		fprintf(ptout,"%e ",coor[j].v[part][dim]);
		fprintf(ptout,"%e ",f1[j][part][dim]);
	      }
	    }
	    fprintf(ptout,"\n");
	    fclose(ptout);*/
	  }
	}
      }
		
      if (cycle%stkfrq == 0) {
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"stacking...\n");
	  fclose(mainout);
	}
	stack();
      }
      
    } else {   /* if polyerr */
      if (verb >= 3) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"polyerr!\n");
	fclose(mainout);
      }
      endtraj = 1;
      polyerr = 0;
      printf("polyerr\n");
      gaexit(2);
    }
         /* end of cycle */

      
    /* compute the rate / update weights / write thist */

    for (j=0;j<per;j++) {
      if (!done[j]) {
	odir = dir[j];
	obead = bead[j];
	
	if (!chkdir(j)) {
	  chkprim(j);
	}
      }
      
    }

    if (me == 0) {
      if (cycle%wrfrq == 0) {
	wrextend(cycle/wrfrq);

	GetFromServer_int(&Coords_weight,0,-1,ti1);
	GetFromServer_int(&Coords_weight,1,-1,ti2);
	mainout = fopen("mdwts.dat","w");    
	for (j=0;j<beads;j++) {
	  fprintf(mainout,"%d %d\n",ti1[j],ti2[j]);
	}
	fclose(mainout);
	
	GetFromServer_int(&Coords_nprob,0,-1,tnprob1);
	GetFromServer_int(&Coords_prob,0,-1,tprob1);
	
	if (tmpswitch) {
	  tmpswitch = 0;
	} else { 
	  tmpswitch = 1;
	}

	sprintf(fname,"%s%d%s","tprob", tmpswitch, ".dat");
	outFile = fopen(fname,"w");
	for (idir=0; idir<=1; idir++) {
	  for (ibead=0; ibead<=beads-1; ibead++) {
	    for (n1=0; n1<2; n1++) {
	      for (n2=0; n2<=beads-1; n2++) {
		it1 = getnprob(idir,ibead,n1,n2);
		if (it1 != 0) {
		  fprintf(outFile,"%d %d %d %d %d\n",idir,ibead,n1,n2,it1);
		  for (i=0; i<2; i++) {
		    for (j=0; j<=beads-1; j++) {
		      it2 = getprob(idir,ibead,n1,n2,i,j);
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
	fclose(outFile);

      }
      /* write weights */

    }
    alldone = 1;
    for (i=0;i<per;i++) {
      if (!done[i]) {
	alldone = 0;
      }
    }
  }
  
  /* ---------------end of dynamics loop---------------------------- */
  
  /* set the restart point: write fluxes, weights, etc. */
  
  if (me == 0) {
    t2 = MPI_Wtime();
    printf("total time of dynamics steps: %f\n",t2-t1); fflush(stdout);
    printf("total time of locking steps: %f\n",globtim); fflush(stdout);
    printf("total number of locking steps: %d\n",nglobtim); fflush(stdout);
    printf("time per locking step: %f\n",globtim/(float)nglobtim); fflush(stdout);

    tcount = per*nproc*T*chkfrq;
    printf("total number of dynamics steps: %d\n",tcount); fflush(stdout);
    printf("time per dynamics step: %f\n",(t2-t1)/(float)tcount); fflush(stdout);
    
  }
  
  /* check point */
//  printf("proc %d: final check\n",me); fflush(stdout);
  GA_Sync();
  if (me ==0) {
    if (cycle%wrfrq == 0) {
      wrextend(0);
    }

    lo[0] = 0;
    lo[1] = 0;
    hi[0] = 0;
    hi[1] = 0;
    NGA_Get(Coords_narr.ga,lo,hi,tw,ld);
    mainout = fopen("mdrate.dat","w");
    fprintf(mainout,"%e \n",(float)nproc/tw[0]);
    fclose(mainout);
 
    GetFromServer_int(&Coords_weight,0,-1,ti1);
    GetFromServer_int(&Coords_weight,1,-1,ti2);
    mainout = fopen("mdwts.dat","w");    
    for (j=0;j<beads;j++) {
      fprintf(mainout,"%d \n",ti1[j]);
    }
    fclose(mainout);
  }

  GA_Sync();

  status = DestroyCoordServer(&Coords_narr);
  status = DestroyCoordServer(&Coords_navg);
  status = DestroyCoordServer(&Coords_prob);
  status = DestroyCoordServer(&Coords_nprob);
  status = DestroyCoordServer(&Coords_stravg);
  status = DestroyCoordServer(&Coords_z);

  status = DestroyCoordServer(&Coords_thist);
  status = DestroyCoordServer(&Coords_nlist);
  status = DestroyCoordServer(&Coords_pts);
  status = DestroyCoordServer(&Coords_wlist);
  status = DestroyCoordServer(&Coords_weight);
  status = DestroyCoordServer(&Coords_nin);
  status = DestroyCoordServer(&Coords_avwalk);
  status = DestroyCoordServer(&Coords_nout);

  if (me == 0) GA_Print_stats();
  
  GA_Terminate();
  MPI_Finalize();
  
  return(0);
}




/*=========================================================
  END OF MAIN
  ============================================================*/

/*-------------------------------------------------------------------------*/
void wrextend(int index) {
    
  /* this function writes out a theta histogram */

  int ind, mcount, mcount2, tsum;
  double x[tres], sum, ext, ootsum;
  FILE * tout;
  char tmpname[30];
  
  for (ind=0;ind<=tres-1;ind++) {
    x[ind] = 0.;
  }

  GetFromServer_dbl(&Coords_thist,0,-1,tthist1);
  GetFromServer_dbl(&Coords_thist,1,-1,tthist2);

  sum = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += tthist1[ind];
  }

  sprintf(tmpname,"%s%d%s","thist",index,"_1.dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin + (tmax-tmin)*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,tthist1[ind]/sum);
  }
  fclose(tout);

  sum = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += tthist2[ind];
  }

  sprintf(tmpname,"%s%d%s","thist",index,"_2.dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin + (tmax-tmin)*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,tthist2[ind]/sum);
  }
  fclose(tout);

  sum = 0.;
  for (ind=0;ind<=tres-1;ind++) {
    sum += tthist1[ind] + tthist2[ind];
  }

  sprintf(tmpname,"%s%d%s","thist",index,"_tot.dat");
  tout = fopen(tmpname,"w");
  for (ind=0;ind<=tres-1;ind++) {
    ext = tmin + (tmax-tmin)*(ind/(float)tres);
    fprintf(tout,"%e %e\n",ext,(tthist1[ind]+tthist2[ind])/sum);
  }
  fclose(tout);
 
  /*  for (ind=0;ind<tres;ind++) {
    tthist1[ind] = 0.;
    tthist2[ind] = 0.;
  }

  PutToServer_dbl(&Coords_thist,0,-1,tthist1);
  PutToServer_dbl(&Coords_thist,1,-1,tthist2);*/

}
/*-------------------------------------------------------------------------*/
int which(int jj) {

  /* this function returns the region that the walker
   is in, presently */
       
   int i, j, off, temp;
  double dist, lodist;

  for (i=0;i<=beads-1;i++) {
    dist = ogtdist(z[dir[jj]][i],ocoor[jj]);
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
int whichx(int xdir, point_t x) {
   
   /* this function returns the region that the point 'x'
      is in, presently */
   
   int i, j, off, temp;
   double dist, lodist;
   opoint_t ox;
   
   ox = projx(x);
   
   for (i=0;i<=beads-1;i++) {
     dist = ogtdist(z[xdir][i],ox);
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
void stack() { 
   
   /* this function adds the current state of the replicas to stravg */
   
   int i, j, index, ind;
   double temps[1], mx, mn;
   point_t rt;

   for (i=0;i<per;i++) {
     if (!done[i]) {

       allind[0] = dir[i];
       allind[1] = bead[i];
       tempi[0] = 1 + NGA_Read_inc(Coords_weight.ga,allind,1);

       if (phase==1) {
   
	 allind[0] = dir[i];
	 allind[1] = bead[i];
	 tempi[0] = 1 + NGA_Read_inc(Coords_navg.ga,allind,1);
	 
	 lockind = dir[i]*beads + bead[i] + 1;
	 GA_Lock(lockind);
	 GA_Init_fence();
	 
	 for (j=0;j<=nop-1;j++) {
	   index = bead[i]*nop + j;
	   
	   GetFromServer_dbl(&Coords_stravg,dir[i],index,temps);
	   temps[0] += ocoor[i].x[j];
	   PutToServer_dbl(&Coords_stravg,dir[i],index,temps);
	 }
	 GA_Fence();
	 GA_Unlock(lockind);
       }
       
       index = tres*(ocoor[i].x[0]-tmin)/(tmax-tmin);
       if ((index < 0) || (index >= tres)) {
	 printf("Error in stack index! %f\n",ocoor[i].x[0]);
       } else {
	 lockind = index + 2*beads + 1;
	 GA_Lock(lockind);
	 GA_Init_fence();
	 GetFromServer_dbl(&Coords_thist,dir[i],index,tw);
	 tw[0] += 1.;
	 PutToServer_dbl(&Coords_thist,dir[i],index,tw);
	 GA_Fence();
	 GA_Unlock(lockind);
       }
     }
   }
}
/*-------------------------------------------------------------------------*/
void wrpath(int ind) {

  int j, k, op;

  sprintf(fname,"%s%d%s","path", ind, ".dat");
  outPath = fopen(fname,"w");
  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      for (op=0;op<=nop-1;op++) {
	fprintf(outPath,"%e ",z[j][k].x[op]);
      }
      fprintf(outPath,"\n");
    }
  }
  fclose(outPath);
}

/*-------------------------------------------------------------------------
--------------------------------------------
the following functions are system-dependent
--------------------------------------------
-------------------------------------------------------------------------*/

int strinit() {

  /* this function initializes the string */

  int op, i;

  /* string 1: forward & backward */

  for (i=0;i<=beads-1;i++) {
    if (i<=4) {
      z[0][i].x[0] = bmin + bwidth + (bmax-bmin-2*bwidth)*(i)/(float)(beads-1);
    } else {
      z[0][i].x[0] = bmin + bwidth + (bmax-bmin-2*bwidth)*(i)/(float)(beads-1);
    }
    z[1][i].x[0] = z[0][i].x[0];
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

int bascheck(int jj) {
  
  /* this function checks if the walker is in a basin */

  int temp, i;
  double dist;

  temp = -1;

  if (ocoor[jj].x[0] < basin[0].x[0]) {
    temp = 0;
  } else if (ocoor[jj].x[0] > basin[1].x[0]) {
    temp = 1;
  } else if (bwidth > 0.) {
    for (i=0;i<=1;i++) {
      dist = ogtdist(ocoor[jj],basin[i]);
      if (dist < bwidth) {
	temp = i;
      }
    }
  }  

  return temp;
}
/*-------------------------------------------------------------------------*/
opoint_t projx(point_t x) {
  
  /* this function returns the projection of coor onto the order parameter space */
  
  int i, j, k, ind1, ind2;
  opoint_t temp;
  double mn, mx, r2, sum, tdub;
  
  /*
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

  temp.x[0] = mx-mn;*/

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
      sum = sum + tdub;
    }
  }

  temp.x[0] = sum;
  
  return temp;
}
/*-------------------------------------------------------------------------*/
void move() {
  
  /* moves coor, updates ocoor */

  int i, j;
  point_t oldx, rt;
  double psi[2],phi[2],fx,fy,temp1,temp2,dvdx,dvdy;
  void dostep(int jj), stream(int jj), dorotate(int jj);
  
  for (j=0;j<per;j++) {
    if (!done[j]) {
      if (verb >= 4) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"dostep...");
	fclose(mainout);
      }
      dostep(j);
      if (!polyerr) {
	if (myclock%strmfrq == 0) {
	  if (verb >= 4) {
	    mainout = fopen(outname,"a");
	    fprintf(mainout,"stream...");
	    fclose(mainout);
	  }
	  stream(j);
	  if (verb >= 4) {
	    mainout = fopen(outname,"a");
	    fprintf(mainout,"dorotate...");
	    fclose(mainout);
	  }
	  dorotate(j);
	}
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"project...\n");
	  fclose(mainout);
	}
	ocoor[j] = projx(coor[j]);  /* update o.p. projection */
      }
    }
  }
}
/*-------------------------------------------------------------------------*/
void dostep(int jj){
  int i, j, k, x;
  double a, f2[npart][ndim],temp,x1,x2,x3;
  point_t f, force();

#pragma omp parallel for default(shared) private(i,x,a) schedule(static)
  for(i=0; i< npart; i++){
    for(x=0; x< 3; x++){	
      /* calculate a(t) = F/m */        
      // forces are pre-divided by polymass to save calcs
      coor[jj].x[i][x] += 0.5*tstepsq*f1[jj][i][x]/polymass + tstep*coor[jj].v[i][x];
    }
  }
  
#pragma omp parallel for default(shared) private(i,j,k,temp) schedule(static)
  for(i=0; i< npart; i++){
    for(j=i; j< npart; j++){
      temp = 0.;
      for (k=0;k<3;k++) {
	temp += (coor[jj].x[i][k] - coor[jj].x[j][k])*(coor[jj].x[i][k] - coor[jj].x[j][k]);
      }
      rdists[i][j] = sqrt(temp);
      rdists[j][i] = rdists[i][j];
      
      for (k=0;k<ndim;k++) {
	rvecs[i][j].x[k] = coor[jj].x[i][k]-coor[jj].x[j][k];
	rvecs[j][i].x[k] = -rvecs[i][j].x[k];
      }
    }
  }
  
  f = force(jj);

#pragma omp parallel for default(shared) private(i,x,x1,x2,x3) schedule(static)
  for(i=0; i< npart; i++){
    for(x=0; x< 3; x++){
      f2[i][x] = f.x[i][x];
      coor[jj].v[i][x] += 0.5*(f1[jj][i][x]+f2[i][x])*tstep/polymass;
      f1[jj][i][x] = f2[i][x];
    }
    x1=coor[jj].x[i][0];
    x2=coor[jj].x[i][1];
    x3=coor[jj].x[i][2];
    if(((x1<0.)||(x1>Lx)) || (((x2<0.)||(x2>Ly)) || ((x3<0.)||(x3>Lz)))) {
      printf("proc %d: polymer bead %d out of bounds! (%f,%f,%f)\n",me,i,x1,x2,x3);
      polyerr = 1;
    }
  }
}
/*-------------------------------------------------------------------------*/
point_t force(int jj){
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
		  trvec.x[j] = coor[jj].x[p1][j] - phan[j];
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
    if (coor[jj].x[p1][1] < wcacut) { // add WCA repulsion if bead is close to bottom wall
      pre = mypow6(wcasig/coor[jj].x[p1][1]);
      f.x[p1][1] += 4.*ep1_6/coor[jj].x[p1][1]*(2.*pre*pre - pre);
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
  int i, j, val2, ti, ind1, ind2;
  
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
}
/*-----------------------------------------------------------------------------*/
void rdsolv(){

  FILE *filein;
  int jj, j, k, good, ind;
  double temp, nudge=1.e-4, rx;
  float ttw;
  filein = fopen(sname,"r");

  for (j=0; j<N; j++) {
    for (k=0; k<3; k++) {
      fscanf(filein,"%f",&ttw);
      rsolv[0][j][k] = ttw;
    }
    for (k=0; k<3; k++) {
      fscanf(filein,"%f",&ttw);
      vsolv[0][j][k] = ttw;
    }
    good = 0;
    while (!good) {
      good = 1;
      if(rsolv[0][j][0]<0){
	rsolv[0][j][0]+=Lx;
	good = 0;
      } else if(rsolv[0][j][0]>Lx){
	rsolv[0][j][0]-=Lx;
	good = 0;
      }
      if(rsolv[0][j][1]<0){
	rsolv[0][j][1]+=Ly;
	good = 0;
      } else if (rsolv[0][j][1]>Ly){
	rsolv[0][j][1]-=Ly;
	good = 0;
      }
      if (rsolv[0][j][1] == 0) {
	rsolv[0][j][1] += nudge;
      } else if(rsolv[0][j][1] == Ly) {
	rsolv[0][j][1] -= nudge;
      }
      if(rsolv[0][j][2]<0){
	rsolv[0][j][2]+=Lz;
	good = 0;
      } else if(rsolv[0][j][2]>Lz){
	rsolv[0][j][2]-=Lz;
	good = 0;
      }
    }
  }

  fclose(filein);

  if (per > 1) {
    for (j=0; j<N; j++) {
      for (k=0; k<3; k++) {
	for (jj=1;jj<per;jj++) {
	  rsolv[jj][j][k] = rsolv[0][j][k];
	  vsolv[jj][j][k] = vsolv[0][j][k];
	}
      }
    }
  }
}
/*-----------------------------------------------------------------------------*/
void rdemerg(){

  FILE *filein;
  int j, part, dim, ind;
  float ttw;
  opoint_t tx;
  point_t px;
  
  filein = fopen(emname,"r");

  for (j=0;j<nemerg;j++) {
    fscanf(filein,"%f",&ttw);
    emerg[j].op = ttw;
    for (part=0;part<npart;part++) {
      for (dim=0;dim<ndim;dim++) {
	fscanf(filein,"%f",&ttw);
	//	if (part < npart) {
	emerg[j].x[part][dim] = ttw;
	px.x[part][dim] = ttw;
	//	}
      }
      for (dim=0;dim<ndim;dim++) {
	//fscanf(filein,"%f",&ttw);
	emerg[j].v[part][dim] = sqrt(-var2*2.*log(ran2(&seed)))*cos(twoPi*ran2(&seed));
	//emerg[j].v[part][dim] = ttw;
      }
      /*      for (dim=0;dim<ndim;dim++) {
	fscanf(filein,"%f",&ttw);
	emerg[j].f1[part][dim] = ttw;
	}*/
    }
    // addpoint(lat,dir,bead,coor,tempw[0],rstart[0],rstart[1]);
    tx = projx(px);
    emerg[j].op = tx.x[0];
  }
  fclose(filein);

  for (j=0;j<per;j++) {
    ind = floor(nemerg*ran2(&seed)) + 1;
    for (part=0;part<npart;part++) {
      for (dim=0;dim<ndim;dim++) {
	coor[j].x[part][dim] = emerg[ind].x[part][dim];
	coor[j].v[part][dim] = emerg[ind].v[part][dim];
	//f1[j][part][dim] = emerg[ind].f1[part][dim];
      }
    }
    ocoor[j] = projx(coor[j]);
  }
}

/*-----------------------------------------------------------------------------*/
void stream(int jj) {
  
  /* This function streams the solvent particles forward */
  
  int i, j, good;
  double temp3, tstar, sum, fac;
  
  sum =0.;
  for(i=0;i<N;i++){
    
    temp3=rsolv[jj][i][1]+vsolv[jj][i][1]*solvstep;         /*first move along y coordinate as it has fixed wall*/

    if (temp3 < 0.) {
      // what time did it cross?; evolve
      tstar = -rsolv[jj][i][1]/vsolv[jj][i][1];

      for(j=0;j<3;j++){
	rsolv[jj][i][j]+=vsolv[jj][i][j]*tstar;
      }

      // pick random velocities
      for (j=0;j<3;j++) {
	//	vsolv[i][j]=sqrt(-var*2.*log(ran2(&seed)))*cos(2*pi*ran2(&seed));
	vsolv[jj][i][j] = -vsolv[jj][i][j];
      }
      if (vsolv[jj][i][1] < 0.) {
	vsolv[jj][i][1] *= -1.;
      }
      for(j=0;j<3;j++){
	rsolv[jj][i][j]+=vsolv[jj][i][j]*(solvstep-tstar);
      }
    } else if (temp3 > Ly) {
      // what time did it cross?; evolve
      tstar = (Ly-rsolv[jj][i][1])/vsolv[jj][i][1];

      for(j=0;j<3;j++){
	rsolv[jj][i][j]+=vsolv[jj][i][j]*tstar;
      }

      // pick random velocities
      for (j=0;j<3;j++) {
	//	vsolv[i][j]=sqrt(-var*2.*log(ran2(&seed)))*cos(2*pi*ran2(&seed));
	vsolv[jj][i][j] = -vsolv[jj][i][j];
      }
      if (vsolv[jj][i][1] > 0.) {
	vsolv[jj][i][1] *= -1.;
      }
      for(j=0;j<3;j++){
	rsolv[jj][i][j]+=vsolv[jj][i][j]*(solvstep-tstar);
      }
    } else {
      for(j=0;j<3;j++){
	rsolv[jj][i][j]+=vsolv[jj][i][j]*solvstep;
      }
    }
    good = 0;
    while (!good) {
      good = 1;
      if(rsolv[jj][i][0]<0){
	rsolv[jj][i][0]+=Lx;
	good = 0;
      } else if(rsolv[jj][i][0]>Lx){
	rsolv[jj][i][0]-=Lx;
	good = 0;
      }
      if((rsolv[jj][i][1]<0) || (rsolv[jj][i][1]>Ly)){
	printf("Error (Ly): i=%d, x,y,z=(%f,%f,%f)!\n",i,rsolv[jj][i][0],rsolv[jj][i][1],rsolv[jj][i][2]);
	gaexit(2);
      }
      if(rsolv[jj][i][2]<0){
	rsolv[jj][i][2]+=Lz;
	good = 0;
      } else if(rsolv[jj][i][2]>Lz){
	rsolv[jj][i][2]-=Lz;
	good = 0;
      }
    }
    for (j=0;j<3;j++) {
      sum += vsolv[jj][i][j]*vsolv[jj][i][j];
    }
  }
  sum = (sum/(float)N)*solvmass;
  fac = sqrt(3./(beta*sum));
  for (i=0;i<N;i++) {
    for (j=0;j<3;j++) {
      vsolv[jj][i][j] *= fac;
    }
  }
}
/*-----------------------------------------------------------------------------*/
void dorotate(int jj) {

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

    if (verb >= 4) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"check0...");
      fclose(mainout);
    }    


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
	  if(blinv*rsolv[jj][i][0]>deltax[0]){
	    pos[i][0]=(int)(blinv*rsolv[jj][i][0]-deltax[0])+1;
	    if (pos[i][0] == nx) {
	      pos[i][0] = 0;
	    }
	  } else {
	    pos[i][0]=0;
	  }
	  j = pos[i][0];
	  if(blinv*rsolv[jj][i][1]>deltax[1]){
	    pos[i][1]=(int)(blinv*rsolv[jj][i][1]-deltax[1])+1;
	  } else {
	    pos[i][1]=0;
	  }
	  k = pos[i][1];
	  if(blinv*rsolv[jj][i][2]>deltax[2]){
	    pos[i][2]=(int)(blinv*rsolv[jj][i][2]-deltax[2])+1;
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
	    if (verb >= 4) {
	      mainout = fopen(outname,"a");
	      fprintf(mainout,"(%d,%d,%d) rsolv[%d][%d][] = (%f,%f,%f)\n ",j,k,l,jj,i,rsolv[jj][i][0],rsolv[jj][i][1],rsolv[jj][i][2]);
	      fclose(mainout);
	    }    
	  }
	}
      } /* end of parallel region */
      if (!good) {
	maxd += 5;
	free(lab);
	//	    printf("proc %d: Increasing maxd to %d\n",me,maxd);
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"%d ",maxd);
	  fclose(mainout);
	}    
	lab = (int ****) alloc4d(sizeof(int),nx,ny,nz,maxd);
	if (verb >= 4) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"! ");
	  fclose(mainout);
	}    
	for (j=0;j<nx;j++) {
	  for (k=0;k<ny;k++) {
	    for (l=0;l<nz;l++) {
	      nlab[j][k][l] = 0;
	    }
	  }
	}
      }
    }

    if (verb >= 4) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"check1...");
      fclose(mainout);
    }    
    /*----------------------------What box is each polymer bead in?-------------------------*/
    
    err = 0;
#pragma omp parallel for default(shared) private(i) schedule(static)
    for(i=0;i<npart;i++){
      if(blinv*coor[jj].x[i][0]>deltax[0]){
	pos_poly[i][0]=(int)(blinv*coor[jj].x[i][0]-deltax[0])+1;
      } else {
	pos_poly[i][0]=0;
      }
      if(blinv*coor[jj].x[i][1]>deltax[1]){
	pos_poly[i][1]=(int)(blinv*coor[jj].x[i][1]-deltax[1])+1;
      } else {    
	pos_poly[i][1]=0;
      }
      if(blinv*coor[jj].x[i][2]>deltax[2]){
	pos_poly[i][2]=(int)(blinv*coor[jj].x[i][2]-deltax[2])+1;
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

    if (verb >= 4) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"check2...");
      fclose(mainout);
    }    

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
	v_cmax[i][j][k]+=solvmass*vsolv[jj][m][0];
	v_cmay[i][j][k]+=solvmass*vsolv[jj][m][1];
	v_cmaz[i][j][k]+=solvmass*vsolv[jj][m][2];
      }
	
#pragma omp parallel for default(shared) private(i,j,k,n) schedule(static)
      for(n=0;n<npart;n++){
	i= pos_poly[n][0];
	j= pos_poly[n][1];
	k= pos_poly[n][2];
	
	v_cmax[i][j][k]+= polymass*coor[jj].v[n][0];
	v_cmay[i][j][k]+= polymass*coor[jj].v[n][1];
	v_cmaz[i][j][k]+= polymass*coor[jj].v[n][2];
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
		  odisp[0] = vsolv[jj][ind][0] - v_cmx[i][j][k];
		  odisp[1] = vsolv[jj][ind][1] - v_cmy[i][j][k];
		  odisp[2] = vsolv[jj][ind][2] - v_cmz[i][j][k];
		  ndisp = distchg(odisp,rndvec,rr);
		  vsolv[jj][ind][0] = v_cmx[i][j][k] + ndisp.x[0];
		  vsolv[jj][ind][1] = v_cmy[i][j][k] + ndisp.x[1];
		  vsolv[jj][ind][2] = v_cmz[i][j][k] + ndisp.x[2];
		  if ((j!=0)&&(j!=ny-1)) {
		    vsolv[jj][ind][0] += solvstep*grav;  
		  }
		}
			    
		// polymer beads
		
#pragma omp for schedule(static)	      
		for (ind=0;ind<npart;ind++) {
		  if (((pos_poly[ind][0]==i)&&(pos_poly[ind][1]==j))&&(pos_poly[ind][2]==k)) {
		    odisp[0] = coor[jj].v[ind][0] - v_cmx[i][j][k];
		    odisp[1] = coor[jj].v[ind][1] - v_cmy[i][j][k];
		    odisp[2] = coor[jj].v[ind][2] - v_cmz[i][j][k];
		    ndisp = distchg(odisp,rndvec,rr);
		    coor[jj].v[ind][0] = v_cmx[i][j][k] + ndisp.x[0];
		    coor[jj].v[ind][1] = v_cmy[i][j][k] + ndisp.x[1];
		    coor[jj].v[ind][2] = v_cmz[i][j][k] + ndisp.x[2];
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
  double temp;
  
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

  int MaxNum = beads;

//  printf("proc %d is creating coord servers...\n",me); fflush(stdout);

  status = CreateCoordServer_dbl(&Coords_extrawt, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_weight, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nout, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nlist, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_avwalk, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_navg, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_rcount, &MaxNum, "test");

  MaxNum = 2*beads*beads;
  status = CreateCoordServer_int(&Coords_nin, &MaxNum, "test");
  MaxNum = 2*beads*2*beads*2*beads;
  status = CreateCoordServer_int(&Coords_prob, &MaxNum, "test");
  MaxNum = 2*beads*2*beads;
  status = CreateCoordServer_int(&Coords_nprob, &MaxNum, "test");

  MaxNum = 1;
  status = CreateCoordServer_dbl(&Coords_narr, &MaxNum, "test");

  MaxNum = per*nproc;
  status = CreateCoordServer_int(&Coords_reg, &MaxNum, "test");
  
  MaxNum = beads*nop;
  status = CreateCoordServer_dbl(&Coords_stravg, &MaxNum, "test");
  status = CreateCoordServer_dbl(&Coords_z, &MaxNum, "test");
  
  MaxNum = tres;
  status = CreateCoordServer_dbl(&Coords_thist, &MaxNum, "test");
  
  MaxNum = beads*mxlist*3*npart*ndim;
  status = CreateCoordServer_dbl(&Coords_pts, &MaxNum, "test");

  MaxNum = beads*mxlist;
  status = CreateCoordServer_dbl(&Coords_wlist, &MaxNum, "test");

  MaxNum = tres + 2*beads + 1;
  status = GA_Create_mutexes(MaxNum);

  /*  lock indices:   one for each region, and each of the 5 variables ->

0: narr
1 to 2*beads:  extrawt

    */

//  printf("all coord servers created by proc %d\n",me); fflush(stdout);
}

void initvar() {

  int i, j, k, ind, pt, op, mcount, dir1, bead1, dir2, bead2;

  /* initialize arrays */

  basin[0].x[0] = bmin;
  basin[1].x[0] = bmax;

  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"Initializing variables...\n");
    fclose(mainout);
  }

  for (k=0;k<per;k++) {
    for(i=0;i<npart;i++){
      for(j=0;j<3;j++){
	f1[k][i][j]=0.;
      }
    }
  }
  
  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      navg[j][k] = 0;
      for (op=0; op<=nop-1;op++) {
	stravg[j][k].x[op] = 0.;     /* super-globals */
      }
      for (op=0; op<=nop-1; op++) {
	z[j][k].x[op] = 0.;
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

  int i, j, mcount;

/* get string from server */

  GetFromServer_dbl(&Coords_z,0,-1,tz1);
  GetFromServer_dbl(&Coords_z,1,-1,tz2);

  mcount = 0;
  for (i=0;i<beads;i++) {
    for (j=0;j<nop;j++) {
      z[0][i].x[j] = tz1[mcount];
      z[1][i].x[j] = tz2[mcount];
      mcount++;
      if (verb >= 3) {
	mainout = fopen(outname,"a");
	fprintf(mainout,"%d %f %f\n", i,tz1[mcount-1],tz2[mcount-1]);
	fclose(mainout);
      }
    }
  }
}
void wrxyz(int step, int jj) {
  
  /* This function writes a point to the xyz function open in 'xyzout' */

  int j, k;

  xyzout = fopen(filename,"a");
  if (step/xyzfrq != 0) {
    fprintf(xyzout,"\n %d \n",npart);
  }
  for(j=0; j< npart; j++){
    if(j<9){
      fprintf(xyzout," atom00%i", j+1);
      for(k=0; k<3; k++){
	if(coor[jj].x[j][k] >= 0 ){
	  fprintf(xyzout,"   %e", coor[jj].x[j][k]);
	}
	if(coor[jj].x[j][k] < 0 ){
	  fprintf(xyzout,"  %e", coor[jj].x[j][k]);
	}
      }
      fprintf(xyzout,"\n");
    }
    if(j<99 && j >= 9){
      fprintf(xyzout," atom0%i", j+1);
      for(k=0; k<3; k++){
	if(coor[jj].x[j][k] >= 0 ){
	  fprintf(xyzout,"   %e", coor[jj].x[j][k]);
	}
	if(coor[jj].x[j][k] < 0 ){
	  fprintf(xyzout,"  %e", coor[jj].x[j][k]);
	}
      }
      fprintf(xyzout,"\n");
    }
    if(j>=99){
      fprintf(xyzout," atom%i", j+1);
      for(k=0; k<3; k++){
	if(coor[jj].x[j][k] >= 0 ){
	  fprintf(xyzout,"   %e", coor[jj].x[j][k]);
	}
	if(coor[jj].x[j][k] < 0 ){
	  fprintf(xyzout,"  %e", coor[jj].x[j][k]);
	}
      }
      fprintf(xyzout,"\n");
    }
  }
  fclose(xyzout);
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

int chkdir(int jj) {
  
  /* returns 1 if there is a basin crossing, 0 otherwise */
  int temp, tbas, odir, n1, n2, olat, newreg, ind, mydir, mybead;

  mydir = dir[jj];
  mybead = bead[jj];

  temp = 0;
  tbas = bascheck(jj);
		    
  if (tbas + mydir == 1) {  /* either tbas = 0, dir = 1 or tbas = 1, dir = 0 */
    ind = mydir*beads + mybead;
    temp = 1;
			
    /* basin crossing */
		
    lo[0] = mydir;
    hi[0] = mydir;
    lo[1] = 0;
    hi[1] = 0;
    GA_Lock(0);
    GA_Init_fence();
    NGA_Get(Coords_narr.ga,lo,hi,tw,ld);
    tw[0] += cycle*chkfrq;
    NGA_Put(Coords_narr.ga,lo,hi,tw,ld);
    GA_Fence();
    GA_Unlock(0);

    if (mydir == 1) {
      odir = 0;
    }else{
      odir = 1;
    }
	
    newreg = whichx(odir,coor[jj]);
    dir[jj] = odir;
    bead[jj] = newreg;
    if (verb >= 3) {
      //  if (((lat==0)&&(dir==0))&&(bead==0)) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Ended traj(%d): dir cross into %d\n",jj,newreg);
      fclose(mainout);
    }
    n1 = rstart[jj][0];
    n2 = rstart[jj][1];
    
    allind[0] = 0;
    allind[1] = 4*beadsp3*dir[jj] + 4*beadsp2*bead[jj] + 2*beadsp2*n1 + 2*beads*n2 + beads*odir + newreg;
    tempi[0] = 1 + NGA_Read_inc(Coords_prob.ga,allind,1);			  
    
    allind[1] = 2*beadsp2*dir[jj] + 2*beads*bead[jj] + beads*n1 + n2;
    tempi[0] = 1 + NGA_Read_inc(Coords_nprob.ga,allind,1);
    done[jj] = 1;
  }
  return temp;
}

void chkprim(int jj) {

  int newreg, ind, n1, n2, mybead;

  mybead = bead[jj];
  
  /* prim crossing */
  newreg = which(jj);
  if (newreg != mybead) {

    if (verb >=3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Ended traj: prim cross into %d\n",newreg);
      fclose(mainout);
    }
    allind[0] = 0;
    n1 = rstart[jj][0];
    n2 = rstart[jj][1];
		    
    allind[1] = 4*beadsp3*dir[jj] + 4*beadsp2*bead[jj] + 2*beadsp2*n1 + 2*beads*n2 + beads*dir[jj] + newreg;
    tempi[0] = 1 + NGA_Read_inc(Coords_prob.ga,allind,1); 
    
    allind[1] = 2*beadsp2*dir[jj] + 2*beads*bead[jj] + beads*n1 + n2;
    tempi[0] = 1 + NGA_Read_inc(Coords_nprob.ga,allind,1); 

    rstart[jj][0] = dir[jj];
    rstart[jj][1] = bead[jj];
    bead[jj] = newreg;
  }
}
void addpoint (int jj) {

  /* this function adds the point 'x' to the (dir,bead) point list
     along with the weight 'wt' */
  
  int tgt, n, i, ind, tfull, wind, part, dim, ptind, mcount, tn;
  double ttw[1], t1,t2;
  
  allind[0] = dir[jj];
  allind[1] = bead[jj];
  
  tnlist = 1 + NGA_Read_inc(Coords_nlist.ga,allind,1);

  wind = mxlist*bead[jj] + (tnlist - 1);

  tw[0] = myweight[jj];
  //  printf("proc %d: %d %d tnlist = %d\n",me,dir,bead,tnlist); fflush(stdout);
  PutToServer_dbl(&Coords_wlist,dir[jj],wind,tw);

  mcount = 0;
  for (part=0; part<npart;part++) {
    for (dim=0; dim<=ndim-1; dim++) {
      tpt[mcount] = coor[jj].x[part][dim];  // position
      mcount++;
      tpt[mcount] = coor[jj].v[part][dim];  // velocity
      mcount++;
      tpt[mcount] = f1[jj][part][dim];  // old force
      mcount++;
    }
  }	    
    
  lo[0] = dir[jj];
  hi[0] = dir[jj];
  lo[1] = 3*npart*ndim*mxlist*bead[jj] + 3*npart*ndim*(tnlist-1);
  hi[1] = 3*npart*ndim*mxlist*bead[jj] + 3*npart*ndim*tnlist -1;
  NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);

}
/*-------------------------------------------------------------------------*/
int getprob(int dir, int bead, int i, int j, int n1, int n2) {
  
  int ind;
  int temp;

  if (((dir > 1) || ((i > 1) || (n1 > 1))) || (((bead >= beads) || (j >= beads)) || (n2 >= beads))) {
    printf("Error in getprob!  Array indices out of bounds\n"); fflush(stdout);
    printf("%d %d %d %d %d %d\n",dir,bead,i,j,n1,n2); fflush(stdout);
    gaexit(98);
  }
  ind = 4*beadsp3*dir + 4*beadsp2*bead + 2*beadsp2*i + 2*beads*j + beads*n1 + n2;
  temp = tprob1[ind];

  return(temp);
}
int getnprob(int dir, int bead, int i, int j) {

  int ind;
  int temp;

  if (((dir > 1) || (i > 1)) || ((bead >= beads) || (j >= beads)))  {
    printf("Error in getnprob!  Array indices out of bounds\n"); fflush(stdout);
    printf("%d %d %d %d\n",dir,bead,i,j); fflush(stdout);
    gaexit(97);
  }
  ind = 2*beadsp2*dir + 2*beads*bead + beads*i + j; 

  temp = tnprob1[ind];

  return(temp);
}
