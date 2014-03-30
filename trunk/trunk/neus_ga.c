#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> /* required for randomize() and random() */
#include <time.h>
#include <unistd.h>
#include "neusglob.h"
#include "rna.h"
#include "../aryeh.h"
#include "../myCoordServer.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

/* external functions -- system dependent */
void initsystem();
int strinit();
void move();
double gethist1(point_t coor);
double gethist2(point_t coor);
double gethist3(point_t coor);
/* end external */

/* super-global variables */

ServedCoords Coords_weight;
ServedCoords Coords_wadj;
ServedCoords Coords_wavg;
ServedCoords Coords_elaps;
ServedCoords Coords_z;
ServedCoords Coords_count;
ServedCoords Coords_narr;
ServedCoords Coords_on;
ServedCoords Coords_navg;
ServedCoords Coords_bdedit;
ServedCoords Coords_nflux;
ServedCoords Coords_hist1;
ServedCoords Coords_hist2;
ServedCoords Coords_hist3;
ServedCoords Coords_full;
ServedCoords Coords_nlist;
ServedCoords Coords_twlist;
ServedCoords Coords_pts;
ServedCoords Coords_wlist;
ServedCoords Coords_from;
ServedCoords Coords_clock;
ServedCoords Coords_mytim;

double myweight=0.;
double **weight, **wavg;
opoint_t **z;
int **count, **on, **navg, **fmat;
int ***thist, ***thist2, ***thist3;
int allind[2];

/* global variables */

double *twavg1, *tweight1;
double *tz1, *tpt;
double globtim=0.;
int nglobtim=0;
int *tint1, *tint2, *ton1, *ton2, *tthist1, *tthist2, *tthist3, *tthist5, ****nin, **nout, *tnflux1, *tnflux2;
int *tmpfrom, **allcnt, **allon;
int verb=2;
int me, nproc, status, beadsp2, beadsp3, lockind, lo[2], hi[2], ld[1];
double tempw[1];

const int mxlist=MXLIST_c;

int bas, mytim=0, dynerr=0;
opoint_t ocoor;
point_t coor;
int myclock=0;
void scalflux();
void wtmat3();
opoint_t projx(point_t x);
opoint_t projxsig(point_t x);
void wrtflux(int find, int cycle);
int bascheck();
int which(int dir);
int whichx(int dir, point_t x);
void gaexit(int sig);
void addpoint (int dir, int bead, point_t x, double wt, int from1, int from2);
double calcsum(int dir);
double ogtdist(opoint_t px, opoint_t py);
void getstring();
void chkprim(int dir, int bead);
int chkdir(int dir, int bead);
int getlockind(int var, int dir, int bead);
void wtread(), fxread(), nfluxread(), wrpath(int ind);
int rstart[2], tempi[1], tempi2[1], endtraj;

void **alloc2d(int varsize, int n, int p) ; 
void ***alloc3d(int varsize, int n, int p, int q) ;
void ****alloc4d(int varsize, int n, int p, int q, int r) ;
void *****alloc5d(int varsize, int n, int p, int q, int r, int s) ;
void *******alloc7d(int varsize, int n, int p, int q, int r, int s, int t, int u) ;
double gssran(long *idum);
double mypow2(double x);
double mypow4(double x);
double mypow6(double x);
double mypow7(double x);
double mypow8(double x);
double mypow14(double x);
double ran2(long *idum) ;
FILE *mainout;
char fname[30];
char outname[30];
char ratname[30];
char onname[30], twname[30];
char wtname[30];
FILE *outPath, *outFile;

/*-------------------------------------------------------------------------*/

int main(int argc,char** argv){

  int i, j, k, l, m, dim, n1, n2, part, lim;
  int dir1, dir2, from, to, limit, pt, bead1, bead2, tnlist;
  int dir, bead, tcount, listind;
  int ind, mcount;
  int cycle, endcycle;
  unsigned int iseed = (unsigned int)time(NULL);
  int back(int dir, int bead);
  void createservers(), initvar();
  void stack(int dir, int bead);
  int getmin(int *dir, int *bead);
  void wtstack(int dir, int bead);
  void updtwts(int ind);
  void wrhist(int index);

  double rate[2], sum1, sum2, t1, t2, telaps;
  FILE * ratFile, *filein, *onfile, *twfile;

  ld[0] = 0;

  if (seed ==0) {
    seed = iseed;
  }

  beadsp2 = beads*beads;
  beadsp3 = beadsp2*beads;

#ifdef _OPENMP
  int desired = MPI_THREAD_MULTIPLE;
  int provided;
  MPI_Init_thread(&argc, &argv, desired, &provided);
#else
  MPI_Init(&argc, &argv);
#endif

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  
  GA_Initialize_ltd(256*1024*1024*sizeof(double));

  createservers();
  initsystem();

  GA_Sync();
  if (me == 0) {
    printf("Starting..\n"); fflush(stdout);
  }

  /* open files for output */

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
  
  /* declare arrays, initialize variables */

  initvar();
  
  if (me == 0) {
    GA_Init_fence();
    
    readtoglobal();
    printf("proc %d is done reading global info\n",me); fflush(stdout);
    GA_Fence();
  }

  /* check point */
  GA_Sync();

  if (me == 0) {
    t1 = MPI_Wtime();
  }
  getstring();

  /* ---------------start of dynamics loop-------------------------- */
  
  for (cycle=startcycle;cycle<=T;cycle++) { 
    
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
    GA_Sync();

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
	ind = dir*beads + bead;
	GetFromServer_int(&Coords_nlist,0,ind,tempi);
	tnlist = tempi[0];
	GetFromServer_int(&Coords_nlist,1,ind,tempi);
	tnlist += tempi[0];

	endtraj = 0;
	if (tnlist) {
	  if (!back(dir,bead)) {;  /* initialize trajectory */
	    endtraj = 1;
	  } else {
	    ind = dir*beads + bead;
	    GetFromServer_int(&Coords_count,0,ind,tempi);
	    if (tempi[0] > every) {
	      endtraj = 1;
	    }
	  }
	} else {
	  endtraj = 1;
	  ind = dir*beads + bead;
	  tempi[0] = 0;
	  PutToServer_int(&Coords_on,0,ind,tempi);
	  printf("%d: turning off (%d,%d)! \n",me,dir,bead);
	}
	
	while (!endtraj) { 

	  allind[0] = 0;
	  allind[1] = dir*beads + bead;
	  tcount = chkfrq + NGA_Read_inc(Coords_count.ga,allind,chkfrq);
	  tempi[0] = chkfrq + NGA_Read_inc(Coords_elaps.ga,allind,chkfrq);

	  /* move */
	  for (i=0;i<chkfrq;i++) {  
	    if (!dynerr) {

	      move();
	      mytim++;
	      myclock++;
	    }
	  }
	  if (!dynerr) {

	    /* check */	
	    if (!chkdir(dir,bead)) {
	      chkprim(dir,bead);
	    }
	
	    if (tcount%stkfrq == 0) {
	      stack(dir,bead);
	    }
		
	    if (tcount >= every) {
	      if (verb >= 3) {
		mainout = fopen(outname,"a");
		fprintf(mainout,"Ended traj: count\n");
		fclose(mainout);
	      }
	      /* write coordinates to (someone's) fluxlist */
	      ind = dir*beads + bead;
		
	      if (!chkdir(dir,bead)) {
		chkprim(dir,bead);
	      }
	      if (!endtraj) {
		addpoint(dir,bead,coor,myweight,rstart[0],rstart[1]);
		endtraj = 1;
	      }
	    }
	  } else {   /* if dynerr */
	    if (verb >= 3) {
	      mainout = fopen(outname,"a");
	      fprintf(mainout,"dynerr!\n");
	      fclose(mainout);
	    }
	    endtraj = 1;
	    dynerr = 0;
	    printf("dynerr: (%d,%d) from (%d,%d)\n",dir,bead,rstart[0],rstart[1]);
	    gaexit(2);
	  }
	}  /* end of trajectory loop */
      }
    }       /* end of cycle */
    
    /* compute the rate / update weights / write hists */
    GA_Sync();
    if (me == 0) {
      printf("proc %d: TCOB\n",me); fflush(stdout);

      /* update weights */

      if (cycle%wtupdt == 0) {
	updtwts(cycle/wtupdt);
	scalflux();
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
	  if (tweight1[mcount] != 0.) {
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

	/* write out hists */

	wrhist(cycle/wrfrq);
            
	/* write nflux */

	GetFromServer_int(&Coords_nflux,0,-1,tnflux1);
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
	fclose(outFile);
	
	printf("proc %d: done TCOB\n",me); fflush(stdout);
      }
    } else if (me <= 2) {

      /* write fluxes, scale */
      if (cycle%wrfrq == 0) {
	if ((cycle/wrfrq)%2 == 0) {
	  wrtflux(me,cycle);
	}
      }
    }

    /* checkpoint */
    GA_Sync();
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
  status = DestroyCoordServer(&Coords_navg);

  status = DestroyCoordServer(&Coords_z);

  status = DestroyCoordServer(&Coords_hist1);
  status = DestroyCoordServer(&Coords_hist2);
  status = DestroyCoordServer(&Coords_hist3);
  status = DestroyCoordServer(&Coords_full);
  status = DestroyCoordServer(&Coords_nlist);
  status = DestroyCoordServer(&Coords_twlist);
  status = DestroyCoordServer(&Coords_pts);
  status = DestroyCoordServer(&Coords_wlist);
  status = DestroyCoordServer(&Coords_from);
  status = DestroyCoordServer(&Coords_clock);
  status = DestroyCoordServer(&Coords_mytim);

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
  
  double wsum, sum, a, ttwlist, oran, t1, t2, tim2, tmpwlist[mxlist];
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
    t1 = tempw[0];
    GetFromServer_dbl(&Coords_twlist,1,ind,tempw);
    t2 = tempw[0];
    oran = (t1+t2)*ran2(&seed);
    if (oran < t1) {
      listind = 0;
    } else {
      listind = 1;
    }
    lockind = getlockind(0,dir,bead);
    GA_Lock(lockind);
    GA_Init_fence();
    GetFromServer_int(&Coords_nlist,listind,ind,tempi);
    tnlist = tempi[0];
    GetFromServer_int(&Coords_full,listind,ind,tempi);
    tfull = tempi[0];
    GA_Fence();
    GA_Unlock(lockind);
    
    ttwlist = 0.;
    if (tnlist != 0) {
      wsum = 0.;
      oran = ran2(&seed);
      a = oran;

      // get complete wlist
      lo[0] = listind;
      hi[0] = listind;
      lo[1] = mxlist*beads*dir + mxlist*bead;
      hi[1] = mxlist*beads*dir + mxlist*(bead + 1) - 1;

      lockind = getlockind(4+listind,dir,bead);
      GA_Lock(lockind);
      GA_Init_fence();
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
	while (wsum < a) {
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
	  wsum += tmpwlist[j-1];
	  j++;
	}

	wind = mxlist*beads*dir + mxlist*bead + j-2;      
	myweight = spltfac*tmpwlist[j-2];
	tempw[0] = (1.-spltfac)*tmpwlist[j-2];
	if (tempw[0] < 1.e-300) {
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
	sum = 0.;
	for (i=0;i<limit;i++) {
	  sum += tmpwlist[i];
	}
	if (fabs(sum - tempw[0]) > 0.001*sum) {
	  tempw[0] = sum;
	  PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	}

	GA_Fence();
	GA_Unlock(lockind);
	
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
	newreg = which(dir);
	if (newreg != bead) {
	  if (verb >= 2) {
	    printf("proc %d:  Error:  Replica reinitialized in wrong region!\n",me);
	    printf("which(%d) = %d\n", dir, newreg);
	    printf("opoint: %e %e %e\n",ocoor.x[0],z[dir][bead].x[0],z[dir][newreg].x[0]);
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
	    lockind = getlockind(4+listind,dir,bead);
	    GA_Lock(lockind);
	    GA_Init_fence();
	    GetFromServer_dbl(&Coords_twlist,listind,ind,tempw);
	    tempw[0] -= tmpwlist[j-2];
	    PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	    GA_Fence();
	    GA_Unlock(lockind);
	  } else {
	    printf("j <=2!!\n");
	    printf("wsum = %e; t1, t2 = %e, %e; a = %e\n",wsum,t1,t2,a);
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
	fprintf(mainout,"in region (%d,%d), no points in list %d (%f,%f,%f)!\n",dir,bead,listind,t1,t2,oran);
	fclose(mainout);
      }
      gotit = 0;
    }
  }
  return gotit;
}
/*-------------------------------------------------------------------------*/
void updtwts(int ind) {
   
  double lim2, sum;
  int tgt,dir,bead,good,lim,i, mcount,j,n1,n2,mcount2, nav;
  FILE * wtout, *tmpfile;
  char tmpname[30];

  GA_Init_fence();
  lim = 2*beads;
  
  if (wtalg == 3) {

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

	/* reset tnflux */

	if (clrtnfx) {
	  tnflux1[mcount] = 0;
	}

	mcount++;
      }
    }
    if (clrtnfx) {
      PutToServer_int(&Coords_nflux,0,-1,tnflux1);
      PutToServer_int(&Coords_nflux,1,-1,tnflux1);
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
	tweight1[mcount] = twavg1[mcount]/sum;
	weight[dir][bead] = tweight1[mcount];
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

  int i, j, k, lim, ind1, n1, n2, ind2, ind3, offset, nreps=5000;
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
void wrhist(int index) {
    
  /* this function writes out the conditional and total theta histograms 
    for the normal order parameter and the sig-modified o.p.*/

  int dir, bead, ind, mcount, mcount2, tsum;
  double x[tres], x2[tres], x3[tres], sum, sum2, sum3, ext, ootsum;
  FILE * tout;
  char tmpname[30];
  
  GetFromServer_dbl(&Coords_weight,0,-1,tweight1);
  GetFromServer_int(&Coords_hist1,0,-1,tthist1);
  GetFromServer_int(&Coords_hist2,0,-1,tthist2);
  GetFromServer_int(&Coords_hist3,0,-1,tthist3);
  
  mcount = 0;
  mcount2 = 0;
  for (dir=0;dir<=1;dir++) {
    for (bead=0;bead<=beads-1;bead++) {
      weight[dir][bead] = tweight1[mcount];
      mcount++;
      for (ind=0;ind<tres;ind++) {
	thist1[dir][bead][ind] = tthist1[mcount2];
	thist2[dir][bead][ind] = tthist2[mcount2];
	thist3[dir][bead][ind] = tthist3[mcount2];
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
      
    /* thist1 */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist1[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist1[dir][bead][ind]) {
	  x[ind] += thist1[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }

    /* thist2 */
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

    /* thist3 */
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
      
    /* thist1 */
    tsum = 0;
    for (ind=0;ind<=tres-1;ind++) {
      tsum += thist1[dir][bead][ind];
    }
    if (tsum) {
      ootsum = 1./(float)tsum;
      for (ind=0;ind<=tres-1;ind++) {
	if (thist1[dir][bead][ind]) {
	  x[ind] += thist1[dir][bead][ind]*weight[dir][bead]*ootsum;
	}
      }
    }

    /* thist2 */
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

    /* thist3 */
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
      
      /* thist1 */
      tsum = 0;
      for (ind=0;ind<=tres-1;ind++) {
	tsum += thist1[dir][bead][ind];
      }
      if (tsum) {
	ootsum = 1./(float)tsum;
	for (ind=0;ind<=tres-1;ind++) {
	  if (thist1[dir][bead][ind]) {
	    x[ind] += thist1[dir][bead][ind]*weight[dir][bead]*ootsum;
	  }
	}
      }
      
      /* thist2 */
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

      /* thist3 */
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

  PutToServer_int(&Coords_hist1,0,-1,tthist1);
  PutToServer_int(&Coords_hist2,0,-1,tthist1);
  PutToServer_int(&Coords_hist3,0,-1,tthist1);

}
/*-------------------------------------------------------------------------*/
int which(int dir) {

  /* this function returns the region that the walker
   in (dir,bead) is in, presently */
       
  int i, j, temp;
  double dist, lodist;

  for (i=0;i<beads;i++) {
    dist = ogtdist(z[dir][i],ocoor);
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
int whichx(int dir, point_t x) {
   
   /* this function returns the region that the point 'x'
      is in, presently */
   
   int i, j, temp;
   double dist, lodist;
   opoint_t ox;
   
   ox = projx(x);
   
   for (i=0;i<beads;i++) {
     dist = ogtdist(z[dir][i],ox);
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
void addpoint (int dir, int bead, point_t x, double wt, int from1, int from2) {

  /* this function adds the point 'x' to the (dir,bead) fluxlist
     along with the weight 'wt' */
  
  int tgt, n, i, ind, tfull, wind, part, dim, ptind, mcount, tn, tfrom2[2], listind;
  double tw[1], ttw[1], t1,t2, extra, newextra, tmpwlist[mxlist], sum, oosum;
  
  ind = dir*beads + bead;

  if ((dir==0)&&(bead==0)) {
    if (from1==1) {
      listind = 0;
    } else {
      listind = 1;
    }
  } else if ((dir==1)&&(bead==beads-1)) {
    if (from1==0) {
      listind = 1;
    } else {
      listind = 0;
    }
  } else {
    if (from2 < bead) {
      listind = 0;
    } else {
      listind = 1;
    }
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

    lockind = getlockind(4+listind,dir,bead);
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

    // flist[dir][bead].twlist += wt - flist[dir][bead].wlist[n];
    GA_Fence();
    GA_Unlock(lockind);

    //flist[dir][bead].wlist[n] = wt;      

  } else {
    
    lockind = getlockind(4+listind,dir,bead);
    GA_Lock(lockind);
    GA_Init_fence();

    tw[0] = wt;
    PutToServer_dbl(&Coords_wlist,listind,wind,tw);
    //flist[dir][bead].wlist[n] = wt;

    GetFromServer_dbl(&Coords_twlist,listind,ind,ttw);
    ttw[0] += wt;
    PutToServer_dbl(&Coords_twlist,listind,ind,ttw);

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
   
   /* this function adds the current state of the replicas to the histograms*/
   
   int i, j, index, ind;
   double temps[1], mx, mn, tw[1];
   point_t rt;
   double tpt;

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

   tpt = gethist1(coor);
   index = tres*(tpt-tmin[0])/(tmax[0]-tmin[0]);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index (1)! %f\n",tpt);
   } else {
     allind[0] = 0;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_hist1.ga,allind,1);
   }

   tpt = gethist2(coor);
   index = tres*(tpt-tmin[1])/(tmax[1]-tmin[1]);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index (2)! %f\n",tpt);
   } else {
     allind[0] = 0;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_thist2.ga,allind,1);
   }

   tpt = gethist3(coor);
   index = tres*(tpt-tmin[2])/(tmax[2]-tmin[2]);
   if ((index < 0) || (index >= tres)) {
     printf("Error in stack index (3)! %f\n",tpt);
   } else {
     allind[0] = 0;
     allind[1] = tres*beads*dir + tres*bead + index;
     tempi[0] = 1 + NGA_Read_inc(Coords_thist3.ga,allind,1);
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
	fprintf(outPath,"%e ",z[j][k].x[op]);
      }
      fprintf(outPath,"\n");
    }
  }
  fclose(outPath);

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
  double temp;
  
  wtin = fopen(wname,"r");

  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      fscanf(wtin,"%le",&temp);
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
/*-------------------------------------------------------------------------*/
void fxread() {

  /* this function reads the initial values of the fluxes
     from the input file 'fxname' */

  int j, k, op, dir, ti, tj, tk, limit, pt, dim, part, ind, mcount, listind;
  int scale;
  double fac;
  FILE * fxin;
  float temp;
  float tw;
  
  fxin = fopen(fxname,"r");

  for (j=0; j<=1; j++) {
    for (k=0; k<=beads-1; k++) {
      fscanf(fxin,"%d %d",&tj,&tk);
      if ((tj != j) || (tk != k)) {
	printf("Error!\n");
	gaexit(2);
      }
      for (listind=0;listind<2;listind++) {
	  
	fscanf(fxin,"%d %f", &tempi[0], &tw);
	tempw[0] = tw;
	printf("(%d,%d,%d) = %e\n",j,k,listind,tempw[0]);
	ind = j*beads + k;
	//LINE_STAMP;
	PutToServer_int(&Coords_nlist,listind,ind,tempi);
	//LINE_STAMP;
	PutToServer_dbl(&Coords_twlist,listind,ind,tempw);
	
	if (tempi[0] == mxlist) {
	  tempi2[0] = 1;
	} else {
	  tempi2[0] = 0;
	}
	//LINE_STAMP;
	PutToServer_int(&Coords_full,listind,ind,tempi2);
	limit = tempi[0];
	if (limit) {
	  for (pt=0; pt<=limit-1; pt++) {
	    fscanf(fxin,"%d %d ", &tempi[0], &tempi2[0]);
	    ind = 2*mxlist*beads*j + 2*mxlist*k + 2*pt;
	    PutToServer_int(&Coords_from,listind,ind,tempi);
	    PutToServer_int(&Coords_from,listind,ind+1,tempi2);
	    
	    fscanf(fxin,"%f %d ",&tw,&tempi[0]);
	    tempw[0] = tw;
	    ind = mxlist*beads*j + mxlist*k + pt;
	    PutToServer_dbl(&Coords_wlist,listind,ind,tempw);
	    PutToServer_int(&Coords_clock,listind,ind,tempi);
	    fscanf(fxin,"%d ",&tempi[0]);
	    PutToServer_int(&Coords_mytim,listind,ind,tempi);
	    
	    mcount = 0;
	    for (part=0; part<npart;part++) {
	      for (dim=0; dim<=ndim-1; dim++) {
		fscanf(fxin,"%f ", &tw);
		tpt[mcount] = tw;
		mcount++;
		fscanf(fxin,"%f ", &tw);
		tpt[mcount] = tw;
		mcount++;
		fscanf(fxin,"%f ", &tw);
		tpt[mcount] = tw;
		mcount++;
	      }
	    }	    
	    
	    lo[0] = listind;
	    hi[0] = listind;
	    lo[1] = 3*npart*ndim*mxlist*beads*j + 3*npart*ndim*mxlist*k + 3*npart*ndim*pt;
	    hi[1] = 3*npart*ndim*mxlist*beads*j + 3*npart*ndim*mxlist*k + 3*npart*ndim*(pt+1) -1;
	    NGA_Put(Coords_pts.ga,lo,hi,tpt,ld);
	  }
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
  
  /* this function checks if the walker is in a basin
     |||| currently set up for a 1D system ||||
 */

  int temp, i;
  double dist;

  temp = -1;

  if (ocoor.x[0] < basin[0].x[0]) {
    temp = 0;
  } else if (ocoor.x[0] > basin[1].x[0]) {
    temp = 1;
  } else if (bwidth > 0.) {
    for (i=0;i<=1;i++) {
      dist = ogtdist(ocoor,basin[i]);
      if (dist < bwidth) {
	temp = i;
      }
    }
  }  

  return temp;
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
  status = CreateCoordServer_dbl(&Coords_twlist, &MaxNum, "test");
  
  status = CreateCoordServer_int(&Coords_elaps, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_count, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_narr, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_on, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_bdedit, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_full, &MaxNum, "test");
  status = CreateCoordServer_int(&Coords_nlist, &MaxNum, "test");
  
  MaxNum = beads;
  status = CreateCoordServer_int(&Coords_navg, &MaxNum, "test");

  MaxNum = 2*beads*2*beads;
  status = CreateCoordServer_int(&Coords_nflux, &MaxNum, "test");
  
  MaxNum = 2*beads*nop;
  status = CreateCoordServer_dbl(&Coords_z, &MaxNum, "test");
  
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

  MaxNum = 6*2*beads;
  status = GA_Create_mutexes(MaxNum);

  /*  lock indices:   one for each region, and each of the 5 variables ->

0: nlist/full
2: wadj / wavg
4/5: twlist    */

//  printf("all coord servers created by proc %d\n",me); fflush(stdout);
}

void initvar() {

  int i, j, k, ind, pt, op, mcount, dir1, bead1, dir2, bead2, dir, bead, lim;

  lim = 2*beads;
  twavg1 = (double *) malloc(2*beads*sizeof(double));
  tweight1 = (double *) malloc(2*beads*sizeof(double));
  tz1 = (double *) malloc(2*beads*nop*sizeof(double));
  tpt = (double *) malloc(3*ndim*npart*sizeof(double));
  tmpfrom = (int *) malloc(2*beads*mxlist*2*sizeof(int));

  tint1 = (int *) malloc(2*beads*sizeof(int));
  ton1 = (int *) malloc(2*beads*sizeof(int));
  tthist1 = (int *) malloc(2*beads*tres*sizeof(int));
  tthist3 = (int *) malloc(2*beads*tres*sizeof(int));
  tthist5 = (int *) malloc(2*beads*tres*sizeof(int));
  tnflux1 = (int *) malloc(2*beads*2*beads*sizeof(int));
  tnflux2 = (int *) malloc(2*beads*2*beads*sizeof(int));

  nout = (int **) alloc2d(sizeof(int), 2, beads);
  nin = (int ****) alloc4d(sizeof(int), 2, beads, 2, beads);
  z = (opoint_t **) alloc2d(sizeof(opoint_t), 2, beads);
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

  /* initialize arrays */

  basin[0].x[0] = bmin;
  basin[1].x[0] = bmax;

  if (verb >= 3) {
    mainout = fopen(outname,"a");
    fprintf(mainout,"Initializing variables...\n");
    fclose(mainout);
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

  int dir, i, j, mcount;

/* get string from server */

  GetFromServer_dbl(&Coords_z,0,-1,tz1);

  mcount = 0;
  for (dir=0;dir<=1;dir++) {
    for (i=0;i<beads;i++) {
      for (j=0;j<nop;j++) {
	z[dir][i].x[j] = tz1[mcount];
	mcount++;
	if (verb >= 3) {
	  mainout = fopen(outname,"a");
	  fprintf(mainout,"%d %d %f\n",dir, i,tz1[mcount-1]);
	  fclose(mainout);
	}
      }
    }
  }
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

  return var*2*beads + dir*beads + bead;
}

int chkdir(int dir, int bead) {
  
  /* returns 1 if there is a basin crossing, 0 otherwise */
  int temp, tbas, odir, n1, n2, newreg, ind, ind1, ind2;

  temp = 0;
  tbas = bascheck();
		    
  if (tbas + dir == 1) {  /* either tbas = 0, dir = 1 or tbas = 1, dir = 0 */
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
  
  /* prim crossing */
  newreg = which(dir);
  if (newreg != bead) {
    ind = dir*beads + bead;
    
    if (verb >=3) {
      mainout = fopen(outname,"a");
      fprintf(mainout,"Ended traj: prim cross into %d\n",newreg);
      fclose(mainout);
    }

    addpoint(dir,newreg,coor,myweight,dir,bead);  /* send point to other flux list */

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
}
void wrtflux(int find, int cycle) {

  int j,k,pt,part,dim,mcount,ind,limit,listind;
  double extra, fac,tsum[2];

  /* write fluxes, divide flux weights */
  
  sprintf(fname,"%s%d%s%d%s","tflux", cycle/wrfrq/2, "_", find, ".dat");
  outFile = fopen(fname,"w");
  if (find==1) {
    j = 0;
  } else {
    j = 1;
  }
  for (k=0; k<=beads-1; k++) {
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
      
      fprintf(outFile,"%d %e\n", limit, tsum[listind]);
      
      for (pt=0; pt<limit; pt++) {
	ind = 2*mxlist*beads*j + 2*mxlist*k + 2*pt;
	GetFromServer_int(&Coords_from,listind,ind,tempi);
	GetFromServer_int(&Coords_from,listind,ind+1,tempi2);
	
	fprintf(outFile,"%d %d ", tempi[0], tempi2[0]);
	
	ind = mxlist*beads*j + mxlist*k + pt;
	
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
	    fprintf(outFile,"%e %e %e ", tpt[mcount],tpt[mcount+1],tpt[mcount+2]);
	    mcount += 3;
	  }
	}
	fprintf(outFile,"\n");
      }
    }
  }
  fclose(outFile);
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
void readtoglobal() {

  int dir, i, j, k, mcount;
  double t1;

  /* initialize string */

  if (verb >= 2) {
    printf("Automatically initializing string\n"); fflush(stdout);
  }
  if (strinit())
    gaexit(2) ;  
  
  wrpath(0);
  
  /* send string to server */
  
  mcount = 0;
  for (dir=0;dir<=1;dir++) {
    for (i=0;i<beads;i++) {
      for (j=0;j<nop;j++) {
	tz1[mcount] = z[dir][i].x[j];
	mcount++;
      }
    }
  }
  PutToServer_dbl(&Coords_z,0,-1,tz1);
  
  if (verb >= 2) {
    printf("Using flux configuration from %s\n", fxname); fflush(stdout);
  }
  fxread();
  
  /* initialize weights */
  
  if (strcmp(wname,"NOREAD") == 0) {
    if (verb >= 2) {
      printf("Automatically initializing weights\n"); fflush(stdout);
    }
    mcount = 0;
    for (j=0;j<2;j++) {
      for (k=0;k<beads;k++) {
	tweight1[mcount] = 1./(float)(2*beads);
	mcount++;
      }
    }
    PutToServer_dbl(&Coords_weight,0,-1,tweight1);
    if (wtalg == 4) {
      mcount = 0;
      for (j=0;j<2;j++) {
	for (k=0;k<beads;k++) {
	  tweight1[mcount] = 0.;
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
  
  /* read in nflux data */
  
  if (strcmp(nfname,"NOREAD") != 0) {
    if (verb >= 2) {
      printf("Reading nflux data from %s\n", nfname); fflush(stdout);
    }
    nfluxread();
    if (begscale) {
      if (strcmp(wtname,"NOREAD") ==0) {
	updtwts(0);
      }
      scalflux();
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
}


