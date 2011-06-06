#define MXLIST_c 250
#define PI_c 3.14159265

const int phase = 2;
int wtalg = 3;       /* 1 = local ; 2 = global ; 3 = VdE ; 4 = ww-avg */
const char flname[30] = "NOREAD";
const char wname[30] = "iwt";
const char fxname[30] = "iflux";
const char pname[30] = "NOREAD";
const char nfname[30] = "influx";
const char tname[30] = "NOREAD";
const int begscale = 0;
const int clrtnfx = 1;
const int startcycle = 3001; 

long seed = 1;
const int beads = 40;
const int wtupdt = 600;     /* measured in cycles */
const int T = 6000;         /* measured in number of cycles */
const int every = 2000;  /* steps */
const int stkfrq = 100;   /* steps */
const int wrfrq = 50;      /* cycles */
const int pzero = 5000;      /* cycles */
const int xyzfrq = 1000;  /* steps */
const int globfrq = 5;   /* steps */
const int chkfrq = 5;    /* steps */
const double bfrac = 1.;
const double spltfac = 0.9;

const double Pi = PI_c;
const double twoPi = 2*PI_c;

typedef struct _fluxlist_t
{
  int full;       /* 0 = not full; 1 = full */
  int nlist;      /* how many points in the list are used up */
  double twlist;  /* total weight of flux list */
  point_t pts[MXLIST_c];  /* the list of flux points for each region */
  double wlist[MXLIST_c];
  int from[MXLIST_c][2];
  int clock[MXLIST_c];
} fluxlist_t ;

