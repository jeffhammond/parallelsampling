#define NDIM_c 3
#define NOP_c 1
#define MXLIST_c 500
#define PI_c 3.14159265
#define NPART_c 262 
#define NEMERG_c 385

const int phase = 2;
const int wtalg = 1;       /* 1 = local ; 2 = global */
const char flname[30] = "NOREAD";
const char wname[30] = "NOREAD";
const char fxname[30] = "NOREAD";

const int beads = 30;
const int wtupdt = 100000;
const int T = 200000;         /* measured in number of steps */
const int every = 200000;
const int stkfrq = 1000;
const int wrfrq = 200000;
const int xyzfrq = 1000;
const int chkfrq = 1000;  
const int globfrq = 100;  
const int tres = 100;
const int nemerg = NEMERG_c;
const double tmin = 100.;
const double tmax = 400.;

const double s = 5.e-2;
const double frac = 0.1;
const double wfrac = 1.0;
const double kappa = 0.1;

const double grav = 0.06;
const double tstep = 5.e-4;
const double beta = 10./6.;

const double Pi = PI_c;
const double twoPi = 2*PI_c;

typedef struct _xpoint_t
{
  double x[NDIM_c];
} xpoint_t;

typedef struct _point_t
{
  double x[NPART_c][NDIM_c];
  double v[NPART_c][NDIM_c];
} point_t;

typedef struct _empoint_t
{
  double op;
  double x[NPART_c][NDIM_c];
  double v[NPART_c][NDIM_c];
} empoint_t;

typedef struct _opoint_t
{
  double x[NOP_c];
} opoint_t;

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

empoint_t emerg[NEMERG_c];

opoint_t basin[2];
double bwidth = 10.;
double bmin = 130.;
double bmax = 245.;
