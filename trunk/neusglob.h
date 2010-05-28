#define NDIM_c 3
#define NOP_c 1
#define MXLIST_c 500
#define PI_c 3.14159265
#define NPART_c 262 
#define NEMERG_c 790

const int phase = 2;
const int wtalg = 2;       /* 1 = local ; 2 = global */
const char flname[30] = "NOREAD";
const char wname[30] = "NOREAD";
const char fxname[30] = "NOREAD";
const char pname[30] = "NOREAD";

long seed = 1;
const int beads = 30;
const int wtupdt = 2000;
const int T = 1;         /* measured in number of cycles */
const int every = 2000;
const int stkfrq = 100;
const int wrfrq = 2000;
const int xyzfrq = 2000;
const int globfrq = 10;
const int chkfrq = 10;  
const int tres = 1000;
const int nemerg = NEMERG_c;
const double tmin = 50.;
const double tmax = 300.;

const double s = 2.e-3;
const double frac = 0.1;
const double wfrac = 1.0;
const double kappa = 0.1;

const double grav = 0.08;
const double tstep = 0.1;
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
double bwidth = 1.;
double bmin = 140.;
double bmax = 245.;
