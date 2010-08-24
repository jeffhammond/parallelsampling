#define NDIM_c 3
#define NOP_c 1
#define PI_c 3.14159265
#define NPART_c 64
#define NEMERG_c 400

const int phase = 2;
int wtalg = 3;       /* 1 = local ; 2 = global */
const char flname[30] = "NOREAD";
const char wname[30] = "NOREAD";
const char fxname[30] = "NOREAD";
const char pname[30] = "NOREAD";
const char tname[30] = "NOREAD";

long seed = 1;
const int beads = 15;
const int wtupdt = 20000;
const int T = 1;         /* measured in number of cycles */
const int every = 20000;
const int stkfrq = 100;
const int wrfrq = 20000;
const int pzero = 500000;
const int xyzfrq = 200;
const int globfrq = 5;
const int chkfrq = 5;  
const int tres = 100;
const int nemerg = NEMERG_c;
const double tmin = 40.;
const double tmax = 200.;
const double divfac = 0.7;
const double bfrac = 1.;

const double s = 2.e-2;
const double frac = 0.1;
const double wfrac = 1.0;
const double kappa = 0.1;

const double grav = 0.001;
const double tstep = 0.01;
const double beta = 0.0002489;

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

typedef struct _fluxpoint_t
{
  double x[NPART_c][NDIM_c];
  double v[NPART_c][NDIM_c];
  double f1[NPART_c][NDIM_c];
  double wt;
  int from[2];
  int clock;
  int mytim;
} fluxpoint_t ;

empoint_t emerg[NEMERG_c];

opoint_t basin[2];
double bwidth = 0.;
double bmin = 87.;
double bmax = 115.;
