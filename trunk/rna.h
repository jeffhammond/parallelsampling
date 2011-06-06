#define NDIM_c 3
#define NOP_c 1
#define NPART_c 262

const int ndim=NDIM_c, nop=NOP_c, npart=NPART_c;

const int tres = 200;  // resolution of the histograms
const double tmin[3] = {700.,800.,0.};  // minimum values of the 3 histograms
const double tmax[3] = {1000.,970.,200.};  // max values

const double grav = 0.01;
const double tstep = 0.01;
const double beta = 0.0002489;

typedef struct _xpoint_t
{
  double x[NDIM_c];
} xpoint_t;

typedef struct _point_t
{
  double x[NPART_c][NDIM_c];
  double v[NPART_c][NDIM_c];
} point_t;

typedef struct _opoint_t
{
  double x[NOP_c];
} opoint_t;

opoint_t basin[2];
double bwidth = 0.;
double bmin = 900.;
double bmax = 960.;

int nsec, ntert, **secind, **tertind;
double tstepsq, wcasig, wcacut;
double sigfacsq=25.;
double rsigfacsq=4.;

char crname[70]="../cutRNAin.dat";      // crystal structure file
char secname[70]="../cutRNAsec.con";    // secondary contacts
char tertname[70]="../cutRNAtert.con";  // tertiary contacts
char sname[30]="isolv";         // solvent file
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
double calpha[2], salpha[2], var, var2, solvstep, phan[NDIM_c], n_av, alpha;
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
int onpart=262;

void rddelta(), rdcryst(), rdsolv();
