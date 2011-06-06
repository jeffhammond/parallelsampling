void initsystem() {

  int i, j;
  
  tstepsq=tstep*tstep;
  
  for(i=0;i<npart;i++){
    for(j=0;j<3;j++){
      f1[i][j]=0.;
    }
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

  f1 = (double **) alloc2d(sizeof(double), npart, ndim);

  rdcryst();
  rddelta();
  rdsolv();


}


int strinit() {

  /* this function initializes the string */

  int dir, bead, op, i;

  /* string 1: forward & backward */

  for (i=0;i<=beads-1;i++) {
    z[0][i].x[0] = bmin + bwidth + (bmax-bmin-2*bwidth)*(i)/(float)(beads-1);
    z[1][i].x[0] = z[0][i].x[0];
  }

  return 0;
}

  
double gethist1(point_t coor) {

  ocoor = projx(coor);
  return ocoor.x[0];
}
double gethist2(point_t coor) {

  opoint_t tpt;

  tpt = projxsig(coor);
  return tpt.x[0];
}
double gethist3(point_t coor) {

  double temps;

  temps = 0.;
  for (i=0;i<ndim;i++) {
    temps += (coor.x[0][i]-coor.x[npart-1][i])*(coor.x[0][i]-coor.x[npart-1][i]);
  }
  temps = sqrt(temps);
  return temps;
}
/*-------------------------------------------------------------------------*/
opoint_t projxsig(point_t x) {
  
  /* this function returns the projection of coor onto the order parameter space, using sigfacsq */
  
  int i, j, k, ind1, ind2;
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
      sum = sum + tdub;
    }
  }

  temp.x[0] = sum;
  
  return temp;
}
/*-------------------------------------------------------------------------*/
opoint_t projx(point_t x) {
  
  /* this function returns the projection of coor onto the order parameter space, using rsigfacsq */
  
  int i, j, k, ind1, ind2;
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
      if (r2 < rsigfacsq*r02[ind1][ind2]) {
	tdub = 1.;
      } else if (r2 < 4.*rsigfacsq*r02[ind1][ind2]) {
	tdub = mypow4(rsigfacsq*r02[ind1][ind2]/r2);
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
      if (r2 < rsigfacsq*r02[ind1][ind2]) {
	tdub = 1.;
      } else if (r2 < 4.*rsigfacsq*r02[ind1][ind2]) {
	tdub = mypow4(rsigfacsq*r02[ind1][ind2]/r2);
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

  int i;
  point_t oldx, rt;
  double psi[2],phi[2],fx,fy,temp1,temp2,dvdx,dvdy;
  void dostep(), stream(), dorotate();
  
  dostep();
  if (!dynerr) {
    if (myclock%strmfrq == 0) {
      stream();
      dorotate();
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
      dynerr = 1;
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
    if (coor.x[p1][0] < wcacut) { // add WCA repulsion if bead is close to x=0 wall
      /* only comes into play for no-flow case (even then, extremely rarely) */
      pre = mypow6(wcasig/coor.x[p1][0]);
      f.x[p1][0] += 4.*ep1_6/coor.x[p1][0]*(2.*pre*pre - pre);
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
	dynerr = 1;
    }

    if (!dynerr) {

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

