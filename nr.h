#include <math.h>
#include"nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void gaussj(double **a, int n,double **b,int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;

  indxc=ivector(0,n-1);
  indxr=ivector(0,n-1);
  ipiv=ivector(0,n-1);
  for (j=0;j<=n-1;j++) ipiv[j]=0;
  for (i=0;i<=n-1;i++) {
    big=0.0;
    for (j=0;j<=n-1;j++)
      if (ipiv[j] != 1)
	for (k=0;k<=n-1;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<=n-1;l++) SWAP(a[irow][l],a[icol][l])
      for (l=0;l<=m-1;l++) SWAP(b[irow][l],b[icol][l])
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<=n-1;l++) a[icol][l] *= pivinv;
    for (l=0;l<=m-1;l++) b[icol][l] *= pivinv;
    for (ll=0;ll<=n-1;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=0;l<=n-1;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=0;l<=m-1;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<=n-1;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,0,n-1);
  free_ivector(indxr,0,n-1);
  free_ivector(indxc,0,n-1);
}

#undef SWAP

