/*************************************************************************/
/*************************** aryeh's functions ***************************/
/*************************************************************************/

/*
 *  Function:  alloc2d
 *  ------------------
 *  Allocate a two-dimensional array.
 */

void **alloc2d(int varsize, int n, int p) {
    int k ;
    void **a ;

    if ((a = (void **) malloc(n*sizeof(void *))) == NULL) { 
        printf ("Memory Error in alloc2d.\n") ;
        exit(1) ;
    }

    for (k = 0 ; k < n ; k++) {
        if ((a[k] = (void *) malloc(p*varsize)) == NULL) { 
            printf ("Memory Error in alloc2d.\n") ;
            exit(1) ;
        }
    }
    return a ;
}



void ***alloc3d(int varsize, int n, int p, int q) {
    int k,j ;
    void ***a ;

    if ((a = (void ***) malloc(n*sizeof(void **))) == NULL) { 
        printf ("Memory Error in alloc2d.\n") ;
        exit(1) ;
    }

    for (k = 0 ; k < n ; k++) {
        if ((a[k] = (void **) malloc(p*sizeof(void *))) == NULL) { 
            printf ("Memory Error in alloc2d.\n") ;
            exit(1) ;
        }
		for( j = 0 ; j < p ; j++){
			if ((a[k][j] = (void *) malloc(q*varsize)) == NULL) { 
				printf ("Memory Error in alloc2d.\n") ;
				exit(1) ;
			}
		}

    }
    return a ;
}

void ****alloc4d(int varsize, int n, int p, int q, int r) {
    int k,j,i ;
    void ****a ;

    if ((a = (void ****) malloc(n*sizeof(void ***))) == NULL) { 
        printf ("Memory Error in alloc2d.\n") ;
        exit(1) ;
    }

    for (k = 0 ; k < n ; k++) {
        if ((a[k] = (void ***) malloc(p*sizeof(void **))) == NULL) { 
            printf ("Memory Error in alloc2d.\n") ;
            exit(1) ;
        }
	for( j = 0 ; j < p ; j++){
		if ((a[k][j] = (void **) malloc(q*sizeof(void *))) == NULL) { 
			printf ("Memory Error in alloc2d.\n") ;
			exit(1) ;
		}
		for( i = 0 ; i < q; i++){	
			if( (a [k][j][i] = (void *) malloc(r*varsize))
				==NULL){
				printf("Memory error in alloc4d.\n");
				exit(1);
			}
		}
	}

    }
    return a;
}

void *****alloc5d(int varsize, int n, int p, int q, int r, int s) {
    int k,j,i,m;
    void *****a ;

    if ((a = (void *****) malloc(n*sizeof(void ****))) == NULL) { 
        printf ("Memory Error in alloc5d.\n") ;
        exit(1) ;
    }

    for (k = 0 ; k < n ; k++) {
        if ((a[k] = (void ****) malloc(p*sizeof(void ***))) == NULL) { 
            printf ("Memory Error in alloc5d.\n") ;
            exit(1) ;
        }
		for( j = 0 ; j < p ; j++){
			if ((a[k][j] = (void ***) malloc(q*sizeof(void **))) == NULL) { 
				printf ("Memory Error in alloc5d.\n") ;
				exit(1) ;
			}
			for( i = 0 ; i < q; i++){	
				if( (a [k][j][i] = (void **) malloc(r*sizeof(void *)))
					==NULL){
					printf("Memory error in alloc5d.\n");
					exit(1);
				}
				for( m = 0 ; m < r; m++){
					if( (a[k][j][i][m] = (void *) malloc(s*varsize)) == NULL){
						printf("Memory error in alloc5d.\n");
						exit(1);
					}
				}
			}
		}
			
	}

    return a;
}

void *******alloc7d(int varsize, int n, int p, int q, int r, int s, int t, int u) {
  int k,j,i,m,g,h;
  void *******a ;

  if ((a = (void *******) malloc(n*sizeof(void ******))) == NULL) { 
    printf ("Memory Error in alloc7d.\n") ;
    exit(1) ;
  }

  for (k = 0 ; k < n ; k++) {
    if ((a[k] = (void ******) malloc(p*sizeof(void *****))) == NULL) { 
      printf ("Memory Error in alloc7d.\n") ;
      exit(1) ;
    }
    for( j = 0 ; j < p ; j++){
      if ((a[k][j] = (void *****) malloc(q*sizeof(void ****))) == NULL) { 
	printf ("Memory Error in alloc7d.\n") ;
	exit(1) ;
      }
      for( i = 0 ; i < q; i++){	
	if( (a [k][j][i] = (void ****) malloc(r*sizeof(void ***))) ==NULL){
	  printf("Memory error in alloc7d.\n");
	  exit(1);
	}
	for( g = 0 ; g < r ; g++){
	  if ((a[k][j][i][g] = (void ***) malloc(s*sizeof(void **))) == NULL) { 
	    printf ("Memory Error in alloc7d.\n") ;
	    exit(1) ;
	  }
	  for( h = 0 ; h < s; h++){	
	    if( (a [k][j][i][g][h] = (void **) malloc(t*sizeof(void *))) ==NULL){
	      printf("Memory error in alloc7d.\n");
	      exit(1);
	    }
	    for( m = 0 ; m < u; m++){
	      if( (a[k][j][i][g][h][m] = (void *) malloc(u*varsize)) == NULL){
		printf("Memory error in alloc7d.\n");
		exit(1);
	      }
	    }
	  }
	}
      }
    }
  }

  return a;
}


/************************************************************************
*  Function: gssran
*  --------------------------------------------------------------------
*  pick a random number from a normal distribution
*  (average zero, unit variance)                    
************************************************************************/
/*
double gssran(long *idum)
{
  double fac, rsq, v1, v2;
  static int iset = 0;
  static double gset;

  if (iset == 0){
    do{
      v1 = 2.0 * ran1(idum) - 1.0;
      v2 = 2.0 * ran1(idum) - 1.0;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1.0) || (rsq ==0.0));

    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  }

  else{
    iset =0;
    return gset;
  }
}
  
*/
