#ifndef DRIVER_H
#define DRIVER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "macdecls.h"
#include "armci.h"
#include "sndrcv.h"
#include "ga.h"
#include "mpi.h"


typedef struct servedcoords
{
    int MaxNum;

    int ndim;
    int dims[2];
    int chunk[2];

    int ga; // Global Array handle
    int pg; // Processor Group handle

} ServedCoords;

int CreateCoordServer_int(ServedCoords* NewCoords,
			  int* MaxNum,
			  char* Name);
int CreateCoordServer_dbl(ServedCoords* NewCoords,
			  int* MaxNum,
			  char* Name);

int DestroyCoordServer(ServedCoords* Coords);

int GetFromServer_int(ServedCoords* Coords, int ind1, int ind2, int* Coord);
int GetFromServer_dbl(ServedCoords* Coords, int ind1, int ind2, double* Coord);
int PutToServer_int(ServedCoords* Coords, int ind1, int ind2, int* Coord);
int PutToServer_dbl(ServedCoords* Coords, int ind1, int ind2, double* Coord);

#endif // DRIVER_H
