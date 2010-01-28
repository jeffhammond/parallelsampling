#include "myCoordServer.h"

int CreateCoordServer_int(ServedCoords* NewCoords,
                         int* MaxNum,
			 char* Name)
{

    int me;
    int status;

    MPI_Comm_rank(MPI_COMM_WORLD,&me);

    NewCoords->MaxNum  = *MaxNum;
    NewCoords->pg = GA_Pgroup_get_world();
    NewCoords->ndim = 2;

    NewCoords->dims[0] = 2;
    NewCoords->dims[1] = NewCoords->MaxNum;

    /*    if(0 == me){
        fprintf(stdout,"TotalSize = %d\n",TotalSize);
        fflush(stdout);
	}*/

    NewCoords->chunk[0] = 1024*3;   // localize coords unless huge
    NewCoords->chunk[1] = 1024*3;   // localize coords unless huge

    NewCoords->ga = GA_Create_handle();
    GA_Set_array_name(NewCoords->ga,Name);
    GA_Set_data(NewCoords->ga,NewCoords->ndim,NewCoords->dims,MT_INT);
    GA_Set_chunk(NewCoords->ga,NewCoords->chunk);
    GA_Set_pgroup(NewCoords->ga,NewCoords->pg);
    status = GA_Allocate(NewCoords->ga);
    if(1 != status){
        fprintf(stderr,"proc %d: GA_Allocate = %d\n",me,status);
        fflush(stderr);
    };
    GA_Zero(NewCoords->ga);

    return(status);
}


int CreateCoordServer_dbl(ServedCoords* NewCoords,
                         int* MaxNum,
			 char* Name)
{

    int me;
    int status;

    MPI_Comm_rank(MPI_COMM_WORLD,&me);

    NewCoords->MaxNum  = *MaxNum;
    NewCoords->pg = GA_Pgroup_get_world();
    NewCoords->ndim = 2;

    NewCoords->dims[0] = 2;
    NewCoords->dims[1] = NewCoords->MaxNum;

    /*    if(0 == me){
        fprintf(stdout,"TotalSize = %d\n",TotalSize);
        fflush(stdout);
	}*/

    NewCoords->chunk[0] = 1024*3;   // localize coords unless huge
    NewCoords->chunk[1] = 1024*3;   // localize coords unless huge

    NewCoords->ga = GA_Create_handle();
    GA_Set_array_name(NewCoords->ga,Name);
    GA_Set_data(NewCoords->ga,NewCoords->ndim,NewCoords->dims,MT_DBL);
    GA_Set_chunk(NewCoords->ga,NewCoords->chunk);
    GA_Set_pgroup(NewCoords->ga,NewCoords->pg);
    status = GA_Allocate(NewCoords->ga);
    if(1 != status){
        fprintf(stderr,"proc %d: GA_Allocate = %d\n",me,status);
        fflush(stderr);
    };
    GA_Zero(NewCoords->ga);

    return(status);
}

int DestroyCoordServer(ServedCoords* Coords)
{
    GA_Destroy(Coords->ga);

    return(0);
}

int GetFromServer_int(ServedCoords* Coords, int ind1, int ind2, int* Coord)
{
    int lo[2];
    int hi[2];
    int ld[1], i;

    if (ind2 == -1) { // get all
      lo[0] = ind1;
      lo[1] = 0;
      hi[0] = ind1;
      hi[1] = Coords->MaxNum-1;
    } else {
      lo[0] = ind1;
      lo[1] = ind2;
      hi[0] = ind1;
      hi[1] = ind2;
    }

    ld[0] = 0;

    NGA_Get(Coords->ga,lo,hi,Coord,ld);
    return(0);
}

int GetFromServer_dbl(ServedCoords* Coords, int ind1, int ind2, double* Coord)
{
    int lo[2];
    int hi[2];
    int ld[1];

    if (ind2 == -1) { // get all
      lo[0] = ind1;
      lo[1] = 0;
      hi[0] = ind1;
      hi[1] = Coords->MaxNum-1;
    } else {
      lo[0] = ind1;
      lo[1] = ind2;
      hi[0] = ind1;
      hi[1] = ind2;
    }

    ld[0] = 0;

    NGA_Get(Coords->ga,lo,hi,Coord,ld);
    return(0);
}

int PutToServer_int(ServedCoords* Coords, int ind1, int ind2, int* Coord)
{
    int lo[2];
    int hi[2];
    int ld[1];

    if (ind2 == -1) { //  put all
      lo[0] = ind1;
      lo[1] = 0;
      hi[0] = ind1;
      hi[1] = Coords->MaxNum-1;
    } else {
      lo[0] = ind1;
      lo[1] = ind2;
      hi[0] = ind1;
      hi[1] = ind2;
    }

    ld[0] = 0;

    NGA_Put(Coords->ga,lo,hi,Coord,ld);
    return(0);
}

int PutToServer_dbl(ServedCoords* Coords, int ind1, int ind2, double* Coord)
{
    int lo[2];
    int hi[2];
    int ld[1];

    if (ind2 == -1) { //  put all
      lo[0] = ind1;
      lo[1] = 0;
      hi[0] = ind1;
      hi[1] = Coords->MaxNum-1;
    } else {
      lo[0] = ind1;
      lo[1] = ind2;
      hi[0] = ind1;
      hi[1] = ind2;
    }

    ld[0] = 0;

    NGA_Put(Coords->ga,lo,hi,Coord,ld);
    return(0);
}
