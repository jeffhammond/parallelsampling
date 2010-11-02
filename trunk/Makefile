   RM = rm
   RMFLAGS = -f

   AR = ar
   ARFLAGS = -r

   #
   # How to build Open MPI 1.3.3 (http://www.open-mpi.org/software/ompi/v1.3/downloads/openmpi-1.3.3.tar.bz2)
   #
   # ./configure --prefix=/software/open-mpi/gnu-build --with-devel-headers --disable-shared --enable-static --enable-mpi-threads
   # make; make install
   #
   #MPI_PREFIX=/bgsys/drivers/ppcfloor/comm
   #MPI_INC=-I$(MPI_PREFIX)/include
   #MPI_LIB=-L$(MPI_PREFIX)/lib -lmpi -lopen-rte -lopen-pal -ldl -lutil

   # How to build GA 4.2 (http://www.emsl.pnl.gov/docs/global/.download/ga-4-2.tgz)
   #
   # export MPI_PREFIX=(whatever you use above)
   # make TARGET=LINUX64 LARGE_FILES=y USE_MPI=y MPI_LIB=/home/adickson/software/open-mpi/gnu-build/lib /home/adickson/software/open-mpi/gnu-build/include LIBMPI="-lmpi -lopen-rte -lopen-pal -ldl -lutil"
   #
   #GA_PREFIX=/home/adickson/software/ga-4-2
   #GA_INC=-I$(GA_PREFIX)/include
   #GA_LIB=-L$(GA_PREFIX)/lib/LINUX -lglobal -lma -larmci -ltcgmsg-mpi -llinalg

   #EXTRAS=-lgfortran -lm -lpthread

   GA_PREFIX=/home/projects/nwchem/ga-alex
   GA_INC=-I$(GA_PREFIX)/include
   GA_LIB=-L$(GA_PREFIX)/lib/BGP -lglobal -lma -larmci -ltcgmsg-mpi -llinalg -lgfortran

   #EXTRAS=-lgfortran -lm -lpthread

##   OMP_FLAGS=-qsmp=omp

   LIB=$(GA_LIB) $(MPI_LIB) $(EXTRAS)
   INC=$(GA_INC) $(MPI_INC)

   CC=tau_cc.sh
   CFLAGS=-O0 $(INC) $(OMP_FLAGS)

   LD=tau_f90.sh
   LDFLAGS=-O0 $(LIB) $(OMP_FLAGS)

   MYPATH=/home/adickson/parallelsampling

#############################################
#
#               End of Targets
#
#############################################

all: flow_ga.x

refresh: realclean all

nrutil.o: $(MYPATH)/nrutil.c $(MYPATH)/nrutil.h
	$(CC) -c $(MYPATH)/nrutil.c $(CFLAGS) -o nrutil.o

flow_ga.x: flow_ga.o $(MYPATH)/myCoordServer.o nrutil.o
	$(LD) flow_ga.o $(MYPATH)/myCoordServer.o nrutil.o $(LDFLAGS) -o flow_ga.x

flow_ga.o: flow_ww2.c neusglob.h
	$(CC) -c flow_ww2.c $(CFLAGS) -o flow_ga.o

myCoordServer.o: $(MYPATH)/myCoordServer.c $(MYPATH)/myCoordServer.h
	$(CC) -c $(MYPATH)/myCoordServer.c $(CFLAGS) -o myCoordServer.o

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.x *.a



