#makefile to build a kdtree2 library


FFLAGS  = $(FCFLAGS) -O3
AR      = ar


.SUFFIXES:      .F90 .o 

all:	kdtree2.a

kdtree2.a : kdtree2.o
	$(AR) cr kdtree2.a kdtree2.o

kdtree2.o : 

%.o :: %.F90
	${FC} ${FFLAGS} -c $<

clean:
	/bin/rm -f *.o *.mod 

superclean:
	/bin/rm -f *.o *.mod 
