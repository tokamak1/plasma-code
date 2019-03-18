ALL : B.out
FPPFLAGS = -c -g -traceback
FLAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
OBJ = module.o construct_spline.o setup.o orbit.o

FC  = ifort
FLINKER  = ifort

module.o : module.F90
	${FC} ${FPPFLAGS} module.F90 -O2
construct_spline.o : construct_spline.F90
	${FC} ${FPPFLAGS} construct_spline.F90 -O2
setup.o : setup.F90
	${FC} ${FPPFLAGS} setup.F90 -O2
orbit.o : orbit.F90
	${FC} ${FPPFLAGS} $^ -O2

B.out : ${OBJ} 
	${FLINKER} ${OBJ} -o B.out ${FLAPACK} 

.PHONY : clean
clean : 
	-rm -r *.mod *.o *.out
	
