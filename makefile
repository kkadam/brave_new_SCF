.SUFFIXES: .F90 .c .o
#OFILE_DIR= obj

F90FILES= main.F90 binary_scf.F90 binary_initialize.F90 binary_output.F90 binary_sum.F90 compute_virial_error.F90 \
potential_solver.F90 bessel.F90 helmadi.F90 tridagr.F90 tridagz.F90 realft.F90 \
output.F90 \
potsetup.F90 tm.F90 sm.F90 elle.F90 ellf.F90 gammln.F90 rd.F90 rf.F90 setup.F90

OFILES= $(F90FILES:.F90=.o) $(F90_SAFE_FILES:.F90=.o)

hydro:$(OFILES)
# normal linking step
#	ifort -traceback -O3 -o scf -fpe0 $(OFILES)
#	ifort -O0 -o scf $(OFILES)
#       ifort -O0 -fPIC -shared-intel -o scf $(OFILES)  #works low res
	ifort -O0 -mcmodel=medium -shared-intel -o scf $(OFILES)

$(OFILES): runscf.h

.F90.o: runscf.h
# normal compilation
#	ifort -traceback -c -O3 -r8 -fpe0 $<
#	ifort -c -O0 -r8 $<
#       ifort -c -O0 -r8 -fPIC -shared-intel $<         #works low res
	ifort -c -O0 -r8 -mcmodel=medium -shared-intel $<

clean:
	/bin/rm -f *.o *.lst scf
