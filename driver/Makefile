# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

include ../user_build_options

OBJS = soil_temperature_output.o soil_temperature_driver.o

all:	$(OBJS)

soil_temperature_driver.o: soil_temperature_driver.f90

	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) -I../src -I./ $(*).f90 

soil_temperature_output.o: soil_temperature_output.f90

	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(*).f90 

#
# This command cleans up object (etc) files:
#

clean:
	rm -f *.o *.mod

#
# Dependencies:
#

soil_temperature_driver.o: soil_temperature_output.o ../src/soil_temperature_routines.o 

