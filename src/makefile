all: fragment.x cal_energy_pes.x runfrag-pesmpi.x 
IF=ifort
MP=mpif90

fragment.o: fragment.f90
	$(IF) -c $<
fragment.x: fragment.o
	$(IF) fragment.o -o $@

cal_energy_pes.o: cal_energy_pes.f90
	$(IF) -c $<
cal_energy_pes.x: cal_energy_pes.o
	$(IF) cal_energy_pes.o -o $@

runfrag-pesmpi.o: runfrag-pesmpi.f90
	$(MP) -c $<
runfrag-pesmpi.x: runfrag-pesmpi.o
	$(MP) runfrag-pesmpi.o -o $@
clean:
	rm -f *.o fragment.x cal_energy_pes.x runfrag-pesmpi.x
