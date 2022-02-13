./src/fragment.x
mpirun -np 10 ./src/runfrag-pesmpi.x >& log
./src/cal_energy_pes.x > energy.txt

