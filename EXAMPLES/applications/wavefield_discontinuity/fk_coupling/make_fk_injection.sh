MPIFC=mpif90
FLAGS='-g -xHost -O3'

${MPIFC} ${FLAGS} -o compute_fk_injection_field coupling_fk.f90 utils.f90 compute_fk_injection_field.f90

${MPIFC} ${FLAGS} -o compute_fk_receiver coupling_fk.f90 utils.f90 compute_fk_receiver.f90
