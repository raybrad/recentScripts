# Usage:  ./mpiLodestar.sh [nNode nCore]
OMP='-env OMP_PROC_BIND true -env OMP_NUM_THREADS' 
MKL='-env MKL_NUM_THREADS' 
nNode=1
nCore=1
if [ $# -gt 0 ]; then
	nNode=$1
fi
if [ $# -gt 1 ]; then
	nCore=$2
fi
echo "mpirun -machinefile hostfile -np $nNode $OMP $nCore $MKL $nCore ./lodestar_mpi in.td"
mpirun -machinefile hostfile -np $nNode $OMP $nCore $MKL $nCore ./lodestar_mpi in.td 2>err
#DEBUG='-verbose -genv I_MPI_HYDRA_DEBUG=1 -check_mpi'
#DEBUG='-verbose -genv I_MPI_HYDRA_DEBUG=1 -check_mpi I_MPI_HYDRA_BRANCH_COUNT=128 I_MPI_DAPL_UD=1'
#mpirun -machinefile hostfile $DEBUG -np $nNode $OMP $nCore $MKL $nCore ./lodestar_mpi in.td 2>err
