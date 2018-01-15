#!/bin/sh
#QSUB2 queue qM
#QSUB2 core 12
#QSUB2 mpi 12
#QSUB2 smp 1
#PBS -N test_300

. /etc/profile.d/modules.sh
module purge
module load g16/A03
module load anaconda/2
module list

export GAUSS_SCRDIR=$HOME/scr
mkdir -p $GAUSS_SCRDIR

cd $PBS_O_WORKDIR
uname -a

#export INPUTFILE=test
#export LINDA=`cat $PBS_NODEFILE | uniq | tr '\n' "," | sed 's/,$//' | sed 's/:$//'`
#export NPPROC=`wc -l $PBS_NODEFILE|awk '{print $1}'`
#export GAUSS_LFLAGS="-nodefile $PBS_NODEFILE"
#echo "LINDA=",$LINDA
#echo "INPUTFILE=",$INPUTFILE
#cat $INPUTFILE.com | sed "s/LINDA/$LINDA/" >temp.com1
#cat temp.com1 | sed "s/NPPROC/$NPPROC/" >temp.com
#g16 < temp.com > test.log 2>&1
#g16 < CheckMolopt.com > CheckMolopt.log
mpiexec -n 4 python mpi_thread_chemts_tree_vl.py  > log
