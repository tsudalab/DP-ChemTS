# DP-ChemTS
A distributed framework based on Monte Carlo tree search for accelerating molecular discovery. DP-ChemTS is the distributed parallel version of our publish ChemTS library. DP-ChemTS implemented three distributed parallel methods: Distributed leaf parallel, distributed tree parallel and distributed tree parallel with virtual loss(is used to reduce the search overhead). 

# Requirements
1.mpi4py

2.MPICH

3.Cluster machines

4. Other requirements for ChemTS.

# How to use DP-ChemTS?
1. Distributed leaf parallel.

cd leaf_parallel_test/simulation1/4core

run: mpiexec -n 4 python mpi_thread_leaf_parallel.py


2. Distributed tree parallel.

cd tree_parallel_test/simulation1/4core

run: mpiexec -n 4 python mpi_thread_chemts_tree_parallel.py


3. Distributed tree parallel with virtual loss.

cd virtual_loss_test/simulation1/4core

run: mpiexec -n 4 python mpi_thread_chemts_tree_vl.py
