'''
	para_agnspec_grid.py

Testing mpi4py

Dependencies:
	mpi4py installed, agnsepc fortran routines compiled

Usage:
	run with
	mpirun -n 4 python para_agnspec_grid.py
	where 4 is number of cores

Documentation:
	See http://mpi4py.scipy.org/docs/usrman/index.html
	Note that you can also communicate information between threads 
	using MPI Pack and Broadcast commands.

	See agnspec/Agnspec.guide for guide on agnspec

	README.md on how to use my scripts

'''

import os, sys
import time
from mpi4py import MPI
from agnspec import run_agnspec, mdot_from_edd
import numpy as np


nproc = MPI.COMM_WORLD.Get_size()   	# number of processes

my_rank = MPI.COMM_WORLD.Get_rank()   	# The number/rank of this process

my_node = MPI.Get_processor_name()    	# Node where this MPI process runs


# prepare arrays to run grid over
FOLDER = "/Users/jmatthews/Documents/runs/QSO_clumped/disk_spec/agnspec_grids/outputs"

spin = [0,0.998]		# spin parameter
masses = 10.0 ** np.arange(7.0,10.0,0.2)  # black hole masses
edd_frac = np.arange(0,1,0.025)			  # eddington fractions
edd_frac[0] = 0.01
incs = np.arange(10,100,10)				  # inclinations in degrees
incs[-1] = 89.0



m_list = []
mdot_list = []
inc_list = []
spin_list = []
fnames = []

for ispin in range(len(spin)):
		for im in range(len(masses)):
			for iedd in range(len(edd_frac)):

				mdot = mdot_from_edd(edd_frac[iedd], masses[im])

				for ii in range(len(incs)):
					logm = np.log10(masses[im])
					if logm - int(logm) < 0.2:
						m2 = 0
					else: m2 = 5


					mbh_string = "%ip%i" % (int(logm), m2)
					fnames.append("sp%i_inc%2i_mbh%s_edd%itminus3" % (int(spin[ispin]+0.5),incs[ii], mbh_string, edd_frac[iedd]/0.001))

					
					m_list.append(masses[im])
					mdot_list.append(mdot)
					inc_list.append(incs[ii])
					spin_list.append(spin[ispin])





N_MODELS_TOTAL = len(fnames)		# total number of models to run	

print "Total models %i" %  N_MODELS_TOTAL


n_models = N_MODELS_TOTAL / nproc		# number of models for each thread

remainder = N_MODELS_TOTAL - ( n_models * nproc )	# the remainder. e.g. your number of models may 


# little trick to spread remainder out among threads
# if say you had 19 total models, and 4 threads
# then n_models = 4, and you have 3 remainder
# this little loop would distribute these three 
if remainder < my_rank + 1:
	my_extra = 0
	extra_below = remainder
else:
	my_extra = 1
	extra_below = my_rank

# where to start and end your loops for each thread
my_nmin = (my_rank * n_models) + extra_below
my_nmax = my_nmin + n_models + my_extra

# total number you actually do
ndo = my_nmax - my_nmin


print "This is thread %i calculating models %i to %i" % (my_rank, my_nmin, my_nmax)	


# set barrier so print output doesn't look muddled
# just waits for other thread
MPI.COMM_WORLD.Barrier()

# start a timer for each thread
time_init = time.time()



# now we can actually do the loop
# we'll open a file for each thread
# then for each file write the number of the loop

for i in range( my_nmin, my_nmax): 

	'''this is where you would do your stuff!!!'''

	run_agnspec(fnames[i], mass=m_list[i], mdot=mdot_list[i], inc=inc_list[i], angm=spin_list[i])

	# print statement which keeps track for the 0th thread only
	if my_rank == 0:
		print "SYS: Run %i of %i complete. Elapsed time %8.4es" % (n, nruns, t2 - tinit)
		print "SYS: Projected finish = %.2f hours" % ( (1.0*nruns) / float(n) * (t2 - tinit) / n / 3600.0)



# get the time taken
time2 = time.time()	
time_tot = time2 - time_init


# set barrier so print output doesn't look muddled
if my_rank == 0: print 'Waiting for threads to finish...'

# another barrier, wait for all to finish
MPI.COMM_WORLD.Barrier()


print "Thread %i took %.2f seconds to calculate %i models" % (my_rank, time_tot, ndo)

# always call this when finishing up
MPI.Finalize()


