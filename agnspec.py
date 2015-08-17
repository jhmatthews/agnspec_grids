#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import pretty
from pretty import *

''' 
subroutines and functions relating to running agnspec grids
'''

D = 100.0 * PARSEC

def Ledd (m):

	'''
	calculates eddington luminosity for a solar mass m
	
	Args:
		m	mass in solar masses
	
	returns:
		Eddington luminosity in ergs s^-1, float
	'''
	
	m *= MSOL
	
	consts = (4.0 * PI * G * C * MPROT ) / THOMPSON
	
	L = consts * m
	
	return L

def run_agnspec(fname, mass=1e9, mdot=5, angm=0, alpha=0.01, inc=80):

	mu = np.cos(inc / 180.0 * np.pi)

	print inc, mu

	param_string = "mass=%8.4e,mdot=%8.4e,angm=%.6f,alpha=%.6f,mu=%8.4e,savename=\"%s\"" % (mass, mdot, angm, alpha, mu, fname)

	print param_string
	os.system("idl -e 'agnspec/agnspec,%s' > idl.out" % param_string)

	return 0

def read_agnspec(fname):

	nu, lnu = np.loadtxt(fname, unpack=True, usecols=(0,1))

	wave = (C / nu) / ANGSTROM

	llambda = nu * lnu / wave

	flambda = llambda / 4.0 / PI / (D*D)

	return wave, flambda

def mdot_from_edd ( edd, m , eta = 0.1):

	''' 
	calculates an accretion rate from an eddington fraction
	Args:
		edd_frac		eddington fraction
		m			mass of central object in solar masses
		eta			accretion efficiency, set to 1 (0.1 more realistic)
	
	returns:
		mdot in solar masses / yr
	'''
	
	L = Ledd (m)		# eddington luminosity
	
	mdot = edd * L / ( (C ** 2) )
	
	mdot *= 1.0 / eta
	
	mdot = mdot * ( YR ) / MSOL	# normalise units
	
	return mdot

def run_grid(masses, edd_frac, spin, incs, folder):

	import time

	nruns = len(spin) * len(masses) * len(edd_frac) * len(incs)

	tinit = time.time()

	n = 0

	for ispin in range(len(spin)):
		for im in range(len(masses)):
			for iedd in range(len(edd_frac)):

				mdot = mdot_from_edd(edd_frac[iedd], masses[im])

				# figure(figsize=(8,8))

				# pretty.set_pretty()

				for ii in range(len(incs)):
					logm = np.log10(masses[im])
					if logm - int(logm) < 0.2:
						m2 = 0
					else: m2 = 5


					mbh_string = "%ip%i" % (int(logm), m2)
					fname = "sp%i_inc%2i_mbh%s_edd%itminus3" % (int(spin[ispin]+0.5),incs[ii], mbh_string, edd_frac[iedd]/0.001)

					run_agnspec(fname, mass=masses[im], mdot=mdot, inc=incs[ii], angm=spin[ispin])

					os.system("cp %s %s" % (fname, folder))
					t2 = time.time()
					n+= 1

					print "SYS: Run %i of %i complete. Elapsed time %8.4es" % (n, nruns, t2 - tinit)
					print "SYS: Projected finish = %.2f hours" % ( (1.0*nruns) / float(n) * (t2 - tinit) / n / 3600.0)



def set_set2():
	import brewer2mpl
	color = "Set1"
	set2 = brewer2mpl.get_map(color, 'qualitative', 8).mpl_colors
	g = set2[2]
	b = set2[1]
	r = set2[0]
	set2[0] = b
	set2[1] = g
	set2[2] = r

	return set2





# Next lines permit one to run the routine from the command line with various options -- see docstring
if __name__ == "__main__":

	#run_grid([1e9],[0.2],[0],[0])
	set_pretty()

	if sys.argv[1] == "grid":
		run_grid()
	elif sys.argv[1] == "rel" or sys.argv[1] == "relmax":

		lll = [1500,500,1000,2000,2500]
		for il in range(len(lll)):
			flambda_look = lll[il]

			set2 = set_set2()

			lstyles=["-", "--"]
			masses = [1e8,1e9]
			edd_frac = [.1, .2]
			spin = [0,0.998]
			#.998
			incs = np.arange(10,90,10)
			for im in range(len(masses)):
				figure(figsize=(8,8))
				for iedd in range(len(edd_frac)):

					mdot = mdot_from_edd(edd_frac[iedd], masses[im])

					

					pretty.set_pretty()

					for ispin in range(len(spin)):
						f2000 = []
						for ii in range(len(incs)):
							fname = "sp%i_inc%2i_mbh%i_edd%itminus2" % (int(spin[ispin]+0.5),incs[ii], np.log10(masses[im]), edd_frac[iedd]/0.01)

							w,f=read_agnspec(fname)

							ff = util.get_flux_at_wavelength(w,f,flambda_look)
							f2000.append(ff)
							#plot(w,f,label=incs[ii])
						if sys.argv[1] == "rel":
							plot(incs, f2000, linestyle=lstyles[ispin], c=set2[iedd], label="$\epsilon=%.2f$" % edd_frac[iedd])

						if spin[ispin] == 0 and edd_frac[iedd] == .2:
							fnorm = np.array(f2000)
				
				if masses[im] == 1e9 and flambda_look > 200 and flambda_look < 30000:
					s = r.read_spectrum('disk_only')
					f2000 = []

					for ii in range(len(incs)):
						f = s["A%iP0.50" % incs[ii]]

						f2000.append(util.get_flux_at_wavelength(s["Lambda"],f,flambda_look))

					if sys.argv[1] == "rel":
						plot(incs, f2000, c="k", label="Standard AD")
					else:
						plot(incs, fnorm/np.array(f2000))



				title("$M_{BH}=10^%i M_\odot$ (dashed=spinning)" % np.log10(masses[im]), fontsize=20)
				xlabel("Inclination", fontsize=16)
				ylabel("$F_{%i}$" % flambda_look, fontsize=16)
				#pretty_axes(
				pretty.float_legend(grey=False)
				#semilogy()

				savefig("normf%i_m%i.png" % (flambda_look, np.log10(masses[im]) ) )
				clf()



	else:
		masses = [1e8,1e9]
		edd_frac = [.1, .2]
		spin = [0,0.998]
		#.998
		incs = np.arange(10,90,10)

		for ispin in range(len(spin)):
			for im in range(len(masses)):
				for iedd in range(len(edd_frac)):

					mdot = mdot_from_edd(edd_frac[iedd], masses[im])

					figure(figsize=(8,8))

					pretty.set_pretty()

					for ii in range(len(incs)):
						fname = "sp%i_inc%2i_mbh%i_edd%itminus2" % (int(spin[ispin]+0.5),incs[ii], np.log10(masses[im]), edd_frac[iedd]/0.01)

						w,f=read_agnspec(fname)

						plot(w,f,label=incs[ii])

					semilogy()
					xlim(200,3000)
					title("$M_{BH}=10^%i M_\odot, \epsilon=%.2f, spin=%.3f$" % (np.log10(masses[im]), edd_frac[iedd], spin[ispin]) )

					savefig("plotsp%i_mbh%i_edd%itminus2.png" % (int(spin[ispin]+0.5),np.log10(masses[im]), edd_frac[iedd]/0.01))
					clf()



