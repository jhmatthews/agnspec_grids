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
from agnspec import *


spins = [0,0.998]
masses = 10.0 ** np.arange(7.0,10.0,0.5)
edd_frac = [0.01,0.1,1]
incs = np.arange(10,110,20)

f = [f for f in os.listdir(".") if 'tminus' in f]

set_pretty()
'''
for fname in f:
	w,flux=read_agnspec(fname)
	plot(w,flux)
	loglog()
	xlabel("Wavelength")
	ylabel("Flux")
	savefig(fname + ".png")
	clf()
'''
for ispin in range(len(spins)):
	for im in range(len(masses)):
		for iedd in range(len(edd_frac)):
			#mdot = mdot_from_edd(edd_frac[iedd], masses[
			for ii in range(len(incs)):
				logm = np.log10(masses[im])
				if logm - int(logm) < 0.2:
					m2 = 0
				else: m2 = 5
				mbh_string = "%ip%i" % (int(logm), m2)

				fname = "sp%i_inc%2i_mbh%s_edd%itminus2" % (int(spins[ispin]+0.5),incs[ii], mbh_string, edd_frac[iedd]/0.01)
				
				w,flux=read_agnspec(fname)
				plot(w,flux,label="$i=%2i^\circ$" % incs[ii])				
				loglog()


			pretty_legend()
			xlabel("Wavelength")
			ylabel("Flux")
			print "$a_* = %.3f, M_{BH}=10^{%i}, \epsilon = %.1f$" % (spins[ispin], np.log10(masses[im]), edd_frac[iedd])
			title("$a_* = %.3f, M_{BH}=10^{%i}, \epsilon = %.1f$" % (spins[ispin], np.log10(masses[im]), edd_frac[iedd]))
			savefig("spec_sp%i_mbh%i_edd%i.png" % (int(spins[ispin]+0.5), np.log10(masses[im]), edd_frac[iedd]/0.01))
			clf()
