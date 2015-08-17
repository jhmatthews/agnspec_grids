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

'''
run an agnspec grid in serial, save to FOLDER
'''

FOLDER = "/Users/jmatthews/Documents/runs/QSO_clumped/disk_spec/agnspec_grids/outputs"

spins = [0,0.998]
masses = 10.0 ** np.arange(7.0,10.0,0.2)
edd_frac = np.arange(0,1,0.025)
edd_frac[0] = 0.01

incs = np.arange(10,100,10)
incs[-1] = 89.0

run_grid(masses,  edd_frac, spins, incs, "/Users/jmatthews/Dropbox/agnspec_grid2")
