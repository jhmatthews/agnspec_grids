AGNSPEC GRID
-------------

Files of format

sp0_inc10_mbh7p0_edd1tminus2

sp0 -- gives you value of spin, either 0 for sp0 or 0.998 for sp1
inc10 -- gives inclination in degrees
mbh 7p0 -- gives log of bh mass. 7p5 means log(m) = 7.5
edd1tminus2 -- gives eddington fraction. 1tminus2 means 0.01, 100tminus2 means 1.00.

-------------

Columns:

Nu (Hz), LNu (erg /s /Hz), 3 other columns I can't remember- maybe convergence things


-------------


Cautionary Notes from Shane Davis:

'''
There are some files describing how to use this, but I just want to add some additional words of caution. 

The code works by generating a full accretion disk spectrum by numerically integrating over the spectra from a range of radii.  The spectra at any single radius is computed by interpolation in a table of annuli spectra parameterized by the effective temperature, surface density (mass/area) and a gravity parameters (essential the square of the Keplerian frequency).  

The problem is that the table of annuli that goes into this angered.save file only covers a limited range.  Since there are not a large number of very hot annuli in the table, it sometimes will extrapolate spectra beyond the table of converged models and it won’t warn you that it is doing this.  This can lead to rather inaccurate spectra — particularly for lower masses (M ~< 10^7 Msun) and in the extreme UV or soft X-rays.  At some point, I had modified the IDL code to warn me when this happened so I could reject models that were wrong but I can’t seem to find this code. However, if you are focussing on spectral features at or long ward of the Lyman edge, I think you should be ok.
'''

-------------

And Note about H/He opacities:

'''
For these models, I only used H and He since the bound-free metal opacities aren’t that important and they made it more difficult to get converged non-LTE models.  Although, I think this also the case for the AGNSPEC models if memory serves — that they are also based on H/He models.  I am sure how closely my table of annuli matches the table of annuli that Ivan used to generate the agngrid.save file, but it is not sampled as densely and uniformly as what I have used for BHSPEC models.  Basically, I just had trouble getting models to converge for certain parameter combinations and I never worked particularly hard at getting a widely converged table of annuli in this range.  The few times I have generated AGN models, I was just computing a very specific set of full disk models so I didn’t attempt a general table. And, I did not make any improvements — they should be virtually identical to the ones in AGNSPEC.  
'''

-------------
