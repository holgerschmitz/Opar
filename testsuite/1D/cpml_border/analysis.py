#!/usr/bin/env python3

import numpy as np
import h5py
import math
from matplotlib import pyplot as plt

# Simulation parameters
clight = 299792458.0
l0 = 1e-6
dx = 0.05*l0
cflFactor = 0.99
dt = cflFactor*dx/clight
omega = 2*math.pi*clight/l0
# todo, we can make the analysis better by using the numerical dispersion relation


time = np.arange(0, 101, 1)
Bz_err = 0.0*time
Ey_err = 0.0*time


for tstep in range(101):
  print(tstep)
  fBz = h5py.File('Bz_{}.h5'.format(tstep), 'r')
  Bz = fBz['data'][:]
  fEy = h5py.File('Ey_{}.h5'.format(tstep), 'r')
  Ey = fEy['data'][:]

  L = Bz.size
  Lx = L*dx

  xs = np.arange(start=0.5*dx, stop=Lx, step=dx)
  x = np.arange(start=0.0*dx, stop=Lx, step=dx)
  t = tstep*dt
  ts = (tstep + 0.5)*dt

  Bz_err[tstep] = clight**2 * np.square(Bz).max()
  Ey_err[tstep] = np.square(Ey).max()
  # if tstep % 5 == 0:
  #   plt.plot(Ey, 'k-', linewidth = 2)
  #   plt.show()

plt.plot(time, Bz_err, 'k-', linewidth = 2)
plt.plot(time, Ey_err, 'b-', linewidth = 2)
plt.yscale('log')
plt.show()
