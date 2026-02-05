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


time = np.arange(0, 201, 1)
Bz_err = 0.0*time
Ey_err = 0.0*time


for tstep in range(201):
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

  Bz_expect = np.sin(2*math.pi*(xs)/l0 - omega*ts)
  Ey_expect = np.sin(2*math.pi*(x)/l0 - omega*t)

  Bz_err[tstep] = (clight*Bz - Bz_expect).max()
  Ey_err[tstep] = (Ey - Ey_expect).max()
  # plt.plot(Ey, 'k-', linewidth = 2)
  # plt.plot(Ey_expect, 'b-')
  # plt.show()

plt.plot(time, Bz_err, 'k-', linewidth = 2)
plt.plot(time, Ey_err, 'b-', linewidth = 2)
plt.show()
