#!/usr/bin/env python3

import numpy as np
import h5py
import math
from matplotlib import pyplot as plt

# Simulation parameters
clight = 299792458.0
l0 = 1e-6
kx = 1
ky = 1
dx = 0.05*l0
cflFactor = 0.99/math.sqrt(2.0)
dt = cflFactor*dx/clight
omega = 2*math.pi*clight*math.sqrt(kx**2 + ky**2)/l0
phi = math.atan2(ky, kx)
# todo, we can make the analysis better by using the numerical dispersion relation


time = np.arange(0, 201, 1)
Bz_err = 0.0*time
Ex_err = 0.0*time
Ey_err = 0.0*time





for tstep in range(201):
  print(tstep)
  fBz = h5py.File('Bz_{}.h5'.format(tstep), 'r')
  Bz = fBz['data'][:,:]
  fEx = h5py.File('Ex_{}.h5'.format(tstep), 'r')
  Ex = fEx['data'][:,:]
  fEy = h5py.File('Ey_{}.h5'.format(tstep), 'r')
  Ey = fEy['data'][:,:]

  Nx = Bz.shape[0]
  Ny = Bz.shape[1]
  Lx = Nx*dx
  Ly = Ny*dx

  xs = np.arange(start=0.5*dx, stop=Lx, step=dx)
  x = np.arange(start=0.0*dx, stop=Lx, step=dx)
  ys = np.arange(start=0.5*dx, stop=Ly, step=dx)
  y = np.arange(start=0.0*dx, stop=Ly, step=dx)

  xxs, yys = np.meshgrid(xs, ys, indexing='ij')
  xx, yy = np.meshgrid(x, y, indexing='ij')


  t = tstep*dt
  ts = (tstep + 0.5)*dt

  Bz_expect = np.sin(2*math.pi*(kx*xxs + ky*yys)/l0 - omega*ts)
  Ex_expect = -math.sin(phi) * np.sin(2*math.pi*(kx*xxs + ky*yy)/l0 - omega*t)
  Ey_expect =  math.cos(phi) * np.sin(2*math.pi*(kx*xx + ky*yys)/l0 - omega*t)

  Bz_err[tstep] = (clight*Bz - Bz_expect).max()
  Ex_err[tstep] = (Ex - Ex_expect).max()
  Ey_err[tstep] = (Ey - Ey_expect).max()
  # plt.plot(Ey, 'k-', linewidth = 2)
  # plt.show()

plt.plot(time, Bz_err, 'k-', linewidth = 2)
plt.plot(time, Ex_err, 'b-', linewidth = 2)
plt.plot(time, Ey_err, 'r-', linewidth = 2)
plt.show()
