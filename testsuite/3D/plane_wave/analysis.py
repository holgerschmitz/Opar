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
kz = 1
k = math.sqrt(kx*kx + ky*ky + kz*kz);
kt = math.sqrt(kx*kx + ky*ky);
dx = 0.05*l0
cflFactor = 0.99/math.sqrt(3.0)
dt = cflFactor*dx/clight
omega = 2*math.pi*clight*k/l0
phi = math.atan2(ky, kx)
theta = math.atan2(kz, kt)
# todo, we can make the analysis better by using the numerical dispersion relation


time = np.arange(0, 201, 1)
Bx_err = 0.0*time
By_err = 0.0*time
Bz_err = 0.0*time
Ex_err = 0.0*time
Ey_err = 0.0*time


for tstep in range(201):
  print(tstep)
  fBx = h5py.File('Bx_{}.h5'.format(tstep), 'r')
  Bx = fBx['data'][:,:]
  fBy = h5py.File('By_{}.h5'.format(tstep), 'r')
  By = fBy['data'][:,:]
  fBz = h5py.File('Bz_{}.h5'.format(tstep), 'r')
  Bz = fBz['data'][:,:]
  fEx = h5py.File('Ex_{}.h5'.format(tstep), 'r')
  Ex = fEx['data'][:,:]
  fEy = h5py.File('Ey_{}.h5'.format(tstep), 'r')
  Ey = fEy['data'][:,:]

  Nx = Bz.shape[0]
  Ny = Bz.shape[1]
  Nz = Bz.shape[2]
  Lx = Nx*dx
  Ly = Ny*dx
  Lz = Nz*dx

  xs = np.arange(start=0.5*dx, stop=Lx, step=dx)
  x = np.arange(start=0.0*dx, stop=Lx, step=dx)
  ys = np.arange(start=0.5*dx, stop=Ly, step=dx)
  y = np.arange(start=0.0*dx, stop=Ly, step=dx)
  zs = np.arange(start=0.5*dx, stop=Lz, step=dx)
  z = np.arange(start=0.0*dx, stop=Lz, step=dx)

  xxs, yys, zzs = np.meshgrid(xs, ys, zs, indexing='ij')
  xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')


  t = tstep*dt
  ts = (tstep + 0.5)*dt

  Bx_expect = -math.sin(theta) * math.cos(phi) * np.sin(2*math.pi*(kx*xx + ky*yys + kz*zzs)/l0 - omega*ts)
  By_expect = -math.sin(theta) * math.sin(phi) * np.sin(2*math.pi*(kx*xxs + ky*yy + kz*zzs)/l0 - omega*ts)
  Bz_expect =  math.cos(theta) * np.sin(2*math.pi*(kx*xxs + ky*yys + kz*zz)/l0 - omega*ts)
  Ex_expect = -math.sin(phi) * np.sin(2*math.pi*(kx*xxs + ky*yy + kz*zz)/l0 - omega*t)
  Ey_expect =  math.cos(phi) * np.sin(2*math.pi*(kx*xx + ky*yys + kz*zz)/l0 - omega*t)

  Bx_err[tstep] = (clight*Bx - Bx_expect).max()
  By_err[tstep] = (clight*By - By_expect).max()
  Bz_err[tstep] = (clight*Bz - Bz_expect).max()
  Ex_err[tstep] = (Ex - Ex_expect).max()
  Ey_err[tstep] = (Ey - Ey_expect).max()
  print(f' Max errors at step {tstep}: Bx: {Bx_err[tstep]}, By: {By_err[tstep]}, Bz: {Bz_err[tstep]}, Ex: {Ex_err[tstep]}, Ey: {Ey_err[tstep]}')
  # plt.plot(Ey, 'k-', linewidth = 2)
  # plt.show()

plt.plot(time, Bx_err, 'k-', linewidth = 2)
plt.plot(time, By_err, 'g-', linewidth = 2)
plt.plot(time, Bz_err, 'y-', linewidth = 2)
plt.plot(time, Ex_err, 'b-', linewidth = 2)
plt.plot(time, Ey_err, 'r-', linewidth = 2)
plt.show()
