/*
 * species.cpp
 *
 * Created on: 16 Nov 2012
 * Author: hschmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2012 Holger Schmitz
 *
 * This file is part of OPar.
 *
 * OPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "defs.hpp"
#include "species.hpp"
#include "currents.hpp"
#include "globals.hpp"
#include "random.hpp"

#include <schnek/grid.hpp>

#include <boost/foreach.hpp>


inline double ipow(double x, int y)
{
  double result = 1.0;
  for (int i = 0; i < y; ++i)
    result *= x;
  return result;
}

template<class Weighting>
Particle &Species<Weighting>::addParticle()
{
  return particles.addParticle();
}

template<class Weighting>
void Species<Weighting>::removeParticle(const ParticleStorage::iterator &p)
{
  particles.removeParticle(p);
}

template<class Weighting>
void Species<Weighting>::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  blockPars.addParameter("charge", &charge);
  blockPars.addParameter("mass", &mass);
  blockPars.addParameter("ppc", &ppc);

  densityParam = blockPars.addArray("density", &density);
  temperatureParam = blockPars.addArrayParameter("temperature", temperature);
  driftParam = blockPars.addArrayParameter("drift", drift);
}

template<class Weighting>
void Species<Weighting>::registerData()
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
  addData("Jx_species", pJx);
  addData("Jy_species", pJy);
  addData("Jz_species", pJz);
}

template<class Weighting>
void Species<Weighting>::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  dynamic_cast<OPar&>(*this->getParent()).addSpecies(this);

  SIntVector low = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  pJx = pDataField(new DataField(low, high, grange, SStagger(true, false), 2));
  pJy = pDataField(new DataField(low, high, grange, SStagger(false, true), 2));
  pJz = pDataField(new DataField(low, high, grange, SStagger(false, false), 2));
  pDensity = pDataField(new DataField(low, high, grange, SStagger(false, false), 2));

  Currents::instance().addCurrent(pJx, pJy, pJz);

  this->getData("Ex", pEx);
  this->getData("Ey", pEy);
  this->getData("Ez", pEz);
  this->getData("Bx", pBx);
  this->getData("By", pBy);
  this->getData("Bz", pBz);


  SDomain d = Weighting::getDomain();
  gx.resize(d);
  hx.resize(d);
  d.grow(2);
  jHelper.resize(d);
}

template<class Weighting>
void Species<Weighting>::initParticles()
{
  SIntVector lo = Globals::instance().getLocalGridMin();
  SIntVector hi = Globals::instance().getLocalGridMax();

  SIntVector i;

  SVector &coords = Globals::instance().getX();
  SVector dx = Globals::instance().getDx();
  pDependencyUpdater updater = Globals::instance().getUpdater(var_space);

  updater.addDependent(densityParam);
  updater.addDependentArray(temperatureParam);
  updater.addDependentArray(driftParam);

  double weight_factor = 1/ppc;

  // TODO this loop should be somewhere in schnek
  for (int i=0; i<dimension; ++i) weight_factor *= dx[i];

  // TODO check that we are only initialising in the interior of the domain
  SPACE_LOOP(pos, lo, hi)
  {
    SVector r;
    for (int n=0; n<ppc; ++n)
    {
      for (int i=0; i<dimension; ++i)
        coords[i] = (pos[i] + Random::uniform())*dx[i];

      Particle &p = particles.addParticle();

      // The updater changes the value of densityInit and temperatureInit by calculating the formulas from the user input
      updater.update();

      p.x = coords;
      p.weight = density*weight_factor;
      for (int i=0; i<3; ++i)
        p.u[i] = drift[i] + Random::gaussian(temperature[i]);
    }
  }
}

template<class Weighting>
void Species<Weighting>::pushParticles(double dt)
{
  pJx = 0.0;
  pJy = 0.0;
  pJz = 0.0;

  // Unvarying multiplication factor
  const SVector dx = Globals::instance().getDx();
  const SVector idx = 1.0 / dx;

  const double idt = 1.0 / dt;
  const double dto2 = dt / 2.0;
  const double dtco2 = clight * dto2;
  // particle weighting multiplication factor
  const double fac = Weighting::particleShapeFactor();
  const double facd = ipow(fac, dimension) * charge;

  const double part_mc = clight * mass;
  const double ipart_mc = 1.0_num / part_mc;

  const double cmratio = charge * dtfac * ipart_mc;
  const double ccmratio = clight * cmratio;

  BOOST_FOREACH(Particle &p_old, particles)
  {
    Particle p(p_old);
    double weight = p.weight;

    double gamma = sqrt(
        p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[1] * p.u[1] + 1.0);

    p.x = p.x + SVector(p.u[0], p.u[1]) * (0.5 * dt / gamma);

    FieldIndex cell1, cell2, dcell;
    SVector cell_frac;

    Weighting::toCellIndex(p.x * idx, cell1, cell_frac);
    Weighting::getShape(cell1, cell_frac, gx);

    Weighting::toCellIndexStagger(p.x * idx, cell2, cell_frac);
    Weighting::getShape(cell2, cell_frac, hx);

    PVector E = Weighting::interpolateE(gx, hx, pEx, pEy, pEz);
    PVector B = Weighting::interpolateB(gx, hx, pBx, pBy, pBz);

    PVector um = p.u + cmratio * E;
    gamma = sqrt(um[0] * um[0] + um[1] * um[1] + um[1] * um[1] + 1.0);

    PVector tau = B * (0.5 * dt / gamma);
    PVector tau2 = tau * tau;
    PVector tau_ifac = 1.0 / (tau2[0] + tau2[1] + tau2[1] + 1.0);

    // This is the EPOCH rotation code translated
    // TODO Check the two algorithms against each other
    PVector
        ur(
            ((1.0 + tau2[0] - tau2[1] - tau2[2]) * um[0] + 2.0
                * ((tau[0] * tau[1] + tau[2]) * um[1] + (tau[0] * tau[2]
                    - tau[1]) * um[2])) * tau_ifac,
            ((1.0 - tau2[0] + tau2[1] - tau2[2]) * um[1] + 2.0
                * ((tau[1] * tau[2] + tau[0]) * um[2] + (tau[1] * tau[0]
                    - tau[2]) * um[0])) * tau_ifac,
            ((1.0 - tau2[0] - tau2[1] + tauz[2]) * um[2] + 2.0
                * ((tau[2] * tau[0] + tau[1]) * um[0] + (tau[2] * tau[1]
                    - tau[0]) * um[1])) * tau_ifac);

    // This is the original OPar rotation code translated
    //    PVector ur(
    //        ((1.0 - tau2[1] - tau2[2]) * um[0] + (tau[0] * tau[1] + tau[2]) * um[1]
    //            + (tau[0] * tau[2] - tau[1]) * um[2]) * tau_ifac,
    //        ((1.0 - tau2[0] - tau2[2]) * um[1] + (tau[0] * tau[1] - tau[2]) * um[0]
    //            + (tau[1] * tau[2] + tau[0]) * um[2]) * tau_ifac,
    //        ((1.0 - tau2[0] - tau2[1]) * um[2] + (tau[0] * tau[2] + tau[1]) * um[0]
    //            + (tau[1] * tau[2] - tau[0]) * um[1]) * tau_ifac);

    p.u = ur + cmratio * E;

    gamma = sqrt(p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[1] * p.u[1] + 1.0);

    SVector delta = SVector(p.u[0], p.u[1]) * (0.5 * dt / gamma);

    p.x = p.x + delta;

    // Calculate the current using the charge conserving algorithm by Esirkepov
    // T.Z. Esirkepov, Comp. Phys. Comm., vol 135, p.144 (2001)

    SVector xplus = p.x + delta;
    Weighting::toCellIndex(p.x * idx, cell2, cell_frac);

    dcell = cell2 - cell1;
    Weighting::getShape(cell2, cell_frac, dcell, hx);

    jHelper = 0.0;
    RecDomain < 3 > d = Weighting::getDomain();

    const PIntVector lo = d.lo() + (dcell - PIntVector::Unity()) / 2;
    const PIntVector hi = d.hi() + (dcell + PIntVector::Unity()) / 2;

    const double sixth = 1.0 / 6.0;
    const double half = 1.0 / 2.0;
    FieldIndex l_ind;

#ifdef ONE_DIMENSIONAL
    const double vy = p.u[1] / gamma;
    const double vz = p.u[2] / gamma;
    const double fjx = idt * facd * weight;
    const double fjy = idx[0] * facd * weight * vy;
    const double fjy = idx[0] * facd * weight * vz;
#endif

#ifdef TWO_DIMENSIONAL
    const double vz = p.u[2] / gamma;
    const double fjx = idt * idx[1] * facd * part_weight;
    const double fjy = idt * idx[0] * facd * part_weight;
    const double fjz = idx[0] * idx[1] * facd * part_weight * vz;
#endif

#ifdef THREE_DIMENSIONAL
    const double fjx = idt * idx[1] * idx[2] * facd * part_weight;
    const double fjy = idt * idx[0] * idx[2] * facd * part_weight;
    const double fjz = idt * idx[0] * idx[1] * facd * part_weight;
#endif

SPACE_LOOP  (l_ind, lo, hi)
  {

#ifdef ONE_DIMENSIONAL
    const double sx = gx[l_ind[0]][0];
    const double rx = hx[l_ind[0]][0];

    const double wx = (rx-sx);
    const double wy = half * (sx+rx);

    jxHelper[l_ind] = jxHelper(i - 1, j, k) - fjx * wx;
    jyHelper[l_ind] = fjy * wy;
    jzHelper[l_ind] = fjz * wy;
#endif

#ifdef TWO_DIMENSIONAL
    const double sx = gx[l_ind[0]][0];
    const double sy = gx[l_ind[1]][1];

    const double rx = hx[l_ind[0]][0];
    const double ry = hx[l_ind[1]][1];

    const double wx = half * (rx-sx) * (sy+ry);
    const double wy = half * (ry-sy) * (sx+rx);
    const double wz = sixth * (rz-sz) * (sx * (2*sy+ry) + rx * (sy+2*ry));

    jxHelper[l_ind] = jxHelper(i - 1, j, k) - fjx * wx;
    jyHelper[l_ind] = jyHelper(i, j - 1, k) - fjy * wy;
    jzHelper[l_ind] = fjz * wz;
#endif

#ifdef THREE_DIMENSIONAL
    const double sx = gx[l_ind[0]][0];
    const double sy = gx[l_ind[1]][1];
    const double sz = gx[l_ind[2]][2];

    const double rx = hx[l_ind[0]][0];
    const double ry = hx[l_ind[1]][1];
    const double rz = hx[l_ind[2]][2];

    const double wx = sixth * (rx-sx) * (sy * (2*sz+rz) + ry * (sz+2*rz));
    const double wy = sixth * (ry-sy) * (sx * (2*sz+rz) + rx * (sz+2*rz));
    const double wz = sixth * (sx * (2*sy+ry) + rx * (sy+2*ry));

    jxHelper[l_ind] = jxHelper(i - 1, j, k) - fjx * wx;
    jyHelper[l_ind] = jyHelper(i, j - 1, k) - fjy * wy;
    jzHelper[l_ind] = jzHelper(i, j, k - 1) - fjz * wz;
#endif

    FieldIndex g_ind = cell1 + l_ind;

    pJx[g_ind] += jxHelper[l_ind];
    pJy[g_ind] += jyHelper[l_ind];
    pJz[g_ind] += jzHelper[l_ind];

  }

}

