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
#include "particle_boundary.hpp"
#include "opar.hpp"
#include "util.hpp"

#include <schnek/grid.hpp>
#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>
#include <boost/functional/factory.hpp>

#undef LOGLEVEL
#define LOGLEVEL 0

int debug_particle_number = 9702;

inline double
ipow(double x, int y)
{
  double result = 1.0;
  for (int i = 0; i < y; ++i)
    result *= x;
  return result;
}

std::map<std::string, Species*> Species::allSpecies;

Particle &Species::addParticle()
{
  return particles.addParticle();
}

void Species::removeParticle(const ParticleStorage::iterator &p)
{
  particles.removeParticle(p);
}

void Species::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  blockPars.addParameter("name", &name);

  blockPars.addParameter("charge", &charge);
  blockPars.addParameter("mass", &mass);
  blockPars.addParameter("ppc", &ppc);

  densityParam = blockPars.addParameter("density", &density);
  temperatureParam = blockPars.addArrayParameter("temperature", temperature);
  driftParam = blockPars.addArrayParameter("drift", drift);

  blockPars.addArrayParameter("boundary_min", bcNamesLo);
  blockPars.addArrayParameter("boundary_max", bcNamesHi);

  particleBCFactories["periodic"] = boost::factory<PeriodicParticleBoundary*>();
  particleBCFactories["open"] = boost::factory<OpenParticleBoundary*>();
  particleBCFactories["reflecting"] = boost::factory<ReflectingParticleBoundary*>();
}

void Species::registerData()
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
  addData("Jx_species", pJx);
  addData("Jy_species", pJy);
  addData("Jz_species", pJz);
}

void Species::preInit()
{
  allSpecies[name] = this;
}

Species* Species::getSpecies(std::string name)
{
  if (allSpecies.count(name)>0)
    return allSpecies[name];

  throw std::string("Species "+name+" could not be found!");
}

void Species::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  dynamic_cast<OPar&>(*this->getParent()).addSpecies(this);

  SIntVector low = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();
  SRange grange = Globals::instance().getDomainRange();

  pJx = pDataField(
      new DataField(low, high, grange, exStaggerYee, 2));
  pJy = pDataField(
      new DataField(low, high, grange, eyStaggerYee, 2));
  pJz = pDataField(
      new DataField(low, high, grange, ezStaggerYee, 2));

  (*pJx) = 0.0;
  (*pJy) = 0.0;
  (*pJz) = 0.0;

  Currents::instance().addCurrent(pJx, pJy, pJz);

  this->retrieveData("Ex", pEx);
  this->retrieveData("Ey", pEy);
  this->retrieveData("Ez", pEz);
  this->retrieveData("Bx", pBx);
  this->retrieveData("By", pBy);
  this->retrieveData("Bz", pBz);

  ScalarDomain ds = Weighting::getScalarDomain();
  gx.resize(ds.getLo(), ds.getHi());
  hx.resize(ds.getLo(), ds.getHi());
  SDomain d = Weighting::getSDomain();
  d.grow(2);
  jxHelper.resize(d.getLo(), d.getHi());
  jyHelper.resize(d.getLo(), d.getHi());
  jzHelper.resize(d.getLo(), d.getHi());

  for (int i = 0; i < dimension; ++i)
  {
    if (particleBCFactories.count(bcNamesLo[i]) == 0) terminateSim(
        "Unknown boundary condition: " + bcNamesLo[i]);

    if (particleBCFactories.count(bcNamesHi[i]) == 0) terminateSim(
        "Unknown boundary condition: " + bcNamesHi[i]);

    boundariesLo[i] = pParticleBoundary(particleBCFactories[bcNamesLo[i]](i,-1));
    boundariesHi[i] = pParticleBoundary(particleBCFactories[bcNamesHi[i]](i, 1));
  }
  particleExchange = pParticleExchange(new ParticleExchange(*this));

  this->initParticles();
}

void Species::initParticles()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  SIntVector lo = Globals::instance().getLocalGridMin();
  SIntVector hi = Globals::instance().getLocalGridMax();

  SIntVector i;

  SVector &coords = Globals::instance().getX();
  SVector dx = Globals::instance().getDx();
  pDependencyUpdater updater = Globals::instance().getUpdater(var_space);

  updater->addDependent(densityParam);
  updater->addDependentArray(temperatureParam);
  updater->addDependentArray(driftParam);

  double weight_factor = 1 / ppc;

  // TODO this loop should be somewhere in schnek
  for (int i = 0; i < dimension; ++i)
  {
    weight_factor *= dx[i];
    --hi[i];
  }

  FieldIndex pos;
  int debug_count = 0;
  SCHNEK_TRACE_LOG(3,"ppc = " << ppc)
  SCHNEK_TRACE_LOG(3,"dx[0] = " << dx[0])
  // TODO check that we are only initialising in the interior of the domain
  SPACE_LOOP(pos,lo,hi)
  {
    SVector r;
    for (int n = 0; n < ppc; ++n)
    {
      ++debug_count;
      for (int i = 0; i < dimension; ++i)
      {
        coords[i] = (pos[i] + Random::uniform()) * dx[i];
        SCHNEK_TRACE_LOG(3,debug_count <<": coords[" << i << "] = " << coords[i])
        SCHNEK_TRACE_LOG(3,"pos[" << i << "] = " << pos[i])
      }

      Particle &p = particles.addParticle();

      // The updater changes the value of densityInit and temperatureInit by calculating the formulas from the user input
      updater->update();

      p.x = coords;
      p.weight = density * weight_factor;
      for (int i = 0; i < 3; ++i)
        p.u[i] = drift[i] + Random::gaussian(temperature[i]);
    }
  }
}

/**
 * Push the particles and calulate the current.
 *
 * This is a relativistic particle pusher based on Boris' scheme
 *
 * @todo check normalisation
 * @todo double check the Esirkepov current calculation and document the maths somewhere
 * @param dt time step
 */
void Species::pushParticles(double dt)
{
  (*pJx) = 0.0;
  (*pJy) = 0.0;
  (*pJz) = 0.0;


  const double clight = 1.0;
  const double dtfac = 1.0;

  // Unvarying multiplication factor
  const SVector dx = Globals::instance().getDx();

#ifdef ONE_DIMENSIONAL
  const SVector idx(1.0/dx[0]);
#endif

#ifdef TWO_DIMENSIONAL
  const SVector idx(1.0/dx[0], 1.0/dx[1]);
#endif

#ifdef THREE_DIMENSIONAL
  const SVector idx(1.0/dx[0], 1.0/dx[1], 1.0/dx[2]);
#endif

  const double idt = 1.0 / dt;
  const double dto2 = dt / 2.0;
  const double dtco2 = clight * dto2;
  // particle weighting multiplication factor
  const double fac = Weighting::particleShapeFactor();
  const double facd = ipow(fac, dimension) * charge;

  const double part_mc = clight * mass;
  const double ipart_mc = 1.0 / part_mc;

  const double cmratio = charge * dtfac * ipart_mc;
  const double ccmratio = clight * cmratio;

  int debug_count = 0;
  for (ParticleStorage::iterator it=particles.begin(); it!=particles.end(); ++it)
  {
    ++debug_count;
    Particle &p_old = *it;
    Particle p(p_old);

    if (debug_count==debug_particle_number)
    {
      std::cout << "Debug Particle " << debug_particle_number << std::endl;
    }

    double weight = p.weight;

    double gamma = sqrt(
        p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[1] * p.u[1] + 1.0);

    for (int i=0; i<dimension; ++i)
      p.x[i] = p.x[i] + p.u[i] * (0.5 * dt / gamma);

    FieldIndex cell1, cell2, dcell;
    SVector cell_frac;
    SVector cell_pos(p.x);
    for (int i=0; i<dimension; ++i) cell_pos[i] *= idx[i];

    Weighting::toCellIndex(cell_pos, cell1, cell_frac);
    Weighting::getShape(cell1, cell_frac, gx);

    Weighting::toCellIndexStagger(cell_pos, cell2, cell_frac);
    Weighting::getShape(cell2, cell_frac, hx);

    PVector E = Weighting::interpolateE(gx, hx, cell1, cell2, *pEx, *pEy, *pEz);
    PVector B = Weighting::interpolateB(gx, hx, cell1, cell2, *pBx, *pBy, *pBz);

    PVector um;
    for (int i=0; i<3; ++i)
      um[i] = p.u[i] + cmratio * E[i];

    gamma = sqrt(um[0] * um[0] + um[1] * um[1] + um[1] * um[1] + 1.0);

    PVector tau;
    for (int i=0; i<3; ++i) tau[i] = B[i] * (0.5 * dt / gamma);

    PVector tau2 = tau;
    for (int i=0; i<3; ++i) tau2[i] *= tau[i];
    double tau_ifac = 1.0 / (tau2[0] + tau2[1] + tau2[2] + 1.0);

    // This is the EPOCH rotation code translated
    // TODO Check the two algorithms against each other
    PVector ur(
        ((1.0 + tau2[0] - tau2[1] - tau2[2]) * um[0]
            + 2.0
                * ((tau[0] * tau[1] + tau[2]) * um[1]
                    + (tau[0] * tau[2] - tau[1]) * um[2])) * tau_ifac,
        ((1.0 - tau2[0] + tau2[1] - tau2[2]) * um[1]
            + 2.0
                * ((tau[1] * tau[2] + tau[0]) * um[2]
                    + (tau[1] * tau[0] - tau[2]) * um[0])) * tau_ifac,
        ((1.0 - tau2[0] - tau2[1] + tau2[2]) * um[2]
            + 2.0
                * ((tau[2] * tau[0] + tau[1]) * um[0]
                    + (tau[2] * tau[1] - tau[0]) * um[1])) * tau_ifac);

    // This is the original OPar rotation code translated
    //    PVector ur(
    //        ((1.0 - tau2[1] - tau2[2]) * um[0] + (tau[0] * tau[1] + tau[2]) * um[1]
    //            + (tau[0] * tau[2] - tau[1]) * um[2]) * tau_ifac,
    //        ((1.0 - tau2[0] - tau2[2]) * um[1] + (tau[0] * tau[1] - tau[2]) * um[0]
    //            + (tau[1] * tau[2] + tau[0]) * um[2]) * tau_ifac,
    //        ((1.0 - tau2[0] - tau2[1]) * um[2] + (tau[0] * tau[2] + tau[1]) * um[0]
    //            + (tau[1] * tau[2] - tau[0]) * um[1]) * tau_ifac);

    for (int i=0; i<3; ++i) p.u[i] = ur[i] + cmratio * E[i];

    gamma = sqrt(p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[1] * p.u[1] + 1.0);

    SVector delta;

    for (int i=0; i<3; ++i) delta[i] = p.u[i] * (0.5 * dt / gamma);
    for (int i=0; i<dimension; ++i) p.x[i] = p.x[i] + delta[i];

    cell_pos = p.x;
    for (int i=0; i<dimension; ++i) cell_pos[i] *= idx[i];

    // Calculate the current using the charge conserving algorithm by Esirkepov
    // T.Z. Esirkepov, Comp. Phys. Comm., vol 135, p.144 (2001)

    // SVector xplus = p.x + delta;
    Weighting::toCellIndex(cell_pos, cell2, cell_frac);

    for (int i=0; i<dimension; ++i)
      dcell[i] = cell2[i] - cell1[i];

    Weighting::getShape(cell2, cell_frac, hx);

    jxHelper = 0.0;
    jyHelper = 0.0;
    jzHelper = 0.0;
    SDomain d = Weighting::getSDomain();

    SIntVector lo, dlo = d.getLo();
    SIntVector hi, dhi = d.getHi();
    for (int i=0; i<dimension; ++i)
    {
      //lo[i] = d.getLo()[i] + (dcell[i] - 1) / 2;
      //hi[i] = d.getHi()[i] + (dcell[i] + 1) / 2;

      lo[i] = d.getLo()[i] + std::min(0,dcell[i]);
      hi[i] = d.getHi()[i] + std::max(0,dcell[i]);
    }

    const double sixth = 1.0 / 6.0;
    const double half = 1.0 / 2.0;
    FieldIndex l_ind;

#ifdef ONE_DIMENSIONAL
    const double vy = p.u[1] / gamma;
    const double vz = p.u[2] / gamma;
    const double fjx = idt * facd * weight;
    const double fjy = idx[0] * facd * weight * vy;
    const double fjz = idx[0] * facd * weight * vz;
#endif

#ifdef TWO_DIMENSIONAL
    const double vz = p.u[2] / gamma;
    const double fjx = idt * idx[1] * facd * weight;
    const double fjy = idt * idx[0] * facd * weight;
    const double fjz = idx[0] * idx[1] * facd * weight * vz;
#endif

#ifdef THREE_DIMENSIONAL
    const double fjx = idt * idx[1] * idx[2] * facd * weight;
    const double fjy = idt * idx[0] * idx[2] * facd * weight;
    const double fjz = idt * idx[0] * idx[1] * facd * weight;
#endif

    SPACE_LOOP (l_ind, lo, hi)
    {
      // Particle weighting according to
      // T. Z. Esirkepov, Computer Physics Communications 135, 144 (2001).
      //
      // The weighting coefficients have been double checked and should be OK
#ifdef ONE_DIMENSIONAL
      int i = l_ind[0];
      int ih = l_ind[0] - dcell[0];
      const double sx = ((i<dlo[0])||(i>dhi[0]))?0.0:gx[i][0];
      const double rx = ((ih<dlo[0])||(ih>dhi[0]))?0.0:hx[ih][0];

      const double wx = (rx-sx);
      const double wy = half * (sx+rx);

      jxHelper[l_ind] = jxHelper(i - 1) - fjx * wx;
      jyHelper[l_ind] = fjy * wy;
      jzHelper[l_ind] = fjz * wy;
#endif

#ifdef TWO_DIMENSIONAL
      int i = l_ind[0];
      int j = l_ind[1];
      int ih = l_ind[0] - dcell[0];
      int jh = l_ind[1] - dcell[1];
      const double sx = ((i<dlo[0])||(i>dhi[0]))?0.0:gx[i][0];
      const double sy = ((j<dlo[1])||(j>dhi[1]))?0.0:gx[j][1];

      const double rx = ((ih<dlo[0])||(ih>dhi[0]))?0.0:hx[ih][0];
      const double ry = ((jh<dlo[1])||(jh>dhi[1]))?0.0:hx[jh][1];

      const double wx = half * (rx-sx) * (sy+ry);
      const double wy = half * (ry-sy) * (sx+rx);
      const double wz = sixth * (sx * (2*sy+ry) + rx * (sy+2*ry));

      jxHelper[l_ind] = jxHelper(i - 1, j) - fjx * wx;
      jyHelper[l_ind] = jyHelper(i, j - 1) - fjy * wy;
      jzHelper[l_ind] = fjz * wz;
#endif

#ifdef THREE_DIMENSIONAL
      int i = l_ind[0];
      int j = l_ind[1];
      int k = l_ind[2];
      int ih = l_ind[0] - dcell[0];
      int jh = l_ind[1] - dcell[1];
      int kh = l_ind[2] - dcell[2];
      const double sx = ((i<dlo[0])||(i>dhi[0]))?0.0:gx[i][0];
      const double sy = ((j<dlo[1])||(j>dhi[1]))?0.0:gx[j][1];
      const double sz = ((k<dlo[2])||(k>dhi[2]))?0.0:gx[k][2];

      const double rx = ((ih<dlo[0])||(ih>dhi[0])?0.0:hx[ih][0];
      const double ry = ((jh<dlo[1])||(jh>dhi[1])?0.0:hx[jh][1];
      const double rz = ((kh<dlo[2])||(kh>dhi[2])?0.0:hx[kh][2];

      const double wx = sixth * (rx - sx)
          * (sy * (2 * sz + rz) + ry * (sz + 2 * rz));
      const double wy = sixth * (ry - sy)
          * (sx * (2 * sz + rz) + rx * (sz + 2 * rz));
      const double wz = sixth * (rz - sz)
          * (sx * (2 * sy + ry) + rx * (sy + 2 * ry));

      jxHelper[l_ind] = jxHelper(i - 1, j, k) - fjx * wx;
      jyHelper[l_ind] = jyHelper(i, j - 1, k) - fjy * wy;
      jzHelper[l_ind] = jzHelper(i, j, k - 1) - fjz * wz;
#endif

      FieldIndex g_ind;
      for (int i=0; i<dimension; ++i)
        g_ind[i] = cell1[i] + l_ind[i];

      (*pJx)[g_ind] += jxHelper[l_ind];
      (*pJy)[g_ind] += jyHelper[l_ind];
      (*pJz)[g_ind] += jzHelper[l_ind];

    }
  }

  // here all the boundary conditions and the MPI happens
  particleExchange->exchange(particles);
}

