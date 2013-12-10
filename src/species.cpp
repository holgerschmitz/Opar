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
#include <schnek/tools/literature.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/factory.hpp>

#include <fstream>

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

  blockPars.addParameter("charge", &charge, 1.0);
  blockPars.addParameter("mass", &mass, 1.0);
  blockPars.addParameter("ppc", &ppc, 100);

  densityParam = blockPars.addParameter("density", &density, 1.0);
  temperatureParam = blockPars.addArrayParameter("temperature", temperature);
  driftParam = blockPars.addArrayParameter("drift", drift, 0.0);

  blockPars.addArrayParameter("boundary_min", bcNamesLo, std::string("periodic"));
  blockPars.addArrayParameter("boundary_max", bcNamesHi, std::string("periodic"));

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
  static schnek::LiteratureArticle Esirkepov2001("Esirkepov2001", "Esirkepov, T Z",
      "Exact charge conservation scheme for Particle-in-Cell simulation with an arbitrary form-factor",
      "Computer Physics Communications", "2001", "135", "144--153");

  schnek::LiteratureManager::instance().addReference(
      "Particle weighting for calculation of currents uses charge conserving scheme by Esirkepov.",
      Esirkepov2001);

  SCHNEK_TRACE_ENTER_FUNCTION(2)
  dynamic_cast<OPar&>(*this->getParent()).addSpecies(this);

  SIntVector low = Globals::instance().getLocalInnerGridMin();
  SIntVector high = Globals::instance().getLocalInnerGridMax();
  SRange grange = Globals::instance().getDomainRange();

  pJx = pDataField(
      new DataField(low, high, grange, exStaggerYee, 3));
  pJy = pDataField(
      new DataField(low, high, grange, eyStaggerYee, 3));
  pJz = pDataField(
      new DataField(low, high, grange, ezStaggerYee, 3));

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
  SIntVector lo = Globals::instance().getLocalInnerGridMin();
  SIntVector hi = Globals::instance().getLocalInnerGridMax();

  SIntVector i;

  SVector &coords = Globals::instance().getX();
  SVector dx = Globals::instance().getDx();
  pDependencyUpdater updater = Globals::instance().getUpdater(var_space);

  updater->addDependent(densityParam);
  updater->addDependentArray(temperatureParam);
  updater->addDependentArray(driftParam);

  double weight_factor = 1.0 / (double)ppc;

  // TODO this loop should be somewhere in schnek
  for (int i = 0; i < dimension; ++i)
  {
    weight_factor *= dx[i];
  }

  FieldIndex pos;
  int debug_count = 0;
  double debug_min_x = hi[0];
  double debug_max_x = lo[0];
  int debug_pro_id = Globals::instance().getSubdivision()->procnum();
  SCHNEK_TRACE_LOG(4,"ppc = " << ppc)
  SCHNEK_TRACE_LOG(4,"dx[0] = " << dx[0])

//  std::string fname = str(boost::format("init%1%.out")%debug_pro_id);
//  std::ofstream debug_out(fname.c_str());

  SPACE_LOOP(pos,lo,hi)
  {
    SVector r;
    for (int n = 0; n < ppc; ++n)
    {
      ++debug_count;
      for (int i = 0; i < dimension; ++i)
      {
        coords[i] = (pos[i] + Random::uniform()) * dx[i];
        SCHNEK_TRACE_LOG(4,debug_count <<": coords[" << i << "] = " << coords[i])
        SCHNEK_TRACE_LOG(4,"pos[" << i << "] = " << pos[i])
//        debug_out << debug_count <<": coords[" << i << "] = " << coords[i] << std::endl;
//        debug_out << "pos[" << i << "] = " << pos[i] << std::endl;
      }

      Particle &p = particles.addParticle();

      // The updater changes the value of densityInit and temperatureInit by calculating the formulas from the user input
      updater->update();

      SCHNEK_TRACE_LOG(5,"density=" << density << "    weight_factor="<<weight_factor)

      p.x = coords;
      p.weight = density * weight_factor;
      for (int i = 0; i < 3; ++i)
        p.u[i] = drift[i] + Random::gaussian(temperature[i]);

      if (p.x[0]<debug_min_x) debug_min_x=p.x[0];
      if (p.x[0]>debug_max_x) debug_max_x=p.x[0];
    }
  }

//  debug_out << " Min Max Particle X "
//      << debug_min_x << " " << debug_max_x << std::endl;
//  debug_out << "Added " << debug_count << " particles\n";
//  debug_out.close();

}

void debug_check_out_of_bounds(std::string checkpoint, Particle p_debug_old, Particle p_old, FieldIndex debug_cell1, FieldIndex debug_cell2)
{
//  if (GridArgCheck<dimension>::getErrorFlag())
//  {
//    std::cerr << "Particle out of bounds at checkpoint " << checkpoint << std::endl;
//    std::cerr << "Before " << p_debug_old.x[0] << std::endl;
//    std::cerr << "After " << p_old.x[0] << std::endl;
//    std::cerr << "cell1 " << debug_cell1[0] << std::endl;
//    std::cerr << "cell2 " << debug_cell2[0] << std::endl;
//    std::cerr << "offending " << GridArgCheck<dimension>::getOffending()[0] << std::endl;
//    std::cerr << "Grid " << Globals::instance().getLocalGridMin()[0] << " "
//        << Globals::instance().getLocalGridMax()[0] << std::endl;
//    std::cerr << "Limits " << Globals::instance().getLocalDomainMin()[0] << " "
//        << Globals::instance().getLocalDomainMax()[0] << std::endl;
//    exit(-1);
//  }

}

/**
 * Push the particles and calulate the current.
 *
 * This is a relativistic particle pusher based on Boris' scheme
 *
 * @param dt time step
 */
void Species::pushParticles(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  (*pJx) = 0.0;
  (*pJy) = 0.0;
  (*pJz) = 0.0;


  const double clight = 1.0;

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

  const double cmratio = charge * dto2 * ipart_mc;
  const double ccmratio = clight * cmratio;

//  int debug_count = 0;
  for (ParticleStorage::iterator it=particles.begin(); it!=particles.end(); ++it)
  {
//    ++debug_count;
    Particle &p_old = *it;
    Particle p_debug_old(p_old);
    Particle p(p_old);

//    if (debug_count==debug_particle_number)
//    {
//      std::cerr << "Debug Particle " << debug_particle_number << std::endl;
//      std::cerr << "dt " << dt << " time " << Globals::instance().getT() << std::endl;
//      std::cerr << "Before " << p.x[0] << " (" <<p.u[0]<<", "<<p.u[1]<<", "<<p.u[2]<<")"<< std::endl;
//    }

    double weight = p.weight;

    double gamma = sqrt(
        p.u[0] * p.u[0] + p.u[1] * p.u[1] + p.u[2] * p.u[2] + 1.0);

    for (int i=0; i<dimension; ++i)
      p.x[i] = p.x[i] + p.u[i] * (0.5 * dt / gamma);

    FieldIndex cell1, cell2, dcell;
    FieldIndex debug_cell1, debug_cell2;
    SVector cell_frac;
    SVector cell_pos(p.x);

    debug_check_out_of_bounds("AA", p_debug_old, p_old, debug_cell1, debug_cell2);

    for (int i=0; i<dimension; ++i) cell_pos[i] *= idx[i];

    Weighting::toCellIndex(cell_pos, cell1, cell_frac);

    Weighting::getShape(cell1, cell_frac, gx);

    Weighting::toCellIndexStagger(cell_pos, cell2, cell_frac);

    Weighting::getShape(cell2, cell_frac, hx);

    debug_cell1 = cell1;
    debug_cell2 = cell2;

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

    p_old = p;

//    if (debug_count==debug_particle_number)
//    {
//      std::cerr << "After " << p.x[0] << " (" <<p.u[0]<<", "<<p.u[1]<<", "<<p.u[2]<<")"<< std::endl;
//    }

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

    debug_check_out_of_bounds("AM", p_debug_old, p_old, debug_cell1, debug_cell2);

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
      SCHNEK_TRACE_LOG(5,"particle current bits "<< l_ind[0] << " " << fjx << " " << wx << " " << idt << " " << facd << " " << weight)
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

      debug_check_out_of_bounds("AS", p_debug_old, p_old, debug_cell1, debug_cell2);

      FieldIndex g_ind;
      for (int i=0; i<dimension; ++i)
        g_ind[i] = cell1[i] + l_ind[i];

      (*pJx)[g_ind] += jxHelper[l_ind];
      (*pJy)[g_ind] += jyHelper[l_ind];
      (*pJz)[g_ind] += jzHelper[l_ind];

      SCHNEK_TRACE_LOG(5,"particle current "<< l_ind[0] << " " << jxHelper[l_ind] << " " << jyHelper[l_ind] << " " << jyHelper[l_ind])

      debug_check_out_of_bounds("AY", p_debug_old, p_old, debug_cell1, debug_cell2);
    }

    debug_check_out_of_bounds("AZ", p_debug_old, p_old, debug_cell1, debug_cell2);
  }

  // here all the boundary conditions and the MPI happens
  particleExchange->exchange(particles);
}


