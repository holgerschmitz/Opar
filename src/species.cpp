/*
 * species.cpp
 *
 * Created on: 16 Nov 2012
 * Author: Holger Schmitz
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
#include "constants.hpp"

#include <schnek/grid.hpp>
#include <schnek/util/logger.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/functions.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/factory.hpp>

#include <fstream>
#include <limits>

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
  blockPars.addParameter("densityCutoff", &densityCutoff, 0.0);

  densityParam = blockPars.addParameter("density", &density, 1.0);
  temperatureParam = blockPars.addArrayParameter("temperature", temperature, 0.0);
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

#undef LOGLEVEL
#define LOGLEVEL 0

void Species::initParticles()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  const SIntVector lo = Globals::instance().getLocalInnerGridMin();
  const SIntVector hi = Globals::instance().getLocalInnerGridMax();
  const SVector dMin = Globals::instance().getDomainMin();

  SIntVector i;

  SVector &coords = Globals::instance().getX();
  SVector dx = Globals::instance().getDx();
  pDependencyUpdater updater = Globals::instance().getUpdater(var_space);

  updater->addDependent(densityParam);
  updater->addDependentArray(temperatureParam);
  updater->addDependentArray(driftParam);

  double weight_factor = dx.product() / (double)ppc;

  SIntVector pos;
  int debug_count = 0;
  double debug_min_x = hi[0];
  double debug_max_x = lo[0];
//  int debug_pro_id = Globals::instance().getSubdivision()->procnum();

  SCHNEK_TRACE_LOG(4,"ppc = " << ppc)
  SCHNEK_TRACE_LOG(4,"dx[0] = " << dx[0])

  SPACE_LOOP(pos,lo,hi)
  {
    SVector r;
    for (int n = 0; n < ppc; ++n)
    {
      ++debug_count;
      for (int i = 0; i < dimension; ++i)
      {
        coords[i] = (pos[i] + Random::uniform()) * dx[i] + dMin[i];
        SCHNEK_TRACE_LOG(4,debug_count <<": coords[" << i << "] = " << coords[i])
        SCHNEK_TRACE_LOG(4,"pos[" << i << "] = " << pos[i])
//        debug_out << debug_count <<": coords[" << i << "] = " << coords[i] << std::endl;
//        debug_out << "pos[" << i << "] = " << pos[i] << std::endl;
      }

      // The updater changes the value of densityInit and temperatureInit by calculating the formulas from the user input
      updater->update();
      if (density < densityCutoff) continue;

      Particle &p = particles.addParticle();

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

#undef LOGLEVEL
#define LOGLEVEL 0

void debug_check_out_of_bounds(std::string checkpoint, Particle p_debug_old, Particle p_old, SIntVector debug_cell1, SIntVector debug_cell2)
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

inline bool debug_particle_nan(std::string checkpoint, const Particle &p)
{
  if (isnan(p.x[0]) || isnan(p.x[1]) || isnan(p.u[0]) || isnan(p.u[1]) || isnan(p.u[2]))
  {
    std::cerr << "Particle NaN at checkpoint " << checkpoint << std::endl;
    std::cerr << "X = " << p.x << "\nU = " << p.u << "\n";
    return true;
  }
  return false;
}

#undef LOGLEVEL
#define LOGLEVEL 0

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

  // Unvarying multiplication factor
  const SVector dx = Globals::instance().getDx();
  const SVector dMin = Globals::instance().getDomainMin();

  const SVector idx(1.0/dx);
  const double idt = 1.0 / dt;

  // particle weighting multiplication factor
  const double fac = ipow(Weighting::particleShapeFactor(), dimension) * charge;

  const double cmratio = fac * dt / mass;


//  double maxJxHelper = 0.0;
//  double maxJyHelper = 0.0;
//  double maxJzHelper = 0.0;

#if BOOST_PP_GREATER_EQUAL( LOGLEVEL, 5 )
  double minVal = 0, maxVal = 0;
  for (ParticleStorage::iterator it=particles.begin(); it!=particles.end(); ++it)
  {
    Particle &p = *it;
    for (int i=0; i<dimension; ++i)
    {
      minVal = std::min(minVal, p.x[i]);
      maxVal = std::max(maxVal, p.x[i]);
    }
    for (int i=0; i<3; ++i)
    {
      minVal = std::min(minVal, p.u[i]);
      maxVal = std::max(maxVal, p.u[i]);
    }
  }

  SCHNEK_TRACE_LOG(5,"Particle minimum val = " << minVal << "; maximum val " << maxVal)
#endif

//  int debug_count = 0;
  for (ParticleStorage::iterator it=particles.begin(); it!=particles.end(); ++it)
  {
//    ++debug_count;
    Particle &particle = *it;
    Particle p_debug_old(particle);
    Particle p(particle);
    debug_particle_nan("A", particle);
    debug_particle_nan("A+", p);

//    if (debug_count==debug_particle_number)
//    {
//      std::cerr << "Debug Particle " << debug_particle_number << std::endl;
//      std::cerr << "dt " << dt << " time " << Globals::instance().getT() << std::endl;
//      std::cerr << "Before " << p.x[0] << " (" <<p.u[0]<<", "<<p.u[1]<<", "<<p.u[2]<<")"<< std::endl;
//    }

    double wfac = p.weight*fac;
    double dt_gamma = 0.5*dt/sqrt(p.u.sqr()/clight2 + 1.0);

    p.x = p.x + p.u.project<dimension>() * dt_gamma;
    debug_particle_nan("B", p);

    SIntVector cell1, cell2, dcell;
    SIntVector debug_cell1, debug_cell2;
    SVector cell_frac;
    SVector cell_pos((p.x - dMin)*idx);

    debug_check_out_of_bounds("AA", p_debug_old, particle, debug_cell1, debug_cell2);

    Weighting::toCellIndex(cell_pos, cell1, cell_frac);
    Weighting::getShape(cell1, cell_frac, gx);
    Weighting::toCellIndexStagger(cell_pos, cell2, cell_frac);
    Weighting::getShape(cell2, cell_frac, hx);

    debug_cell1 = cell1;
    debug_cell2 = cell2;

    PVector E = Weighting::interpolateE(gx, hx, cell1, cell2, *pEx, *pEy, *pEz);
    PVector B = Weighting::interpolateB(gx, hx, cell1, cell2, *pBx, *pBy, *pBz);

    PVector um = p.u + 0.5*cmratio * E;

    dt_gamma = 1.0/sqrt(um.sqr()/clight2 + 1.0);

    const PVector tau = B * (cmratio*dt_gamma);
    const PVector tau2 = tau*tau;
    PVector ud, urot;

    const double tau_ifac = 2.0 / (tau2[0] + tau2[1] + tau2[2] + 1.0);

    schnek::crossProduct(ud, um, tau);
    ud += um;
    schnek::crossProduct(urot, ud,tau);

    for (int i=0; i<3; ++i)
      p.u[i] = um[i] + tau_ifac*urot[i] + 0.5*cmratio * E[i];

#ifdef TWO_DIMENSIONAL

    if (debug_particle_nan("C", p))
    {
      std::cerr << "E = "<< E << "\nB = "<<B<< "\num = "<<um<<"\ntau_ifac = "<<tau_ifac<<"\nurot = "<<urot<<"\ncmratio = "<<cmratio<<"\n";
      std::cerr << "gx = "<< gx<< "\nhx = "<< hx<< "\ncell1 = "<< cell1<< "\ncell2 = "<< cell2<<"\n";
      SDomain d = Weighting::getSDomain();

      for (int j=d.getLo()[1]; j<d.getHi()[1]; ++j)
      {
        for (int i=d.getLo()[0]; i<d.getHi()[0]; ++i)
        {
          std::cerr << "Ex("<< cell2[0]+i <<", "<< cell1[1]+j<<") = " << (*pEx)(cell2[0]+i, cell1[1]+j)<<"\n";
        }
      }
      exit(-1);
    }
#endif

    const double igamma = 1.0/sqrt(p.u.sqr()/clight2 + 1.0);

    for (int i=0; i<dimension; ++i)
      p.x[i] = p.x[i] + p.u[i] * (0.5 * dt * igamma);

    debug_particle_nan("D", p);

    particle = p;

//    if (debug_count==debug_particle_number)
//    {
//      std::cerr << "After " << p.x[0] << " (" <<p.u[0]<<", "<<p.u[1]<<", "<<p.u[2]<<")"<< std::endl;
//    }

    cell_pos = (p.x - dMin)*idx;

    // Calculate the current using the charge conserving algorithm by Esirkepov
    // T.Z. Esirkepov, Comp. Phys. Comm., vol 135, p.144 (2001)

    // SVector xplus = p.x + delta;
    Weighting::toCellIndex(cell_pos, cell2, cell_frac);

    for (int i=0; i<dimension; ++i)
    {
      dcell[i] = cell2[i] - cell1[i];
      if ((dcell[i]>1) || (dcell[i]<-1))
      {
        std::cerr << "Particle is moving more that one grid cell per time step.\nStopping now!\n";
        std::cerr << "X = " << p.x << "\nU = " << p.u << "\ndt = " << dt << "\nigamma = "<<igamma<< "\n";
        exit(-1);
      }
    }

    Weighting::getShape(cell2, cell_frac, hx);

    jxHelper = 0.0;
    jyHelper = 0.0;
    jzHelper = 0.0;
    SDomain d = Weighting::getSDomain();

    SIntVector lo, dlo = d.getLo();
    SIntVector hi, dhi = d.getHi();
    for (int i=0; i<dimension; ++i)
    {
      lo[i] = d.getLo()[i] + std::min(0,dcell[i]);
      hi[i] = d.getHi()[i] + std::max(0,dcell[i]);
    }

    const double sixth = 1.0 / 6.0;
#ifndef THREE_DIMENSIONAL
    const double half = 1.0 / 2.0;
#endif
    SIntVector l_ind;

    debug_check_out_of_bounds("AM", p_debug_old, particle, debug_cell1, debug_cell2);

#ifdef ONE_DIMENSIONAL
    const double vy = p.u[1] * igamma;
    const double vz = p.u[2] * igamma;
    const double fjx = idt * wfac;
    const double fjy = idx[0] * wfac * vy;
    const double fjz = idx[0] * wfac * vz;
#endif

#ifdef TWO_DIMENSIONAL
    const double vz = p.u[2] * igamma;
    const double fjx = idt * idx[1] * wfac;
    const double fjy = idt * idx[0] * wfac;
    const double fjz = idx[0] * idx[1] * wfac * vz;
#endif

#ifdef THREE_DIMENSIONAL
    const double fjx = idt * idx[1] * idx[2] * wfac;
    const double fjy = idt * idx[0] * idx[2] * wfac;
    const double fjz = idt * idx[0] * idx[1] * wfac;
#endif

    SCHNEK_TRACE_LOG(5,"loop "<< lo[0] << " " << hi[0] << " " << d.getLo()[0] << " " << d.getHi()[0] << " " << dcell[0])
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
      SCHNEK_TRACE_LOG(5,"particle current bits "<< l_ind[0] << " " << fjx << " " << wx << " " << idt << " " << wfac << " " << weight)
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

      const double rx = ((ih<dlo[0])||(ih>dhi[0]))?0.0:hx[ih][0];
      const double ry = ((jh<dlo[1])||(jh>dhi[1]))?0.0:hx[jh][1];
      const double rz = ((kh<dlo[2])||(kh>dhi[2]))?0.0:hx[kh][2];

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

      debug_check_out_of_bounds("AS", p_debug_old, particle, debug_cell1, debug_cell2);

      SIntVector g_ind;
      for (int i=0; i<dimension; ++i)
        g_ind[i] = cell1[i] + l_ind[i];

//      if (fabs(jxHelper[l_ind])>maxJxHelper) maxJxHelper=fabs(jxHelper[l_ind]);
//      if (fabs(jyHelper[l_ind])>maxJyHelper) maxJyHelper=fabs(jyHelper[l_ind]);
//      if (fabs(jzHelper[l_ind])>maxJzHelper) maxJzHelper=fabs(jzHelper[l_ind]);

//      if ((jxHelper[l_ind]>0) || (jyHelper[l_ind]>0) || (jzHelper[l_ind]>0))
//      if ((Globals::instance().getTCount() == 358) || (Globals::instance().getTCount() == 300))
//        SCHNEK_TRACE_LOG(0,"pos "<< Globals::instance().getTCount() << " "
//                         << cell1[0] << " "  << cell1[1]  << " "
//                         << l_ind[0] << " " << l_ind[1] << " "
//                         << g_ind[0] << " " << g_ind[1] << " "
//                         << jxHelper.getLo()[0] << " " << jxHelper.getHi()[0] << " " << pJx->getLo()[0] << " " << pJx->getHi()[0])

      (*pJx)[g_ind] += jxHelper[l_ind];
      (*pJy)[g_ind] += jyHelper[l_ind];
      (*pJz)[g_ind] += jzHelper[l_ind];

//      if ((jxHelper[l_ind]>0) || (jyHelper[l_ind]>0) || (jzHelper[l_ind]>0))
      SCHNEK_TRACE_LOG(5,"particle current "<< l_ind[0] << " " << jxHelper[l_ind] << " " << jyHelper[l_ind] << " " << jzHelper[l_ind])

      debug_check_out_of_bounds("AY", p_debug_old, particle, debug_cell1, debug_cell2);
    }

    debug_check_out_of_bounds("AZ", p_debug_old, particle, debug_cell1, debug_cell2);
  }

//  std::cerr << "maxJxHelper = " << maxJxHelper << std::endl;
//  std::cerr << "maxJyHelper = " << maxJyHelper << std::endl;
//  std::cerr << "maxJzHelper = " << maxJzHelper << std::endl;

  // here all the boundary conditions and the MPI happens
  particleExchange->exchange(particles);
}


#undef LOGLEVEL
#define LOGLEVEL 0
