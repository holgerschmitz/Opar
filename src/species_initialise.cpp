/*
 * species_initialise.cpp
 *
 * Created on: 20 May 2015
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


#include "species_initialise.hpp"
#include "defs.hpp"
#include "random.hpp"

void SpeciesInitialiser::preInit()
{
  dynamic_cast<Species&>(*this->getParent()).addInitialiser(this);
}

void MaxwellianInitialiser::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  blockPars.addParameter("ppc", &ppc, 100);
  blockPars.addParameter("densityCutoff", &densityCutoff, 0.0);

  densityParam = blockPars.addParameter("density", &density, 1.0);
  temperatureParam = blockPars.addArrayParameter("temperature", temperature, 0.0);
  driftParam = blockPars.addArrayParameter("drift", drift, 0.0);
}

void MaxwellianInitialiser::initialiseSpecies(Species &species, ParticleStorage &particles)
{
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


void TestParticle::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  blockPars.addParameter("weight", &weight, 1.0);

  blockPars.addArrayParameter("p", pos);
  blockPars.addArrayParameter("u", u, 0.0);
}

void TestParticle::initialiseSpecies(Species &species, ParticleStorage &particles)
{

  const SVector dx = Globals::instance().getDx();
  const SVector dMin = Globals::instance().getDomainMin();
  const SVector dMax = Globals::instance().getDomainMax();

  if ((dMin<=pos) && (pos<dMax))
  {
    Particle &p = particles.addParticle();
    p.x = pos;
    p.u = u;
    p.weight = weight;
  }
}
