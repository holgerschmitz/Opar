/*
 * particle_boundary.cpp
 *
 * Created on: 16 Apr 2013
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2013 Holger Schmitz
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

#include "particle_boundary.hpp"

#include <schnek/util/logger.hpp>

#undef LOGLEVEL
#define LOGLEVEL 0

ParticleBoundary::ParticleBoundary(int dim_, int direction_) :
    dim(dim_), direction(direction_)
{
  if (direction > 0)
    limit = 2*(Globals::instance().getGlobalGridSize())[dim];
  else
    limit = 0;
}

PeriodicParticleBoundary::PeriodicParticleBoundary(int dim_, int direction_)
  : ParticleBoundary(dim_, direction_)
{
  shift = 2*(Globals::instance().getGlobalGridSize())[dim];
}

void PeriodicParticleBoundary::apply(ParticleExchange::PartList &particles)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  // wrapping the coordinates. We do not remove the particle from the list
  // because it still needs to be transferred to the right node.

  if (direction > 0)
  {
    for (ParticleExchange::PartList::iterator it = particles.begin();
        it != particles.end(); ++it)
    {
      SCHNEK_TRACE_LOG(5,"Shifting particle from "<<(*it)->xi[dim]<<" by "<< -shift)
      (*it)->xi[dim] -= shift;
    }
  }
  else
  {
    for (ParticleExchange::PartList::iterator it = particles.begin();
        it != particles.end(); ++it)
    {
      SCHNEK_TRACE_LOG(5,"Shifting particle from "<<(*it)->xi[dim]<<" by "<< shift)
      (*it)->xi[dim] += shift;
    }
  }

}

void ReflectingParticleBoundary::apply(ParticleExchange::PartList &particles)
{
  while (!particles.empty())
  {
    Particle &p = *(particles.front());

    p.xi[dim] = 2 * limit - p.xi[dim] - 1;
    p.xf[dim] = 1.0 - p.xf[dim];
    p.u[dim] = -p.u[dim];

    particles.pop_front();
  }
}

