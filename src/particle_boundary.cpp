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

ParticleBoundary::ParticleBoundary(int dim_, int direction_, SimulationContext &context) :
    dim(dim_), direction(direction_), context(context)
{
  if (direction > 0)
    limit = 0.0;
  else
    limit = context.getSize()[dim];
}

PeriodicParticleBoundary::PeriodicParticleBoundary(int dim_, int direction_, SimulationContext &context)
  : ParticleBoundary(dim_, direction_, context)
{
  shift = context.getSize()[dim];
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
      SCHNEK_TRACE_LOG(5,"Shifting particle from "<<(*it)->x[dim]<<" by "<< -shift)
      (*it)->x[dim] -= shift;
    }
  }
  else
  {
    for (ParticleExchange::PartList::iterator it = particles.begin();
        it != particles.end(); ++it)
    {
      SCHNEK_TRACE_LOG(5,"Shifting particle from "<<(*it)->x[dim]<<" by "<< shift)
      (*it)->x[dim] += shift;
    }
  }

}

void ReflectingParticleBoundary::apply(ParticleExchange::PartList &particles)
{
  while (!particles.empty())
  {
    Particle &p = *(particles.front());

    p.x[dim] = 2.0 * limit - p.x[dim];
    p.u[dim] = -p.u[dim];

    particles.pop_front();
  }
}

