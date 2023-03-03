/*
 * particle_exchange.cpp
 *
 * Created on: 15 Apr 2013
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

#include "particle_exchange.hpp"
#include "species.hpp"

#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>

#undef LOGLEVEL
#define LOGLEVEL 0

ParticleExchange::ParticleExchange(Species& species_) : species(species_) {}

void ParticleExchange::exchange(ParticleStorage &particles)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  auto &subdivision = species.getContext().getSubdivision();

  SRange localExtent = subdivision.getInnerExtent(species.getContext().getSize());

  SVector locMin = localExtent.getLo();
  SVector locMax = localExtent.getHi();
  SVector dx = species.getContext().getDx();

  locMin = locMin - 0.5*dx;
  locMax = locMax + 0.5*dx;

  for (int d=0; d<dimension; ++d)
  {
    //
    // Low boundary
    //

    for (ParticleStorage::iterator it = particles.begin();
        it != particles.end();
        ++it)
    {
      if (it->x[d] <  locMin[d]) listSend.push_front(it);
    }

    SCHNEK_TRACE_LOG(2,"Pushing " << listSend.size() << " particles < " << locMin[d])
    if (subdivision.isBoundLo(d))
      species.getBoundaryLo(d).apply(listSend);

    doExchange(particles, d, -1);

    //
    // High boundary
    //

    for (ParticleStorage::iterator it = particles.begin();
        it != particles.end();
        ++it)
    {
      if (it->x[d] >= locMax[d]) listSend.push_front(it);
    }

    SCHNEK_TRACE_LOG(2,"Pushing " << listSend.size() << " particles > " << locMax[d])
    if (subdivision.isBoundHi(d))
      species.getBoundaryHi(d).apply(listSend);

    doExchange(particles, d, +1);
  }
}


void ParticleExchange::doExchange(ParticleStorage &particles, int dim, int direction)
{
  auto &subdivision = species.getContext().getSubdivision();
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  bufferSend.makeBuffer(listSend);

  while (!listSend.empty())
  {
    particles.removeElement(listSend.front());
    listSend.pop_front();
  }

  SCHNEK_TRACE_LOG(3,"subdivision->exchangeData")
  subdivision.exchangeData(dim, direction, bufferSend.getBuffer(), bufferReceive.getBuffer());

  for (ParticleBuffer::iterator it = bufferReceive.begin(); it != bufferReceive.end(); ++it)
  {
    SCHNEK_TRACE_LOG(5,"particles.count "<<particles.getCount())
    SCHNEK_TRACE_LOG(5,"setting "<< it->x[0])
    Particle &p = particles.addElement();
    p = *it;
  }
  SCHNEK_TRACE_LOG(3,"doExchange received")
}


#undef LOGLEVEL
#define LOGLEVEL 0


