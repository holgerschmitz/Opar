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
  subdivision = Globals::instance().getSubdivision();

  SVector locMin = Globals::instance().getLocalDomainMin();
  SVector locMax = Globals::instance().getLocalDomainMax();

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
    if (Globals::instance().getSubdivision()->isBoundLo(d))
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
    if (Globals::instance().getSubdivision()->isBoundHi(d))
      species.getBoundaryHi(d).apply(listSend);

    doExchange(particles, d, +1);
  }
}

void ParticleExchange::doExchange(ParticleStorage &particles, int dim, int direction)
{

  bufferSend.makeBuffer(listSend);

  while (!listSend.empty())
  {

    SCHNEK_TRACE_LOG(2,"doExchange send " << listSend.front()->x[0] << " : "
        << Globals::instance().getLocalDomainMin()[0] << " "
        << Globals::instance().getLocalDomainMax()[0])
    particles.removeParticle(listSend.front());
    listSend.pop_front();
  }

  subdivision->exchangeData(dim, direction, bufferSend.getBuffer(), bufferReceive.getBuffer());

  for (ParticleBuffer::iterator it = bufferReceive.begin(); it != bufferReceive.end(); ++it)
  {
    particles.addParticle() = *it;
    SCHNEK_TRACE_LOG(2,"doExchange receive " << it->x[0] << " : "
        << Globals::instance().getLocalDomainMin()[0] << " "
        << Globals::instance().getLocalDomainMax()[0])
  }
}




