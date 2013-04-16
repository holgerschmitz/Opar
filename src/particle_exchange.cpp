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

#include <boost/foreach.hpp>


ParticleExchange::ParticleExchange(Species& species_) : species(species_) {}

void ParticleExchange::exchange(ParticleStorage &particles)
{
  subdivision = Globals.instance().getSubdivision();

  SVector locMin = Globals.instance().getLocalDomainMin();
  SVector locMax = Globals.instance().getLocalDomainMax();

  for (int d=0; d<dimension; ++d)
  {
    //
    // Low boundary
    //

    for (ParticleStorage::iterator it = ParticleStorage.begin();
        it != ParticleStorage.end();
        ++it)
    {
      if ((*it).x[d] <  locMin[d]) listSend.push_front(it);
    }

    // species.getBoundaryLo(d).apply(listSend);

    doExchange(particles, -1);

    //
    // High boundary
    //

    for (ParticleStorage::iterator it = ParticleStorage.begin();
        it != ParticleStorage.end();
        ++it)
    {
      if ((*it).x[d] >= locMax[d]) listSend.push_front(it);
    }

    // species.getBoundaryHi(d).apply(listSend);

    doExchange(particles, +1);
  }
}

void ParticleExchange::doExchange(ParticleStorage &particles, int direction)
{

  bufferSend.makeBuffer(listSend);

  while (!listSend.empty())
  {
    particles.removeParticle(listSend.front());
    listSend.pop_front();
  }

  subdivision->exchangeData(dim, direction, bufferSend, bufferReceive);

  BOOST_FOREACH(Particle &p, bufferReceive)
  {
    particles.addParticle() = p;
  }
}




