/*
 * particle_exchange.hpp
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

#include "particles.hpp"

#include <schnek/grid/grid.hpp>
#include <schnek/grid/domainsubdivision.hpp>
#include <schnek/util/databuffer.hpp>

#include <list>

#ifndef PARTICLE_EXCHANGE_HPP_
#define PARTICLE_EXCHANGE_HPP_

class ParticleExchange
{
  public:
    typedef std::list<ParticleStorage::iterator> PartList;

  private:
    typedef schnek::DataBuffer<Particle> ParticleBuffer;

    PartList listSend;

    ParticleBuffer bufferSend;
    ParticleBuffer bufferReceive;

    Species& species;
    pSubdivision subdivision;

  public:
    ParticleExchange(Species& species_);

    void exchange(ParticleStorage &particles);
  private:
    void doExchange(ParticleStorage &particles, int direction);

};


#endif /* PARTICLE_EXCHANGE_HPP_ */
