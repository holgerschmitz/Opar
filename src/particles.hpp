/*
 * particles.hpp
 *
 * Created on: 13 Nov 2012
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

#ifndef PARTICLES_HPP_
#define PARTICLES_HPP_

#include "types.hpp"

#include "../huerto/storage/small_object_storage.hpp"

/** @Particle
 *  A low impact particle class for holding all the particle data
 */
class Particle
{
  private:
    static long nextIndex;
    long index;
  public:
    /// This is the position
    SVector x;
    /// This is the relativistic four velocity
    PVector u;
    real weight;

    Particle() : index(0), weight(1.0)
    {
    }
    Particle(const SVector &x_, const PVector &u_, real weight_) :
      index(++nextIndex), x(x_), u(u_), weight(weight_)
    {
    }
    Particle(const Particle &P) :
      index(P.index), x(P.x), u(P.u), weight(P.weight)
    {
    }

    Particle &operator=(const Particle &P)
    {
      SVector tmpx = P.x;
      x = tmpx; // P.x;
      u = P.u;
      weight = P.weight;
      index = P.index;
      return *this;
    }

    void copyFrom(const Particle &P)
    {
      x = P.x;
      u = P.u;
      weight = P.weight;
      index = ++nextIndex;
    }

    long getIndex()
    {
      return index;
    }

    void setValues(SVector &x_, PVector &u_, real weight_)
    {
      x = x_;
      u = u_;
      weight = weight_;
    }

    void setValues(SVector &x_, PVector &u_)
    {
      x = x_;
      u = u_;
    }
};

typedef SmallObjectStorage<Particle> ParticleStorage;

#endif /* PARTICLES_HPP_ */
