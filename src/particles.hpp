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

#include <list>
#include <iterator>

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

    Particle()
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
      x = P.x;
      u = P.u;
      weight = P.weight;
      index = P.index;
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
      weight = weight;
    }

    void setValues(SVector &x_, PVector &u_)
    {
      x = x_;
      u = u_;
    }
};

class ParticleStorage
{
  private:
    struct DataBlock
    {
        Particle *data;
        long count;
        DataBlock();
        DataBlock(const DataBlock &block);
        void free();
        Particle &addParticle();
    };
    typedef std::list<DataBlock> BlockList;
    typedef BlockList::iterator BlockIterator;
    BlockList blocks;
    BlockIterator freeBlock;
  public:
    ParticleStorage()
    {
      freeBlock = blocks.end();
    }

    class iterator : public std::iterator<std::forward_iterator_tag, Particle>
    {
      private:
        friend class ParticleStorage;
        BlockIterator blockIter;
        long pos;

        iterator(BlockIterator blockIter_, long pos_ = 0) :
          blockIter(blockIter_), pos(pos_)
        {
        }
      public:
        iterator()
        {
        }
        iterator(const iterator &it)
        {
          blockIter = it.blockIter;
          pos = it.pos;
        }

        iterator& operator++()
        {
          if (++pos >= blockIter->count)
          {
            ++blockIter;
            pos = 0;
          }
          return *this;
        }

        iterator operator++(int)
        {
          iterator tmp(*this);
          operator++();
          return tmp;
        }

        bool operator==(const iterator& rhs)
        {
          return (blockIter == rhs.blockIter) && (pos == rhs.pos);
        }

        bool operator!=(const iterator& rhs)
        {
          return (blockIter != rhs.blockIter) || (pos != rhs.pos);
        }

        Particle& operator*()
        {
          return blockIter->block[pos];
        }
    };

    iterator begin()
    {
      return iterator(BlockList.begin());
    }
    iterator end()
    {
      return iterator(BlockList.end());
    }

    Particle &addParticle();

    /** Remove a particle from the storage.
     *
     * After deletion the iterator will point to the position after the deleted particle.
     * This is in line with STL behaviour
     */
    iterator removeParticle(const iterator&);
};

#endif /* PARTICLES_HPP_ */
