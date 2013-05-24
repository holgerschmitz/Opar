/*
 * particles.cpp
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
 */

#include "particles.hpp"

#include <algorithm>

static const long STORAGE_BLOCK_SIZE = 1000000;

ParticleStorage::DataBlock::DataBlock()
{
  data = new Particle[STORAGE_BLOCK_SIZE];
}

ParticleStorage::DataBlock::DataBlock(const DataBlock &block) : data(block.data), count(block.count) {}

void ParticleStorage::DataBlock::free()
{
  delete[] data;
}

Particle &ParticleStorage::DataBlock::addParticle()
{
  return data[count++];
}

Particle &ParticleStorage::addParticle()
{
  // find a free block or create one
  if ((STORAGE_BLOCK_SIZE == freeBlock->count) || (freeBlock == blocks.end()))
  {
    bool foundBlock = false;
    freeBlock = blocks.begin();
    while(!foundBlock && (freeBlock != blocks.end()))
    {
      foundBlock = (STORAGE_BLOCK_SIZE > freeBlock->count);
      ++freeBlock;
    }
    if (!foundBlock)
    {
      freeBlock = blocks.insert(blocks.begin(), DataBlock());
    }
  }
  return freeBlock->addParticle();

}

ParticleStorage::iterator ParticleStorage::removeParticle(const iterator &it_)
{
  iterator it = it_;
  DataBlock &dblock = *(it.blockIter);
  Particle *data = dblock.data;
  long pos = it.pos;
  long last = dblock.count - 1;

  if (last>pos) std::swap(data[pos], data[last]);

  // if no more particles in the block then delete the block
  if (0 == --dblock.count)
  {
    // manually free the data
    dblock.free();
    it.blockIter = blocks.erase(it.blockIter);
    it.pos = 0;
  }

  // for performance, see if there are more free spaces in this block than in the
  // freeBlock, and update freeBlock accordingly
  if ((it.blockIter != blocks.end()) && (it.blockIter->count < freeBlock->count))
  {
    freeBlock = it.blockIter;
  }
  return it;
}

long ParticleStorage::getCount() const
{
  long count = 0;
  for (BlockConstIterator b = blocks.begin(); b!=blocks.end(); ++b)
    count += b->count;
  return count;
}
