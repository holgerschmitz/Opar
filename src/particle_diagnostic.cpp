/*
 * particle_diagnostic.cpp
 *
 * Created on: 7 May 2013
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

#include "particle_diagnostic.hpp"

#include "particles.hpp"
#include "globals.hpp"

#include <schnek/util/logger.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#undef LOGLEVEL
#define LOGLEVEL 0

void ParticleDiagnostic::write()
{
  localCount = species->getStorage().getCount();
#ifdef HAVE_MPI
  int rank;
  int procCount = Globals::instance().getSubdivision()->procCount();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Allgather(&localCount, 1, MPI_LONG, sizes.getRawData(), 1, MPI_LONG, MPI_COMM_WORLD);
  localStart = 0;
  for (int i=0; i<rank; ++i) localStart += sizes(i);
  totalCount = 0;
  for (int i=0; i<procCount; ++i) totalCount += sizes(i);
#else
  sizes(0) = localCount;
  totalCount = localCount;
  localStart = 0;
#endif
  data->resize(IndexType(localStart), IndexType(localStart+localCount-1));
  container.grid = &(*data);
  container.global_min = IndexType(0);
  container.global_max = IndexType(totalCount-1);
  container.local_min = IndexType(localStart);
  container.local_max = IndexType(localStart + localCount -1);

  SCHNEK_TRACE_LOG(2,"ParticleDiagnostic::write "<< 0<<" " << totalCount-1<<" "<< localStart<<" "<<localStart + localCount -1)

  if (totalCount==0) return;

  std::string coordNames[] = {"x","y","z"};
  ParticleStorage particles = species->getStorage();

  for (int d=0; d<dimension; ++d)
  {
    long pos = localStart;
    for (ParticleStorage::iterator it=particles.begin(); it!=particles.end(); ++it)
    {
      (*data)(pos++) = it->x[d];
    }
    output.setBlockName(coordNames[d]);
    output.writeGrid(container);
  }

  for (int d=0; d<3; ++d)
  {
    long pos = localStart;
    for (ParticleStorage::iterator it=particles.begin(); it!=particles.end(); ++it)
    {
      (*data)(pos++) = it->u[d];
    }
    output.setBlockName("u"+coordNames[d]);
    output.writeGrid(container);
  }

}

void ParticleDiagnostic::registerData()
{
  data = pDataField1d(new DataField1d());
  addData(getLongFieldName(), data);
}

void ParticleDiagnostic::init()
{
  Super::init();
  species = Species::getSpecies(getFieldName());
  int procCount = Globals::instance().getSubdivision()->procCount();
  sizes.resize(procCount);
}

ParticleDiagnostic::IndexType ParticleDiagnostic::getGlobalMin()
{
  return IndexType(0);
}

ParticleDiagnostic::IndexType ParticleDiagnostic::getGlobalMax()
{
  return IndexType(totalCount-1);
}

std::string ParticleDiagnostic::getLongFieldName()
{
  return "particle_data_" + getFieldName();
}


