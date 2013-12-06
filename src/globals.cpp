/*
 * globals.cpp
 *
 * Created on: 6 Sep 2012
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


#include "globals.hpp"
#include "random.hpp"

#include <schnek/grid/mpisubdivision.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/util/logger.hpp>

#undef LOGLEVEL
#define LOGLEVEL 5

void Globals::initCommonParameters(BlockParameters &blockPars)
{
  std::cout << "Setting common parameters\n";

  blockPars.addParameter("end_time", &endTime, 1.0);
  blockPars.addParameter("dt", &dt, 0.005);
  blockPars.addArrayParameter("min_", domainMin, 0.0);
  blockPars.addArrayParameter("max_", domainMax, 1.0);
  blockPars.addArrayParameter("N_", globalGridSize, 100);
}

void Globals::initGlobalParameters(BlockParameters &blockPars)
{
  x_parameters = blockPars.addArrayParameter("", x, BlockParameters::readonly);
  t_parameter  = blockPars.addParameter("t", &t, BlockParameters::readonly);

  spaceVars->addArray(x_parameters);
  timeVars->add(t_parameter);
}

void Globals::setup(VariableStorage &vars)
{
  blockVars = vars.getRootBlock();
  //depMap = pDependencyMap(new DependencyMap(vars.getRootBlock()));
  spaceVars = pParametersGroup(new ParametersGroup());
  timeVars  = pParametersGroup(new ParametersGroup());
}

void Globals::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  //for (int i=0; i<dimension; ++i) globalGridSize[i] += 3;
#ifdef HAVE_MPI
  subdivision = pSubdivision(new MPICartSubdivision<DataField>());
#else
  subdivision = pSubdivision(new SerialSubdivision<DataField>());
#endif
  subdivision->init(globalGridSize, 2);
  localGridMin = subdivision->getLo();
  localGridMax = subdivision->getHi();

  localInnerGridMin = subdivision->getInnerLo();
  localInnerGridMax = subdivision->getInnerHi();

  for (int i=0; i<dimension; ++i)
  {
    dx[i] = (domainMax[i]-domainMin[i]) / (real)globalGridSize[i]; // account for ghost cells
    localDomainMin[i] = domainMin[i] + localInnerGridMin[i]*dx[i];
    localDomainMax[i] = domainMin[i] + (localInnerGridMax[i]+1)*dx[i];

    std::cout << "II Domain: "<< subdivision->getUniqueId() << " " << i << " " << localGridMin[i] << " " << localGridMax[i] << std::endl;
    std::cout << "II Inner Domain: "<< subdivision->getUniqueId() << " " << i << " " << localInnerGridMin[i] << " " << localInnerGridMax[i] << std::endl;
    std::cout << "II Extent: "<< subdivision->getUniqueId() << " " << i << " " << localDomainMin[i] << " " << localDomainMax[i] << std::endl;
  }

  DiagnosticManager::instance().setTimeCounter(&t_count);
  DiagnosticManager::instance().setMaster(subdivision->master());

  SCHNEK_TRACE_LOG(2,"Setting Rank: calling getUniqueId " << subdivision->getUniqueId())
  DiagnosticManager::instance().setRank(subdivision->getUniqueId());

  Random::seed(subdivision->getUniqueId());
}

bool Globals::stepTime()
{
  ++t_count;
  return (t += dt) <= endTime;
}

pDependencyUpdater Globals::getUpdater(VarGroup gr)
{
  pDependencyMap depMap(new DependencyMap(blockVars));
  pDependencyUpdater updater(new DependencyUpdater(depMap));
  switch (gr)
  {
    case var_space:
      updater->addIndependentArray(x_parameters);
      break;
    case var_time:
      updater->addIndependent(t_parameter);
      break;
    case var_spacetime:
      updater->addIndependentArray(x_parameters);
      updater->addIndependent(t_parameter);
      break;
    case var_none:
    default:
      break;
  }
  return updater;
}
