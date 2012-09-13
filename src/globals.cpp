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

void Globals::initCommonParameters(BlockParameters &blockPars)
{
  std::cout << "Setting common parameters\n";

  blockPars.addParameter("end_time", &endTime);
  blockPars.addParameter("dt", &dt);
  blockPars.addArrayParameter("min_", domainMin);
  blockPars.addArrayParameter("max_", domainMax);
  blockPars.addArrayParameter("N_", globalGridSize);
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
  localGridMin = SIntVector(0,0);
  localGridMax = globalGridSize;
  for (int i=0; i<dimension; ++i) dx[i] = (domainMax[i]-domainMin[i]) / (real)globalGridSize[i];
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