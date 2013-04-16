/*
 * globals.hpp
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

#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include "types.hpp"

#include <schnek/util.hpp>
#include <schnek/variables.hpp>

using namespace schnek;

class CommonBlock;

enum VarGroup { var_none, var_space, var_time, var_spacetime };

class Globals : public Singleton<Globals>
{
  public:
    typedef boost::shared_ptr<schnek::DomainSubdivision<DataField> > pSubdivision;
  private:
    /// The total duration of the simulation
    double endTime;
    /// The time step
    double dt;

    /// The physical size of the simulation domain
    SVector domainMin;
    SVector domainMax;

    /// The physical size of the local simulation domain
    SVector localDomainMin;
    SVector localDomainMax;

    /// The global grid size of the simulation
    SIntVector globalGridSize;

    // derived quantities
    SIntVector localGridMin;
    SIntVector localGridMax;
    SVector dx;

    /// Global space variable
    SVector x;
    /// The parameters associated with x
    SParameterVector x_parameters;

    /// Global time variable
    double t;
    /// Global time counter
    int t_count;

    /// The parameter associated with t
    pParameter t_parameter;

    pSubdivision subdivision;

    pParametersGroup spaceVars;
    pParametersGroup timeVars;

    pBlockVariables blockVars;
    pDependencyMap depMap;

    friend class Singleton<Globals>;
    friend class CreateUsingNew<Globals>;

    Globals() : t(0.0), t_count(0), spaceVars(new ParametersGroup()), timeVars(new ParametersGroup()) {}
  public:
    void setup(VariableStorage &vars);
    void initGlobalParameters(BlockParameters &blockPars);
    void initCommonParameters(BlockParameters &blockPars);
    void init();

    double getEndTime() const { return endTime; }
    double getDt() const { return dt; }
    bool stepTime();

    const SVector &getDomainMin() const { return domainMin; }
    const SVector &getDomainMax() const { return domainMax; }

    const SVector &getLocalDomainMin() const { return localDomainMin; }
    const SVector &getLocalDomainMax() const { return localDomainMax; }

    SRange getDomainRange() const { return  SRange(domainMin, domainMax); }
    SRange getLocalDomainRange() const { return  SRange(localDomainMin, localDomainMax); }

    const SIntVector &getGlobalGridSize() const { return globalGridSize; }

    const SIntVector &getLocalGridMin() const { return localGridMin; }
    const SIntVector &getLocalGridMax() const { return localGridMax; }
    const SVector &getDx() const { return dx; }

    const pParametersGroup getSpaceVars() const { return spaceVars; }
    const pParametersGroup getTimeVars() const { return timeVars; }

    SVector &getX() { return x; }
    double &getT() { return t; }

    pDependencyUpdater getUpdater(VarGroup gr=var_none);
    pSubdivision getSubdivision() { return subdivision; }
};

#endif // GLOBALS_HPP_ 
