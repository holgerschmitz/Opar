/*
 * particle_diagnostic.hpp
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

#ifndef PARTICLE_DIAGNOSTIC_HPP_
#define PARTICLE_DIAGNOSTIC_HPP_

#include "types.hpp"
#include "species.hpp"

#include <schnek/grid.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class ParticleDiagnostic 
  : public schnek::HDFGridDiagnostic<DataGrid1d, schnek::DeltaTimeDiagnostic>, 
    public SimulationEntity 
{
  private:
    typedef schnek::HDFGridDiagnostic<DataGrid1d, schnek::DeltaTimeDiagnostic> Super;
    Species *species;
    DataGrid1d data;

    schnek::Grid<long, 1, GridArgCheck> sizes;
    long totalCount;
    long localCount;
    long localStart;
  protected:
    typedef Super::IndexType IndexType;
    void write() override;
    void registerData() override;
    void init() override;
    IndexType getGlobalMin() override;
    IndexType getGlobalMax() override;
    std::string getLongFieldName();
    bool isDerived() override { return true; }
  public:
    virtual ~ParticleDiagnostic() {}
};


#endif /* PARTICLE_DIAGNOSTIC_HPP_ */
