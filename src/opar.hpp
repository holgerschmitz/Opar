/*
 * opar.hpp
 *
 * Created on: 30 Jul 2012
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2012-2023 Holger Schmitz
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

#ifndef OPAR_HPP_
#define OPAR_HPP_

#include "../huerto/simulation/simulation_context.hpp"
#include "../huerto/constants.hpp"
#include "../huerto/electromagnetics/fieldsolver.hpp"

#include <schnek/variables.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <list>

class Species;

enum VarGroup { var_none, var_space, var_time, var_spacetime };

class OPar : public schnek::Block,
             public SimulationContext,
             public schnek::BlockContainer<FieldSolver>
{
  private:
    /**
     * The CLF factor \f$f_{\mathrm{CLF}}\f$
     *
     * A factor that is multiplied into the timestep #dt of the simulation.
     */
    double cflFactor;

    /**
     * If true, then don't do an initial half time step to bring E and B field in sync.
     */
    bool ignore_initial_time_stagger;

    std::list<Species*> species;

    schnek::pParametersGroup spaceVars;
    schnek::pParametersGroup timeVars;

    schnek::pBlockVariables blockVars;

    /// The parameter associated with t
    schnek::pParameter t_parameter;

  protected:
    void initParameters(schnek::BlockParameters &blockPars) override;
  public:
    void execute();
    void registerData() override;
    void init() override;
    schnek::pDependencyUpdater getUpdater(VarGroup gr);
    void addSpecies(Species *s);
    void addDiagnostic(schnek::DiagnosticInterface *d);

};


#endif /* OPAR_HPP_ */
