/*
 * species_initialise.hpp
 *
 * Created on: 20 May 2015
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

#ifndef SRC_SPECIES_INITIALISE_HPP_
#define SRC_SPECIES_INITIALISE_HPP_

#include "types.hpp"
#include "particles.hpp"
#include "species.hpp"

#include <schnek/variables.hpp>

class SpeciesInitialiser : public schnek::Block
{
  protected:
    void preInit();
  public:
    virtual void initialiseSpecies(Species &species, ParticleStorage &particles) = 0;
};

class MaxwellianInitialiser : public SpeciesInitialiser
{
  private:
    double density;
    PVector temperature;
    PVector drift;

    double densityCutoff;
    int ppc;

    pParameter densityParam;
    PParameterVector temperatureParam;
    PParameterVector driftParam;
  protected:
    void initParameters(BlockParameters &blockPars);
  public:
    void initialiseSpecies(Species &species, ParticleStorage &particles);
};

class TestParticle : public SpeciesInitialiser
{
  private:
    SVector pos;
    PVector u;
    double weight;
  protected:
    void initParameters(BlockParameters &blockPars);
  public:
    void initialiseSpecies(Species &species, ParticleStorage &particles);
};


#endif /* SRC_SPECIES_INITIALISE_HPP_ */
