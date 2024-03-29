/*
 * species.hpp
 *
 * Created on: 15 Nov 2012
 * Author: hschmitz
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

#ifndef SPECIES_HPP_
#define SPECIES_HPP_

#include "types.hpp"
#include "particles.hpp"
#include "particle_boundary.hpp"
#include "weighting.hpp"

#include "../huerto/simulation/simulation_context.hpp"

#include <schnek/variables.hpp>

class Species : public schnek::Block, public SimulationEntity
{
  private:
    std::string name;
    ParticleStorage particles;

    pDataField pJx;
    pDataField pJy;
    pDataField pJz;

    DataGrid jxHelper;
    DataGrid jyHelper;
    DataGrid jzHelper;
    pDataField pEx;
    pDataField pEy;
    pDataField pEz;
    pDataField pBx;
    pDataField pBy;
    pDataField pBz;


    double charge;
    double mass;
    double density;
    // double minDensity;
    PVector temperature;
    PVector drift;

    double densityCutoff;

    int ppc;

    schnek::pParameter densityParam;
    PParameterVector temperatureParam;
    PParameterVector driftParam;

    typedef Weighting::WeightingCoefficients WeightingCoefficients;
    WeightingCoefficients gx, hx;

    typedef boost::function<ParticleBoundary*(int&, int&, SimulationContext&)> particleBCFactoryFunction;
    std::map<std::string, particleBCFactoryFunction> particleBCFactories;

    schnek::Array<std::string, dimension> bcNamesLo, bcNamesHi;
    schnek::Array<pParticleBoundary, dimension> boundariesLo, boundariesHi;

    pParticleExchange particleExchange;

    void initParticles();

    static std::map<std::string, Species*> allSpecies;

  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void registerData();
    void preInit();
    void init();
    // void postInit();
  public:
    void pushParticles(double dt);
    Particle &addParticle();
    void removeParticle(const ParticleStorage::iterator &p);

    ParticleBoundary &getBoundaryLo(int dim) { return *(boundariesLo[dim]); }
    ParticleBoundary &getBoundaryHi(int dim) { return *(boundariesHi[dim]); }

    static Species* getSpecies(std::string name);
    ParticleStorage& getStorage() { return particles; }

};

#endif // SPECIES_HPP_ 
