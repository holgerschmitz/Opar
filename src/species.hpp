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

#include <schnek/variables.hpp>

class Species : public Block
{
  private:
    ParticleStorage particles;

    pDataField pJx;
    pDataField pJy;
    pDataField pJz;

    DataField jxHelper;
    DataField jyHelper;
    DataField jzHelper;
    pDataField pEx;
    pDataField pEy;
    pDataField pEz;
    pDataField pBx;
    pDataField pBy;
    pDataField pBz;


    double charge;
    double mass;
    double density;
    PVector temperature;
    PVector drift;

    int ppc;

    pParameter densityParam;
    PParameterVector temperatureParam;
    PParameterVector driftParam;

    typedef Weighting::WeightingCoefficients WeightingCoefficients;
    WeightingCoefficients gx, hx;

    typedef boost::function<ParticleBoundary*(int, int)> particleBCFactoryFunction;
    std::map<std::string, particleBCFactoryFunction> particleBCFactories;

    schnek::Array<std::string, dimension> bcNamesLo, bcNamesHi;
    schnek::Array<pParticleBoundary, dimension> boundariesLo, boundariesHi;

    pParticleExchange particleExchange;

    void initParticles();

  protected:
    void initParameters(BlockParameters &blockPars);
    void registerData();
    void init();
    // void postInit();
  public:
    void pushParticles(double dt);
    Particle &addParticle();
    void removeParticle(const ParticleStorage::iterator &p);

    ParticleBoundary &getBoundaryLo(int dim) { return *(boundariesLo[dim]); }
    ParticleBoundary &getBoundaryHi(int dim) { return *(boundariesHi[dim]); }

};

#endif // SPECIES_HPP_ 
