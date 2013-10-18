/*
 * particle_boundary.hpp
 *
 * Created on: 16 Apr 2013
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2013 Holger Schmitz
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

#ifndef PARTICLE_BOUNDARY_HPP_
#define PARTICLE_BOUNDARY_HPP_

#include "particle_exchange.hpp"

/** An abstract particle boundary class.
 *
 * Any particle boundary must implement the apply method of this abstract class.
 */
class ParticleBoundary
{
  protected:
    /// The dimension in which the boundary condition is applied
    int dim;
    /// The direction either -1 for the lower bound or +1 for the higher bound
    int direction;
    /// The limiting value of the coordinate
    double limit;
  public:
    /** Constructor of the boundary class
     *
     * @param dim_ The dimension in which the boundary condition is applied
     * @param direction_ The direction either -1 for the lower bound or +1 for the higher bound
     */
    ParticleBoundary(int dim_, int direction_);

    /// Virtual destructor needed because of virtual methods.
    virtual ~ParticleBoundary() {}

    /** Apply the boundary condition
     *
     * @param particles A list of iterators over particles. The list may be modified.
     * Any particle removed will not be exchanged via MPI.
     */
    virtual void apply(ParticleExchange::PartList &particles) = 0;
};

typedef boost::shared_ptr<ParticleBoundary> pParticleBoundary;

/** Periodic boundary conditions for particles.
 */
class PeriodicParticleBoundary : public ParticleBoundary
{
  private:
    double shift;
  public:
    /** Constructor of the boundary class
     *
     * @param dim_ The dimension in which the boundary condition is applied
     * @param direction_ The direction either -1 for the lower bound or +1 for the higher bound
     */
    PeriodicParticleBoundary(int dim_, int direction_);

    /** Apply periodic boundary conditions to particles.
     *
     * Periodic boundary conditions are already implemented in ParticleExchange.
     * No need to do anything here.
     *
     * @param particles A list of iterators over particles. The list is not modified.
     */
    void apply(ParticleExchange::PartList &particles);
};


/** Open boundary conditions for particles.
 */
class OpenParticleBoundary : public ParticleBoundary
{
  public:
    /** Constructor of the boundary class
     *
     * @param dim_ The dimension in which the boundary condition is applied
     * @param direction_ The direction either -1 for the lower bound or +1 for the higher bound
     */
    OpenParticleBoundary(int dim_, int direction_) : ParticleBoundary(dim_, direction_) {}

    /** Apply open boundary conditions to particles.
     *
     *  At this point the particles are simply removed from the system by clearing
     *  the particle list.
     *
     *  @todo Correct the current of the removed particles in open boudary conditions.
     *  In the particle pusher the particles move by an undefined amount into the boundary
     *  region. The particles should be removed exactly on the boundary and no current
     *  contribution should be added when the particle is outside the box.
     *
     * @param particles A list of iterators over particles. The list is cleared.
     */
    void apply(ParticleExchange::PartList &particles) { particles.clear(); }
};


/** Reflecting boundary conditions for particles.
 */
class ReflectingParticleBoundary : public ParticleBoundary
{
  public:
    /** Constructor of the boundary class
     *
     * @param dim_ The dimension in which the boundary condition is applied
     * @param direction_ The direction either -1 for the lower bound or +1 for the higher bound
     */
    ReflectingParticleBoundary(int dim_, int direction_) : ParticleBoundary(dim_, direction_) {}

    /** Apply reflecting boundary conditions to particles.
     *
     *  At this point the particles are simply reflected and placed back inside the box.
     *
     *  @todo Correct the current of the reflected particles in reflecting boudary conditions.
     *  In the particle pusher the particles move by an undefined amount into the boundary
     *  region. The particles should be reflected exactly on the boundary and no current
     *  contribution should be added when the particle is outside the box.
     *
     * @param particles A list of iterators over particles. The particles are reflected
     * and the list is cleared.
     */
    void apply(ParticleExchange::PartList &particles);
};

#endif /* PARTICLE_BOUNDARY_HPP_ */
