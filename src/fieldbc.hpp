/*
 * fieldbc.hpp
 *
 * Created on: 8 Oct 2012
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

#ifndef SCHNEK_FIELDBC_HPP_
#define SCHNEK_FIELDBC_HPP_

#include "types.hpp"

namespace schnek {

class FieldBC
{
  public:
    FieldBC() {}
    virtual ~FieldBC() {}
    virtual void applyEx(DataField &grid, int dim) = 0;
    virtual void applyEy(DataField &grid, int dim) = 0;
    virtual void applyEz(DataField &grid, int dim) = 0;
    virtual void applyBx(DataField &grid, int dim) = 0;
    virtual void applyBy(DataField &grid, int dim) = 0;
    virtual void applyBz(DataField &grid, int dim) = 0;
};
typedef boost::shared_ptr<FieldBC> pFieldBC;

class FieldPeriodicBC : public FieldBC
{
  public:
    FieldPeriodicBC() {}

    /** No need to do anything. This has been taken care of by the
     *  DomainSubdivision object
     */
    void applyEx(DataField &grid, int) {}
    void applyEy(DataField &grid, int) {}
    void applyEz(DataField &grid, int) {}
    void applyBx(DataField &grid, int) {}
    void applyBy(DataField &grid, int) {}
    void applyBz(DataField &grid, int) {}
};

} // namespace schnek

#endif // SCHNEK_FIELDBC_HPP_
