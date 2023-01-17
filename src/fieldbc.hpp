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

#ifndef OPAR_FIELDBC_HPP_
#define OPAR_FIELDBC_HPP_

#include "types.hpp"
#include "numboundary.hpp"

class FieldBC
{
  public:
    typedef enum {lo = 1, hi = -1} Direction;
    FieldBC() {}
    virtual ~FieldBC() {}
    virtual void applyEx(DataField &grid, int dim, Direction dir) = 0;
    virtual void applyEy(DataField &grid, int dim, Direction dir) = 0;
    virtual void applyEz(DataField &grid, int dim, Direction dir) = 0;
    virtual void applyBx(DataField &grid, int dim, Direction dir) = 0;
    virtual void applyBy(DataField &grid, int dim, Direction dir) = 0;
    virtual void applyBz(DataField &grid, int dim, Direction dir) = 0;
};
typedef boost::shared_ptr<FieldBC> pFieldBC;

class FieldPeriodicBC : public FieldBC
{
  public:
    FieldPeriodicBC() {}

    /** No need to do anything. This has been taken care of by the
     *  DomainSubdivision object
     */
    void applyEx(DataField &, int, FieldBC::Direction) {}
    void applyEy(DataField &, int, FieldBC::Direction) {}
    void applyEz(DataField &, int, FieldBC::Direction) {}
    void applyBx(DataField &, int, FieldBC::Direction) {}
    void applyBy(DataField &, int, FieldBC::Direction) {}
    void applyBz(DataField &, int, FieldBC::Direction) {}
};

class FieldConductingBC : public FieldBC
{
  private:
    FieldDirichletBC dirichlet;
    FieldDirichletBC neumann;
  public:
    FieldConductingBC() {}

    void applyEx(DataField &grid, int dim, Direction dir)
    {
      if (dim==0) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyEy(DataField &grid, int dim, Direction dir)
    {
      if (dim==1) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyEz(DataField &grid, int dim, Direction dir)
    {
      if (dim==2) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyBx(DataField &grid, int dim, Direction dir)
    {
      if (dim==0) dirichlet.apply(grid, dim, dir);
      else neumann.apply(grid, dim, dir);
    }

    void applyBy(DataField &grid, int dim, Direction dir)
    {
      if (dim==1) dirichlet.apply(grid, dim, dir);
      else neumann.apply(grid, dim, dir);
    }

    void applyBz(DataField &grid, int dim, Direction dir)
    {
      if (dim==2) dirichlet.apply(grid, dim, dir);
      else neumann.apply(grid, dim, dir);
    }
};


class FieldSymmetryBC : public FieldBC
{
  private:
    FieldDirichletBC dirichlet;
    FieldDirichletBC neumann;
  public:
    FieldSymmetryBC() {}

    void applyEx(DataField &grid, int dim, Direction dir)
    {
      if (dim==0) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyEy(DataField &grid, int dim, Direction dir)
    {
      if (dim==1) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyEz(DataField &grid, int dim, Direction dir)
    {
      if (dim==2) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyBx(DataField &grid, int dim, Direction dir)
    {
      if (dim==0) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyBy(DataField &grid, int dim, Direction dir)
    {
      if (dim==1) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }

    void applyBz(DataField &grid, int dim, Direction dir)
    {
      if (dim==2) neumann.apply(grid, dim, dir);
      else dirichlet.apply(grid, dim, dir);
    }
};
#endif // OPAR_FIELDBC_HPP_
