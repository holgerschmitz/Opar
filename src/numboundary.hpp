/*
 * numboundary.hpp
 *
 * Created on: 15 Oct 2012
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2012 Holger Schmitz
 *
 * This file is part of Schnek.
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

#ifndef OPAR_NUMBOUNDARY_HPP_
#define OPAR_NUMBOUNDARY_HPP_

#include "types.hpp"

class NumBoundary
{
  private:
    int lo, dim, dir;

  protected:
    real& get(int pos, DataField &field, SIntVector &gpos)
    {
      gpos[dim] = lo + dir*pos;
      return field[gpos];
    }

    void setParameters(int lo_, int dim_, int dir_)
    {
      lo = lo_;
      dim = dim_;
      dir = dir_;
    }

    int getDim() { return dim; }
    int getDir() { return dir; }
};

template<int ghostCells>
class BoundDirichletImpl : public NumBoundary
{
  protected:
    void applyRow(DataField &field, SIntVector &gpos)
    {
      for (int i=0; i<ghostCells; ++i)
        get(-1-i,field,gpos) = -get(i,field,gpos);
    }

    void applyRowStaggerLow(DataField &field, SIntVector &gpos)
    {
      get(-1,field,gpos) = 0.0;
      for (int i=1; i<ghostCells; ++i)
        get(-1-i,field,gpos) = -get(-1+i,field,gpos);
    }

    void applyRowStaggerHigh(DataField &field, SIntVector &gpos)
    {
      get(0,field,gpos) = 0.0;
      for (int i=0; i<ghostCells; ++i)
        get(-1-i,field,gpos) = -get(1+i,field,gpos);
    }
};

template<int ghostCells>
class BoundNeumannImpl : public NumBoundary
{
  protected:
    void applyRow(DataField &field, SIntVector &gpos)
    {
      for (int i=0; i<ghostCells; ++i)
        get(-1-i,field,gpos) = get(i,field,gpos);
    }

    void applyRowStaggerLow(DataField &field, SIntVector &gpos)
    {
      for (int i=1; i<ghostCells; ++i)
        get(-1-i,field,gpos) = get(-1+i,field,gpos);
    }

    void applyRowStaggerHigh(DataField &field, SIntVector &gpos)
    {
      for (int i=0; i<ghostCells; ++i)
        get(-1-i,field,gpos) = get(1+i,field,gpos);
    }
};

template<template<int> class Implementation, int ghostCells>
class FieldBoundary : public Implementation<ghostCells>
{
  private:
    bool advancePos(const SIntVector &Lo, const SIntVector &Hi, const SIntVector &itdims, SIntVector &pos)
    {
      int ind=dimension-2;
      while (ind>=0)
      {
        int dim = itdims[ind];
        if (++pos[dim] <= Hi[dim]) return true; // done
        pos[dim] = Lo[dim];
        --ind;
      }
      return false;
    }

  public:
    void apply(DataField &field, int dim, int dir)
    {
      int lo;
      if (dir>0) lo = ghostCells;
      else lo = field.getHi()[dim]-ghostCells;

      this->setParameters(lo,dim,dir);

      SIntVector gpos;
      SIntVector itdims;
      SIntVector Lo, Hi;

      for (int i=0, k=0; i<dimension; ++i)
        if (i != this->getDim()) itdims[k++] = i;

      for (int i=0; i<dimension-1; ++i)
      {
        Lo[i] = field.getLo()[itdims[i]];
        Hi[i] = field.getHi()[itdims[i]];
        gpos[itdims[i]] = Lo[i];
      }

      bool loop=true;
      if (!field.getStagger(this->getDim()))
      {
        while (loop)
        {
          this->applyRow(field, gpos);
          loop = advancePos(Lo, Hi, itdims, gpos);
        }
      }
      else if (dir>0)
      {
        while (loop)
        {
          this->applyRowStaggerLow(field, gpos);
          loop = advancePos(Lo, Hi, itdims, gpos);
        }
      }
      else
      {
        while (loop)
        {
          this->applyRowStaggerHigh(field, gpos);
          loop = advancePos(Lo, Hi, itdims, gpos);
        }
      }
    }
};

typedef FieldBoundary<BoundDirichletImpl, ghostCells> FieldDirichletBC;
typedef FieldBoundary<BoundNeumannImpl, ghostCells> FieldNeumannBC;


#endif // OPAR_NUMBOUNDARY_HPP_
