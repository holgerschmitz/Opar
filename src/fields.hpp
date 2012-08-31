/*
 * fields.hpp
 *
 * Created on: 8 Aug 2012
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

#ifndef FIELDS_HPP_
#define FIELDS_HPP_

#include "config.hpp"

#include <schnek/variables.hpp>

using namespace schnek;

class Fields;

class FieldsImplementation
{
  private:
    pDataGrid pEx;
    pDataGrid pEy;
    pDataGrid pEz;
    pDataGrid pBx;
    pDataGrid pBy;
    pDataGrid pBz;

//    typedef std::list<Current*> CurrentList;
//    CurrentList currents;
    friend class Fields;

    OParImplementation *sim;

  public:
    void stepSchemeInit(double dt);
    void stepScheme(double dt);

  private:
    void stepD(double dt);
    void stepB(double dt);
    void fdtdStepD(double dt, int i, int j, int k, double dx, double dy, double dz, double Jx, double Jy, double Jz);
    void fdtdStepB(double dt, int i, int j, int k, double dx, double dy, double dz, double Jx, double Jy, double Jz);
};

class Fields : public Block
{
  private:
    FieldsImplementation impl;
    double ExInit;
    double EyInit;
    double EzInit;
    double BxInit;
    double ByInit;
    double BzInit;
  protected:
    void initParameters(BlockParameters &blockPars);
    void registerData();
    void init();
};

#endif /* FIELDS_HPP_ */