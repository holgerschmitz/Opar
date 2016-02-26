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
 *
 */

#ifndef FIELDS_HPP_
#define FIELDS_HPP_

#include "types.hpp"
#include "fieldbc.hpp"

#include <schnek/variables.hpp>

class FieldInitialiser;

using namespace schnek;

class Fields : public Block
{
  protected:
    pDataField pEx;
    pDataField pEy;
    pDataField pEz;
    pDataField pBx;
    pDataField pBy;
    pDataField pBz;

    std::list<FieldInitialiser*> initialisers;


    void initParameters(BlockParameters &blockPars);
    void registerData();
    void init();
  public:
    void addInitialiser(FieldInitialiser *init) { initialisers.push_back(init); }

    virtual void stepSchemeInit(double dt) = 0;
    virtual void stepScheme(double dt) = 0;


    pDataField getEx() { return pEx; }
    pDataField getEy() { return pEy; }
    pDataField getEz() { return pEz; }
    pDataField getBx() { return pBx; }
    pDataField getBy() { return pBy; }
    pDataField getBz() { return pBz; }
};

class EMFields : public Fields
{
  private:
    pDataField pJx;
    pDataField pJy;
    pDataField pJz;

    typedef boost::function<FieldBC*()> fieldBCFactoryFunction;
    std::map<std::string, fieldBCFactoryFunction> fieldBCFactories;

    schnek::Array<std::string, dimension> bcNamesLo, bcNamesHi;
    schnek::Array<pFieldBC, dimension> boundariesLo, boundariesHi;

  protected:
    void initParameters(BlockParameters &blockPars);
    void registerData();
    void init();
    // void postInit();

  public:
    void stepSchemeInit(double dt);
    void stepScheme(double dt);

    void writeAsTextFiles(int n);

  private:
    void stepD(double dt);
    void stepB(double dt);
#ifdef THREE_DIMENSIONAL
    void fdtdStepD(double dt, int i, int j, int k, SVector dx, double Jx, double Jy, double Jz);
    void fdtdStepB(double dt, int i, int j, int k, SVector dx, double Jx, double Jy, double Jz);
#endif
#ifdef TWO_DIMENSIONAL
    void fdtdStepD(double dt, int i, int j, SVector dx, double Jx, double Jy, double Jz);
    void fdtdStepB(double dt, int i, int j, SVector dx, double Jx, double Jy, double Jz);
#endif
#ifdef ONE_DIMENSIONAL
    void fdtdStepD(double dt, int i, SVector dx, double Jx, double Jy, double Jz);
    void fdtdStepB(double dt, int i, SVector dx, double Jx, double Jy, double Jz);
#endif
};


class ConstantFields : public Fields
{
  protected:
    void init();

  public:
    void stepSchemeInit(double dt) {}
    void stepScheme(double dt) {}
};

#endif /* FIELDS_HPP_ */
