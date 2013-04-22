/*
 * fields.cpp
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

#include "currents.hpp"
#include "fields.hpp"
#include "globals.hpp"
#include "opar.hpp"
#include "util.hpp"

#include <schnek/tools/fieldtools.hpp>
#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/functional/factory.hpp>
#include <fstream>


#undef LOGLEVEL
#define LOGLEVEL 1


#ifdef THREE_DIMENSIONAL

inline void Fields::fdtdStepD(double dt,
                              int i, int j, int k,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaEdx = 1.0;
  double kappaEdy = 1.0;
  double kappaEdz = 1.0;

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;

  (*pEx)(i,j,k) =
    dt*(
        ((*pBz)(i,j,k) - (*pBz)(i,j-1,k))/kappaEdy
      - ((*pBy)(i,j,k) - (*pBy)(i,j,k-1))/kappaEdz
      + Jx
    );

  (*pEy)(i,j,k) =
    dt*(
        ((*pBx)(i,j,k) - (*pBx)(i,j,k-1))/kappaEdz
      - ((*pBz)(i,j,k) - (*pBz)(i-1,j,k))/kappaEdx
      + Jy
    );

  (*pEz)(i,j,k) =
    dt*(
        ((*pBy)(i,j,k) - (*pBy)(i-1,j,k))/kappaEdx
      - ((*pBx)(i,j,k) - (*pBx)(i,j-1,k))/kappaEdy
      + Jz
    );
}

inline void Fields::fdtdStepB(double dt,
                              int i, int j, int k,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaHdx = 1.0;
  double kappaHdy = 1.0;
  double kappaHdz = 1.0;

//  double kappaHdx = (*pKappaHdx)(i)*dx;
//  double kappaHdy = (*pKappaHdy)(j)*dy;
//  double kappaHdz = (*pKappaHdz)(k)*dz;

  (*pBx)(i,j,k) = (*pBx)(i,j,k)
    + dt*(
        ((*pEy)(i,j,k+1) - (*pEy)(i,j,k))/kappaHdz
      - ((*pEz)(i,j+1,k) - (*pEz)(i,j,k))/kappaHdy
     + Jx
    );

  (*pBy)(i,j,k) = (*pBy)(i,j,k)
    + dt*(
        ((*pEz)(i+1,j,k) - (*pEz)(i,j,k))/kappaHdx
      - ((*pEx)(i,j,k+1) - (*pEx)(i,j,k))/kappaHdz
     + Jy
    );

  (*pBz)(i,j,k) = (*pBz)(i,j,k)
    + dt*(
        ((*pEx)(i,j+1,k) - (*pEx)(i,j,k))/kappaHdy
      - ((*pEy)(i+1,j,k) - (*pEy)(i,j,k))/kappaHdx
     + Jz
    );
}


void Fields::stepD(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
  Currents::instance().update();

  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      for (int k=low[2]+1; k<high[2]; ++k)
      {
        jx = (*this->pJx)(i,j,k);
        jy = (*this->pJy)(i,j,k);
        jz = (*this->pJz)(i,j,k);

        this->fdtdStepD(dt, i, j, k, dx, jx, jy, jz);
      }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pEx,i);
    sub->exchange(*pEy,i);
    sub->exchange(*pEz,i);

    boundariesLo[i]->applyEx(*pEx,i,FieldBC::lo);
    boundariesLo[i]->applyEy(*pEy,i,FieldBC::lo);
    boundariesLo[i]->applyEy(*pEz,i,FieldBC::lo);

    boundariesHi[i]->applyEx(*pEx,i,FieldBC::hi);
    boundariesHi[i]->applyEy(*pEy,i,FieldBC::hi);
    boundariesHi[i]->applyEy(*pEz,i,FieldBC::hi);
  }
}

void Fields::stepB(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
//  if (this->pMx != 0) sumMagCurrents();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      for (int k=low[2]; k<high[2]; ++k)
      {
//        if (this->pMx != 0)
//        {
//          jx = (*this->pMx)(i,j,k);
//          jy = (*this->pMy)(i,j,k);
//          jz = (*this->pMz)(i,j,k);
//        }

        this->fdtdStepB(dt, i, j, k, dx, dy, dz, jx, jy, jz);
      }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pBx,i);
    sub->exchange(*pBy,i);
    sub->exchange(*pBz,i);

    boundariesLo[i]->applyBx(*pBx,i,FieldBC::lo);
    boundariesLo[i]->applyBy(*pBy,i,FieldBC::lo);
    boundariesLo[i]->applyBy(*pBz,i,FieldBC::lo);

    boundariesHi[i]->applyBx(*pBx,i,FieldBC::hi);
    boundariesHi[i]->applyBy(*pBy,i,FieldBC::hi);
    boundariesHi[i]->applyBy(*pBz,i,FieldBC::hi);
  }
}

#endif

#ifdef TWO_DIMENSIONAL

inline void Fields::fdtdStepD(double dt,
                              int i, int j,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaEdx = dx[0];
  double kappaEdy = dx[1];

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;

  (*pEx)(i,j) +=
    dt*(
        ((*pBz)(i,j) - (*pBz)(i,j-1))/kappaEdy
      + Jx
    );

  (*pEy)(i,j) +=
    dt*(
      - ((*pBz)(i,j) - (*pBz)(i-1,j))/kappaEdx
      + Jy
    );

  (*pEz)(i,j) +=
    dt*(
        ((*pBy)(i,j) - (*pBy)(i-1,j))/kappaEdx
      - ((*pBx)(i,j) - (*pBx)(i,j-1))/kappaEdy
      + Jz
    );
}

inline void Fields::fdtdStepB(double dt,
                              int i, int j,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaHdx = dx[0];
  double kappaHdy = dx[1];

//  double kappaHdx = (*pKappaHdx)(i)*dx;
//  double kappaHdy = (*pKappaHdy)(j)*dy;
//  double kappaHdz = (*pKappaHdz)(k)*dz;

  (*pBx)(i,j) +=
    + dt*(
      - ((*pEz)(i,j+1) - (*pEz)(i,j))/kappaHdy
     + Jx
    );

  (*pBy)(i,j) +=
    + dt*(
        ((*pEz)(i+1,j) - (*pEz)(i,j))/kappaHdx
     + Jy
    );
  SCHNEK_TRACE_LOG(6, i << " " << j << " " << (*pBz)(i,j))

  (*pBz)(i,j) +=
    + dt*(
        ((*pEx)(i,j+1) - (*pEx)(i,j))/kappaHdy
      - ((*pEy)(i+1,j) - (*pEy)(i,j))/kappaHdx
     + Jz
    );
  SCHNEK_TRACE_LOG(6,(*pBz)(i,j))

}


void Fields::stepD(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
  Currents::instance().update();

  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      {
        jx = (*this->pJx)(i,j);
        jy = (*this->pJy)(i,j);
        jz = (*this->pJz)(i,j);

        this->fdtdStepD(dt, i, j, dx, jx, jy, jz);
      }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pEx,i);
    sub->exchange(*pEy,i);
    sub->exchange(*pEz,i);

    boundariesLo[i]->applyEx(*pEx,i,FieldBC::lo);
    boundariesLo[i]->applyEy(*pEy,i,FieldBC::lo);
    boundariesLo[i]->applyEy(*pEz,i,FieldBC::lo);

    boundariesHi[i]->applyEx(*pEx,i,FieldBC::hi);
    boundariesHi[i]->applyEy(*pEy,i,FieldBC::hi);
    boundariesHi[i]->applyEy(*pEz,i,FieldBC::hi);
  }
}

void Fields::stepB(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
//  if (this->pMx != 0) sumMagCurrents();

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      {
//        if (this->pMx != 0)
//        {
//          jx = (*this->pMx)(i,j,k);
//          jy = (*this->pMy)(i,j,k);
//          jz = (*this->pMz)(i,j,k);
//        }

        this->fdtdStepB(dt, i, j, dx, jx, jy, jz);
      }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pBx,i);
    sub->exchange(*pBy,i);
    sub->exchange(*pBz,i);

    boundariesLo[i]->applyBx(*pBx,i,FieldBC::lo);
    boundariesLo[i]->applyBy(*pBy,i,FieldBC::lo);
    boundariesLo[i]->applyBy(*pBz,i,FieldBC::lo);

    boundariesHi[i]->applyBx(*pBx,i,FieldBC::hi);
    boundariesHi[i]->applyBy(*pBy,i,FieldBC::hi);
    boundariesHi[i]->applyBy(*pBz,i,FieldBC::hi);
  }
}

#endif

#ifdef ONE_DIMENSIONAL

inline void Fields::fdtdStepD(double dt,
                              int i,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaEdx = dx[0];

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;

  (*pEx)(i) += dt*Jx;

  (*pEy)(i) +=
    dt*(
      - ((*pBz)(i) - (*pBz)(i-1))/kappaEdx
      + Jy
    );

  (*pEz)(i) +=
    dt*(
        ((*pBy)(i) - (*pBy)(i-1))/kappaEdx
      + Jz
    );
}

inline void Fields::fdtdStepB(double dt,
                              int i,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaHdx = dx[0];

//  double kappaHdx = (*pKappaHdx)(i)*dx;
//  double kappaHdy = (*pKappaHdy)(j)*dy;
//  double kappaHdz = (*pKappaHdz)(k)*dz;

  (*pBx)(i) += dt*Jx;

  (*pBy)(i) +=
    + dt*(
        ((*pEz)(i+1) - (*pEz)(i))/kappaHdx
     + Jy
    );
  SCHNEK_TRACE_LOG(6, i << " " << j << " " << (*pBz)(i))

  (*pBz)(i) +=
    + dt*(
      - ((*pEy)(i+1) - (*pEy)(i))/kappaHdx
     + Jz
    );
  SCHNEK_TRACE_LOG(6,(*pBz)(i))

}


void Fields::stepD(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
  Currents::instance().update();

  for (int i=low[0]+1; i<high[0]; ++i)
  {
    jx = (*this->pJx)(i);
    jy = (*this->pJy)(i);
    jz = (*this->pJz)(i);

    this->fdtdStepD(dt, i, dx, jx, jy, jz);
  }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pEx,i);
    sub->exchange(*pEy,i);
    sub->exchange(*pEz,i);

    boundariesLo[i]->applyEx(*pEx,i,FieldBC::lo);
    boundariesLo[i]->applyEy(*pEy,i,FieldBC::lo);
    boundariesLo[i]->applyEy(*pEz,i,FieldBC::lo);

    boundariesHi[i]->applyEx(*pEx,i,FieldBC::hi);
    boundariesHi[i]->applyEy(*pEy,i,FieldBC::hi);
    boundariesHi[i]->applyEy(*pEz,i,FieldBC::hi);
  }
}

void Fields::stepB(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
//  if (this->pMx != 0) sumMagCurrents();

  for (int i=low[0]; i<high[0]; ++i)
  {
//        if (this->pMx != 0)
//        {
//          jx = (*this->pMx)(i,j,k);
//          jy = (*this->pMy)(i,j,k);
//          jz = (*this->pMz)(i,j,k);
//        }

    this->fdtdStepB(dt, i, dx, jx, jy, jz);
  }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pBx,i);
    sub->exchange(*pBy,i);
    sub->exchange(*pBz,i);

    boundariesLo[i]->applyBx(*pBx,i,FieldBC::lo);
    boundariesLo[i]->applyBy(*pBy,i,FieldBC::lo);
    boundariesLo[i]->applyBy(*pBz,i,FieldBC::lo);

    boundariesHi[i]->applyBx(*pBx,i,FieldBC::hi);
    boundariesHi[i]->applyBy(*pBy,i,FieldBC::hi);
    boundariesHi[i]->applyBy(*pBz,i,FieldBC::hi);
  }
}

#endif


void Fields::stepSchemeInit(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  stepB(0.5*dt);

//  BOOST_FOREACH(Current *cur, this->currents)
//  {
//    cur->stepSchemeInit(dt);
//  }
//
//  BOOST_FOREACH(Current *cur, this->magCurrents)
//  {
//    cur->stepSchemeInit(dt);
//  }
}
void Fields::stepScheme(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  Currents::instance().update();

  stepD(dt);

//  BOOST_FOREACH(Current *cur, this->magCurrents)
//  {
//    cur->stepScheme(dt);
//  }

  stepB(dt);
}



void Fields::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  std::cout << "Fields::initParameters()" << std::endl;
  EParam = blockPars.addArrayParameter("E", EInit);
  BParam = blockPars.addArrayParameter("B", BInit);

  blockPars.addArrayParameter("boundary_min", bcNamesLo);
  blockPars.addArrayParameter("boundary_max", bcNamesHi);

  fieldBCFactories["periodic"] = boost::factory<FieldPeriodicBC*>();
  fieldBCFactories["conducting"] = boost::factory<FieldConductingBC*>();
  fieldBCFactories["symmetric"] = boost::factory<FieldSymmetryBC*>();
}

void Fields::registerData()
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
  std::cout << "Fields::registerData()" << std::endl;
  addData("Ex", pEx);
  addData("Ey", pEy);
  addData("Ez", pEz);
  addData("Bx", pBx);
  addData("By", pBy);
  addData("Bz", pBz);
  addData("Jx", pJx);
  addData("Jy", pJy);
  addData("Jz", pJz);
}

void Fields::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  dynamic_cast<OPar&>(*this->getParent()).addField(this);

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();
  SRange grange = Globals::instance().getLocalDomainRange();

  pEx = pDataField(new DataField(low, high, grange, SStagger(true,  false),2));
  pEy = pDataField(new DataField(low, high, grange, SStagger(false, true),2));
  pEz = pDataField(new DataField(low, high, grange, SStagger(false, false),2));

  pBx = pDataField(new DataField(low, high, grange, SStagger(false, true),2));
  pBy = pDataField(new DataField(low, high, grange, SStagger(true,  false),2));
  pBz = pDataField(new DataField(low, high, grange, SStagger(true,  true),2));

  pJx = pDataField(new DataField(low, high, grange, SStagger(true,  false),2));
  pJy = pDataField(new DataField(low, high, grange, SStagger(false, true),2));
  pJz = pDataField(new DataField(low, high, grange, SStagger(false, false),2));

  Currents::instance().setGlobalCurrent(pJx, pJy, pJz);
//  for (int i=0; i<dimension; ++i)
//  {
//    std::cout << "Field: "<< Globals::instance().getSubdivision()->getUniqueId() << " " << i <<
//        " (" << low[i] << " " << high[i] <<
//        ") (" << grange.getLo()[i] << " " << grange.getHi()[i] << ")" << std::endl;
//  }

  SVector &coords = Globals::instance().getX();
  pDependencyUpdater updater = Globals::instance().getUpdater(var_space);

  fill_field(*pEx, coords, EInit[0], *updater, EParam[0]);
  fill_field(*pEy, coords, EInit[1], *updater, EParam[1]);
  fill_field(*pEz, coords, EInit[2], *updater, EParam[2]);

  fill_field(*pBx, coords, BInit[0], *updater, BParam[0]);
  fill_field(*pBy, coords, BInit[1], *updater, BParam[1]);
  fill_field(*pBz, coords, BInit[2], *updater, BParam[2]);

  for (int i=0; i<dimension; ++i)
  {
    if (fieldBCFactories.count(bcNamesLo[i]) == 0)
      terminateSim("Unknown boundary condition: "+bcNamesLo[i]);

    if (fieldBCFactories.count(bcNamesHi[i]) == 0)
      terminateSim("Unknown boundary condition: "+bcNamesHi[i]);

    boundariesLo[i] = pFieldBC(fieldBCFactories[bcNamesLo[i]]());
    boundariesHi[i] = pFieldBC(fieldBCFactories[bcNamesHi[i]]());
  }

}

void Fields::writeAsTextFiles(int n)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  std::string num = boost::lexical_cast<std::string>(n);
  std::string rank = boost::lexical_cast<std::string>(
      Globals::instance().getSubdivision()->procnum());
  SIntVector min(pEx->getLo());
  SIntVector max(pEx->getHi());;

  std::ofstream exFile(("ex"+rank+"_"+num+".dat").c_str());
  for (int i=min[0]; i<=max[0]; ++i)
  {
    for (int j=min[1]; j<=max[1]; ++j)
      exFile << i << " " << j << " " << (*pEx)(i,j) << std::endl;
    exFile << std::endl;
  }
  exFile.close();

  std::ofstream eyFile(("ey"+rank+"_"+num+".dat").c_str());
  for (int i=min[0]; i<=max[0]; ++i)
  {
    for (int j=min[1]; j<=max[1]; ++j)
      eyFile << i << " " << j << " " << (*pEy)(i,j) << std::endl;
    eyFile << std::endl;
  }
  eyFile.close();

  std::ofstream ezFile(("ez"+rank+"_"+num+".dat").c_str());
  for (int i=min[0]; i<=max[0]; ++i)
  {
    for (int j=min[1]; j<=max[1]; ++j)
      ezFile << i << " " << j << " " << (*pEz)(i,j) << std::endl;
    ezFile << std::endl;
  }
  ezFile.close();

  std::ofstream bxFile(("bx"+rank+"_"+num+".dat").c_str());
  for (int i=min[0]; i<=max[0]; ++i)
  {
    for (int j=min[1]; j<=max[1]; ++j)
      bxFile << i << " " << j << " " << (*pBx)(i,j) << std::endl;
    bxFile << std::endl;
  }
  bxFile.close();

  std::ofstream byFile(("by"+rank+"_"+num+".dat").c_str());
  for (int i=min[0]; i<=max[0]; ++i)
  {
    for (int j=min[1]; j<=max[1]; ++j)
      byFile << i << " " << j << " " << (*pBy)(i,j) << std::endl;
    byFile << std::endl;
  }
  byFile.close();

  std::ofstream bzFile(("bz"+rank+"_"+num+".dat").c_str());
  for (int i=min[0]; i<=max[0]; ++i)
  {
    for (int j=min[1]; j<=max[1]; ++j)
      bzFile << i << " " << j << " " << (*pBz)(i,j) << std::endl;
    bzFile << std::endl;
  }
  bzFile.close();
}



