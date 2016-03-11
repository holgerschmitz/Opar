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

#include "config.hpp"
#include "currents.hpp"
#include "fields.hpp"
#include "globals.hpp"
#include "opar.hpp"
#include "util.hpp"

#include "constants.hpp"
#include "field_initialise.hpp"

#include <schnek/tools/fieldtools.hpp>
#include <schnek/util/logger.hpp>
#include <schnek/tools/literature.hpp>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/functional/factory.hpp>
#include <fstream>


#undef LOGLEVEL
#define LOGLEVEL 0

void Fields::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
}


void Fields::registerData()
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
  addData("Ex", pEx);
  addData("Ey", pEy);
  addData("Ez", pEz);
  addData("Bx", pBx);
  addData("By", pBy);
  addData("Bz", pBz);
}



void Fields::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  BOOST_FOREACH(FieldInitialiser *init, initialisers)
  {
    init->initialiseFields(*this);
  }
}

#ifdef THREE_DIMENSIONAL

inline void EMFields::fdtdStepD(double dt,
                              int i, int j, int k,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaEdx = dx[0]/clight2;
  double kappaEdy = dx[1]/clight2;
  double kappaEdz = dx[2]/clight2;

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;

  (*pEx)(i,j,k) +=
    dt*(
        ((*pBz)(i,j,k) - (*pBz)(i,j-1,k))/kappaEdy
      - ((*pBy)(i,j,k) - (*pBy)(i,j,k-1))/kappaEdz
      - eps_0_inv*Jx
    );

  (*pEy)(i,j,k) +=
    dt*(
        ((*pBx)(i,j,k) - (*pBx)(i,j,k-1))/kappaEdz
      - ((*pBz)(i,j,k) - (*pBz)(i-1,j,k))/kappaEdx
      - eps_0_inv*Jy
    );

  (*pEz)(i,j,k) +=
    dt*(
        ((*pBy)(i,j,k) - (*pBy)(i-1,j,k))/kappaEdx
      - ((*pBx)(i,j,k) - (*pBx)(i,j-1,k))/kappaEdy
      - eps_0_inv*Jz
    );
}

inline void EMFields::fdtdStepB(double dt,
                              int i, int j, int k,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaHdx = dx[0];
  double kappaHdy = dx[1];
  double kappaHdz = dx[2];

//  double kappaHdx = (*pKappaHdx)(i)*dx;
//  double kappaHdy = (*pKappaHdy)(j)*dy;
//  double kappaHdz = (*pKappaHdz)(k)*dz;

  (*pBx)(i,j,k) +=
    + dt*(
        ((*pEy)(i,j,k+1) - (*pEy)(i,j,k))/kappaHdz
      - ((*pEz)(i,j+1,k) - (*pEz)(i,j,k))/kappaHdy
     - Jx
    );

  (*pBy)(i,j,k) +=
    + dt*(
        ((*pEz)(i+1,j,k) - (*pEz)(i,j,k))/kappaHdx
      - ((*pEx)(i,j,k+1) - (*pEx)(i,j,k))/kappaHdz
     - Jy
    );

  (*pBz)(i,j,k) +=
    + dt*(
        ((*pEx)(i,j+1,k) - (*pEx)(i,j,k))/kappaHdy
      - ((*pEy)(i+1,j,k) - (*pEy)(i,j,k))/kappaHdx
     - Jz
    );
}


void EMFields::stepD(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);

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

void EMFields::stepB(double dt)
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

        this->fdtdStepB(dt, i, j, k, dx, jx, jy, jz);
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

inline void EMFields::fdtdStepD(double dt,
                              int i, int j,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaEdx = dx[0]/clight2;
  double kappaEdy = dx[1]/clight2;

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;

  (*pEx)(i,j) +=
    dt*(
        ((*pBz)(i,j) - (*pBz)(i,j-1))/kappaEdy
      - eps_0_inv*Jx
    );

  (*pEy)(i,j) +=
    dt*(
      - ((*pBz)(i,j) - (*pBz)(i-1,j))/kappaEdx
      - eps_0_inv*Jy
    );

  (*pEz)(i,j) +=
    dt*(
        ((*pBy)(i,j) - (*pBy)(i-1,j))/kappaEdx
      - ((*pBx)(i,j) - (*pBx)(i,j-1))/kappaEdy
      - eps_0_inv*Jz
    );
}

inline void EMFields::fdtdStepB(double dt,
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
     - Jx
    );

  (*pBy)(i,j) +=
    + dt*(
        ((*pEz)(i+1,j) - (*pEz)(i,j))/kappaHdx
     - Jy
    );
  SCHNEK_TRACE_LOG(6, i << " " << j << " " << (*pBz)(i,j))

  (*pBz)(i,j) +=
    + dt*(
        ((*pEx)(i,j+1) - (*pEx)(i,j))/kappaHdy
      - ((*pEy)(i+1,j) - (*pEy)(i,j))/kappaHdx
     - Jz
    );
  SCHNEK_TRACE_LOG(6,(*pBz)(i,j))

}


void EMFields::stepD(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);

//  double maxEx  = 0.0;
//  double maxEy  = 0.0;
//  double maxEz  = 0.0;
//  double maxExa = 0.0;
//  double maxEya = 0.0;
//  double maxEza = 0.0;
//  double maxBx  = 0.0;
//  double maxBy  = 0.0;
//  double maxBz  = 0.0;
//  double maxJx  = 0.0;
//  double maxJy  = 0.0;
//  double maxJz  = 0.0;
  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      {
        jx = (*this->pJx)(i,j);
        jy = (*this->pJy)(i,j);
        jz = (*this->pJz)(i,j);

//        if (fabs((*pEx)(i,j))>maxEx) maxEx=(*pEx)(i,j);
//        if (fabs((*pEy)(i,j))>maxEy) maxEy=(*pEy)(i,j);
//        if (fabs((*pEz)(i,j))>maxEz) maxEz=(*pEz)(i,j);
//        if (fabs((*pBx)(i,j))>maxBx) maxBx=(*pBx)(i,j);
//        if (fabs((*pBy)(i,j))>maxBy) maxBy=(*pBy)(i,j);
//        if (fabs((*pBz)(i,j))>maxBz) maxBz=(*pBz)(i,j);
//
//        if (fabs(jx)>maxJx) maxJx=jx;
//        if (fabs(jy)>maxJy) maxJy=jy;
//        if (fabs(jz)>maxJz) maxJz=jz;

        this->fdtdStepD(dt, i, j, dx, jx, jy, jz);

//        if (fabs((*pEx)(i,j))>maxExa) maxExa=(*pEx)(i,j);
//        if (fabs((*pEy)(i,j))>maxEya) maxEya=(*pEy)(i,j);
//        if (fabs((*pEz)(i,j))>maxEza) maxEza=(*pEz)(i,j);
      }

//  std::cerr << "====== stepD ======\n";
//  std::cerr << "Ex  = " << maxEx  << std::endl;
//  std::cerr << "Ey  = " << maxEy  << std::endl;
//  std::cerr << "Ez  = " << maxEz  << std::endl;
//  std::cerr << "Exa = " << maxExa << std::endl;
//  std::cerr << "Eya = " << maxEya << std::endl;
//  std::cerr << "Eza = " << maxEza << std::endl;
//  std::cerr << "Bx  = " << maxBx  << std::endl;
//  std::cerr << "By  = " << maxBy  << std::endl;
//  std::cerr << "Bz  = " << maxBz  << std::endl;
//  std::cerr << "Jx  = " << maxJx  << std::endl;
//  std::cerr << "Jy  = " << maxJy  << std::endl;
//  std::cerr << "Jz  = " << maxJz  << std::endl;

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pEx,i);
    sub->exchange(*pEy,i);
    sub->exchange(*pEz,i);

//    boundariesLo[i]->applyEx(*pEx,i,FieldBC::lo);
//    boundariesLo[i]->applyEy(*pEy,i,FieldBC::lo);
//    boundariesLo[i]->applyEz(*pEz,i,FieldBC::lo);
//
//    boundariesHi[i]->applyEx(*pEx,i,FieldBC::hi);
//    boundariesHi[i]->applyEy(*pEy,i,FieldBC::hi);
//    boundariesHi[i]->applyEz(*pEz,i,FieldBC::hi);
  }
}

void EMFields::stepB(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);
//  if (this->pMx != 0) sumMagCurrents();

//  double maxEx  = 0.0;
//  double maxEy  = 0.0;
//  double maxEz  = 0.0;
//  double maxBxa = 0.0;
//  double maxBya = 0.0;
//  double maxBza = 0.0;
//  double maxBx  = 0.0;
//  double maxBy  = 0.0;
//  double maxBz  = 0.0;
//  double maxJx  = 0.0;
//  double maxJy  = 0.0;
//  double maxJz  = 0.0;

  for (int i=low[0]; i<high[0]; ++i)
    for (int j=low[1]; j<high[1]; ++j)
      {
//        if (this->pMx != 0)
//        {
//          jx = (*this->pMx)(i,j,k);
//          jy = (*this->pMy)(i,j,k);
//          jz = (*this->pMz)(i,j,k);
//        }


//        if (fabs((*pEx)(i,j))>maxEx) maxEx=(*pEx)(i,j);
//        if (fabs((*pEy)(i,j))>maxEy) maxEy=(*pEy)(i,j);
//        if (fabs((*pEz)(i,j))>maxEz) maxEz=(*pEz)(i,j);
//        if (fabs((*pBx)(i,j))>maxBx) maxBx=(*pBx)(i,j);
//        if (fabs((*pBy)(i,j))>maxBy) maxBy=(*pBy)(i,j);
//        if (fabs((*pBz)(i,j))>maxBz) maxBz=(*pBz)(i,j);
//
//        if (fabs(jx)>maxJx) maxJx=jx;
//        if (fabs(jy)>maxJy) maxJy=jy;
//        if (fabs(jz)>maxJz) maxJz=jz;

        this->fdtdStepB(dt, i, j, dx, jx, jy, jz);

//        if (fabs((*pBx)(i,j))>maxBxa) maxBxa=(*pBx)(i,j);
//        if (fabs((*pBy)(i,j))>maxBya) maxBya=(*pBy)(i,j);
//        if (fabs((*pBz)(i,j))>maxBza) maxBza=(*pBz)(i,j);
      }

//  std::cerr << "====== stepB ======\n";
//  std::cerr << "Ex  = " << maxEx  << std::endl;
//  std::cerr << "Ey  = " << maxEy  << std::endl;
//  std::cerr << "Ez  = " << maxEz  << std::endl;
//  std::cerr << "Bxa = " << maxBxa << std::endl;
//  std::cerr << "Bya = " << maxBya << std::endl;
//  std::cerr << "Bza = " << maxBza << std::endl;
//  std::cerr << "Bx  = " << maxBx  << std::endl;
//  std::cerr << "By  = " << maxBy  << std::endl;
//  std::cerr << "Bz  = " << maxBz  << std::endl;
//  std::cerr << "Jx  = " << maxJx  << std::endl;
//  std::cerr << "Jy  = " << maxJy  << std::endl;
//  std::cerr << "Jz  = " << maxJz  << std::endl;

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pBx,i);
    sub->exchange(*pBy,i);
    sub->exchange(*pBz,i);

//    boundariesLo[i]->applyBx(*pBx,i,FieldBC::lo);
//    boundariesLo[i]->applyBy(*pBy,i,FieldBC::lo);
//    boundariesLo[i]->applyBy(*pBz,i,FieldBC::lo);
//
//    boundariesHi[i]->applyBx(*pBx,i,FieldBC::hi);
//    boundariesHi[i]->applyBy(*pBy,i,FieldBC::hi);
//    boundariesHi[i]->applyBy(*pBz,i,FieldBC::hi);
  }
}

#endif

#ifdef ONE_DIMENSIONAL

inline void EMFields::fdtdStepD(double dt,
                              int i,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaEdx = dx[0]/clight2;

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;

  (*pEx)(i) -= eps_0_inv*dt*Jx;

  (*pEy)(i) +=
    dt*(
      - ((*pBz)(i) - (*pBz)(i-1))/kappaEdx
      - eps_0_inv*Jy
    );

  (*pEz)(i) +=
    dt*(
        ((*pBy)(i) - (*pBy)(i-1))/kappaEdx
      - eps_0_inv*Jz
    );
}

inline void EMFields::fdtdStepB(double dt,
                              int i,
                              SVector dx,
                              double Jx, double Jy, double Jz)
{
  SCHNEK_TRACE_ENTER_FUNCTION(5)
  double kappaHdx = dx[0];

//  double kappaHdx = (*pKappaHdx)(i)*dx;
//  double kappaHdy = (*pKappaHdy)(j)*dy;
//  double kappaHdz = (*pKappaHdz)(k)*dz;

  (*pBx)(i) -= dt*Jx;

  (*pBy)(i) +=
    + dt*(
        ((*pEz)(i+1) - (*pEz)(i))/kappaHdx
     - Jy
    );
  SCHNEK_TRACE_LOG(6, i << " " << j << " " << (*pBz)(i))

  (*pBz)(i) +=
    + dt*(
      - ((*pEy)(i+1) - (*pEy)(i))/kappaHdx
     - Jz
    );
  SCHNEK_TRACE_LOG(6,(*pBz)(i))

}


void EMFields::stepD(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  SIntVector low  = Globals::instance().getLocalGridMin();
  SIntVector high = Globals::instance().getLocalGridMax();

  SVector dx = Globals::instance().getDx();

  double jx(0), jy(0), jz(0);

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

void EMFields::stepB(double dt)
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


void EMFields::stepSchemeBefore(double dt)
{
  debug_check_out_of_bounds("Fields B");
  stepD(0.5*dt);

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

void EMFields::stepSchemeAfter(double dt)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  debug_check_out_of_bounds("Fields A");
  Currents::instance().update();

  stepB(0.5*dt);
  debug_check_out_of_bounds("Fields C");

  debug_check_out_of_bounds("Fields B");
  stepD(0.5*dt);

//  BOOST_FOREACH(Current *cur, this->magCurrents)
//  {
//    cur->stepScheme(dt);
//  }

}



void EMFields::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)

  std::cout << "Fields::initParameters()" << std::endl;
  Fields::initParameters(blockPars);

  blockPars.addArrayParameter("boundary_min", bcNamesLo, std::string("periodic"));
  blockPars.addArrayParameter("boundary_max", bcNamesHi, std::string("periodic"));

  fieldBCFactories["periodic"] = boost::factory<FieldPeriodicBC*>();
  fieldBCFactories["conducting"] = boost::factory<FieldConductingBC*>();
  fieldBCFactories["symmetric"] = boost::factory<FieldSymmetryBC*>();
}

void EMFields::registerData()
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
  Fields::registerData();
  addData("Jx", pJx);
  addData("Jy", pJy);
  addData("Jz", pJz);
}

void EMFields::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  static schnek::LiteratureArticle Yee1966("Yee:1966", "Yee, K",
      "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media.",
      "IEEE Transactions on Antennas and Propagation", "1966", "AP-14", "302--307");

  schnek::LiteratureManager::instance().addReference(
      "Integration of electrodynamic fields by Finite Difference Time Domain method.",
      Yee1966);

  dynamic_cast<OPar&>(*this->getParent()).addField(this);

  SIntVector low  = Globals::instance().getLocalInnerGridMin();
  SIntVector high = Globals::instance().getLocalInnerGridMax();
  SRange grange = Globals::instance().getLocalDomainRange();

  pEx = pDataField(new DataField(low, high, grange, exStaggerYee, ghostCells));
  pEy = pDataField(new DataField(low, high, grange, eyStaggerYee, ghostCells));
  pEz = pDataField(new DataField(low, high, grange, ezStaggerYee, ghostCells));

  pBx = pDataField(new DataField(low, high, grange, bxStaggerYee, ghostCells));
  pBy = pDataField(new DataField(low, high, grange, byStaggerYee, ghostCells));
  pBz = pDataField(new DataField(low, high, grange, bzStaggerYee, ghostCells));

  pJx = pDataField(new DataField(low, high, grange, exStaggerYee, ghostCells));
  pJy = pDataField(new DataField(low, high, grange, eyStaggerYee, ghostCells));
  pJz = pDataField(new DataField(low, high, grange, ezStaggerYee, ghostCells));

  Currents::instance().setGlobalCurrent(pJx, pJy, pJz);

  Fields::init();
//  for (int i=0; i<dimension; ++i)
//  {
//    std::cout << "Field: "<< Globals::instance().getSubdivision()->getUniqueId() << " " << i <<
//        " (" << low[i] << " " << high[i] <<
//        ") (" << grange.getLo()[i] << " " << grange.getHi()[i] << ")" << std::endl;
//  }


  for (int i=0; i<dimension; ++i)
  {
    if (fieldBCFactories.count(bcNamesLo[i]) == 0)
      terminateSim("Unknown boundary condition: "+bcNamesLo[i]);

    if (fieldBCFactories.count(bcNamesHi[i]) == 0)
      terminateSim("Unknown boundary condition: "+bcNamesHi[i]);

    boundariesLo[i] = pFieldBC(fieldBCFactories[bcNamesLo[i]]());
    boundariesHi[i] = pFieldBC(fieldBCFactories[bcNamesHi[i]]());
  }

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pEx,i);
    sub->exchange(*pEy,i);
    sub->exchange(*pEz,i);

    sub->exchange(*pBx,i);
    sub->exchange(*pBy,i);
    sub->exchange(*pBz,i);
  }

}

void ConstantFields::init()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  dynamic_cast<OPar&>(*this->getParent()).addField(this);

  SIntVector low  = Globals::instance().getLocalInnerGridMin();
  SIntVector high = Globals::instance().getLocalInnerGridMax();
  SRange grange = Globals::instance().getLocalDomainRange();

  pEx = pDataField(new DataField(low, high, grange, exStaggerYee, ghostCells));
  pEy = pDataField(new DataField(low, high, grange, eyStaggerYee, ghostCells));
  pEz = pDataField(new DataField(low, high, grange, ezStaggerYee, ghostCells));

  pBx = pDataField(new DataField(low, high, grange, bxStaggerYee, ghostCells));
  pBy = pDataField(new DataField(low, high, grange, byStaggerYee, ghostCells));
  pBz = pDataField(new DataField(low, high, grange, bzStaggerYee, ghostCells));

  Fields::init();

  Globals::pSubdivision sub = Globals::instance().getSubdivision();

  for (int i=0; i<dimension; ++i)
  {
    sub->exchange(*pEx,i);
    sub->exchange(*pEy,i);
    sub->exchange(*pEz,i);

    sub->exchange(*pBx,i);
    sub->exchange(*pBy,i);
    sub->exchange(*pBz,i);
  }

}

