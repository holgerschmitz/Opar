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

#include "fields.hpp"


inline void FieldImplementation::fdtdStepD(double dt,
                                           int i, int j, int k,
                                           double dx, double dy, double dz,
                                           double Jx, double Jy, double Jz)
{
  double &ex = (*pEx)(i,j,k);
  double &ey = (*pEy)(i,j,k);
  double &ez = (*pEz)(i,j,k);

#ifndef NDEBUG
  {
      double test = ex*ey*ez;
      if ( !((test>0) || (test<1)) )
      {
        std::cerr << "NaN in Field (1) " << ex << " " << ey << " " << ez << "\n";
        exit(-1);
      }
  }
#endif

  double kappaEdx = 1.0;
  double kappaEdy = 1.0;
  double kappaEdz = 1.0;

//  double kappaEdx = (*pKappaEdx)(i)*dx;
//  double kappaEdy = (*pKappaEdy)(j)*dy;
//  double kappaEdz = (*pKappaEdz)(k)*dz;
//
////  std::cerr << i << " " << j << " " << k << " " << kappaEdx << " " << kappaEdy << " " << kappaEdz << " " << std::endl;
//
//#ifndef NDEBUG
//  {
//      double test = kappaEdx*kappaEdy*kappaEdz;
//      if ( !((test>0) || (test<1)) || (test==0) )
//      {
//        std::cerr << "Error in Kappa " << kappaEdx << " " << kappaEdy << " " << kappaEdz << "\n";
//        exit(-1);
//      }
//  }
//#endif

  double exn =
    dt*(
        ((*pBz)(i,j,k) - (*pBz)(i,j-1,k))/kappaEdy
      - ((*pBy)(i,j,k) - (*pBy)(i,j,k-1))/kappaEdz
      + Jx
    );

  double eyn =
    dt*(
        ((*pBx)(i,j,k) - (*pBx)(i,j,k-1))/kappaEdz
      - ((*pBz)(i,j,k) - (*pBz)(i-1,j,k))/kappaEdx
      + Jy
    );

  double ezn =
    dt*(
        ((*pBy)(i,j,k) - (*pBy)(i-1,j,k))/kappaEdx
      - ((*pBx)(i,j,k) - (*pBx)(i,j-1,k))/kappaEdy
      + Jz
    );

  ex = exn;
  ey = eyn;
  ez = ezn;

//#ifndef NDEBUG
// {
//     double test = ex*ey*ez;
//     if ( !((test>0) || (test<1)) )
//     {
//         std::cerr << "NaN in Field (2) " << ex << " " << ey << " " << ez << "\n";
//         exit(-1);
//     }
// }
//#endif

}

inline void FieldImplementation::fdtdStepB(double dt,
                                           int i, int j, int k,
                                           double dx, double dy, double dz,
                                           double Jx, double Jy, double Jz)
{

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


void FieldsImplementation::stepSchemeInit(double dt)
{
  stepB(0.5*dt);

  BOOST_FOREACH(Current *cur, this->currents)
  {
    cur->stepSchemeInit(dt);
  }

  BOOST_FOREACH(Current *cur, this->magCurrents)
  {
    cur->stepSchemeInit(dt);
  }
}
void FieldsImplementation::stepScheme(double dt)
{


  BOOST_FOREACH(Current *cur, this->currents)
  {
    cur->stepScheme(dt);
  }

  stepD(dt);

  BOOST_FOREACH(Current *cur, this->magCurrents)
  {
    cur->stepScheme(dt);
  }

  stepB(dt);
}

void FieldsImplementation::stepD(double dt)
{
  GridIndex low = this->storage->getLow();
  GridIndex high = this->storage->getHigh();

  double dx = this->storage->getDx();
  double dy = this->storage->getDy();
  double dz = this->storage->getDz();


//  double jx(0), jy(0), jz(0);
//  if (this->pJx != 0) sumCurrents();

  for (int i=low[0]+1; i<high[0]; ++i)
    for (int j=low[1]+1; j<high[1]; ++j)
      for (int k=low[2]+1; k<high[2]; ++k)
      {
//        if (this->pJx != 0)
//        {
//          jx = (*this->pJx)(i,j,k);
//          jy = (*this->pJy)(i,j,k);
//          jz = (*this->pJz)(i,j,k);
//        }

        this->fdtdStepD(dt, i, j, k, dx, dy, dz, jx, jy, jz);
      }


  this->storage->applyBoundary("E");
}

void FieldsImplementation::stepB(double dt)
{
  GridIndex low = this->storage->getLow();
  GridIndex high = this->storage->getHigh();

  double dx = this->storage->getDx();
  double dy = this->storage->getDy();
  double dz = this->storage->getDz();

//  double jx(0), jy(0), jz(0);
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

  this->storage->applyBoundary("B");
}



void Fields::initParameters(BlockParameters &blockPars)
{
  blockPars.addParameter(&ExInit, "Ex");
  blockPars.addParameter(&EyInit, "Ey");
  blockPars.addParameter(&EzInit, "Ez");
  blockPars.addParameter(&BxInit, "Bx");
  blockPars.addParameter(&ByInit, "By");
  blockPars.addParameter(&BzInit, "Bz");
}

void Fields::registerData()
{
  addData("Ex", impl.pEx);
  addData("Ey", impl.pEy);
  addData("Ez", impl.pEz);
  addData("Bx", impl.pBx);
  addData("By", impl.pBy);
  addData("Bz", impl.pBz);
}

void Fields::init()
{}




