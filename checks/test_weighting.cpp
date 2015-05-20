/*
 * test_weighting.cpp
 *
 * Created on: 25 Mar 2015
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

#include "../src/weighting.hpp"


bool test_weighting()
{
#ifdef ONE_DIMENSIONAL
  DataField::RangeLimit rlo(-2.0), rhi(12.0);
  DataField::RangeType range(rlo, rhi);
  DataField::IndexType lo(-2), hi(12);
  DataField::Stagger noStagger(false);
  DataField::Stagger withStagger(true);
  DataField Ex(lo, hi, range, withStagger, 2);
  DataField Ey(lo, hi, range, noStagger, 2);
  DataField Ez(lo, hi, range, noStagger, 2);
  DataField Bx(lo, hi, range, noStagger, 2);
  DataField By(lo, hi, range, withStagger, 2);
  DataField Bz(lo, hi, range, noStagger, 2);

  Weighting::WeightingCoefficients gx, hx;

  ScalarDomain ds = Weighting::getScalarDomain();
  gx.resize(ds.getLo(), ds.getHi());
  hx.resize(ds.getLo(), ds.getHi());

  Ex = 0.0;
  Ey = 0.0;
  Ez = 0.0;
  Bx = 0.0;
  By = 0.0;
  Bz = 0.0;

  Ex(5) = 1.0;
  Ey(5) = 1.0;
  Ez(5) = 1.0;
  Bx(5) = 1.0;
  By(5) = 1.0;
  Bz(5) = 1.0;

  SVector x;
  SIntVector cell1, cell2;
  SVector cell_frac1, cell_frac2;
  for (x[0]=0; x[0]<=10.0; x[0]+=0.01)
  {
    Weighting::toCellIndex(x, cell1, cell_frac1);
    Weighting::getShape(cell1, cell_frac1, gx);
    Weighting::toCellIndexStagger(x, cell2, cell_frac2);
    Weighting::getShape(cell2, cell_frac2, hx);

    PVector E = Weighting::interpolateE(gx, hx, cell1, cell2, Ex, Ey, Ez);
    PVector B = Weighting::interpolateB(gx, hx, cell1, cell2, Bx, By, Bz);
    std::cerr << x[0] << " " << cell1[0] << " " << cell_frac1[0] <<
        " " << (cell1[0]+cell_frac1[0]) <<
        " " << cell2[0] << " " << cell_frac2[0] <<
        " " << (cell2[0]+cell_frac2[0]) <<
        " -- " << gx(-1) << " " << gx(0) << " " << gx(1) << " " << (gx(-1)[0] + gx(0)[0] + gx(1)[0]) <<
        " -- " << hx(-1) << " " << hx(0) << " " << hx(1) << " " << (hx(-1)[0] + hx(0)[0] + hx(1)[0]) <<
        " -- " << E[0] << " " << E[1] << " " << E[2] <<
        " -- " << B[0] << " " << B[1] << " " << B[2] << std::endl;
  }
  return true;
#else
  return true;
#endif
}

int main()
{
  return test_weighting()?0:-1;
}
