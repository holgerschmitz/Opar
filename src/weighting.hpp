/*
 * weighting.hpp
 *
 * Created on: 30 Nov 2012
 * Author: hschmitz
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

#ifndef WEIGHTING_HPP_
#define WEIGHTING_HPP_

#include "types.hpp"

class TriangularWeighting
{
  public:
    typedef schnek::Grid<SVector, 1, GridArgCheck> WeightingCoefficients;

    static ScalarDomain getScalarDomain()
    {
      ScalarLimit lo(-1);
      ScalarLimit hi(1);
      return ScalarDomain(lo, hi);
    }

    static SDomain getSDomain()
    {
#ifdef ONE_DIMENSIONAL
      SLimit lo(-1);
      SLimit hi(1);
#endif
#ifdef TWO_DIMENSIONAL
      SLimit lo(-1,-1);
      SLimit hi(1,1);
#endif
#ifdef THREE_DIMENSIONAL
      SLimit lo(-1,-1,-1);
      SLimit hi(1,1,1);
#endif
      return SDomain(lo, hi);
    }

    static PDomain getPDomain()
    {
      PLimit lo(-1,-1,-1);
      PLimit hi(1,1,1);
      return PDomain(lo, hi);
    }

    static double particleShapeFactor()
    {
      return 1.0;
    }

    static void toCellIndex(const SVector pos, SIntVector &cell, SVector &cell_frac)
    {
      for (int d = 0; d < dimension; ++d)
      {
        cell[d] = int(pos[d] + 0.5);
        cell_frac[d] = pos[d] - cell[d];
      }
    }

    static void toCellIndexStagger(const SVector pos, SIntVector &cell,
        SVector &cell_frac)
    {
      for (int d = 0; d < dimension; ++d)
      {
        cell[d] = int(pos[d]);
        cell_frac[d] = pos[d] - cell[d] - 0.5;
      }
    }

    static void getShape(const SIntVector cell, const SVector cell_frac,
        WeightingCoefficients &gx)
    {
      for (int d = 0; d < dimension; ++d)
      {
        double cf = -cell_frac[d];
        double cf2 = cf * cf;
        gx(-1)[d] = 0.5*(0.25 + cf2 + cf);
        gx(+0)[d] = 0.5*(1.5 - 2.0 * cf2);
        gx(+1)[d] = 0.5*(0.25 + cf2 - cf);
      }
    }

//    static void getShape(const SIntVector cell, const SVector cell_frac,
//        const SIntVector shift, WeightingCoefficients &gx)
//    {
//      for (int d = 0; d < dimension; ++d)
//      {
//        double cf = cell_frac[d];
//        double cf2 = cf * cf;
//        gx(shift[d] - 1)[d] = 0.25 + cf2 + cf;
//        gx(shift[d]    )[d] = 1.5 - 2.0 * cf2;
//        gx(shift[d] + 1)[d] = 0.25 + cf2 - cf;
//      }
//    }

    // TODO check speed of loops versus unrolled
    static PVector interpolateE(const WeightingCoefficients &wc,
        const WeightingCoefficients &ws, const SIntVector &ic,
        const SIntVector &is, const DataField &Ex, const DataField &Ey,
        const DataField &Ez)
    {
#ifdef ONE_DIMENSIONAL
      PVector result(0.0,0.0,0.0);
      SDomain d = getSDomain();
      for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
      {
        result[0] += ws(i)[0]*Ex(is[0]+i);
        result[1] += wc(i)[0]*Ey(ic[0]+i);
        result[2] += wc(i)[0]*Ez(ic[0]+i);
      }
      return result;

      //      return PVector(
      //          ws(-1)[0] * Ex(is[0]-1)
      //          + ws( 0)[0] * Ex(is[0] )
      //          + ws( 1)[0] * Ex(is[0]+1),
      //          wc(-1)[0] * Ey(ic[0]-1)
      //          + wc( 0)[0] * Ey(ic[0] )
      //          + wc( 1)[0] * Ey(ic[0]+1),
      //          wc(-1)[0] * Ez(ic[0]-1)
      //          + wc( 0)[0] * Ez(ic[0] )
      //          + wc( 1)[0] * Ez(ic[0]+1)
      //      );
#endif
#ifdef TWO_DIMENSIONAL
      PVector result(0.0,0.0,0.0);
      SDomain d = getSDomain();

      for (int j=d.getLo()[1]; j<=d.getHi()[1]; ++j)
      {
        const double Wc = wc(j)[1];
        const double Ws = ws(j)[1];
        for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
        {
          result[0] += Wc*ws(i)[0]*Ex(is[0]+i, ic[1]+j);
          result[1] += Ws*wc(i)[0]*Ey(ic[0]+i, is[1]+j);
          result[2] += Wc*wc(i)[0]*Ez(ic[0]+i, ic[1]+j);
        }
      }
      return result;
      //
      //      return PVector(
      //          wc(-1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]-1)
      //              + ws( 0)[0] * ex(is[0] ,ic[1]-1)
      //              + ws( 1)[0] * ex(is[0]+1,ic[1]-1))
      //          + wc( 0)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1] )
      //              + ws( 0)[0] * ex(is[0] ,ic[1] )
      //              + ws( 1)[0] * ex(is[0]+1,ic[1] ))
      //          + wc( 1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]+1)
      //              + ws( 0)[0] * ex(is[0] ,ic[1]+1)
      //              + ws( 1)[0] * ex(is[0]+1,ic[1]+1)),
      //          ws(-1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]-1)
      //              + wc( 0)[0] * ey(ic[0] ,is[1]-1)
      //              + wc( 1)[0] * ey(ic[0]+1,is[1]-1))
      //          + ws( 0)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1] )
      //              + wc( 0)[0] * ey(ic[0] ,is[1] )
      //              + wc( 1)[0] * ey(ic[0]+1,is[1] ))
      //          + ws( 1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]+1)
      //              + wc( 0)[0] * ey(ic[0] ,is[1]+1)
      //              + wc( 1)[0] * ey(ic[0]+1,is[1]+1)),
      //          wc(-1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]-1)
      //              + wc( 0)[0] * ez(ic[0] ,ic[1]-1)
      //              + wc( 1)[0] * ez(ic[0]+1,ic[1]-1))
      //          + wc( 0)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1] )
      //              + wc( 0)[0] * ez(ic[0] ,ic[1] )
      //              + wc( 1)[0] * ez(ic[0]+1,ic[1] ))
      //          + wc( 1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]+1)
      //              + wc( 0)[0] * ez(ic[0] ,ic[1]+1)
      //              + wc( 1)[0] * ez(ic[0]+1,ic[1]+1)) );
#endif
#ifdef THREE_DIMENSIONAL
      PVector result(0.0,0.0,0.0);
      SDomain d = getSDomain();

      for (int k=d.getLo()[2]; k<=d.getHi()[2]; ++k)
      {
        const double Wzc = wc(k)[2];
        const double Wzs = ws(k)[2];
        for (int j=d.getLo()[1]; j<=d.getHi()[1]; ++j)
        {
          const double Wyc = wc(j)[1];
          const double Wys = ws(j)[1];
          for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
          {
            result[0] += Wzc*Wyc*ws(i)[0]*Ex(is[0]+i, ic[1]+j, ic[2]+k);
            result[1] += Wzc*Wys*wc(i)[0]*Ey(ic[0]+i, is[1]+j, ic[2]+k);
            result[2] += Wzs*Wyc*wc(i)[0]*Ez(ic[0]+i, ic[1]+j, is[2]+k);
          }
        }
      }
      return result;

        //        return PVector(
        //            wc(-1)[2] * (wc(-1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]-1,ic[2]-1)
        //                    + ws( 0)[0] * ex(is[0] ,ic[1]-1,ic[2]-1)
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1]-1,ic[2]-1))
        //                + wc( 0)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1] ,ic[2]-1)
        //                    + ws( 0)[0] * ex(is[0] ,ic[1] ,ic[2]-1)
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1] ,ic[2]-1))
        //                + wc( 1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]+1,ic[2]-1)
        //                    + ws( 0)[0] * ex(is[0] ,ic[1]+1,ic[2]-1)
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1]+1,ic[2]-1)))
        //            + wc( 0)[2] * (wc(-1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]-1,ic[2] )
        //                    + ws( 0)[0] * ex(is[0] ,ic[1]-1,ic[2] )
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1]-1,ic[2] ))
        //                + wc( 0)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1] ,ic[2] )
        //                    + ws( 0)[0] * ex(is[0] ,ic[1] ,ic[2] )
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1] ,ic[2] ))
        //                + wc( 1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]+1,ic[2] )
        //                    + ws( 0)[0] * ex(is[0] ,ic[1]+1,ic[2] )
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1]+1,ic[2] )))
        //            + wc( 1)[2] * (wc(-1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]-1,ic[2]+1)
        //                    + ws( 0)[0] * ex(is[0] ,ic[1]-1,ic[2]+1)
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1]-1,ic[2]+1))
        //                + wc( 0)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1] ,ic[2]+1)
        //                    + ws( 0)[0] * ex(is[0] ,ic[1] ,ic[2]+1)
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1] ,ic[2]+1))
        //                + wc( 1)[1] * (ws(-1)[0] * ex(is[0]-1,ic[1]+1,ic[2]+1)
        //                    + ws( 0)[0] * ex(is[0] ,ic[1]+1,ic[2]+1)
        //                    + ws( 1)[0] * ex(is[0]+1,ic[1]+1,ic[2]+1))),
        //            wc(-1)[2] * (ws(-1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]-1,ic[2]-1)
        //                    + wc( 0)[0] * ey(ic[0] ,is[1]-1,ic[2]-1)
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1]-1,ic[2]-1))
        //                + ws( 0)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1] ,ic[2]-1)
        //                    + wc( 0)[0] * ey(ic[0] ,is[1] ,ic[2]-1)
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1] ,ic[2]-1))
        //                + ws( 1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]+1,ic[2]-1)
        //                    + wc( 0)[0] * ey(ic[0] ,is[1]+1,ic[2]-1)
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1]+1,ic[2]-1)))
        //            + wc( 0)[2] * (ws(-1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]-1,ic[2] )
        //                    + wc( 0)[0] * ey(ic[0] ,is[1]-1,ic[2] )
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1]-1,ic[2] ))
        //                + ws( 0)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1] ,ic[2] )
        //                    + wc( 0)[0] * ey(ic[0] ,is[1] ,ic[2] )
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1] ,ic[2] ))
        //                + ws( 1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]+1,ic[2] )
        //                    + wc( 0)[0] * ey(ic[0] ,is[1]+1,ic[2] )
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1]+1,ic[2] )))
        //            + wc( 1)[2] * (ws(-1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]-1,ic[2]+1)
        //                    + wc( 0)[0] * ey(ic[0] ,is[1]-1,ic[2]+1)
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1]-1,ic[2]+1))
        //                + ws( 0)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1] ,ic[2]+1)
        //                    + wc( 0)[0] * ey(ic[0] ,is[1] ,ic[2]+1)
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1] ,ic[2]+1))
        //                + ws( 1)[1] * (wc(-1)[0] * ey(ic[0]-1,is[1]+1,ic[2]+1)
        //                    + wc( 0)[0] * ey(ic[0] ,is[1]+1,ic[2]+1)
        //                    + wc( 1)[0] * ey(ic[0]+1,is[1]+1,ic[2]+1))),
        //            ws(-1)[2] * (wc(-1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]-1,is[2]-1)
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1]-1,is[2]-1)
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1]-1,is[2]-1))
        //                + wc( 0)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1] ,is[2]-1)
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1] ,is[2]-1)
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1] ,is[2]-1))
        //                + wc( 1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]+1,is[2]-1)
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1]+1,is[2]-1)
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1]+1,is[2]-1)))
        //            + ws( 0)[2] * (wc(-1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]-1,is[2] )
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1]-1,is[2] )
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1]-1,is[2] ))
        //                + wc( 0)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1] ,is[2] )
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1] ,is[2] )
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1] ,is[2] ))
        //                + wc( 1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]+1,is[2] )
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1]+1,is[2] )
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1]+1,is[2] )))
        //            + ws( 1)[2] * (wc(-1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]-1,is[2]+1)
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1]-1,is[2]+1)
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1]-1,is[2]+1))
        //                + wc( 0)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1] ,is[2]+1)
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1] ,is[2]+1)
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1] ,is[2]+1))
        //                + wc( 1)[1] * (wc(-1)[0] * ez(ic[0]-1,ic[1]+1,is[2]+1)
        //                    + wc( 0)[0] * ez(ic[0] ,ic[1]+1,is[2]+1)
        //                    + wc( 1)[0] * ez(ic[0]+1,ic[1]+1,is[2]+1)))
        //
        //        );
#endif

    }

    static PVector interpolateB(const WeightingCoefficients &wc,
        const WeightingCoefficients &ws, const SIntVector &ic,
        const SIntVector &is, const DataField &Ex, const DataField &Ey,
        const DataField &Ez)
    {
#ifdef ONE_DIMENSIONAL
      PVector result(0.0,0.0,0.0);
      SDomain d = getSDomain();
      for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
      {
        result[0] += wc(i)[0]*Ex(ic[0]+i);
        result[1] += ws(i)[0]*Ey(is[0]+i);
        result[2] += ws(i)[0]*Ez(is[0]+i);
      }
      return result;
#endif

#ifdef TWO_DIMENSIONAL
      PVector result(0.0,0.0,0.0);
      SDomain d = getSDomain();

      for (int j=d.getLo()[1]; j<=d.getHi()[1]; ++j)
      {
        const double Wc = wc(j)[1];
        const double Ws = ws(j)[1];
        for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
        {
          result[0] += Ws*wc(i)[0]*Ex(ic[0]+i, is[1]+j);
          result[1] += Wc*ws(i)[0]*Ey(is[0]+i, ic[1]+j);
          result[2] += Ws*ws(i)[0]*Ez(is[0]+i, is[1]+j);
        }
      }
      return result;
#endif

#ifdef THREE_DIMENSIONAL
      PVector result(0.0,0.0,0.0);
      SDomain d = getSDomain();

      for (int k=d.getLo()[2]; k<=d.getHi()[2]; ++k)
      {
        const double Wzc = wc(k)[2];
        const double Wzs = ws(k)[2];
        for (int j=d.getLo()[1]; j<=d.getHi()[1]; ++j)
        {
          const double Wyc = wc(j)[1];
          const double Wys = ws(j)[1];
          for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
          {
            result[0] += Wzs*Wys*wc(i)[0]*Ex(ic[0]+i, is[1]+j, is[2]+k);
            result[1] += Wzs*Wyc*ws(i)[0]*Ey(is[0]+i, ic[1]+j, is[2]+k);
            result[2] += Wzc*Wys*wc(i)[0]*Ez(is[0]+i, is[1]+j, ic[2]+k);
          }
        }
      }
      return result;
#endif

    }

    /**
     * Weighting a cell centered scalar field
     * @param wc the cell centered weighting coefficients
     * @param ic the cell centered grid index
     * @param F the field to interpolate
     * @return the interpolated result
     */
    static double interpolateScalar(const WeightingCoefficients &wc,
        const SIntVector &ic, const DataField &F)
    {
#ifdef ONE_DIMENSIONAL
      double result = 0.0;
      SDomain d = getSDomain();
      for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
      {
        result += wc(i)[0]*F(ic[0]+i);
      }
      return result;
#endif

#ifdef TWO_DIMENSIONAL
      double result = 0.0;
      SDomain d = getSDomain();

      for (int j=d.getLo()[1]; j<=d.getHi()[1]; ++j)
      {
        const double Wc = wc(j)[1];
        for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
        {
          result += Wc*wc(i)[0]*F(ic[0]+i, ic[1]+j);
        }
      }
      return result;
#endif

#ifdef THREE_DIMENSIONAL
      double result = 0.0;
      SDomain d = getSDomain();

      for (int k=d.getLo()[2]; k<=d.getHi()[2]; ++k)
      {
        const double Wzc = wc(k)[2];
        for (int j=d.getLo()[1]; j<=d.getHi()[1]; ++j)
        {
          const double Wyc = wc(j)[1];
          for (int i=d.getLo()[0]; i<=d.getHi()[0]; ++i)
          {
            result += Wzc*Wyc*wc(i)[0]*F(ic[0]+i, ic[1]+j, ic[2]+k);
          }
        }
      }
      return result;
#endif

    }
};



#endif // WEIGHTING_HPP_ 
