/*
 * functions.cpp
 *
 * Created on: 11 Feb 2015
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

#include "functions.hpp"

#include <cmath>
#include <limits>

double step(double x)
{
  return (x<0.0)?0.0:1.0;
}

double logistic(double x, double w, double x0)
{
  static const double lim = M_LN2*std::numeric_limits<double>::digits + 1;

  double z = (x-x0)/w;
  if (z>lim)
    return 1.0;
  else if (z<-lim)
    return exp(-z);
  else
    return 1/(1+exp(-z));
}
