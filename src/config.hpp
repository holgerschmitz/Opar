/*
 * config.hpp
 *
 *  Created on: 1 Aug 2012
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

#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include <schnek/array.hpp>
#include <schnek/grid.hpp>

#ifdef ONE_DIMENSIONAL
static const int dimension = 1;
#endif

#ifdef TWO_DIMENSIONAL
static const int dimension = 2;
#endif

#ifdef THREE_DIMENSIONAL
static const int dimension = 3;
#endif

#ifdef SINGLE_PRECISION
typedef float real;
#else
typedef double real;
#endif

#ifdef DEBUG
#define ArrayArgCheck schnek::ArrayNoArgCheck
#else
#define ArrayArgCheck schnek::ArrayNoArgCheck
#endif

typedef schnek::Array<real, dimension, ArrayArgCheck> SVector;
typedef schnek::Array<real, 3, ArrayArgCheck> PVector;

typedef schnek::Array<int, dimension, ArrayArgCheck> SIntVector;
typedef schnek::Array<int, 3, ArrayArgCheck> PIntVector;


#endif /* CONFIG_HPP_ */
