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

#include <schnek/grid.hpp>
#include <schnek/variables.hpp>

#include <boost/shared_ptr.hpp>

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

#ifdef NDEBUG
#define ArrayArgCheck schnek::ArrayNoArgCheck
#define GridArgCheck schnek::GridNoArgCheck
#else
#define ArrayArgCheck schnek::ArrayNoArgCheck
#define GridArgCheck schnek::GridAssertCheck
#endif


typedef schnek::Field<double, 1, GridArgCheck> DataField1d;
typedef schnek::Field<double, 2, GridArgCheck> DataField2d;
typedef schnek::Field<double, 3, GridArgCheck> DataField3d;
typedef DataField1d DataLine;

typedef schnek::Field<double, dimension, GridArgCheck> DataField;
typedef DataField::IndexType FieldIndex;

typedef boost::shared_ptr<DataField1d> pDataField1d;
typedef boost::shared_ptr<DataField2d> pDataField2d;
typedef boost::shared_ptr<DataField3d> pDataField3d;
typedef boost::shared_ptr<DataField> pDataField;

typedef schnek::Array<real, dimension, ArrayArgCheck> SVector;
typedef schnek::Range<real, dimension, ArrayArgCheck> SRange;
typedef schnek::Array<bool, dimension> SStagger;

typedef schnek::Array<real, 3, ArrayArgCheck> PVector;

typedef schnek::Array<int, dimension, ArrayArgCheck> SIntVector;
typedef schnek::Array<int, 3, ArrayArgCheck> PIntVector;

typedef schnek::Array<schnek::pParameter, dimension, ArrayArgCheck> SParameterVector;
typedef schnek::Array<schnek::pParameter, 3, ArrayArgCheck> PParameterVector;


enum Direction {north, south, west, east, up, down};

#endif /* CONFIG_HPP_ */
