/*
 * types.hpp
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

#include <iostream>
#include <schnek/config.hpp>

#undef LOGLEVEL
#define LOGLEVEL 0
#include <schnek/grid.hpp>
#include <schnek/variables.hpp>

#include <boost/shared_ptr.hpp>

#undef LOGLEVEL
#define LOGLEVEL 0

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

//#define ArrayArgCheck schnek::ArrayNoArgCheck
//#define GridArgCheck schnek::GridNoArgCheck
#define ArrayArgCheck schnek::ArrayAssertArgCheck
#define GridArgCheck schnek::GridAssertCheck

class TriangularWeighting;
typedef TriangularWeighting Weighting;

typedef schnek::Grid<real, 1, GridArgCheck> DataGrid1d;
typedef boost::shared_ptr<DataGrid1d> pDataGrid1d;

typedef schnek::Grid<real, dimension, GridArgCheck> DataGrid;
typedef boost::shared_ptr<DataGrid> pDataGrid;

typedef schnek::Field<real, 1, GridArgCheck> DataField1d;
typedef schnek::Field<real, 2, GridArgCheck> DataField2d;
typedef schnek::Field<real, 3, GridArgCheck> DataField3d;
typedef DataField1d DataLine;

typedef schnek::Field<real, dimension, GridArgCheck> DataField;

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

typedef schnek::Range<int, 1> ScalarDomain;
typedef schnek::Range<int, dimension> SDomain;
typedef schnek::Range<int, 3> PDomain;

typedef ScalarDomain::LimitType ScalarLimit;
typedef SDomain::LimitType SLimit;
typedef PDomain::LimitType PLimit;

enum Direction {north, south, west, east, up, down};

#ifdef THREE_DIMENSIONAL
static const SStagger exStaggerYee(true,  false, false);
static const SStagger eyStaggerYee(false, true,  false);
static const SStagger ezStaggerYee(false, false, true );

static const SStagger bxStaggerYee(false, true,  true );
static const SStagger byStaggerYee(true,  false, true );
static const SStagger bzStaggerYee(true,  true,  false);
#endif


#ifdef TWO_DIMENSIONAL
static const SStagger exStaggerYee(true,  false);
static const SStagger eyStaggerYee(false, true );
static const SStagger ezStaggerYee(false, false);

static const SStagger bxStaggerYee(false, true );
static const SStagger byStaggerYee(true,  false);
static const SStagger bzStaggerYee(true,  true );
#endif


#ifdef ONE_DIMENSIONAL
static const SStagger exStaggerYee(true );
static const SStagger eyStaggerYee(false);
static const SStagger ezStaggerYee(false);

static const SStagger bxStaggerYee(false);
static const SStagger byStaggerYee(true );
static const SStagger bzStaggerYee(true );
#endif

void debug_check_out_of_bounds(std::string checkpoint);

#endif /* CONFIG_HPP_ */
