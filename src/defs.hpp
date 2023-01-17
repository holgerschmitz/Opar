/*
 * defs.hpp
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

#ifndef SCHNEK_DEFS_HPP_
#define SCHNEK_DEFS_HPP_

#include <schnek/config.hpp>

#ifdef ONE_DIMENSIONAL

#define SPACE_LOOP(index, lo, hi)                          \
  for (index[0]=lo[0]; index[0]<=hi[0]; ++index[0])

#endif

#ifdef TWO_DIMENSIONAL

#define SPACE_LOOP(index, lo, hi)                          \
  for (index[0]=lo[0]; index[0]<=hi[0]; ++index[0])         \
    for (index[1]=lo[1]; index[1]<=hi[1]; ++index[1])

#endif

#ifdef THREE_DIMENSIONAL

#define SPACE_LOOP(index, lo, hi)                          \
  for (index[0]=lo[0]; index[0]<=hi[0]; ++index[0])         \
    for (index[1]=lo[1]; index[1]<=hi[1]; ++index[1])       \
      for (index[2]=lo[2]; index[2]<=hi[2]; ++index[2])

#endif

#endif // SCHNEK_DEFS_HPP_
