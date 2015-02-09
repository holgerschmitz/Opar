/*
 * vtkdisplay.cpp
 *
 * Created on: 30 Jan 2015
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

#ifdef OPAR_HAVE_VTK

#include "vtkdisplay.hpp"

template<>
const int TypeToVtk<double>::value = VTK_DOUBLE;

template<>
const int TypeToVtk<float>::value = VTK_FLOAT;

template<>
const int TypeToVtk<int>::value = VTK_INT;

template<>
const int TypeToVtk<short>::value = VTK_SHORT;

template<>
const int TypeToVtk<unsigned short>::value = VTK_UNSIGNED_SHORT;

template<>
const int TypeToVtk<unsigned char>::value = VTK_UNSIGNED_CHAR;

#endif // OPAR_HAVE_VTK
