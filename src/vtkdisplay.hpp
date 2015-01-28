/*
 * vtkdisplay.hpp
 *
 * Created on: 20 Jan 2015
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

#ifndef SRC_VTKDISPLAY_HPP_
#define SRC_VTKDISPLAY_HPP_


#include <schnek/diagnostic/diagnostic.hpp>

#include <boost/shared_ptr.hpp>

template<typename Type, typename PointerType = boost::shared_ptr<Type> >
class VTKGridDisplay : public schnek::SimpleDiagnostic<Type, PointerType> {
  protected:

  protected:
    typedef typename Type::IndexType IndexType;
    void write();
    void close();
    void init();
  public:
    virtual ~VTKGridDisplay() {}
};





#endif /* SRC_VTKDIAGNOSTIC_HPP_ */
