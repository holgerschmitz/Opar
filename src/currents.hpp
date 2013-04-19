/*
 * currents.hpp
 *
 * Created on: 19 Nov 2012
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

#ifndef CURRENTS_HPP_
#define CURRENTS_HPP_

#include "types.hpp"

#include <schnek/util/singleton.hpp>

using namespace schnek;

class Currents : public Singleton<Currents>
{
  private:
    pDataField jxg, jyg, jzg;

    std::list<pDataField> jxl, jyl, jzl;

    friend class Singleton<Currents>;
    friend class CreateUsingNew<Currents>;

    Currents() {}
    void updateCurrent(pDataField j, const std::list<pDataField> &jl);
  public:
    void setGlobalCurrent(pDataField jx, pDataField jy, pDataField jz);
    void addCurrent(pDataField jx, pDataField jy, pDataField jz);
    void update();
};

#endif // CURRENTS_HPP_ 
