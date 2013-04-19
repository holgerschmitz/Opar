/*
 * currents.cpp
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

#include <boost/variant.hpp>
#include "currents.hpp"

#include <schnek/grid/domain.hpp>

#include <boost/foreach.hpp>

void Currents::setGlobalCurrent(pDataField jx, pDataField jy, pDataField jz)
{
  jxg = jx;
  jyg = jy;
  jzg = jz;
}

void Currents::addCurrent(pDataField jx, pDataField jy, pDataField jz)
{
  jxl.push_back(jx);
  jyl.push_back(jy);
  jzl.push_back(jz);
}

void Currents::updateCurrent(pDataField j, const std::list<pDataField> &jl)
{
  (*j) = 0.0;

  BOOST_FOREACH(pDataField jc, jl)
  {
    RecDomain<dimension> domain(jc->getLo(), jc->getHi());
    for (RecDomain<dimension>::iterator it = domain.begin(); it != domain.end(); ++it)
    //BOOST_FOREACH(RecDomain<dimension>::LimitType pos, domain)
    {
      (*j)[*it] += (*jc)[*it];
    }
  }
}

void Currents::update()
{
  updateCurrent(jxg, jxl);
  updateCurrent(jyg, jyl);
  updateCurrent(jzg, jzl);
}

