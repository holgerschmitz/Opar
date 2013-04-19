/*
 * common.hpp
 *
 * Created on: 1 Aug 2012
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

#ifndef OPAR_COMMON_HPP_
#define OPAR_COMMON_HPP_

#include <schnek/variables.hpp>

//using namespace schnek;

class CommonBlock : public schnek::Block
{
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void preInit();
};

#endif /* OPAR_COMMON_HPP_ */
