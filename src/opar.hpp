/*
 * opar.hpp
 *
 * Created on: 30 Jul 2012
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

#ifndef OPAR_HPP_
#define OPAR_HPP_

#include <schnek/variables.hpp>

using namespace schnek;


class FieldsImplementation;

class OPar;

class OParImplementation
{
  private:
    friend class OPar;
    FieldsImplementation *fields;
  public:
    OParImplementation() : fields(0) {}
    void execute();
};

class OPar : public Block
{
  private:
    OParImplementation impl;
  protected:
    void initParameters(BlockParameters &blockPars);
  public:
    void execute();
//    void init();
    void addField(FieldsImplementation *fields);
};


#endif /* OPAR_HPP_ */
