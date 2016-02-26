/*
 * field_initialise.hpp
 *
 * Created on: 20 May 2015
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
 *
 */

#ifndef SRC_FIELD_INITIALISE_HPP_
#define SRC_FIELD_INITIALISE_HPP_


#include "types.hpp"
#include "fields.hpp"

#include <schnek/variables.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class FieldInitialiser : public schnek::Block
{
  protected:
    void preInit();
    void checkField(std::string name, const DataField &field);
  public:
    virtual void initialiseFields(Fields &species) = 0;
};

class ExpressionFieldInitialiser : public FieldInitialiser
{
  private:
    PVector EInit;
    PVector BInit;

    PParameterVector EParam;
    PParameterVector BParam;
  protected:
    void initParameters(BlockParameters &blockPars);
  public:
    void initialiseFields(Fields &species);
};

class HdfFieldInitialiser : public FieldInitialiser
{
  private:
    std::string ExFileName;
    std::string EyFileName;
    std::string EzFileName;
    std::string BxFileName;
    std::string ByFileName;
    std::string BzFileName;

    std::string ExBlockName;
    std::string EyBlockName;
    std::string EzBlockName;
    std::string BxBlockName;
    std::string ByBlockName;
    std::string BzBlockName;

    void readField(schnek::GridContainer<DataField> &container,
                   schnek::HdfIStream hdfInput,
                   pDataField field,
                   std::string filename,
                   std::string blockname);

  protected:
    void initParameters(BlockParameters &blockPars);
  public:
    void initialiseFields(Fields &species);
};

#endif /* SRC_FIELD_INITIALISE_HPP_ */
