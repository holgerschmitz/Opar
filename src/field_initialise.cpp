/*
 * field_initialise.cpp
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


#include "field_initialise.hpp"
#include "defs.hpp"
#include "globals.hpp"
#include "random.hpp"

#include <schnek/tools/fieldtools.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

void FieldInitialiser::preInit()
{
  dynamic_cast<Fields&>(*this->getParent()).addInitialiser(this);
}

void FieldInitialiser::checkField(std::string name, const DataField &field)
{
  SIntVector lo = field.getLo();
  SIntVector hi = field.getHi();
  SIntVector pos;
  for (pos[0]=lo[0]; pos[0]<=hi[0]; ++pos[0])
#ifndef ONE_DIMENSIONAL
    for (pos[1]=lo[1]; pos[1]<=hi[1]; ++pos[1])
#ifdef THREE_DIMENSIONAL
      for (pos[2]=lo[2]; pos[2]<=hi[2]; ++pos[2])
#endif
#endif
      {
        if (isnan(field[pos]) || isinf(field[pos]))
        {
          std::string val = isnan(field[pos])?"NaN":"Inf";
          std::cerr << "We have detected a "<<val<<" value in field '" << name
              << "' at position " << pos << "\n"
              << "You probably have an error in the formula initialising this field\n";
          exit(-1);
        }
      }
}



void ExpressionFieldInitialiser::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(3)
  EParam = blockPars.addArrayParameter("E", EInit, 0.0);
  BParam = blockPars.addArrayParameter("B", BInit, 0.0);
}

void ExpressionFieldInitialiser::initialiseFields(Fields &species)
{

  pDataField pEx = species.getEx();
  pDataField pEy = species.getEy();
  pDataField pEz = species.getEz();
  pDataField pBx = species.getBx();
  pDataField pBy = species.getBy();
  pDataField pBz = species.getBz();

  SVector &coords = Globals::instance().getX();
  pDependencyUpdater updater = Globals::instance().getUpdater(var_space);

  fill_field(*pEx, coords, EInit[0], *updater, EParam[0]);
  fill_field(*pEy, coords, EInit[1], *updater, EParam[1]);
  fill_field(*pEz, coords, EInit[2], *updater, EParam[2]);

  fill_field(*pBx, coords, BInit[0], *updater, BParam[0]);
  fill_field(*pBy, coords, BInit[1], *updater, BParam[1]);
  fill_field(*pBz, coords, BInit[2], *updater, BParam[2]);

  checkField("Ex",*pEx);
  checkField("Ey",*pEy);
  checkField("Ez",*pEz);

  checkField("Bx",*pBx);
  checkField("By",*pBy);
  checkField("Bz",*pBz);
}


void HdfFieldInitialiser::initParameters(BlockParameters &blockPars)
{
  blockPars.addParameter("Ex", &ExFileName, std::string(""));
  blockPars.addParameter("Ey", &EyFileName, std::string(""));
  blockPars.addParameter("Ez", &EzFileName, std::string(""));

  blockPars.addParameter("Bx", &BxFileName, std::string(""));
  blockPars.addParameter("By", &ByFileName, std::string(""));
  blockPars.addParameter("Bz", &BzFileName, std::string(""));

  blockPars.addParameter("ExBlock", &ExBlockName, std::string("data"));
  blockPars.addParameter("EyBlock", &EyBlockName, std::string("data"));
  blockPars.addParameter("EzBlock", &EzBlockName, std::string("data"));

  blockPars.addParameter("BxBlock", &BxBlockName, std::string("data"));
  blockPars.addParameter("ByBlock", &ByBlockName, std::string("data"));
  blockPars.addParameter("BzBlock", &BzBlockName, std::string("data"));
}

void HdfFieldInitialiser::initialiseFields(Fields &species)
{
  pDataField pEx = species.getEx();
  pDataField pEy = species.getEy();
  pDataField pEz = species.getEz();
  pDataField pBx = species.getBx();
  pDataField pBy = species.getBy();
  pDataField pBz = species.getBz();

  schnek::GridContainer<DataField> container;
  schnek::HdfIStream hdfInput;

  container.global_min = SIntVector(-2);
  container.global_max = Globals::instance().getGlobalGridSize();
  container.global_max += 1;

  if (ExFileName!="") {
    readField(container, hdfInput, pEx, ExFileName, ExBlockName);
  }
  else {
    *pEx = 0.0;
  }

  if (EyFileName!="") {
    readField(container, hdfInput, pEy, EyFileName, EyBlockName);
  }
  else {
    *pEy = 0.0;
  }

  if (EzFileName!="") {
    readField(container, hdfInput, pEz, EzFileName, EzBlockName);
  }
  else {
    *pEz = 0.0;
  }


  if (BxFileName!="") {
    readField(container, hdfInput, pBx, BxFileName, BxBlockName);
  }
  else {
    *pBx = 0.0;
  }

  if (ByFileName!="") {
    readField(container, hdfInput, pBy, ByFileName, ByBlockName);
  }
  else {
    *pBy = 0.0;
  }

  if (BzFileName!="") {
    readField(container, hdfInput, pBz, BzFileName, BzBlockName);
  }
  else {
    *pBz = 0.0;
  }

}

void HdfFieldInitialiser::readField(schnek::GridContainer<DataField> &container,
                                    schnek::HdfIStream hdfInput,
                                    pDataField field,
                                    std::string filename,
                                    std::string blockname)
{
  container.grid = &(*field);
  container.local_min = field->getInnerLo();
  container.local_max = field->getInnerHi();

  hdfInput.open(filename.c_str());
  hdfInput.setBlockName(blockname);
  hdfInput.readGrid(container);
  hdfInput.close();
}
