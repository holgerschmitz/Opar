/*
 * opar.cpp
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

// TODO Create JavaDoc comments everywhere

#include "config.hpp"
#include "opar.hpp"
#include "common.hpp"
#include "fields.hpp"
#include "globals.hpp"
#include "species.hpp"

#include <schnek/parser.hpp>
#include <schnek/util/logger.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

#include <iostream>
#include <fstream>
#include <cmath>

#undef LOGLEVEL
#define LOGLEVEL 1


class FieldDiagnostic : public schnek::HDFGridDiagnostic<DataField, pDataField>
{
  protected:
    FieldIndex getGlobalMin() { return FieldIndex(0); }
    FieldIndex getGlobalMax() { return Globals::instance().getGlobalGridSize(); }
};

void OPar::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  //schnek::Logger::instance().out() << "Entering " << BOOST_CURRENT_FUNCTION << std::endl;
  Globals::instance().initGlobalParameters(blockPars);
}

void OPar::execute()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  double dt = Globals::instance().getDt();

  BOOST_FOREACH(Fields *f, fields) f->stepSchemeInit(dt);

  do
  {
    // Advance electromagnetic fields
    BOOST_FOREACH(Fields *f, fields) f->stepScheme(dt);

    // run diagnostics
    DiagnosticManager::instance().execute();

    //if ((++n % 10) == 0) { BOOST_FOREACH(Fields *f, fields) f->writeAsTextFiles(n); }
  } while (Globals::instance().stepTime());


}

void OPar::addField(Fields *f)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  fields.push_back(f);
}

void OPar::addSpecies(Species *s);
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  species.push_back(s);
}

void initBlockLayout(BlockClasses &blocks)
{

  SCHNEK_TRACE_ENTER_FUNCTION(2)
  blocks.addBlockClass("opar");
  blocks("opar").addChildren("Common")("Fields")("Species")("FieldDiagnostic");
  blocks("opar").setBlockClass<OPar>();

  blocks("Common").setBlockClass<CommonBlock>();

  blocks("Fields").addChildren("FieldBC")("FieldInit");
  blocks("Fields").setBlockClass<Fields>();

  blocks("Species").addChildren("SpeciesBC")("SpeciesInit");

  blocks("FieldDiagnostic").setBlockClass<FieldDiagnostic>();

  //blocks.addBlockClass("Collection").addChildren("Values")("Constants");
}

void initFunctions(FunctionRegistry &freg)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  registerCMath(freg);
}

void test_species()
{
  Species sp;

}

int main(int argc, char **argv)
{
  test_species();
  exit(0);

  SCHNEK_TRACE_ENTER_FUNCTION(2)

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

  VariableStorage vars("opar", "opar");
  FunctionRegistry freg;
  BlockClasses blocks;

  initFunctions(freg);
  initBlockLayout(blocks);

  Parser P(vars, freg, blocks);
  pBlock application;

  Globals::instance().setup(vars);

  std::ifstream in("opar.config");
  if (!in) {
    std::cerr << "Could not open file\n";
    exit(-1);
  }
  try
  {
    application = P.parse(in, "test_parser_sample.txt");
  }
  catch (ParserError &e)
  {
    std::cerr << "Parse error, " << e.atomToken.getFilename() << "(" << e.atomToken.getLine() << "): "<< e.message << "\n";
    exit(-1);
  }

  OPar &opar = dynamic_cast<OPar&>(*application);
  try
  {
    opar.initAll();
  }
  catch (VariableNotInitialisedException e)
  {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    exit(-1);
  }
  opar.execute();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

