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
#include "particle_diagnostic.hpp"
#include "species.hpp"

#include "constants.hpp"

#include <schnek/parser.hpp>
#include <schnek/util/logger.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

#include <iostream>
#include <fstream>
#include <cmath>

#undef LOGLEVEL
#define LOGLEVEL 1

void debug_check_out_of_bounds(std::string checkpoint)
{
//  if (GridArgCheck<dimension>::getErrorFlag())
//  {
//    std::cerr << "Index out of bounds at checkpoint " << checkpoint << std::endl;
//    std::cerr << "offending " << GridArgCheck<dimension>::getOffending()[0] << std::endl;
//    std::cerr << "Grid " << Globals::instance().getLocalGridMin()[0] << " "
//        << Globals::instance().getLocalGridMax()[0] << std::endl;
//    std::cerr << "Limits " << Globals::instance().getLocalDomainMin()[0] << " "
//        << Globals::instance().getLocalDomainMax()[0] << std::endl;
//    exit(-1);
//  }
}


class FieldDiagnostic : public schnek::HDFGridDiagnostic<
    typename DataField::BaseType, pDataField>
{
  protected:
    typedef HDFGridDiagnostic<typename DataField::BaseType, pDataField>::IndexType IndexType;
    IndexType getGlobalMin()
    {
      // return IndexType(0);
      return IndexType(-2); // we want to write out the ghost cells
    }
    IndexType getGlobalMax()
    {
      IndexType max = Globals::instance().getGlobalGridSize();
//      max -= 1;
      max += 1;  // we want to write out the ghost cells
      return max;
    }
};

void OPar::initParameters(BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  blockPars.addConstant("pi", M_PI);
  blockPars.addConstant("clight", clight);
  blockPars.addConstant("me", mass_e);
  blockPars.addConstant("mp", mass_p);
  blockPars.addConstant("e", unit_charge);
  blockPars.addConstant("mu0", mu_0);
  blockPars.addConstant("eps0", eps_0);

  Globals::instance().initGlobalParameters(blockPars);
}

void OPar::execute()
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  double dt = Globals::instance().getDt();

  BOOST_FOREACH(Fields *f, fields) f->stepSchemeInit(dt);

  do
  {
    //std::cerr << "Time = " << Globals::instance().getT() << std::endl;
    debug_check_out_of_bounds("A");
    // Advance electromagnetic fields
    BOOST_FOREACH(Fields *f, fields) f->stepScheme(dt);
    debug_check_out_of_bounds("B");

    // Advance particle species
    //std::cerr << "push" << std::endl;
    BOOST_FOREACH(Species *s, species) s->pushParticles(dt);
    debug_check_out_of_bounds("C");

    // run diagnostics
    //std::cerr << "diagnostic" << std::endl;
    DiagnosticManager::instance().execute();
    debug_check_out_of_bounds("D");

    //if ((++n % 10) == 0) { BOOST_FOREACH(Fields *f, fields) f->writeAsTextFiles(n); }
  } while (Globals::instance().stepTime());


}

void OPar::addField(Fields *f)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  fields.push_back(f);
}

void OPar::addSpecies(Species *s)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  species.push_back(s);
}

void initBlockLayout(BlockClasses &blocks)
{

  SCHNEK_TRACE_ENTER_FUNCTION(2)
  blocks.registerBlock("opar");
  blocks("opar").addChildren("Common")("Fields")
      ("Species")("FieldDiagnostic")("ParticleDiagnostic");

  blocks("opar").setClass<OPar>();
  blocks("Common").setClass<CommonBlock>();
  blocks("Fields").setClass<Fields>();
  blocks("Species").setClass<Species>();
  blocks("FieldDiagnostic").setClass<FieldDiagnostic>();
  blocks("ParticleDiagnostic").setClass<ParticleDiagnostic>();

  //blocks("Fields").addChildren("FieldBC")("FieldInit");
  //blocks("Species").addChildren("SpeciesBC")("SpeciesInit");
  //blocks.addBlockClass("Collection").addChildren("Values")("Constants");
}

void initFunctions(FunctionRegistry &freg)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  registerCMath(freg);
}

int main(int argc, char **argv)
{
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
    SCHNEK_TRACE_ERR(1,"Parsing opar.config");
    application = P.parse(in, "opar.config");
    SCHNEK_TRACE_ERR(1,"Done parsing opar.config");
  }
  catch (ParserError &e)
  {
    std::cerr << "Parse error, " << e.atomToken.getFilename() << "(" << e.atomToken.getLine() << "): "<< e.message << "\n";
    exit(-1);
  }

  OPar &opar = dynamic_cast<OPar&>(*application);
  try
  {
    SCHNEK_TRACE_ERR(1,"Initialising Variables");
    opar.initAll();
  }
  catch (VariableNotInitialisedException &e)
  {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    exit(-1);
  }
  catch (std::string &err)
  {
    std::cerr << "FATAL ERROR: >>" << err << "<<" << std::endl;
    exit(-1);
  }

  if (Globals::instance().getSubdivision()->master())
  {
    std::ofstream referencesText("information.tex");
    std::ofstream referencesBib("references.bib");

    schnek::LiteratureManager::instance().writeInformation(referencesText,"references.bib");
    schnek::LiteratureManager::instance().writeBibTex(referencesBib);
    referencesText.close();
    referencesBib.close();
  }
  opar.execute();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

