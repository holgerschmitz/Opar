/*
 * opar.cpp
 *
 * Created on: 30 Jul 2012
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2012-2023 Holger Schmitz
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

#include <schnek/config.hpp>
#include "opar.hpp"
// #include "particle_diagnostic.hpp"
// #include "species.hpp"
#include "currents.hpp"

#include "functions.hpp"
#include "random.hpp"
#include "../huerto/types.hpp"
#include "../huerto/constants.hpp"
#include "../huerto/maths/functions/core.hpp"
// #include "../huerto/diagnostic/field_diagnostic.hpp"
#include "../huerto/electromagnetics/em_fields.hpp"
#include "../huerto/electromagnetics/fdtd/fdtd_plain.hpp"
// #include "../huerto/electromagnetics/pml/cpml_border.hpp"

#include <schnek/parser.hpp>
#include <schnek/util/logger.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

#include <iostream>
#include <fstream>
#include <cmath>
#include <fenv.h>


#undef LOGLEVEL
#define LOGLEVEL 1

void debug_check_out_of_bounds(std::string)
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


// typedef FieldDiagnostic<DataField, schnek::DeltaTimeDiagnostic> OparFieldDiagnostic;

void OPar::initParameters(schnek::BlockParameters &blockPars)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  SimulationContext::initParameters(blockPars);
  initConstantParameters(blockPars);
  t_parameter  = blockPars.addParameter("t", &time, schnek::BlockParameters::readonly);

  blockPars.addArrayParameter("N", gridSize, 100);
  blockPars.addArrayParameter("L", size, 1.0);
  blockPars.addParameter("tMax", &tMax, 1.0);
  blockPars.addParameter("cflFactor", &cflFactor, 0.99);
  blockPars.addParameter("ignore_initial_time_stagger", (int*)(&ignore_initial_time_stagger), 0);
}

void OPar::init()
{
  // Currents::instance().setContext(this);

  double minDx = std::numeric_limits<double>::max();
  for (std::size_t i=0; i < DIMENSION; ++i) {
    minDx = std::min(size[i] / gridSize[i], minDx);
  }
  dt = cflFactor*minDx/clight;
}

void OPar::registerData()
{
  for (std::size_t i=0; i < DIMENSION; ++i) {
    dx[i] = size[i] / gridSize[i];
  }

  decomposition = std::make_shared<HuertoDecomposition>();
  decomposition->setGlobalRange({{0}, getGridSize()});
  decomposition->setGlobalDomain({{0.0} , getSize()});
  decomposition->init();
}

schnek::pDependencyUpdater OPar::getUpdater(VarGroup gr)
{
  schnek::pDependencyMap depMap(new schnek::DependencyMap(getVariables()));
  schnek::pDependencyUpdater updater(new schnek::DependencyUpdater(depMap));
  switch (gr)
  {
    case var_space:
      updater->addIndependentArray(x_parameters);
      break;
    case var_time:
      updater->addIndependent(t_parameter);
      break;
    case var_spacetime:
      updater->addIndependentArray(x_parameters);
      updater->addIndependent(t_parameter);
      break;
    case var_none:
    default:
      break;
  }
  return updater;
}

void OPar::execute()
{
  time = 0.0;
  timeStep = 0;

  Random::seed(decomposition->getUniqueId());

  SCHNEK_TRACE_ENTER_FUNCTION(2)
  schnek::DiagnosticManager::instance().setPhysicalTime(&time);
  schnek::DiagnosticManager::instance().setMaster(decomposition->master());

  double dt = getDt();

  if (!ignore_initial_time_stagger) {
    // first half time-step for the EM fields
    for(pFieldSolver f: schnek::BlockContainer<FieldSolver>::childBlocks())
    {
      f->stepSchemeInit(dt);
    }
  }

  schnek::DiagnosticManager::instance().execute();

  while (time<=tMax)
  {
    // run diagnostics
    debug_check_out_of_bounds("A");

    if (getDecomposition().master())
      schnek::Logger::instance().out() << "Time " << getTime() << std::endl;

    //std::cerr << "Time = " << getTime() << std::endl;
    debug_check_out_of_bounds("B");

    // Advance electromagnetic fields
    for(pFieldSolver f: schnek::BlockContainer<FieldSolver>::childBlocks())
    {
      f->stepScheme(dt);
    }

    debug_check_out_of_bounds("C");

    // Advance particle species
    //std::cerr << "push" << std::endl;
    // for (Species *s: species) s->pushParticles(dt);
    debug_check_out_of_bounds("D");

    time += dt;
    ++timeStep;
    // schnek::DiagnosticManager::instance().execute();
  }

  // run diagnostics
  //std::cerr << "diagnostic" << std::endl;
  schnek::DiagnosticManager::instance().execute();
  debug_check_out_of_bounds("E");

}

// void OPar::addSpecies(Species *s)
// {
//   SCHNEK_TRACE_ENTER_FUNCTION(2)
//   species.push_back(s);
// }

void initBlockLayout(schnek::BlockClasses &blocks)
{

  SCHNEK_TRACE_ENTER_FUNCTION(2)
  blocks.registerBlock("opar");

  blocks("opar").setClass<OPar>();
  blocks("EMFields").setClass<EMFields>();
  blocks("FDTD").setClass<FDTD_Plain>();
  // blocks("CPMLBorder").setClass<CPMLBorder>();

  // blocks("Species").setClass<Species>();
  // blocks("FieldDiagnostic").setClass<OparFieldDiagnostic>();
  // blocks("ParticleDiagnostic").setClass<ParticleDiagnostic>();

  blocks("opar").addChildren("EMFields")("FDTD");
  //     ("Species")("FieldDiagnostic")("ParticleDiagnostic");

  // blocks("FDTD").addChildren("CPMLBorder");
  //blocks("Fields").addChildren("FieldBC")("FieldInit");
  //blocks("Species").addChildren("SpeciesBC")("SpeciesInit");
  //blocks.addBlockClass("Collection").addChildren("Values")("Constants");
}

void initFunctions(schnek::FunctionRegistry &freg)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)
  registerCMath(freg);
  freg.registerFunction("logistic", logistic);
  freg.registerFunction("pulse1d", pulse1d);
}

/** Runs the OPar simulation code
 *
 *  This is placed outside the main function so that we can return on error and still
 *  be sure that MPI is properly closed down.
 */
int runOpar(int, char **)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

  schnek::VariableStorage vars("opar", "opar");
  schnek::FunctionRegistry freg;
  schnek::BlockClasses blocks;

  initFunctions(freg);
  registerCoreFunctions(freg);
  initBlockLayout(blocks);

  schnek::Parser P(vars, freg, blocks);
  schnek::pBlock application;

  std::ifstream in("opar.setup");
  if (!in) {
    std::cerr << "Could not open file\n";
    return -1;
  }
  try
  {
    SCHNEK_TRACE_ERR(1,"Parsing opar.setup");
    application = P.parse(in, "opar.setup");
    SCHNEK_TRACE_ERR(1,"Done parsing opar.setup");
  }
  catch (schnek::ParserError &e)
  {
    std::cerr << "Parse error, " << e.atomToken.getFilename() << "(" << e.atomToken.getLine() << "): "<< e.message << "\n";
    return -1;
  }

  OPar &opar = dynamic_cast<OPar&>(*application);
  try
  {
    SCHNEK_TRACE_ERR(1,"Initialising Variables");
    opar.initAll();
  }
  catch (schnek::VariableNotInitialisedException &e)
  {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    return -1;
  }
  catch (schnek::VariableNotFoundException &e)
  {
    std::cerr << "Variable was not found: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (schnek::EvaluationException &e)
  {
    std::cerr << "Error in evaluation: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (std::string &err)
  {
    std::cerr << "FATAL ERROR: >>" << err << "<<" << std::endl;
    return -1;
  }

  if (opar.getDecomposition().master())
  {
    std::ofstream referencesText("information.tex");
    std::ofstream referencesBib("references.bib");

    schnek::LiteratureManager::instance().writeInformation(referencesText,"references.bib");
    schnek::LiteratureManager::instance().writeBibTex(referencesBib);
    referencesText.close();
    referencesBib.close();
  }
  opar.execute();
  return 0;
}

int main(int argc, char **argv)
{
  SCHNEK_TRACE_ENTER_FUNCTION(2)

#ifdef SCHNEK_HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

  int returnCode = runOpar(argc, argv);

#ifdef SCHNEK_HAVE_MPI
    MPI_Finalize();
#endif

  return returnCode;
}

