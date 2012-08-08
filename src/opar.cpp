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

#include "opar.hpp"

#include <schnek/parser.hpp>

#include <iostream>
#include <fstream>
#include <cmath>


void OPar::initParameters(BlockParameters &blockPars)
{}

void OPar::execute()
{
  this->evaluateParameters();

}

void initBlockLayout(BlockClasses &blocks)
{

  blocks.addBlockClass("opar");
  blocks("opar").addChildren("Common")("Fields")("Species");
  blocks("opar").setBlockClass<OPar>();

  blocks("Fields").addChildren("FieldBC")("FieldInit");
  blocks("Species").addChildren("SpeciesBC")("SpeciesInit");

  //blocks.addBlockClass("Collection").addChildren("Values")("Constants");
}

void initFunctions(FunctionRegistry &freg)
{
  registerCMath(freg);
}

int main()
{
  VariableStorage vars("opar", "opar");
  FunctionRegistry freg;
  BlockClasses blocks;

  initFunctions(freg);
  initBlockLayout(blocks);

  Parser P(vars, freg, blocks);
  pBlock application;


  std::ifstream in("test_parser_sample.txt");
  if (!in) {
    std::cerr << "Could not open file\n";
  }
  try
  {
    application = P.parse(in, "test_parser_sample.txt");
  }
  catch (ParserError &e)
  {
    std::cerr << "Parse error, " << e.atomToken.getFilename() << "(" << e.atomToken.getLine() << "): "<< e.message << "\n";
    throw -1;
  }

  application->execute();

}
