/*
 * random.hpp
 *
 * Created on: 3 Dec 2012
 * Author: hschmitz
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

#ifndef RANDOM_HPP_
#define RANDOM_HPP_

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class Random
{
  private:
    static boost::random::mt11213b generator_;
    static boost::random::uniform_01<real> uniform_;
    static boost::random::normal_distribution<real> normal_;

  public:
    static double uniform()
    {
      return uniform_(generator_);
    }

    static double uniform(double range)
    {
      return range*uniform_(generator_);
    }

    static double gaussian()
    {
      return normal_(generator_);
    }

    static double gaussian(double stddev)
    {
      return stddev*normal_(generator_);
    }
};

#endif // RANDOM_HPP_ 
