/*
 * test_particle_storage.cpp
 *
 * Created on: 24 Feb 2015
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

#include "../src/particles.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <iostream>

// safe comparison of two floating point numbers
bool compare_double(double a, double b)
{
  return ((a==0.0) && (b==0.0)) ||
      ( fabs(a-b)/(fabs(a)+fabs(b)) < 100*std::numeric_limits<double>::epsilon() );
}

ParticleStorage storage;

boost::random::mt19937 rGen;

const int NRepeat = 1; //100;
const int PartSize = 100; //5000001;

//void diagStorage()
//{
//  for (ParticleStorage::BlockIterator it=storage.blocks.begin();
//      it != storage.blocks.end();
//      ++it)
//  {
//    for (int i=0; i<it->count; ++i) std::cout << '#';
//    for (int i=it->count; i<10; ++i) std::cout << '-';
//    std::cout << ' ';
//  }
//  std::cout<< std::endl;
//}
//
//void diagIter(const ParticleStorage::iterator &it)
//{
//  ParticleStorage::BlockIterator it2=storage.blocks.begin();
//
//  while ((it.blockIter != it2) && (it2 != storage.blocks.end()))
//  {
//    ++it2;
//    std::cout << "           ";
//  }
//
//  if (it2==storage.blocks.end())
//    std::cout << "X";
//  else
//  {
//    for (int i=0; i<it.pos-1; ++i) std::cout << ' ';
//    std::cout << 'i';
//  }
//
//  std::cout<< std::endl;
//}

class FactorSeries
{
  private:
    double f;
  public:
    FactorSeries() : f(1.0) {}
    double operator()()
    {
      return 1.0;
//      f = (f*M_PI);
//      f =  f - int(f);
//      return f;
    }
};

void addParticles()
{
  boost::random::uniform_int_distribution<> randSize(1,PartSize);
  FactorSeries f;
  int numParticles = randSize(rGen);
//  std::cout<< "\nAdding "<< numParticles << " particles\n\n";
//  diagStorage();
  for (int n = 0; n<numParticles; ++n)
  {
    storage.addParticle();
//    diagStorage();
  }
}

void removeParticles()
{
  boost::random::uniform_real_distribution<> mc(0.0,1.0);
  ParticleStorage::iterator it = storage.begin();
//  std::cout<< "\nRemoving\n\n";
//  diagIter(it);
//  diagStorage();
  while (it!=storage.end())
  {
    if (mc(rGen) < 0.5)
    {
//      std::cout << "del\n";
      it = storage.removeParticle(it);
    } else
      ++it;
//    diagIter(it);
//    diagStorage();
  }
}

double updateParticles()
{
  boost::random::uniform_real_distribution<> randValue(-1.0,1.0);
  double sum = 0.0;
  FactorSeries f;

  for (ParticleStorage::iterator it = storage.begin(); it!=storage.end(); ++it)
  {
    for (int i=0; i<dimension; ++i)
    {
      double x = randValue(rGen);
      sum += f()*x;
      it->x[i] = x;
    }

    for (int i=0; i<3; ++i)
    {
      double u = randValue(rGen);
      sum += f()*u;
      it->u[i] = u;
    }

    double w = randValue(rGen);
    sum += f()*w;
    it->weight = w;
  }

  return sum;
}

double readParticles()
{
  double sum = 0.0;
  FactorSeries f;

  for (ParticleStorage::iterator it = storage.begin(); it!=storage.end(); ++it)
  {
    for (int i=0; i<dimension; ++i)
      sum += f()* it->x[i];

    for (int i=0; i<3; ++i)
      sum += f()* it->u[i];

    sum += f()*it->weight;
  }

  return sum;
}


int main()
{
  std::cout << "Testing Particle Storage" << std::endl;
  for (int i=0; i<NRepeat; ++i)
  {
    //std::cout << '.' << std::flush;
    addParticles();
    double sumWrite = updateParticles();
    double sumRead = readParticles();

    if (!compare_double(sumWrite, sumRead))
    {
      std::cerr << "Sums don't match" << std::endl;
      exit(-1);
    }

    removeParticles();
  }
  std::cout << std::endl << "SUCESS" << std::endl;

  return 0;
}

