/*
 * functions.hpp
 *
 * Created on: 11 Feb 2015
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
 *
 */

/** @file functions.hpp
 *
 *  Contains additional mathematical functions used in the input deck and
 *  across the rest of the code.
 */

#ifndef SRC_FUNCTIONS_HPP_
#define SRC_FUNCTIONS_HPP_

/** The unit step function
 *
 * \f$\text{step}(x) = \begin{cases}1&:\quad x\ge0 \\ 0 &:\quad x<0 \end{cases} \f$
 */
double step(double x);

/** The logistic function
 *
 * \f$\text{logistix}(x) = \frac{1}{1+\exp\left(-\frac{x-x0}{w}\right) \f$
 */
double logistic(double x, double w, double x0);


#endif /* SRC_FUNCTIONS_HPP_ */
