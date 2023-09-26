/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "JpowerIntegrals.h"


namespace mrcpp {


JpowerIntegrals::JpowerIntegrals(double a, int scaling, int M, double treshold)
{
    this->scaling = scaling;
    int N = 1 << scaling;
    for(int l = 0; l < N; l++  )
        integrals.push_back( calculate_J_power_inetgarls(l, a, M, treshold) );
    for(int l = 1 - N; l < 0; l++  )
        integrals.push_back( calculate_J_power_inetgarls(l, a, M, treshold) );
}

JpowerIntegrals::JpowerIntegrals(const JpowerIntegrals &other)
{
    scaling = other.scaling;

    // Copy integrals (deep copy)
    integrals.clear();
    for (const auto& row : other.integrals)
    {
        std::vector<std::complex<double>> newRow(row);
        integrals.push_back(newRow);
    }    
}

JpowerIntegrals::~JpowerIntegrals()
{
    //for (auto& integral : integrals) 
    //{
    //    integral.clear();  
    //}
    // The loop is probably not needed for the complete memory deallocation
    integrals.clear();
}

std::vector<std::complex<double>> & JpowerIntegrals::operator[](int index)
{
    if( index < 0 ) index += integrals.size();
    return integrals[index];
}

std::vector<std::complex<double>> JpowerIntegrals::calculate_J_power_inetgarls(int l, double a, int M, double treshold)
{
    const double pi = 3.1415926535897932384626433832795;
    const std::complex<double> I(0.0, 1.0);  // Imaginary unit
    const std::complex<double> J_0 = 0.25 * std::exp(-0.25 * I * pi) / std::sqrt(pi * a) * std::exp(0.25 * I * static_cast<double>(l * l) / a);
    const std::complex<double> beta(0, 0.5 / a);
    const std::complex<double> alpha = static_cast<double>(l) * beta;

    std::vector<std::complex<double>> J = {0.0, J_0};

    for (int m = 0; m < M; m++)
    {
        std::complex<double> term1 = J[J.size() - 1] * alpha;
        std::complex<double> term2
        = J[J.size() - 2] * beta * static_cast<double>(m) / static_cast<double>(m + 2);
        std::complex<double> last = (term1 + term2) / static_cast<double>(m + 3);
        if ( last.real() <  treshold && last.imag() <  treshold && last.real() > -treshold && last.imag() > -treshold ) break;
        J.push_back(last);
    }

    J.erase(J.begin());
    return J;
}

} // namespace mrcpp
