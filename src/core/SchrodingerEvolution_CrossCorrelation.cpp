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

/*
 *
 *
 *  \date Jul 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "SchrodingerEvolution_CrossCorrelation.h"

#include <fstream>

#include "MRCPP/config.h"
#include "MRCPP/constants.h"

#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {


/** @brief SchrodingerEvolution_CrossCorrelation constructor.
 *
 * @param[in] amount: the integer specifying the maximum amount of matrices \f$ C^k \f$
 *                    to be used in calculations
 * @param[in] k: the integer specifying the polynomial order
 * @param[in] t: the integer specifying the scaling basis type
 *
 * @details It checks if the order and type are meaningful and then reads matrices from a file.
 * By default the file has some information about the data stored,
 * so the first interger to read is describing the size of the documentation text.
 * 
 * 
 */
SchrodingerEvolution_CrossCorrelation::SchrodingerEvolution_CrossCorrelation(int amount, int k, int t)
    : type(t), order(k), amount(amount)
{
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order: " << this->order);
    switch (this->type) {
        case (Interpol):
            MSG_ERROR("Not implemented yet filter type: " << this->type);
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    setCCCPath(details::find_filters());

    readCCCBin();
}


void SchrodingerEvolution_CrossCorrelation::setCCCPath(const std::string &lib) {
    switch (this->type) {
        case (Interpol):
            MSG_ERROR("Not implemented yet filter type: " << this->type);
            break;
        case (Legendre):
            this->path = lib + "/Schrodinger_evolution_cross_correlation_coefficients_Legendre_scaling_type.bin";
            break;
        default:
            MSG_ERROR("Invalid CrossCorrelation type");
    }
}

void SchrodingerEvolution_CrossCorrelation::readCCCBin()
{
    std::ifstream input_file(this->path.c_str(), std::ios::binary);

    if (not input_file) MSG_ABORT("Could not open cross correlation: " << this->path);

    // Read the text length
    int text_length;
    input_file.read(reinterpret_cast<char*>(&text_length), sizeof(text_length));

    // Read the Unicode characters
    std::vector<char32_t> unicode_chars(text_length);
    input_file.read(reinterpret_cast<char*>(unicode_chars.data()), sizeof(char32_t) * text_length);

    // Read the amount of matrices
    int K;
    input_file.read(reinterpret_cast<char*>(&K), sizeof(K));

    // Read the size/order of each matrix
    int order;
    input_file.read(reinterpret_cast<char*>(&order), sizeof(order));

    // Read the matrices
    std::vector<Eigen::MatrixXd> C_even(K, Eigen::MatrixXd(order, order));
    auto data_amount = order * order * sizeof(double);
    for (auto& matrix : C_even) input_file.read(reinterpret_cast<char*>(matrix.data()), data_amount);
/*

    // Print the text length
    std::cout << text_length << std::endl;
    
    // Print the text
    for (char32_t c : unicode_chars) {
        std::wcout << static_cast<wchar_t>(c);
    }
    // Print the matrices
    std::cout << std::endl;
    std::cout << "----------------------------------" << std::endl;
    for (auto& matrix : C_even)
    {
        std::cout << matrix  << std::endl;
        std::cout << "----------------------------------" << std::endl;
    }
*/
    // Create matrix containing the appropriate amount of coefficients
    int Order = this->order + 1;
    for (int k = 0; k < this->amount; k++)
    {
        this->Matrix.push_back( C_even[k].block(0, 0, Order, Order).eval() );
    }

    C_even.clear();
    input_file.close();
}

} // namespace mrcpp
