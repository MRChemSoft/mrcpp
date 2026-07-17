#pragma once

// #include "mpi_utils.h"
#include "trees/FunctionTreeVector.h"
#include "utils/CompFunction.h"
// using namespace std::complex_literals;

namespace mrcpp {
    //NOTE: These functions assume that the input functions are spinors, i.e. they have either 2 or 4 components and are complex valued.

    /** @brief in-place application of Pauli/Gamma matrices to a CompFunction 
     *  @param inp: input (and output) CompFunction, passed by reference, will be modified
     *  @param index: integer index of the gamma matrix to be applied. 0 is the identity. 5 is not yet implemented
     *  
     * WARNING: 4C/Gamma matrices are not implemented. Only 2C/Pauli matrices are, currently.
     */
    template<int D> void apply_gamma(CompFunction<D> &inp, int index);

    // orginial implementation of apply_alpha, may be more efficient because everything is done by reference, but it doesn't work in mrchem
    void apply_Pauli(CompFunction<3> &out, CompFunction<3> &inp, int pauli, double prec = -1.0, bool conjugate = false);

    /** @brief Compute the complete overlap matrix of a set of spinors subject to Kramers time-reversal symmetry.
     *    
     *  S = ( A  B )
     *      (-B* A*)
     *  with * denoting the complex conjugate and
     *  A_ij = <bra_i | ket_j> = <K bra_i | K ket_j>^* 
     *  B_ij = <bra_i | K ket_j> = -<K bra_i | ket_j>^* 
     *  and K the anti-unitary time-reversal operator (2C spinor space here) K = -i σ_y Conj 
     *  (Conj being the complex conjugation operator) 
     */
    ComplexMatrix calc_kramers_overlap_matrix(CompFunctionVector &bra, CompFunctionVector &ket);

    /* @brief Normalization of spinor functions. 
     * The function is normalized in place, i.e. the input function is modified. 
     * The norm is computed as the sum of the norms of the components, i.e. ||Psi||^2 = ||Psi_1||^2 + ||Psi_2||^2 for a 2-component spinor. 
     * For a 4-component spinor, the norm is computed as ||Psi||^2 = ||Psi_1||^2 + ||Psi_2||^2 + ||Psi_3||^2 + ||Psi_4||^2. 
     * The function is normalized by rescaling each component by the same factor, i.e. Psi_i -> Psi_i / ||Psi||. 
     * @param inp The input spinor function to be normalized. It is modified in place.
     * @param prec The precision for the normalization. If the norm is smaller than prec, the norm is considered to be zero.
     * NOTE: This function rescales the components of the input function rather than its overall scaling coefficient.
     *       This makes it slower but allows to reset the scaling coefficient to a reasonable value.
     */
    void normalize_spinor(CompFunction<3> &inp, double prec = -1.0);

}