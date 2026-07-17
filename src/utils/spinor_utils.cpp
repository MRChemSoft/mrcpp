#include "spinor_utils.h"
#include "utils/CompFunction.h"
#include "utils/mpi_utils.h"
#include "utils/parallel.h"
// #include "FunctionTreeVector.h"

#include <complex>
#include <iostream>

// remove below
#include "Printer.h"
// till here

using namespace std::complex_literals;

namespace mrcpp {

    /** @brief in-place application of Pauli/Gamma matrices to a CompFunction 
     *  @param inp: input (and output) CompFunction, passed by reference, will be modified
     *  @param index: integer index of the gamma matrix to be applied. 0 is the identity. 5 is not yet implemented
     *  
     * WARNING: 4C/Gamma matrices are not implemented. Only 2C/Pauli matrices are, currently.
     */
    template<int D> void apply_gamma(CompFunction<D> &inp, int index) {
        ComplexDouble comp_i(0.0, 1.0); // Define the imaginary unit for convenience later
        switch (index) {
        case 0:
            //Identity, base case, nothing to apply
            break;
        case 1:
            // Apply sigma_X matrix
            // Basically amounts to swapping first and second components
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                //swapping trees
                std::swap(inp.CompD[i], inp.CompD[i+1]);
                std::swap(inp.CompC[i], inp.CompC[i+1]);
                //swapping tree metadata
                std::swap(inp.func_ptr->data.Nchunks[i], inp.func_ptr->data.Nchunks[i+1]);
                //swapping prefactors
                std::swap(inp.func_ptr->data.c1[i], inp.func_ptr->data.c1[i+1]);
            }
            break;
        case 2:
            // Apply sigma_Y matrix
            // Amounts to swapping first and second components, and multiplying former by i and the latter by -i
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                //swapping trees
                std::swap(inp.CompD[i], inp.CompD[i+1]);
                std::swap(inp.CompC[i], inp.CompC[i+1]);
                //swapping tree metadata
                std::swap(inp.func_ptr->data.Nchunks[i], inp.func_ptr->data.Nchunks[i+1]);
                //swapping prefactors
                std::swap(inp.func_ptr->data.c1[i], inp.func_ptr->data.c1[i+1]);
                //multiplying the 1st and 2nd prefactors by the complex unit ±i
                inp.func_ptr->data.c1[i] *= (-1.0*comp_i);
                inp.func_ptr->data.c1[i+1] *= comp_i;
            }
            break;
        case 3:
            // Apply Alpha-Z matrix
            // nothing to do except multiply the second element's prefactor by -1
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                inp.func_ptr->data.c1[i+1] *= -1.0;
            }
            break;
        default:
            //identity case again. Nothing to apply.
            break;
        }
    }
    
    /** DEPRECATED -
     * @brief: shuffles the indices of a spinor, simulating the application of a Dirac matrix to it
     * pauli represents the index of the Dirac matrices.
     * For scalar operators, it is unused.
     * For 2 component (Weyl/Pauli) spinors, pauli = 0,1,2,3 corresponds to indentiy, sigma_x, y and z respectively.
     * 
    */
    //TODO: Add a provision in case inp and out are identical 
    void apply_Pauli(CompFunction<3> &out, CompFunction<3> &inp, int pauli, double prec, bool conjugate) { //NOTE: assumes 2-component spinors for now
        MSG_WARN("DEPRECATED METHOD - USE apply_gamma INSTEAD");
        // Implementation of applying Pauli matrices to spinor functions
        // This function will modify 'out' based on the Pauli matrix specified by 'pauli'
        // and the input function 'inp'.
        // The 'prec' parameter is used for precision control.
        // The 'conjugate' parameter indicates whether to apply conjugation.
        ComplexDouble comp_i = {0.0, 1.0}; // Define the imaginary unit

        switch (pauli) {
        case 0:
            //Identity, base case, nothing to apply, just copy the input to the output if they are not the same function
            if (&out != &inp) {
                // out.deep_copy(inp);
                out = inp;
            }
            break;
        case 1:
            // Apply Pauli-X matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    inp.CompD[i]->deep_copy(out.CompD[i+1]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i]);
                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    inp.CompC[i]->deep_copy(out.CompC[i+1]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i]);;
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient. 
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1];
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i];
            }
            break;
        case 2: //WARNING: We may need to enforce out to be complex in this case, rather than just multiplying the whole phase by i
            // Apply Pauli-Y matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    inp.CompD[i]->deep_copy(out.CompD[i+1]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i]);


                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    inp.CompC[i]->deep_copy(out.CompC[i+1]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i]);
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient.
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1] * (-1.0)*comp_i;
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i] * comp_i;
            }
            break;
        case 3:
            // Apply Pauli-Z matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    inp.CompD[i]->deep_copy(out.CompD[i]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i+1]);
                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    inp.CompC[i]->deep_copy(out.CompC[i]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i+1]);
                }
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i];
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i+1] * (-1.0);
            }
            break;
        default:
            // std::cerr << "Invalid Pauli matrix index, values must be 0,1,2,3. Current value: " << pauli << std::endl;
            MSG_ABORT("Invalid Pauli matrix index, values must be 0,1,2,3. Current value: " << pauli);
        }
    }

    
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
    ComplexMatrix calc_kramers_overlap_matrix(CompFunctionVector &bra, CompFunctionVector &ket){
        //Useful variables
        int N = bra.size();
        int M = ket.size();

        // ==== Compute the A block 
        ComplexMatrix A = calc_overlap_matrix(bra, ket); //N x N block 

        // ==== Compute the B block
        // -- first compute the time-reversed ket
        // Build the time-reversed ket vector Kket, keeping the same
        // per-index rank ownership as `ket` (mpi::my_func(j)) so that
        // Kket has exactly the "owned vs. freed" shape calc_overlap_matrix
        // already expects and handles under MPI.
        CompFunctionVector Kket(M);
        for (int j=0; j < M; j++){
            //temp variable to hold the time-reversed ket
            CompFunction<3> Kket_j; 
            deep_copy(Kket_j, ket[j]);
            // Complex conjugate the copy of ket[j]
            Kket_j.conj(); //note: only changes a flag that will affect the dot product. Would be problematic if the rest of the time-reversal was imaginary, as it would be conjugated too during the dot. 
            // apply σ_y (2C ONLY) to it
            if (ket[j].Ncomp()>2) MSG_WARN("ONLY IMPLEMENTED FOR 2C!");
            apply_gamma(Kket_j, 2); //σ_y
            //multiply by -i
            ComplexDouble cplx_i = {0.0, 1.0}; 
            Kket_j.func_ptr->data.c1[0] *= -cplx_i;
            Kket_j.func_ptr->data.c1[1] *= -cplx_i;
            Kket[j] = Kket_j;
        }
        ComplexMatrix B = calc_overlap_matrix(bra, Kket);
        // == fill the total overlap matrix
        ComplexMatrix S(2*N, 2*M); //2N x 2N matrix
        for (int i=0; i < N; i++){
            for (int j=0; j < M; j++){
                //fill the upper left block A
                S(i, j) = A(i,j);
                //fill the upper right block B
                S(i, M+j) = B(i,j);
                //fill the lower left block -B^* 
                S(N+i, j) = -1.0 * std::conj(B(i,j));
                //fill the lower right block A^*
                S(N+i, M+j) = std::conj(A(i,j));
            }
        }
        return S;
    }

    
    void normalize_spinor(CompFunction<3> &inp, double prec) {
        // Implementation of normalization for spinor functions
        // This function normalizes the input spinor function 'inp' in place.
 
        double norm = inp.norm();
        if (norm < prec) {
            std::cerr << "normalize_spinor: Norm is too small, cannot normalize." << std::endl;
            return;
        }
        // inp.rescale(1.0 / norm);
        if (inp.isreal() == 1) {
            for (int i = 0; i < inp.Ncomp(); i++) {
                inp.CompD[i]->rescale(1.0 / norm); //Rescaling each component in place
                inp.func_ptr->data.c1[i] = 1.0; // Resetting the overall multiplicative factor to 1 after normalization, since the components have already been rescaled.
                // std::cout << "normalize_spinor: Component " << i << " normalized." << std::endl;
            }
        } else {
            for (int i = 0; i < inp.Ncomp(); i++) {
                inp.CompC[i]->rescale(1.0 / norm); //Rescaling each component in place
                inp.func_ptr->data.c1[i] = ComplexDouble(1.0,0.0); // Resetting the overall multiplicative factor to 1 after normalization, since the components have already been rescaled.
                // std::cout << "normalize_spinor: Component " << i << " normalized." << std::endl;
            }
        }
    }

    template void apply_gamma(CompFunction<3> &inp, int index);
}
