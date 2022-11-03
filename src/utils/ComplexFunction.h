
#pragma once

#include "mpi_utils.h"
#include "trees/FunctionTree.h"

namespace mrcpp {

struct FunctionData {
    int type{0};
    int order{1};
    int scale{0};
    int depth{0};
    int boxes[3] = {0, 0, 0};
    int corner[3] = {0, 0, 0};
    int real_size{0};
    int imag_size{0};
    bool is_shared{false};
};

class ComplexFunction final {
public:
    explicit ComplexFunction(bool share)
            : shared_mem_re(nullptr)
            , shared_mem_im(nullptr)
            , re(nullptr)
            , im(nullptr) {
        this->func_data.is_shared = share;
        if (this->func_data.is_shared and mpi::share_size > 1) {
            // Memory size in MB defined in input. Virtual memory, does not cost anything if not used.
#ifdef MRCPP_HAS_MPI
            this->shared_mem_re = new mrcpp::SharedMemory(mpi::comm_share, mpi::shared_memory_size);
            this->shared_mem_im = new mrcpp::SharedMemory(mpi::comm_share, mpi::shared_memory_size);
#endif
        }
    }

    ~ComplexFunction() {
        if (this->shared_mem_re != nullptr) delete this->shared_mem_re;
        if (this->shared_mem_im != nullptr) delete this->shared_mem_im;
        if (this->re != nullptr) delete this->re;
        if (this->im != nullptr) delete this->im;
    }

    friend class MPI_Func;

private:
    FunctionData func_data;
    mrcpp::SharedMemory *shared_mem_re;
    mrcpp::SharedMemory *shared_mem_im;
    mrcpp::FunctionTree<3> *re; ///< Real part of function
    mrcpp::FunctionTree<3> *im; ///< Imaginary part of function

    void flushFuncData() {
        this->func_data.real_size = 0;
        this->func_data.imag_size = 0;
        if (this->re != nullptr) {
            this->func_data.real_size = this->re->getNChunksUsed();
            flushMRAData(this->re->getMRA());
        }
        if (this->im != nullptr) {
            this->func_data.imag_size = this->im->getNChunksUsed();
            flushMRAData(this->im->getMRA());
        }
    }

    void flushMRAData(const mrcpp::MultiResolutionAnalysis<3> &mra) {
        const auto &box = mra.getWorldBox();
        this->func_data.type = mra.getScalingBasis().getScalingType();
        this->func_data.order = mra.getOrder();
        this->func_data.depth = mra.getMaxDepth();
        this->func_data.scale = box.getScale();
        this->func_data.boxes[0] = box.size(0);
        this->func_data.boxes[1] = box.size(1);
        this->func_data.boxes[2] = box.size(2);
        this->func_data.corner[0] = box.getCornerIndex().getTranslation(0);
        this->func_data.corner[1] = box.getCornerIndex().getTranslation(1);
        this->func_data.corner[2] = box.getCornerIndex().getTranslation(2);
    }
};

} // namespace mrcpp