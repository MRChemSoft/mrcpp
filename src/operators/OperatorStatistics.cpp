#include "OperatorStatistics.h"
#include "trees/MWNode.h"

using namespace Eigen;

namespace mrcpp {

template<int D>
OperatorStatistics<D>::OperatorStatistics()
        : nThreads(omp_get_max_threads()),
          totFCount(0),
          totGCount(0),
          totGenCount(0),
          fCount(0),
          gCount(0),
          genCount(0),
          totCompCount(0),
          compCount(0) {

    this->totCompCount = new Matrix<int, 8, 8>;
    this->totCompCount->setZero();

    this->fCount = new int[this->nThreads];
    this->gCount = new int[this->nThreads];
    this->genCount = new int[this->nThreads];
    this->compCount = new Matrix<int, 8, 8> *[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->compCount[i] = new Matrix<int, 8, 8>;
        this->compCount[i]->setZero();
        this->fCount[i] = 0;
        this->gCount[i] = 0;
        this->genCount[i] = 0;
    }
}

template<int D>
OperatorStatistics<D>::~OperatorStatistics() {
    for (int i = 0; i < this->nThreads; i++) {
        delete this->compCount[i];
    }
    delete[] this->compCount;
    delete[] this->fCount;
    delete[] this->gCount;
    delete[] this->genCount;
    delete totCompCount;
}

/** Sum all node counters from all threads. */
template<int D>
void OperatorStatistics<D>::flushNodeCounters() {
    for (int i = 0; i < this->nThreads; i++) {
        this->totFCount += this->fCount[i];
        this->totGCount += this->gCount[i];
        this->totGenCount += this->genCount[i];
        *this->totCompCount += *this->compCount[i];
        this->fCount[i] = 0;
        this->gCount[i] = 0;
        this->genCount[i] = 0;
        this->compCount[i]->setZero();
    }
}

/** Increment g-node usage counter. Needed for load balancing. */
template<int D>
void OperatorStatistics<D>::incrementGNodeCounters(const MWNode<D> &gNode) {
    int thread = omp_get_thread_num();
    this->gCount[thread]++;
}

/** Increment operator application counter. */
template<int D>
void OperatorStatistics<D>::incrementFNodeCounters(const MWNode<D> &fNode,
                                                   int ft, int gt) {
    int thread = omp_get_thread_num();
    this->fCount[thread]++;
    (*this->compCount[thread])(ft, gt) += 1;
    if (fNode.isGenNode()) {
        this->genCount[thread]++;
    }
}

template<int D>
std::ostream& OperatorStatistics<D>::print(std::ostream &o) const {
    o << std::setw(8);
    o << "*OperatorFunc statistics: " << std::endl << std::endl;
    o << "  Total calculated gNodes      : " << this->totGCount << std::endl;
    o << "  Total applied fNodes         : " << this->totFCount << std::endl;
    o << "  Total applied genNodes       : " << this->totGenCount << std::endl << std::endl;
    o << "  By components:" << std::endl << *this->totCompCount << std::endl;
    return o;
}

template class OperatorStatistics<1>;
template class OperatorStatistics<2>;
template class OperatorStatistics<3>;

} // namespace mrcpp
