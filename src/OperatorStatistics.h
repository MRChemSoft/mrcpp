#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <iomanip>

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class OperatorStatistics {
public:
    OperatorStatistics();
    virtual ~OperatorStatistics();

    void flushNodeCounters();
    void incrementFNodeCounters(const MWNode<D> &fNode, int ft, int gt);
    void incrementGNodeCounters(const MWNode<D> &gNode);

    friend std::ostream& operator<<(std::ostream &o, const OperatorStatistics &os) { return os.print(o); }

protected:
    int nThreads;
    int totFCount;
    int totGCount;
    int totGenCount;
    int *fCount;
    int *gCount;
    int *genCount;
    Eigen::Matrix<int, 8, 8> *totCompCount;
    Eigen::Matrix<int, 8, 8> **compCount;

    std::ostream& print(std::ostream &o) const;
};

}
