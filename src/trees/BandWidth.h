/*
 * BandWidth.h
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <iomanip>

namespace mrcpp {

class BandWidth final {
public:
    BandWidth(int depth = 0) : widths(depth + 1, 5) { this->clear(); }
    BandWidth(const BandWidth &bw) : widths(bw.widths) { }
    BandWidth &operator=(const BandWidth &bw);

    void clear() { this->widths.setConstant(-1); }

    bool isEmpty(int depth) const;
    int getDepth() const { return this->widths.rows() - 1; }
    int getMaxWidth(int depth) const;
    int getWidth(int depth, int index) const;
    void setWidth(int depth, int index, int wd);

    friend std::ostream& operator<<(std::ostream &o, const BandWidth &bw) { bw.print(o); }

private:
    Eigen::MatrixXi widths; /// column 5 stores max width at depth

    std::ostream& print(std::ostream &o) const;
};

}
