/*
 * BandWidth.h
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <iomanip>

class BandWidth {
public:
    BandWidth(int depth = 0) : widths(depth + 1, 5) { this->clear(); }
    BandWidth &operator=(const BandWidth &bw);
    virtual ~BandWidth() { }

    void clear() { this->widths.setConstant(-1); }

    bool isEmpty(int depth) const;
    int getDepth() const { return this->widths.rows() - 1; }
    int getMaxWidth(int depth) const;
    int getWidth(int depth, int index) const;
    void setWidth(int depth, int index, int wd);

    friend std::ostream& operator<<(std::ostream &o, const BandWidth &bw) {
        o << "  *BandWidths:" << std::endl;
        o << "   n      T   C   B   A  |  max " << std::endl;
        o << " -------------------------------" << std::endl;
        for (int depth = 0; depth <= bw.getDepth(); depth++) {
            o << std::setw(4) << depth << " | ";
            o << std::setw(4) << bw.widths(depth, 0);
            o << std::setw(4) << bw.widths(depth, 1);
            o << std::setw(4) << bw.widths(depth, 2);
            o << std::setw(4) << bw.widths(depth, 3) << "  | ";
            o << std::setw(4) << bw.widths(depth, 4) << std::endl;
        }
        o << std::endl;
        return o;
    }

private:
    Eigen::MatrixXi widths; /// column 5 stores max width at depth
};

