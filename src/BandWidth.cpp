#include "BandWidth.h"
#include "Printer.h"

BandWidth& BandWidth::operator=(const BandWidth &bw) {
    this->widths = bw.widths;
    return *this;
}

bool BandWidth::isEmpty(int depth) const {
    if (depth > getDepth()) {
        return true;
    }
    if (this->widths(depth, 4) < 0) {
        return true;
    }
    return false;
}

int BandWidth::getWidth(int depth, int index) const {
    assert(depth >= 0);
    assert(index >= 0 and index < 4);
    if (depth > getDepth()) { // No bw at requested depth
        return -1;
    }
    return this->widths(depth, index);
}

int BandWidth::getMaxWidth(int depth) const {
    assert(depth >= 0);
    if (depth > getDepth()) { // No bw at requested depth
        return -1;
    }
    return this->widths(depth, 4);
}

void BandWidth::setWidth(int depth, int index, int wd) {
    assert(depth >= 0 and depth < getDepth());
    assert(index >= 0 and index < 4);
    assert(wd >= 0);
    this->widths(depth, index) = wd;
    if (wd > this->widths(depth, 4)) {
        this->widths(depth, 4) = wd;
    }
}
