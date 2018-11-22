#pragma once

namespace mrcpp {

template<int D>
class HilbertPath final {
public:
    HilbertPath() = default;
    HilbertPath(const HilbertPath<D> &p) : path(p.path) { }
    HilbertPath(const HilbertPath<D> &p, int cIdx) {
        int hIdx = p.getHIndex(cIdx);
        this->path = p.getChildPath(hIdx);
    }
    HilbertPath &operator=(const HilbertPath<D> &p) {
        this->path = p.path;
        return *this;
    }

    short int getPath() const { return this->path; }
    short int getChildPath(int hIdx) const { return this->pTable[this->path][hIdx]; }

    int getZIndex(int hIdx) const { return this->zTable[this->path][hIdx]; }
    int getHIndex(int zIdx) const { return this->hTable[this->path][zIdx]; }

private:
    short int path{0};
    static const short int pTable[][8];
    static const int zTable[][8];
    static const int hTable[][8];
};

}
