#include "CopyAdaptor.h"

#include <tuple>

namespace mrcpp {

template<int D>
CopyAdaptor<D>::CopyAdaptor(FunctionTree<D> &t, int ms, int *bw)
        : TreeAdaptor<D>(ms) {
    setBandWidth(bw);
    tree_vec.push_back(std::make_tuple(1.0, &t));
}

template<int D>
CopyAdaptor<D>::CopyAdaptor(FunctionTreeVector<D> &t, int ms, int *bw)
        : TreeAdaptor<D>(ms),
          tree_vec(t) {
    setBandWidth(bw);
}

template<int D>
void CopyAdaptor<D>::setBandWidth(int *bw) {
    for (int d = 0; d < D; d++) {
        if (bw != 0) {
            this->bandWidth[d] = bw[d];
        } else {
            this->bandWidth[d] = 0;
        }
    }
}

template<int D>
bool CopyAdaptor<D>::splitNode(const MWNode<D> &node) const {
    const NodeIndex<D> &idx = node.getNodeIndex();
    for (int c = 0; c < node.getTDim(); c++) {
        NodeIndex<D> cIdx(idx, c);
        for (int d = 0; d < D; d++) {
            for (int bw = -this->bandWidth[d]; bw <= this->bandWidth[d]; bw++) {
                NodeIndex<D> bwIdx(cIdx);
                bwIdx.getTranslation()[d] += bw;
                for (int i = 0; i < this->tree_vec.size(); i++) {
                    const FunctionTree<D> &func_i = get_func(tree_vec, i);
                    const MWNode<D> *node_i = func_i.findNode(bwIdx);
                    if (node_i != nullptr) return true;
                }
            }
        }
    }
    return false;
}

template class CopyAdaptor<1>;
template class CopyAdaptor<2>;
template class CopyAdaptor<3>;

} // namespace mrcpp
