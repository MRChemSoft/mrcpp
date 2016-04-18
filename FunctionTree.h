/**
*
*
*  \date Aug 14, 2009
*  \author Jonas Juselius <jonas.juselius@uit.no> \n
*  CTCC, University of Tromsø
*
*  Basic class for representing functions in a multiwavelet
* representation.
*/

#ifndef FUNCTIONTREE_H_
#define FUNCTIONTREE_H_

#include "MWTree.h"

template<int D>
class FunctionTree: public MWTree<D> {
public:
    virtual ~FunctionTree();

    void clear();

    double integrate() const;
    virtual double dot(const FunctionTree<D> &ket);
    virtual double evalf(const double *r) const;

    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);

    // In place operations
    void square();
    void power(double d);
    void normalize();
    void orthogonalize(const FunctionTree<D> &tree);
    void map(const RepresentableFunction<1> &func);

    FunctionTree<D>& operator *=(double c);
    FunctionTree<D>& operator *=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator +=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator -=(const FunctionTree<D> &tree);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->getEndMWNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->rootBox.getNode(i)); }

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->getEndMWNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->rootBox.getNode(i)); }

    template<int T>
    friend std::ostream& operator <<(std::ostream &o, FunctionTree<T> &tree);

    friend class GridGenerator<D>;
    friend class GridCleaner<D>;
    friend class MWProjector<D>;
    friend class MWAdder<D>;
    friend class MWMultiplier<D>;
    friend class MWOperator<D>;

protected:
    FunctionTree(const MultiResolutionAnalysis<D> &mra);
    FunctionTree(const MWTree<D> &tree);
    FunctionTree(const FunctionTree<D> &tree);
    FunctionTree<D> &operator=(const FunctionTree<D> &tree);

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT;
    }
};

template<int D>
std::ostream& operator<<(std::ostream &o, FunctionTree<D> &tree) {
    o << std::endl << "*FunctionTree: " << tree.name << std::endl;
    o << "  square norm: " << tree.squareNorm << std::endl;
    o << "  root scale: " << tree.getRootScale() << std::endl;
    o << "  order: " << tree.order << std::endl;
    o << "  nodes: " << tree.getNNodes() << std::endl;
    o << "  local end nodes: " << tree.endNodeTable.size() << std::endl;
    o << "  genNodes: " << tree.getNGenNodes() <<
            " (" << tree.getNAllocGenNodes() << ")" <<std::endl;
    o << "  nodes per scale: " << std::endl;
    for (int i = 0; i < tree.nodesAtDepth.size(); i++) {
        o << "    scale=" << i + tree.getRootScale() << "  nodes="
          << tree.nodesAtDepth[i] << std::endl;
    }
    return o;
}


#endif /* FUNCTIONTREE_H_*/
