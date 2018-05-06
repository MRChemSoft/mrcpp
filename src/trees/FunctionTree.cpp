/**
 *  \date 2016
 *          CTCC, University of Tromsø
 */

#include <fstream>

#include "FunctionTree.h"
#include "SerialFunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"
#include "HilbertIterator.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

using namespace std;
using namespace Eigen;

namespace mrcpp {

/** FunctionTree constructor for Serial Tree.
  * */
template<int D>
FunctionTree<D>::FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem)
        : MWTree<D> (mra) {
    this->serialTree_p = new SerialFunctionTree<D>(this, sh_mem);
    this->serialTree_p->allocRoots(*this);
    this->resetEndNodeTable();
}

//** FunctionTree destructor. */
template<int D>
FunctionTree<D>::~FunctionTree() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
    delete this->serialTree_p;
}

/** Leaves the tree inn the same state as after construction, e.i.
  * undefined function containing only root nodes without coefficients.
  * The assigned memory (nodeChunks in SerialTree) is NOT released,
  * but is immediately available to the new function. */
template<int D>
void FunctionTree<D>::clear() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.clearHasCoefs();
        root.clearNorms();
    }
    this->resetEndNodeTable();
    this->clearSquareNorm();
}

/** Write the tree structure to disk, for later use.
  * Argument file name will get a ".tree" file extension. */
template<int D>
void FunctionTree<D>::saveTree(const string &file) {
    // This is basically a copy of MPI send_tree
    Timer t1;
    stringstream fname;
    fname << file << ".tree";

    fstream f;
    f.open(fname.str(), ios::out | ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");

    this->deleteGenerated();
    SerialFunctionTree<D> &sTree = *this->getSerialFunctionTree();

    // Write size of tree
    int nChunks = sTree.getNChunksUsed();
    f.write((char*) &nChunks, sizeof(int));

    // Write tree data, chunk by chunk
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        //set serialIx of the unused nodes to -1
        int iShift = iChunk*sTree.maxNodesPerChunk;
        for (int i = 0; i < sTree.maxNodesPerChunk; i++) {
            if (sTree.nodeStackStatus[iShift+i] != 1) {
                sTree.nodeChunks[iChunk][i].setSerialIx(-1);
            }
        }
        f.write((char *) sTree.nodeChunks[iChunk], count);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        f.write((char *) sTree.nodeCoeffChunks[iChunk], count*sizeof(double));
    }
    f.close();

    t1.stop();
    Printer::printTime(10, "Time write", t1);
}

/** Read a previously stored tree structure from disk.
  * Argument file name will get a ".tree" file extension. */
template<int D>
void FunctionTree<D>::loadTree(const string &file) {
    // This is basically a copy of MPI recv_tree
    Timer t1;
    stringstream fname;
    fname << file << ".tree";

    fstream f;
    f.open(fname.str(), ios::in | ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");

    // Read size of tree
    int nChunks;
    f.read((char*) &nChunks, sizeof(int));
    SerialFunctionTree<D> &sTree = *this->getSerialFunctionTree();

    // Read tree data, chunk by chunk
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        if (iChunk < sTree.nodeChunks.size()) {
            sTree.sNodes = sTree.nodeChunks[iChunk];
        } else {
            double *sNodesCoeff;
            sNodesCoeff = new double[sTree.sizeNodeCoeff*sTree.maxNodesPerChunk];
            sTree.nodeCoeffChunks.push_back(sNodesCoeff);
            sTree.sNodes = (ProjectedNode<D>*) new char[sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>)];
            sTree.nodeChunks.push_back(sTree.sNodes);
        }
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        f.read((char*) sTree.nodeChunks[iChunk], count);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        f.read((char*) sTree.nodeCoeffChunks[iChunk], count*sizeof(double));
    }
    f.close();

    t1.stop();
    Printer::printTime(10, "Time read tree", t1);

    Timer t2;
    sTree.rewritePointers(nChunks);
    t2.stop();
    Printer::printTime(10, "Time rewrite pointers", t2);
}

template<int D>
double FunctionTree<D>::integrate() const {
    double result = 0.0;
    for (int i = 0; i < this->rootBox.size(); i++) {
        const FunctionNode<D> &fNode = getRootFuncNode(i);
        result += fNode.integrate();
    }
    return result;
}

template<int D>
double FunctionTree<D>::evalf(const double *r) {
    MWNode<D> &mr_node = this->getNodeOrEndNode(r);
    FunctionNode<D> &f_node = static_cast<FunctionNode<D> &>(mr_node);
    double result = f_node.evalf(r);
    this->deleteGenerated();
    return result;
}

template<int D>
void FunctionTree<D>::square() {
    this->deleteGenerated();
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim()*this->getKp1_d();
    for (int n = 0; n < nNodes; n++) {
        MWNode<D> &node = *this->endNodeTable[n];
        node.mwTransform(Reconstruction);
        node.cvTransform(Forward);
        double *coefs = node.getCoefs();
        for (int i = 0; i < nCoefs; i++) {
            coefs[i] *= coefs[i];
        }
        node.cvTransform(Backward);
        node.mwTransform(Compression);
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template<int D>
void FunctionTree<D>::rescale(double c) {
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim()*this->getKp1_d();
    for (int i = 0; i < nNodes; i++) {
        MWNode<D> &node = *this->endNodeTable[i];
        if (not node.hasCoefs()) MSG_FATAL("No coefs");
        double *coefs = node.getCoefs();
        for (int j = 0; j < nCoefs; j++) {
            coefs[j] *= c;
        }
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template<int D>
void FunctionTree<D>::normalize() {
    double sq_norm = this->getSquareNorm();
    if (sq_norm < 0.0) MSG_ERROR("Normalizing uninitialized function");
    this->rescale(1.0/sqrt(sq_norm));
}


template<int D>
int FunctionTree<D>::getNChunks() {
    return this->getSerialFunctionTree()->getNChunks();
}

template<int D>
int FunctionTree<D>::getNChunksUsed() {
    return this->getSerialFunctionTree()->getNChunksUsed();
}

template<int D>
void FunctionTree<D>::getEndValues(VectorXd &data) {
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim()*this->getKp1_d();
    data = VectorXd::Zero(nNodes*nCoefs);
    for (int n = 0; n < nNodes; n++) {
        MWNode<D> &node = getEndFuncNode(n);
        node.mwTransform(Reconstruction);
        node.cvTransform(Forward);
        const double *c = node.getCoefs();
        for (int i = 0; i < nCoefs; i++) {
            data(n*nCoefs + i) = c[i];
        }
        node.cvTransform(Backward);
        node.mwTransform(Compression);
    }
}

template<int D>
void FunctionTree<D>::setEndValues(VectorXd &data) {
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim()*this->getKp1_d();
    for (int i = 0; i < nNodes; i++) {
        MWNode<D> &node = getEndFuncNode(i);
        const double *c = data.segment(i*nCoefs, nCoefs).data();
        node.setCoefBlock(0, nCoefs, c);
        node.cvTransform(Backward);
        node.mwTransform(Compression);
        node.setHasCoefs();
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template<int D>
std::ostream& FunctionTree<D>::print(std::ostream &o) {
    o << std::endl << "*FunctionTree: " << this->name << std::endl;
    return MWTree<D>::print(o);
}

template<int D>
double dot(FunctionTree<D> &bra, FunctionTree<D> &ket) {
    if (bra.getMRA() != ket.getMRA()){
        MSG_FATAL("Trees not compatible");
    }
    MWNodeVector nodeTable;
    HilbertIterator<D> it(&bra);
    it.setReturnGenNodes(false);
    while(it.next()) {
        MWNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    double result = 0.0;
    double locResult = 0.0;
//OMP is disabled in order to get EXACT results (to the very last digit), the
//order of summation makes the result different beyond the 14th digit or so.
//OMP does improve the performace, but its not worth it for the time being.
//#pragma omp parallel firstprivate(n_nodes, locResult)
//		shared(nodeTable,rhs,result)
//    {
//#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        const FunctionNode<D> &braNode = static_cast<const FunctionNode<D> &>(*nodeTable[n]);
        const MWNode<D> *mwNode = ket.findNode(braNode.getNodeIndex());
        if (mwNode == 0) continue;

        const FunctionNode<D> &ketNode = static_cast<const FunctionNode<D> &>(*mwNode);
        if (braNode.isRootNode()) {
            locResult += dotScaling(braNode, ketNode);
        }
        locResult += dotWavelet(braNode, ketNode);
    }
//#pragma omp critical
    result += locResult;
//    }
    return result;
}

template double dot(FunctionTree<1> &bra, FunctionTree<1> &ket);
template double dot(FunctionTree<2> &bra, FunctionTree<2> &ket);
template double dot(FunctionTree<3> &bra, FunctionTree<3> &ket);

template class FunctionTree<1>;
template class FunctionTree<2>;
template class FunctionTree<3>;

} //namespace mrcpp
