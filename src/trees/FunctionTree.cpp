/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "FunctionTree.h"

#include <fstream>

#include "FunctionNode.h"
#include "NodeAllocator.h"

#include "treebuilders/grid.h"
#include "utils/Bank.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/mpi_utils.h"
#include "utils/periodic_utils.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

/** @returns New FunctionTree object
 *
 *  @param[in] mra: Which MRA the function is defined
 *  @param[in] sh_mem: Pointer to MPI shared memory block
 *
 *  @details Constructs an uninitialized tree, containing only empty root nodes.
 *  If a shared memory pointer is provided the tree will be allocated in this
 *  shared memory window, otherwise it will be local to each MPI process.
 */
template <int D, typename T>
FunctionTree<D, T>::FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory<T> *sh_mem, const std::string &name)
        : MWTree<D, T>(mra, name)
        , RepresentableFunction<D, T>(mra.getWorldBox().getLowerBounds().data(), mra.getWorldBox().getUpperBounds().data()) {
    int nodesPerChunk = 2048; // Large chunks are required for not leading to memory fragmentation (32 MB on "Betzy" 2023)
    // nodesPerChunk is same for real and complex trees: the size (in MB) of the complex chunks are twice as large
    int coefsGenNodes = this->getKp1_d();
    int coefsRegNodes = this->getTDim() * this->getKp1_d();
    this->nodeAllocator_p = std::make_unique<NodeAllocator<D, T>>(this, sh_mem, coefsRegNodes, nodesPerChunk);
    this->genNodeAllocator_p = std::make_unique<NodeAllocator<D, T>>(this, nullptr, coefsGenNodes, nodesPerChunk);
    this->allocRootNodes();
    this->resetEndNodeTable();
}

template <int D, typename T> void FunctionTree<D, T>::allocRootNodes() {
    auto &allocator = this->getNodeAllocator();
    auto &rootbox = this->getRootBox();

    int nRoots = rootbox.size();
    int sIdx = allocator.alloc(nRoots);

    auto n_coefs = allocator.getNCoefs();
    auto *coef_p = allocator.getCoef_p(sIdx);
    auto *root_p = allocator.getNode_p(sIdx);

    MWNode<D, T> **roots = rootbox.getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        // construct into allocator memory
        new (root_p) FunctionNode<D, T>(this, rIdx);
        roots[rIdx] = root_p;

        root_p->serialIx = sIdx;
        root_p->parentSerialIx = -1;
        root_p->childSerialIx = -1;

        root_p->n_coefs = n_coefs;
        root_p->coefs = coef_p;
        root_p->setIsAllocated();

        root_p->setIsRootNode();
        root_p->setIsLeafNode();
        root_p->setIsEndNode();
        root_p->clearHasCoefs();

        this->incrementNodeCount(root_p->getScale());
        sIdx++;
        root_p++;
        coef_p += n_coefs;
    }
}

// FunctionTree destructor
template <int D, typename T> FunctionTree<D, T>::~FunctionTree() {
    if (this->getNNodes() > 0) this->deleteRootNodes();
}

/** @brief Read a previously stored tree assuming text/ASCII format,
 *   in a representation using MADNESS conventions for n, l and index order.
 * @param[in] file: File name
 * @note This tree must have the exact same MRA the one that was saved(?)
 */
template <int D, typename T> void FunctionTree<D, T>::loadTreeTXT(const std::string &file) {
    std::ifstream in(file);
    int NDIM, k;
    in >> NDIM;
    if (NDIM != D) NOT_IMPLEMENTED_ABORT;
    double coord[D][2];
    for (int d = 0; d < D; d++) in >> coord[d][0] >> coord[d][1];

    int p = 1;
    int rscale = this->getRootScale(); // root scale of target MRA (MRChem) . NB: negative
    for (int i = rscale; i < 0; i++) p *= 2;
    int L = p; // NB for now we assume the world as a cube going from -L to +L and L is a power of 2
    // We require that the world box size is identical and a power of 2
    double TXT_thres = 1.0e-14; // threshold for differences in scaling factors
    for (int d = 0; d < D; d++) {
        if (std::abs(coord[d][0] + L) > TXT_thres) std::cout << coord[d][0] << " " << L << std::endl;
        if (std::abs(coord[d][0] + L) > TXT_thres) NOT_IMPLEMENTED_ABORT;
        if (std::abs(coord[d][1] - L) > TXT_thres) std::cout << coord[d][1] << " " << L << std::endl;
        if (std::abs(coord[d][1] - L) > TXT_thres) NOT_IMPLEMENTED_ABORT;
    }

    int nChildren = 1;
    for (int d = 0; d < D; d++) nChildren *= 2;

    int nmax = 0; // deppeset scale in TXT
    in >> k;
    if (k != this->getKp1()) NOT_IMPLEMENTED_ABORT;
    k--; // MRChem defines k as highest polynomial order. MADNESS as number of polynomials

    int ncoefs = 1; // number of coefficents for one single node (not a full MRChem MWnode which stores 2**D of them)
    for (int i = 0; i < D; i++) ncoefs *= k + 1;

    std::vector<std::vector<MWNode<D, T> *>> NodeTable(50); // to store all the nodes pointers
    std::map<int, int> mp;                                  // to store the number of children stored in each parent node
    // MRChem and MADNESS do not use the same indices order for the qudrature points
    // We read MADNESS convention (note that mapMRC[mapMRC[i]]=i for all i)
    std::vector<int> mapMRC; // mapping vector
    int kx = k;
    int ky = k;
    int kz = k;
    if (D < 3) kz = 0;
    if (D < 2) ky = 0;
    int kp1 = k + 1;
    // MADNESS: zyx and i=k,k-1,k-2... MRChem: xyz, i=0,1,2,3 ...
    for (int x = kx; x >= 0; x--) {
        for (int y = ky; y >= 0; y--) {
            for (int z = kz; z >= 0; z--) { mapMRC.push_back(z * kp1 * kp1 + y * kp1 + x); }
        }
    }

    MWNode<D, T> **roots = this->getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nChildren; rIdx++) {
        roots[rIdx]->deleteChildren();
        roots[rIdx]->zeroCoefs();
    }
    this->clearEndNodeTable();

    int nread; // number of nodes to read
    in >> nread;
    while (nread-- > 0) {
        // NB: MRChem stores quadrature points values in the PARENT node. 2**D nodes are stored in the same parent
        int n;    // TXT scale
        int n_in; // MRChem scale
        in >> n_in;
        n = n_in + rscale - 1; // MRChem does not define root scale as zero.

        std::array<int, D> l_in; // translation index TXT
        std::array<int, D> l;    // translation index MRChem
        std::array<int, D> lp;   // translation index MRChem, parent

        for (int i = 0; i < D; i++) in >> l_in[i];

        // MRChem defines smallest l as -(2**n)*L , where -L is smallest world coordinate.
        // note that root scale has 2**D nodes (if range is -L,L)
        for (int i = 0; i < D; i++) {
            l[i] = l_in[i] - std::pow(2, n) * L;
            lp[i] = l_in[i] / 2 - std::pow(2, n - 1) * L; // for parent
        }
        NodeIndex<D> idx_p(n - 1, lp); // index of parent node
        MWNode<D, T> *node = &this->getNode(idx_p, true);
        // note that node is not necesssarily an endnode, but they children are always endnodes
        // must find to which child of the parent node it corresponds
        int c_ix = 0; // child index in the parent
        int p = 1;
        for (int i = 0; i < D; i++) {
            if (abs(l[i]) % 2 == 1) c_ix += p;
            p *= 2;
        }
        T *values = node->getCoefs();
        if (mp[node->getSerialIx()] == 0) {
            // init to zero
            node->zeroCoefs();
            if (not node->isRootNode()) {
                // also set siblings to zero if not set yet
                MWNode<D, T> *parent = &node->getMWParent();
                for (int cIdx = 0; cIdx < nChildren; cIdx++) {
                    if (mp[parent->getMWChild(cIdx).getSerialIx()] == 0) parent->getMWChild(cIdx).zeroCoefs();
                }
            }
        }
        values += c_ix * ncoefs;                                  // repoint to the right child position (ncoefs is for one child only)
        for (int i = 0; i < ncoefs; i++) in >> values[mapMRC[i]]; // the indice i is mapped
        mp[node->getSerialIx()]++;                                // counts the number of children included
        nmax = std::max(nmax, n_in);                              // deepest scale in TXT
        if (mp[node->getSerialIx()] == 1) NodeTable[n_in].push_back(node);
    }
    in.close();
    // transform all nodes from quadrature point values to scaling coefficients
    for (int n = nmax; n > -1; n--) {
        for (int i = 0; i < NodeTable[n].size(); i++) {
            MWNode<D, T> *node = NodeTable[n][i];
            node->cvTransform(Backward);
            node->calcNorms();
        }
    }
    // now tree has only scaling coefficients or zeros on end nodes

    // Transform into scaling and wavelets, starting by leaf nodes and copying scaling into parents
    for (int n = nmax; n > -1; n--) {
        for (int i = 0; i < NodeTable[n].size(); i++) {
            MWNode<D, T> *node = NodeTable[n][i];
            if (mp[node->getSerialIx()] == nChildren) {
                // node complete: transform into scaling and wavelets
                if (node->isEndNode()) {
                    node->mwTransform(Compression);
                    node->setHasCoefs();
                    node->calcNorms();
                    this->endNodeTable.push_back(node);
                } else {
                    // MRCPP requires that all nodes that have no children are end nodes
                    // and all nodes are groups of 2**D siblings
                    T *pcoefs = node->getCoefs(); // parent coefficients
                    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
                        MWNode<D, T> *cnode = &node->getMWChild(cIdx);
                        if (mp[cnode->getSerialIx()] != nChildren) {
                            // This child is not defined. must take scaling from parent
                            if (mp[cnode->getSerialIx()] > 0) std::cout << "accounting error " << std::endl;
                            T *ccoefs = cnode->getCoefs(); // child coefficients
                            for (int j = 0; j < ncoefs; j++) ccoefs[j] = pcoefs[j + cIdx * ncoefs];
                            for (int j = ncoefs; j < ncoefs * nChildren; j++) ccoefs[j] = 0.0; // the remainder are set to zero
                            this->endNodeTable.push_back(cnode);                               // add to the list of nodes
                            cnode->setHasCoefs();
                            cnode->calcNorms();
                        }
                    }
                    node->mwTransform(Compression);
                    node->setHasCoefs();
                    node->calcNorms();
                }
                if (not node->isRootNode()) {
                    // and copy the new scaling parts into parent
                    MWNode<D, T> *parent = &node->getMWParent();
                    // check if parent exist already, and put in the list if not.
                    if (mp[parent->getSerialIx()] == 0) NodeTable[n - 1].push_back(parent);
                    int my_ix = -1;
                    // find index among siblings
                    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
                        if (&parent->getMWChild(cIdx) == node) my_ix = cIdx;
                    }
                    if (my_ix < 0) std::cout << " DID NOT FIND INDEX" << std::endl;
                    T *ccoefs = node->getCoefs();
                    T *pcoefs = parent->getCoefs();
                    for (int j = 0; j < ncoefs; j++) pcoefs[j + my_ix * ncoefs] = ccoefs[j];
                    mp[parent->getSerialIx()]++;
                }
            } else {
                std::cout << " WARNING: found incomplete node " << std::endl;
            }
        }
    }
    this->calcSquareNorm();
}

/** @brief Write the tree to disk in text/ASCII format in a representation
 *   using MADNESS conventions for n, l and index order.
 * @param[in] file: File name
 */
template <int D, typename T> void FunctionTree<D, T>::saveTreeTXT(const std::string &fname) {
    std::ofstream out(fname);
    out << std::setprecision(14);
    out << D << std::endl;
    int rscale = this->getRootScale();
    std::array<double, D> sf = this->getMRA().getWorldBox().getScalingFactors();
    double LMRChem = 1.0;
    for (int i = 0; i > rscale; i--) LMRChem *= 2; // we assume world is from -L to L, and a cube with 2 root nodes in each direction
    for (int d = 0; d < D; d++) { out << -sf[d] * LMRChem << " " << sf[d] * LMRChem << std::endl; }
    int kp1 = this->getKp1();
    out << kp1 << std::endl;
    int ncoefs = 1;
    for (int d = 0; d < D; d++) ncoefs *= kp1;
    int Tdim = std::pow(2, D);
    std::vector<T> values(ncoefs * Tdim);

    int nout = this->endNodeTable.size();
    out << Tdim * nout << std::endl; // could output only scaling coeff?

    // MRChem and MADNESS do not use the same indices order for the qudrature points
    // We write into MADNESS convention (note that mapMRC[mapMRC[i]]=i for all i)
    std::vector<int> mapMRC; // mapping vector
    int kx = kp1 - 1;
    int ky = kp1 - 1;
    int kz = kp1 - 1;
    if (D < 3) kz = 0;
    if (D < 2) ky = 0;
    // MADNESS: zyx and i=k,k-1,k-2... MRChem: xyz, i=0,1,2,3 ...
    for (int x = kx; x >= 0; x--) {
        for (int y = ky; y >= 0; y--) {
            for (int z = kz; z >= 0; z--) { mapMRC.push_back(z * kp1 * kp1 + y * kp1 + x); }
        }
    }

    int L = std::pow(2, -rscale);
    int count = -1;
    while (++count < nout) {
        std::array<int, D> l;
        NodeIndex<D> idx = this->endNodeTable[count]->getNodeIndex();
        MWNode<D, T> *node = &(this->getNode(idx, false));
	T *coefs = node->getCoefs();
	for (int i = 0; i < ncoefs * Tdim; i++) values[i] = coefs[i];
	    node->attachCoefs(values.data());
        int n = idx.getScale();
        node->mwTransform(Reconstruction);
        node->cvTransform(Forward);
        // we write for each children nodes separately
        for (int i = 0; i < D; i++) {
            // l in interval [0, max], while in MRCPP it is defined in [-max/2, max/2-1]
            l[i] = 2 * (idx.getTranslation(i) + std::pow(2, n) * L); // first child
        }
        for (int cix = 0; cix < Tdim; cix++) {
            out << n - rscale + 2 << " "; // scales start at zero. NB: children are one scale larger than node
            for (int i = 0; i < D; i++) {
                int p = (cix >> i) & 1; // shift by one for odd child indices
                out << l[i] + p << " ";
            }
            out << std::endl;
            for (int i = 0; i < ncoefs; i++) out << values[cix * ncoefs + mapMRC[i]] << " ";
            out << std::endl;
        }
	node->attachCoefs(coefs); // put back original coeff
   }
    out.close();
}
/** @brief Write the tree structure to disk, for later use
 * @param[in] file: File name, will get ".tree" extension
 */
template <int D, typename T> void FunctionTree<D, T>::saveTree(const std::string &file) {
    Timer t1;

    this->deleteGenerated();
    auto &allocator = this->getNodeAllocator();
    std::stringstream fname;
    fname << file << ".tree";
    std::fstream f;
    f.open(fname.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");

    // Write size of tree
    int nChunks = allocator.getNChunksUsed();
    f.write((char *)&nChunks, sizeof(int));
    // Write tree data, chunk by chunk
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
       f.write((char *)allocator.getNodeChunk(iChunk), allocator.getNodeChunkSize());
       f.write((char *)allocator.getCoefChunk(iChunk), allocator.getCoefChunkSize());
    }
    f.close();
    print::time(10, "Time write", t1);
}

/** @brief Read a previously stored tree structure from disk
 * @param[in] file: File name, will get ".tree" extension
 * @note This tree must have the exact same MRA the one that was saved
 */
template <int D, typename T> void FunctionTree<D, T>::loadTree(const std::string &file) {
    Timer t1;

    std::stringstream fname;
    fname << file << ".tree";

    std::fstream f;
    f.open(fname.str(), std::ios::in | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");

    // Read size of tree
    int nChunks;
    f.read((char *)&nChunks, sizeof(int));

    // Read tree data, chunk by chunk
    this->deleteRootNodes();
    auto &allocator = this->getNodeAllocator();
    allocator.init(nChunks);
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        f.read((char *)allocator.getNodeChunk(iChunk), allocator.getNodeChunkSize());
        f.read((char *)allocator.getCoefChunk(iChunk), allocator.getCoefChunkSize());
    }
    f.close();
    print::time(10, "Time read tree", t1);

    Timer t2;
    allocator.reassemble();
    this->resetEndNodeTable();
    this->calcSquareNorm(true);
    print::time(10, "Time rewrite pointers", t2);
}

/** @returns Integral of the function over the entire computational domain */
template <int D, typename T> T FunctionTree<D, T>::integrate() const {

    T result = 0.0;
    for (int i = 0; i < this->rootBox.size(); i++) {
        const FunctionNode<D, T> &fNode = getRootFuncNode(i);
        result += fNode.integrate();
    }

    // Handle potential scaling
    auto scaling_factor = this->getMRA().getWorldBox().getScalingFactors();
    auto jacobian = 1.0;
    for (const auto &sf_i : scaling_factor) jacobian *= std::sqrt(sf_i);
    // Square root of scaling factor in each diection. The seemingly missing
    // multiplication by the square root of sf_i is included in the basis

    return jacobian * result;
}

/** @returns Integral of a representable function over the grid given by the tree */
template <> double FunctionTree<3, double>::integrateEndNodes(RepresentableFunction_M &f) {
    // traverse tree, and treat end nodes only
    std::vector<FunctionNode<3> *> stack; // node from this
    for (int i = 0; i < this->getRootBox().size(); i++) stack.push_back(&(this->getRootFuncNode(i)));
    double result = 0.0;
    while (stack.size() > 0) {
        FunctionNode<3> *Node = stack.back();
        stack.pop_back();
        if (Node->getNChildren() > 0) {
            for (int i = 0; i < Node->getNChildren(); i++) stack.push_back(&(Node->getFuncChild(i)));
        } else {
            // end nodes
            Eigen::MatrixXd fmat = f.evalf(Node->nodeIndex);
            double *coefs = Node->getCoefs(); // save position of coeff, but do not use them!
            // The data in fmat is not organized so that two consecutive points are stored after each other in memory, so needs to copy before mwtransform, cannot use memory adress directly.
            int nc = fmat.cols();
            std::vector<double> cc(nc);
            for (int i = 0; i < nc; i++) cc[i] = fmat(0, i);
            Node->attachCoefs(cc.data());
            result += Node->integrateValues();
            Node->attachCoefs(coefs); // put back original coeff
        }
    }

    // Handle potential scaling
    auto scaling_factor = this->getMRA().getWorldBox().getScalingFactors();
    auto jacobian = 1.0;
    for (const auto &sf_i : scaling_factor) jacobian *= std::sqrt(sf_i);
    // Square root of scaling factor in each diection. The seemingly missing
    // multiplication by the square root of sf_i is included in the basis
    return jacobian * result;
}

/** @returns Function value in a point, out of bounds returns zero
 *
 * @param[in] r: Cartesian coordinate
 *
 * @note This will only evaluate the _scaling_ part of the
 *       leaf nodes in the tree, which means that the function
 *       values will not be fully accurate.
 *       This is done to allow a fast and const function evaluation
 *       that can be done in OMP parallel. If you want to include
 *       also the _final_ wavelet part you can call the corresponding
 *       evalf_precise function, _or_ you can manually extend
 *       the MW grid by one level before evaluating, using
 *       `mrcpp::refine_grid(tree, 1)`
 */
template <int D, typename T> T FunctionTree<D, T>::evalf(const Coord<D> &r) const {
    // Handle potential scaling
    const auto scaling_factor = this->getMRA().getWorldBox().getScalingFactors();
    auto arg = r;
    for (auto i = 0; i < D; i++) arg[i] = arg[i] / scaling_factor[i];

    // The 1.0 appearing in the if tests comes from the period is
    // always 1.0 from the point of view of this function.
    if (this->getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(arg, this->getRootBox().getPeriodic()); }

    // Function is zero outside the domain for non-periodic functions
    if (this->outOfBounds(arg) and not this->getRootBox().isPeriodic()) return 0.0;

    const MWNode<D, T> &mw_node = this->getNodeOrEndNode(arg);
    auto &f_node = static_cast<const FunctionNode<D, T> &>(mw_node);
    auto result = f_node.evalScaling(arg);

    // Adjust for scaling factor included in basis
    auto coef = 1.0;
    for (const auto &fac : scaling_factor) coef /= std::sqrt(fac);

    return coef * result;
}

/** @returns Function value in a point, out of bounds returns zero
 *
 * @param[in] r: Cartesian coordinate
 *
 * @note This will evaluate the _true_ value (scaling + wavelet) of the
 *       leaf nodes in the tree. This requires an on-the-fly MW transform
 *       on the node which makes this function slow and non-const. If you
 *       need fast evaluation, use refine_grid(tree, 1) first, and then
 *       evalf.
 */
template <int D, typename T> T FunctionTree<D, T>::evalf_precise(const Coord<D> &r) {
    // Handle potential scaling
    const auto scaling_factor = this->getMRA().getWorldBox().getScalingFactors();
    auto arg = r;
    for (auto i = 0; i < D; i++) arg[i] = arg[i] / scaling_factor[i];

    // The 1.0 appearing in the if tests comes from the period is
    // always 1.0 from the point of view of this function.
    if (this->getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(arg, this->getRootBox().getPeriodic()); }

    // Function is zero outside the domain for non-periodic functions
    if (this->outOfBounds(arg) and not this->getRootBox().isPeriodic()) return 0.0;

    MWNode<D, T> &mw_node = this->getNodeOrEndNode(arg);
    auto &f_node = static_cast<FunctionNode<D, T> &>(mw_node);
    auto result = f_node.evalf(arg);
    this->deleteGenerated();

    // Adjust for scaling factor included in basis
    auto coef = 1.0;
    for (const auto &fac : scaling_factor) coef /= std::sqrt(fac);

    return coef * result;
}

/** @brief In-place square of MW function representations, fixed grid
 *
 * @details The leaf node point values of the function will be in-place
 * squared, no grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::square() {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");

#pragma omp parallel num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &node = *this->endNodeTable[n];
            node.mwTransform(Reconstruction);
            node.cvTransform(Forward);
            T *coefs = node.getCoefs();
            for (int i = 0; i < nCoefs; i++) { coefs[i] *= coefs[i]; }
            node.cvTransform(Backward);
            node.mwTransform(Compression);
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

/** @brief In-place power of MW function representations, fixed grid
 *
 * @param[in] p: Numerical power
 *
 * @details The leaf node point values of the function will be in-place raised
 * to the given power, no grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::power(double p) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");

#pragma omp parallel num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &node = *this->endNodeTable[n];
            node.mwTransform(Reconstruction);
            node.cvTransform(Forward);
            T *coefs = node.getCoefs();
            for (int i = 0; i < nCoefs; i++) { coefs[i] = std::pow(coefs[i], p); }
            node.cvTransform(Backward);
            node.mwTransform(Compression);
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

/** @brief In-place multiplication by a scalar, fixed grid
 *
 * @param[in] c: Scalar coefficient
 *
 * @details The leaf node point values of the function will be
 * in-place multiplied by the given coefficient, no grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::rescale(T c) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int i = 0; i < nNodes; i++) {
            MWNode<D, T> &node = *this->endNodeTable[i];
            if (not node.hasCoefs()) MSG_ABORT("No coefs");
            T *coefs = node.getCoefs();
            for (int j = 0; j < nCoefs; j++) { coefs[j] *= c; }
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

/** @brief In-place rescaling by a function norm \f$ ||f||^{-1} \f$, fixed grid */
template <int D, typename T> void FunctionTree<D, T>::normalize() {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    double sq_norm = this->getSquareNorm();
    if (sq_norm < 0.0) MSG_ERROR("Normalizing uninitialized function");
    this->rescale(1.0 / std::sqrt(sq_norm));
}

/** @brief In-place addition with MW function representations, fixed grid
 *
 * @param[in] c: Numerical coefficient of input function
 * @param[in] inp: Input function to add
 *
 * @details The input function will be added in-place on the current grid of
 * the function, i.e. no further grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::add(T c, FunctionTree<D, T> &inp) {
    if (this->getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &out_node = *this->endNodeTable[n];
            MWNode<D, T> &inp_node = inp.getNode(out_node.getNodeIndex());
            T *out_coefs = out_node.getCoefs();
            const T *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] += c * inp_coefs[i]; }
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}
/** @brief In-place addition with MW function representations, fixed grid
 *
 * @param[in] c: Numerical coefficient of input function
 * @param[in] inp: Input function to add
 *
 * @details The input function will be added to the union of the current grid of
 * and input the function grid.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::add_inplace(T c, FunctionTree<D, T> &inp) {
    if (this->getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    while (refine_grid(*this, inp)) {};
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &out_node = *this->endNodeTable[n];
            MWNode<D, T> &inp_node = inp.getNode(out_node.getNodeIndex());
            T *out_coefs = out_node.getCoefs();
            const T *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] += c * inp_coefs[i]; }
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}

/** @brief In-place addition of absolute values of MW function representations
 *
 * @param[in] c Numerical coefficient of input function
 * @param[in] inp Input function to add
 *
 * The absolute value of input function will be added in-place on the current grid of the output
 * function, i.e. no further grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::absadd(T c, FunctionTree<D, T> &inp) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &out_node = *this->endNodeTable[n];
            MWNode<D, T> inp_node = inp.getNode(out_node.getNodeIndex()); // Full copy
            out_node.mwTransform(Reconstruction);
            out_node.cvTransform(Forward);
            inp_node.mwTransform(Reconstruction);
            inp_node.cvTransform(Forward);
            T *out_coefs = out_node.getCoefs();
            const T *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] = std::norm(out_coefs[i]) + std::norm(c * inp_coefs[i]); }
            out_node.cvTransform(Backward);
            out_node.mwTransform(Compression);
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}

/** @brief In-place multiplication with MW function representations, fixed grid
 *
 * @param[in] c: Numerical coefficient of input function
 * @param[in] inp: Input function to multiply
 *
 * @details The input function will be multiplied in-place on the current grid
 * of the function, i.e. no further grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::multiply(T c, FunctionTree<D, T> &inp) {
    if (this->getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &out_node = *this->endNodeTable[n];
            MWNode<D, T> inp_node = inp.getNode(out_node.getNodeIndex()); // Full copy
            out_node.mwTransform(Reconstruction);
            out_node.cvTransform(Forward);
            inp_node.mwTransform(Reconstruction);
            inp_node.cvTransform(Forward);
            T *out_coefs = out_node.getCoefs();
            const T *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] *= c * inp_coefs[i]; }
            out_node.cvTransform(Backward);
            out_node.mwTransform(Compression);
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}

/** @brief In-place mapping with a predefined function f(x), fixed grid
 *
 * @param[in] fmap: mapping function
 *
 * @details The input function will be mapped in-place on the current grid
 * of the function, i.e. no further grid refinement.
 *
 */
template <int D, typename T> void FunctionTree<D, T>::map(FMap<T, T> fmap) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    {
        int nNodes = this->getNEndNodes();
#pragma omp parallel for schedule(guided) num_threads(mrcpp_get_num_threads())
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &node = *this->endNodeTable[n];
            node.mwTransform(Reconstruction);
            node.cvTransform(Forward);
            T *coefs = node.getCoefs();
            for (int i = 0; i < node.getNCoefs(); i++) { coefs[i] = fmap(coefs[i]); }
            node.cvTransform(Backward);
            node.mwTransform(Compression);
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template <int D, typename T> void FunctionTree<D, T>::getEndValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &data) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim() * this->getKp1_d();
    data = VectorXd::Zero(nNodes * nCoefs);
    for (int n = 0; n < nNodes; n++) {
        MWNode<D, T> &node = getEndFuncNode(n);
        node.mwTransform(Reconstruction);
        node.cvTransform(Forward);
        const T *c = node.getCoefs();
        for (int i = 0; i < nCoefs; i++) { data(n * nCoefs + i) = c[i]; }
        node.cvTransform(Backward);
        node.mwTransform(Compression);
    }
}

template <int D, typename T> void FunctionTree<D, T>::setEndValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &data) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim() * this->getKp1_d();
    for (int i = 0; i < nNodes; i++) {
        MWNode<D, T> &node = getEndFuncNode(i);
        const T *c = data.segment(i * nCoefs, nCoefs).data();
        node.setCoefBlock(0, nCoefs, c);
        node.cvTransform(Backward);
        node.mwTransform(Compression);
        node.setHasCoefs();
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template <int D, typename T> std::ostream &FunctionTree<D, T>::print(std::ostream &o) const {
    o << std::endl << "*FunctionTree: " << this->name << std::endl;
    o << "  genNodes: " << getNGenNodes() << std::endl;
    return MWTree<D, T>::print(o);
}

/** @brief Reduce the precision of the tree by deleting nodes
 *
 * @param prec: New precision criterion
 * @param splitFac: Splitting factor: 1, 2 or 3
 * @param absPrec: Use absolute precision
 *
 * @details This will run the tree building algorithm in "reverse", starting
 * from the leaf nodes, and perform split checks on each node based on the given
 * precision and the local wavelet norm.
 *
 * @note The splitting factor appears in the threshold for the wavelet norm as
 * \f$ ||w|| < 2^{-sn/2} ||f|| \epsilon \f$. In principal, `s` should be equal
 * to the dimension; in practice, it is set to `s=1`.
 */
template <int D, typename T> int FunctionTree<D, T>::crop(double prec, double splitFac, bool absPrec) {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D, T> &root = this->getRootMWNode(i);
        root.crop(prec, splitFac, absPrec);
    }
    int nChunks = this->getNodeAllocator().compress();
    this->resetEndNodeTable();
    this->calcSquareNorm();
    return nChunks;
}

/** Traverse tree using BFS and returns an array with the address of the coefs.
 * Also returns an array with the corresponding indices defined as the
 * values of serialIx in refTree, and an array with the indices of the parent.
 * Set index -1 for nodes that are not present in refTree */
template <int D, typename T>
void FunctionTree<D, T>::makeCoeffVector(std::vector<T *> &coefs,
                                         std::vector<int> &indices,
                                         std::vector<int> &parent_indices,
                                         std::vector<double> &scalefac,
                                         int &max_index,
                                         MWTree<D, double> &refTree,
                                         std::vector<MWNode<D, double> *> *refNodes) {
    coefs.clear();
    indices.clear();
    parent_indices.clear();
    max_index = 0;
    std::vector<MWNode<D, double> *> refstack; // nodes from refTree
    std::vector<MWNode<D, T> *> thisstack;     // nodes from this Tree
    for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
        refstack.push_back(refTree.getRootBox().getNodes()[rIdx]);
        thisstack.push_back(this->getRootBox().getNodes()[rIdx]);
    }
    int stack_p = 0;
    while (thisstack.size() > stack_p) {
        // refNode and thisNode are the same node in space, but on different trees
        MWNode<D, T> *thisNode = thisstack[stack_p];
        MWNode<D, double> *refNode = refstack[stack_p++];
        coefs.push_back(thisNode->getCoefs());
        if (refNodes != nullptr) refNodes->push_back(refNode);
        if (refNode != nullptr) {
            indices.push_back(refNode->getSerialIx());
            max_index = std::max(max_index, refNode->getSerialIx());
            parent_indices.push_back(refNode->parentSerialIx); // is -1 for root nodes
            double sfac = std::pow(2.0, -0.5 * (refNode->getScale() + 1.0));
            scalefac.push_back(sfac); // could be faster: essentially inverse of powers of 2
        } else {
            indices.push_back(-1); // indicates that the node is not in the refTree
            parent_indices.push_back(-2);
            scalefac.push_back(1.0);
        }
        if (thisNode->getNChildren() > 0) {
            for (int i = 0; i < thisNode->getNChildren(); i++) {
                if (refNode != nullptr and refNode->getNChildren() > 0)
                    refstack.push_back(refNode->children[i]);
                else
                    refstack.push_back(nullptr);
                thisstack.push_back(thisNode->children[i]);
            }
        }
    }
}

/** Traverse tree using DFS and reconstruct it using node info from the
 * reference tree and a list of coefficients.
 * It is the reference tree (refTree) which is traversed, but one does not descend
 * into children if the norm of the tree is smaller than absPrec. */
template <int D, typename T> void FunctionTree<D, T>::makeTreefromCoeff(MWTree<D, double> &refTree, std::vector<T *> coefpVec, std::map<int, int> &ix2coef, double absPrec, const std::string &mode) {
    std::vector<MWNode<D, double> *> stack;
    std::map<int, MWNode<D, T> *> ix2node; // gives the nodes in this tree for a given ix
    int sizecoef = (1 << this->getDim()) * this->getKp1_d();
    int sizecoefW = ((1 << this->getDim()) - 1) * this->getKp1_d();
    this->squareNorm = 0.0;
    this->clearEndNodeTable();
    for (int rIdx = 0; rIdx < refTree.getRootBox().size(); rIdx++) {
        MWNode<D, double> *refNode = refTree.getRootBox().getNodes()[rIdx];
        stack.push_back(refNode);
        int ix = ix2coef[refNode->getSerialIx()];
        ix2node[ix] = this->getRootBox().getNodes()[rIdx];
    }

    while (stack.size() > 0) {
        MWNode<D, double> *refNode = stack.back(); // node in the reference tree refTree
        stack.pop_back();
        assert(ix2coef.count(refNode->getSerialIx()) > 0);
        int ix = ix2coef[refNode->getSerialIx()];
        MWNode<D, T> *node = ix2node[ix]; // corresponding node in this tree
        // copy coefficients into this tree
        int size = sizecoefW;
        if (refNode->isRootNode() or mode == "copy") {
            size = sizecoef; // all coeff
            if ((mode != "copy" and refNode->isRootNode()) or refNode->getNChildren() == 0 or ix2coef.count(refNode->children[0]->getSerialIx()) == 0) {
                for (int k = 0; k < size; k++) node->getCoefs()[k] = coefpVec[ix][k];

            } else {
                for (int k = 0; k < size; k++) node->getCoefs()[k] = 0.0; // coefpVec[ix][k];
            }
        } else {
            // only wavelets are defined in coefVec. Scaling part set below, when creating children
            if (ix < coefpVec.size() and ix >= 0) {
                for (int k = 0; k < size; k++) node->getCoefs()[k + this->getKp1_d()] = coefpVec[ix][k];
            } else {
                // we do not have W coefficients. Set them to zero
                for (int k = 0; k < size; k++) node->getCoefs()[k + this->getKp1_d()] = 0.0;
            }
        }

        node->setHasCoefs();
        node->calcNorms();
        if (mode == "copy") {
            // add children if and only if they exist in the tree
            // refNode->getNChildren() == 0 -> no children in ref tree
            // ix2coef.count(refNode->children[0]->getSerialIx()) == 0 -> the first child is not found in ix2coef
            if (refNode->getNChildren() == 0 or ix2coef.count(refNode->children[0]->getSerialIx()) == 0) {
                // no children-> it is an EndNode
                this->endNodeTable.push_back(node);
                this->squareNorm += node->getSquareNorm();
            } else {
                node->createChildren(true);
                for (int i = 0; i < refNode->getNChildren(); i++) {
                    int ixc = ix2coef[refNode->children[i]->getSerialIx()];
                    ix2node[ixc] = node->children[i];      // corresponding child node in this tree
                    stack.push_back(refNode->children[i]); // means we continue to traverse the reference tree
                }
            }
        } else if ((absPrec < 0 or tree_utils::split_check(*node, absPrec, 1.0, true)) and refNode->getNChildren() > 0) {
            // include children in tree
            node->createChildren(true);
            T *inp = node->getCoefs();
            T *out = node->getMWChild(0).getCoefs();
            tree_utils::mw_transform(*this, inp, out, false, sizecoef, true); // make the scaling part
            for (int i = 0; i < refNode->getNChildren(); i++) {
                stack.push_back(refNode->children[i]); // means we continue to traverse the reference tree
                int ixc = ix2coef[refNode->children[i]->getSerialIx()];
                ix2node[ixc] = node->children[i]; // corresponding child node in this tree
            }
        } else {
            this->endNodeTable.push_back(node);
            this->squareNorm += node->getSquareNorm();
        }
    }
}

/** Traverse tree using DFS and append same nodes as another tree, without coefficients
 *  Note that we do not use coefficients, so it does not matter what is real or complex
 */
template <int D, typename T> void FunctionTree<D, T>::appendTreeNoCoeff(MWTree<D, double> &inTree) {
    std::vector<MWNode<D, double> *> instack; // node from inTree
    std::vector<MWNode<D, T> *> thisstack;    // node from this Tree
    this->clearEndNodeTable();
    for (int rIdx = 0; rIdx < inTree.getRootBox().size(); rIdx++) {
        instack.push_back(inTree.getRootBox().getNodes()[rIdx]);
        thisstack.push_back(this->getRootBox().getNodes()[rIdx]);
    }
    while (thisstack.size() > 0) {
        // inNode and thisNode are the same node in space, but on different trees
        MWNode<D, T> *thisNode = thisstack.back();
        thisstack.pop_back();
        MWNode<D, double> *inNode = instack.back();
        instack.pop_back();
        if (inNode->getNChildren() > 0) {
            thisNode->clearIsEndNode();
            if (thisNode->getNChildren() < inNode->getNChildren()) thisNode->createChildren(false);
            for (int i = 0; i < inNode->getNChildren(); i++) {
                instack.push_back(inNode->children[i]);
                thisstack.push_back(thisNode->children[i]);
            }
        } else {
            // construct EndNodeTable for "This", starting from this branch
            // This could be done more efficiently, if it proves to be time consuming
            std::vector<MWNode<D, T> *> branchstack; // local stack starting from this branch
            branchstack.push_back(thisNode);
            while (branchstack.size() > 0) {
                MWNode<D, T> *branchNode = branchstack.back();
                branchstack.pop_back();
                if (branchNode->getNChildren() > 0) {
                    for (int i = 0; i < branchNode->getNChildren(); i++) { branchstack.push_back(branchNode->children[i]); }
                } else
                    this->endNodeTable.push_back(branchNode);
            }
        }
    }
}

/** Traverse tree using DFS and append same nodes as another tree, without coefficients */
template <int D, typename T> void FunctionTree<D, T>::appendTreeNoCoeff(MWTree<D, ComplexDouble> &inTree) {
    std::vector<MWNode<D, ComplexDouble> *> instack; // node from inTree
    std::vector<MWNode<D, T> *> thisstack;           // node from this Tree
    this->clearEndNodeTable();
    for (int rIdx = 0; rIdx < inTree.getRootBox().size(); rIdx++) {
        instack.push_back(inTree.getRootBox().getNodes()[rIdx]);
        thisstack.push_back(this->getRootBox().getNodes()[rIdx]);
    }
    while (thisstack.size() > 0) {
        // inNode and thisNode are the same node in space, but on different trees
        MWNode<D, T> *thisNode = thisstack.back();
        thisstack.pop_back();
        MWNode<D, ComplexDouble> *inNode = instack.back();
        instack.pop_back();
        if (inNode->getNChildren() > 0) {
            thisNode->clearIsEndNode();
            if (thisNode->getNChildren() < inNode->getNChildren()) thisNode->createChildren(false);
            for (int i = 0; i < inNode->getNChildren(); i++) {
                instack.push_back(inNode->children[i]);
                thisstack.push_back(thisNode->children[i]);
            }
        } else {
            // construct EndNodeTable for "This", starting from this branch
            // This could be done more efficiently, if it proves to be time consuming
            std::vector<MWNode<D, T> *> branchstack; // local stack starting from this branch
            branchstack.push_back(thisNode);
            while (branchstack.size() > 0) {
                MWNode<D, T> *branchNode = branchstack.back();
                branchstack.pop_back();
                if (branchNode->getNChildren() > 0) {
                    for (int i = 0; i < branchNode->getNChildren(); i++) { branchstack.push_back(branchNode->children[i]); }
                } else
                    this->endNodeTable.push_back(branchNode);
            }
        }
    }
}

template <int D, typename T> void FunctionTree<D, T>::deleteGenerated() {
    for (int n = 0; n < this->getNEndNodes(); n++) this->getEndMWNode(n).deleteGenerated();
}

template <int D, typename T> void FunctionTree<D, T>::deleteGeneratedParents() {
    for (int n = 0; n < this->getRootBox().size(); n++) this->getRootMWNode(n).deleteParent();
}

template <> int FunctionTree<3, double>::saveNodesAndRmCoeff() {
    if (this->isLocal) MSG_INFO("Tree is already in local representation");
    NodesCoeff = new BankAccount; // NB: must be a collective call!
    int stack_p = 0;
    if (mpi::wrk_rank == 0) {
        int sizecoeff = (1 << 3) * this->getKp1_d();
        std::vector<MWNode<3, double> *> stack; // nodes from this Tree
        for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) { stack.push_back(this->getRootBox().getNodes()[rIdx]); }
        while (stack.size() > stack_p) {
            MWNode<3, double> *Node = stack[stack_p++];
            NodesCoeff->put_data(Node->getNodeIndex(), sizecoeff, Node->getCoefs());
            for (int i = 0; i < Node->getNChildren(); i++) { stack.push_back(Node->children[i]); }
        }
    }
    this->nodeAllocator_p->deallocAllCoeff();
    mpi::broadcast_Tree_noCoeff(*this, mpi::comm_wrk);
    this->isLocal = true;
    assert(this->NodeIndex2serialIx.size() == getNNodes());
    return this->NodeIndex2serialIx.size();
}

template <> int FunctionTree<3, ComplexDouble>::saveNodesAndRmCoeff() {
    if (this->isLocal) MSG_INFO("Tree is already in local representation");
    NodesCoeff = new BankAccount; // NB: must be a collective call!
    int stack_p = 0;
    if (mpi::wrk_rank == 0) {
        int sizecoeff = (1 << 3) * this->getKp1_d();
        sizecoeff *= 2;                                // double->ComplexDouble. Saved as twice as many doubles
        std::vector<MWNode<3, ComplexDouble> *> stack; // nodes from this Tree
        for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) { stack.push_back(this->getRootBox().getNodes()[rIdx]); }
        while (stack.size() > stack_p) {
            MWNode<3, ComplexDouble> *Node = stack[stack_p++];
            NodesCoeff->put_data(Node->getNodeIndex(), sizecoeff, Node->getCoefs());
            for (int i = 0; i < Node->getNChildren(); i++) { stack.push_back(Node->children[i]); }
        }
    }
    this->nodeAllocator_p->deallocAllCoeff();
    mpi::broadcast_Tree_noCoeff(*this, mpi::comm_wrk);
    this->isLocal = true;
    assert(this->NodeIndex2serialIx.size() == getNNodes());
    return this->NodeIndex2serialIx.size();
}

/**  @brief Deep copy of tree
 *
 * @details Exact copy without any binding between old and new tree
 */
template <int D, typename T> void FunctionTree<D, T>::deep_copy(FunctionTree<D, T> *out) {
    copy_grid(*out, *this);
    copy_func(*out, *this);
}

/**  @brief New tree with only real part
 */
template <int D, typename T> FunctionTree<D, double> *FunctionTree<D, T>::Real() {
    FunctionTree<D, double> *out = new FunctionTree<D, double>(this->getMRA(), this->getName());
    out->setZero();
    // TODO: why does the omp parallelism not always work here?
    //#pragma omp parallel num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
	//#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &inp_node = *this->endNodeTable[n];
            MWNode<D, double> &out_node = out->getNode(inp_node.getNodeIndex(), true);
            double *out_coefs = out_node.getCoefs();
            const T *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] = std::real(inp_coefs[i]); }
            out_node.calcNorms();
        }
    }
    out->resetEndNodeTable();
    out->mwTransform(BottomUp);
    out->calcSquareNorm();
    return out;
}

/**  @brief New tree with only imaginary part
 */
template <int D, typename T> FunctionTree<D, double> *FunctionTree<D, T>::Imag() {
    FunctionTree<D, double> *out = new FunctionTree<D, double>(this->getMRA(), this->getName());
    out->setZero();
    //#pragma omp parallel num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
	//#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D, T> &inp_node = *this->endNodeTable[n];
            MWNode<D, double> &out_node = out->getNode(inp_node.getNodeIndex(), true);
            double *out_coefs = out_node.getCoefs();
            const T *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] = std::imag(inp_coefs[i]); }
            out_node.calcNorms();
        }
    }
    out->mwTransform(BottomUp);
    out->calcSquareNorm();
    return out;
}

/*
 * From real to complex tree. Copy everything, and convert double to ComplexDouble for the coefficents.
 * Should use a deep_copy if generalized in the future.
 */

template <> void FunctionTree<3, double>::CopyTreeToComplex(FunctionTree<3, ComplexDouble> *&outTree) {
    delete outTree;
    outTree = new FunctionTree<3, ComplexDouble>(this->getMRA());
    std::vector<MWNode<3, double> *> instack;         // node from this
    std::vector<MWNode<3, ComplexDouble> *> outstack; // node from outTree
    outTree->clearEndNodeTable();
    for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
        instack.push_back(this->getRootBox().getNodes()[rIdx]);
        outstack.push_back(outTree->getRootBox().getNodes()[rIdx]);
    }
    int ncoefs = this->getNodeAllocator().getNCoefs();
    while (instack.size() > 0) {
        // inNode and outNode are the same node in space, but on different trees
        MWNode<3, ComplexDouble> *outNode = outstack.back();
        outstack.pop_back();
        MWNode<3, double> *inNode = instack.back();
        instack.pop_back();
        // copy coefficients:
        double *incoefs = inNode->getCoefs();
        ComplexDouble *outcoefs = outNode->getCoefs();
        for (int i = 0; i < ncoefs; i++) outcoefs[i] = incoefs[i];
        outNode->setHasCoefs();
        outNode->calcNorms();

        if (inNode->getNChildren() > 0) {
            if (outNode->getNChildren() < inNode->getNChildren()) outNode->createChildren(true);
            for (int i = 0; i < inNode->getNChildren(); i++) {
                instack.push_back(inNode->children[i]);
                outstack.push_back(outNode->children[i]);
            }
        } else {
            outTree->endNodeTable.push_back(outNode);
        }
    }
    outTree->calcSquareNorm();
    outTree->calcSquareNorm(true);
}

template <> void FunctionTree<2, double>::CopyTreeToComplex(FunctionTree<2, ComplexDouble> *&outTree) {
    delete outTree;
    outTree = new FunctionTree<2, ComplexDouble>(this->getMRA());
    std::vector<MWNode<2, double> *> instack;         // node from this
    std::vector<MWNode<2, ComplexDouble> *> outstack; // node from outTree
    outTree->clearEndNodeTable();
    for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
        instack.push_back(this->getRootBox().getNodes()[rIdx]);
        outstack.push_back(outTree->getRootBox().getNodes()[rIdx]);
    }
    int ncoefs = this->getNodeAllocator().getNCoefs();
    while (instack.size() > 0) {
        // inNode and outNode are the same node in space, but on different trees
        MWNode<2, ComplexDouble> *outNode = outstack.back();
        outstack.pop_back();
        MWNode<2, double> *inNode = instack.back();
        instack.pop_back();
        // copy coefficients:
        double *incoefs = inNode->getCoefs();
        ComplexDouble *outcoefs = outNode->getCoefs();
        for (int i = 0; i < ncoefs; i++) outcoefs[i] = incoefs[i];
        outNode->setHasCoefs();
        outNode->calcNorms();

        if (inNode->getNChildren() > 0) {
            if (outNode->getNChildren() < inNode->getNChildren()) outNode->createChildren(true);
            for (int i = 0; i < inNode->getNChildren(); i++) {
                instack.push_back(inNode->children[i]);
                outstack.push_back(outNode->children[i]);
            }
        } else {
            outTree->endNodeTable.push_back(outNode);
        }
    }
    outTree->calcSquareNorm();
    outTree->calcSquareNorm(true);
}

template <> void FunctionTree<1, double>::CopyTreeToComplex(FunctionTree<1, ComplexDouble> *&outTree) {
    delete outTree;
    outTree = new FunctionTree<1, ComplexDouble>(this->getMRA());
    std::vector<MWNode<1, double> *> instack;         // node from this
    std::vector<MWNode<1, ComplexDouble> *> outstack; // node from outTree
    outTree->clearEndNodeTable();
    for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
        instack.push_back(this->getRootBox().getNodes()[rIdx]);
        outstack.push_back(outTree->getRootBox().getNodes()[rIdx]);
    }
    int ncoefs = this->getNodeAllocator().getNCoefs();
    while (instack.size() > 0) {
        // inNode and outNode are the same node in space, but on different trees
        MWNode<1, ComplexDouble> *outNode = outstack.back();
        outstack.pop_back();
        MWNode<1, double> *inNode = instack.back();
        instack.pop_back();
        // copy coefficients:
        double *incoefs = inNode->getCoefs();
        ComplexDouble *outcoefs = outNode->getCoefs();
        for (int i = 0; i < ncoefs; i++) outcoefs[i] = incoefs[i];
        outNode->setHasCoefs();
        outNode->calcNorms();

        if (inNode->getNChildren() > 0) {
            if (outNode->getNChildren() < inNode->getNChildren()) outNode->createChildren(true);
            for (int i = 0; i < inNode->getNChildren(); i++) {
                instack.push_back(inNode->children[i]);
                outstack.push_back(outNode->children[i]);
            }
        } else {
            outTree->endNodeTable.push_back(outNode);
        }
    }
    outTree->calcSquareNorm();
    outTree->calcSquareNorm(true);
}

// for testing
template <> void FunctionTree<3, double>::CopyTreeToReal(FunctionTree<3, double> *&outTree) {
    delete outTree;
    // FunctionTree<3, double>* inTree = this;
    outTree = new FunctionTree<3, double>(this->getMRA());
    std::vector<MWNode<3, double> *> instack;  // node from this
    std::vector<MWNode<3, double> *> outstack; // node from outTree
    outTree->clearEndNodeTable();
    for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
        instack.push_back(this->getRootBox().getNodes()[rIdx]);
        outstack.push_back(outTree->getRootBox().getNodes()[rIdx]);
    }
    int ncoefs = this->getNodeAllocator().getNCoefs();
    while (instack.size() > 0) {
        // inNode and outNode are the same node in space, but on different trees
        MWNode<3, double> *outNode = outstack.back();
        outstack.pop_back();
        MWNode<3, double> *inNode = instack.back();
        instack.pop_back();
        // copy coefficients:
        double *incoefs = inNode->getCoefs();
        double *outcoefs = outNode->getCoefs();
        for (int i = 0; i < ncoefs; i++) outcoefs[i] = incoefs[i];
        outNode->setHasCoefs();
        outNode->calcNorms();

        if (inNode->getNChildren() > 0) {
            outNode->clearIsEndNode();
            if (outNode->getNChildren() < inNode->getNChildren()) outNode->createChildren(true);
            for (int i = 0; i < inNode->getNChildren(); i++) {
                instack.push_back(inNode->children[i]);
                outstack.push_back(outNode->children[i]);
            }
        } else {
            outTree->endNodeTable.push_back(outNode);
        }
    }
}

template class FunctionTree<1, double>;
template class FunctionTree<2, double>;
template class FunctionTree<3, double>;

template class FunctionTree<1, ComplexDouble>;
template class FunctionTree<2, ComplexDouble>;
template class FunctionTree<3, ComplexDouble>;

} // namespace mrcpp
