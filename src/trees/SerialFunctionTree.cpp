#include "SerialFunctionTree.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "GenNode.h"
#include "utils/Printer.h"
#include "utils/mpi_utils.h"

namespace mrcpp {

/** SerialTree class constructor.
 * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
 * Initializes rootNodes to represent the zero function and allocate their nodes.
 * NOTES:
 * Serial trees are made of projected nodes, and include gennodes and loose nodes separately.
 * All created (using class creator) Projected nodes or GenNodes are loose nodes.
 * Loose nodes have their coeff in serial Tree, but not the node part.
 * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~ProjectedNode and ~GenNode.
 * Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
 * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
 */
template<int D>
SerialFunctionTree<D>::SerialFunctionTree(FunctionTree<D> *tree, SharedMemory *mem)
    : SerialTree<D>(tree, mem),
      nGenNodes(0),
      lastNode(0),
      lastGenNode(0) {

    //Size for GenNodes chunks. ProjectedNodes will be 8 times larger
    this->sizeGenNodeCoeff = this->tree_p->getKp1_d();//One block
    this->sizeNodeCoeff = (1<<D)*this->sizeGenNodeCoeff;//TDim  blocks
    println(10, "SizeNode Coeff (kB) " << this->sizeNodeCoeff*sizeof(double)/1024);
    println(10, "SizeGenNode Coeff (kB) " << this->sizeGenNodeCoeff*sizeof(double)/1024);

    int sizePerChunk = 2*1024*1024;// 2 MB small for no waisting place, but large enough so that latency and overhead work is negligible
    if(D<3){
        //define rather from number of nodes per chunk
        this->maxNodesPerChunk = 64;
        sizePerChunk = this->maxNodesPerChunk*this->sizeNodeCoeff;
    }else{
        this->maxNodesPerChunk = (sizePerChunk/this->sizeNodeCoeff/sizeof(double)/8)*8;
    }

    this->lastNode = (ProjectedNode<D>*) this->sNodes;//position of last allocated node
    this->lastGenNode = this->sGenNodes;//position of last allocated Gen node

    //make virtual table pointers
    auto* tmpNode = new ProjectedNode<D>();
    this->cvptr_ProjectedNode =  *(char**)(tmpNode);
    delete tmpNode;

    auto* tmpGenNode = new GenNode<D>();
    this->cvptr_GenNode =  *(char**)(tmpGenNode);
    delete tmpGenNode;

#ifdef HAVE_OPENMP
    omp_init_lock(&Sfunc_tree_lock);
#endif
}

/** SerialTree destructor. */
template<int D>
SerialFunctionTree<D>::~SerialFunctionTree() {
    for (int i = 0; i < this->genNodeCoeffChunks.size(); i++) delete[] this->genNodeCoeffChunks[i];
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[] (char*)(this->nodeChunks[i]);
    if(not this->isShared())//if the data is shared, it must be freed by MPI_Win_free
        for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    for (int i = 0; i < this->genNodeChunks.size(); i++) delete[] (char*)(this->genNodeChunks[i]);

    this->nodeStackStatus.clear();
    this->genNodeStackStatus.clear();

#ifdef HAVE_OPENMP
    omp_destroy_lock(&Sfunc_tree_lock);
#endif
}

template<int D>
void SerialFunctionTree<D>::allocRoots(MWTree<D> &tree) {
    int sIx;
    double *coefs_p;
    //reserve place for nRoots
    int nRoots = tree.getRootBox().size();
    ProjectedNode<D> *root_p = this->allocNodes(nRoots, &sIx, &coefs_p);

    MWNode<D> **roots = tree.getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        roots[rIdx] = root_p;

        *(char**)(root_p) = this->cvptr_ProjectedNode;

        root_p->tree = &tree;
        root_p->parent = 0;
        for (int i = 0; i < root_p->getTDim(); i++) {
            root_p->children[i] = 0;
        }

        root_p->nodeIndex = tree.getRootBox().getNodeIndex(rIdx);
        root_p->hilbertPath = HilbertPath<D>();

        root_p->n_coefs = this->sizeNodeCoeff;
        root_p->coefs = coefs_p;

        root_p->lockX = 0;
        root_p->serialIx = sIx;
        root_p->parentSerialIx = -1;//to indicate rootnode
        root_p->childSerialIx = -1;

        root_p->status = 0;

        root_p->clearNorms();
        root_p->setIsLeafNode();
        root_p->setIsAllocated();
        root_p->clearHasCoefs();
        root_p->setIsEndNode();
        root_p->setIsRootNode();

        tree.incrementNodeCount(root_p->getScale());

        sIx++;
        root_p++;
        coefs_p += this->sizeNodeCoeff;
    }
}

template<int D>
void SerialFunctionTree<D>::allocChildren(MWNode<D> &parent) {
    int sIx;
    double *coefs_p;
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    ProjectedNode<D> *child_p = this->allocNodes(nChildren, &sIx, &coefs_p);

    //position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char**)(child_p) = this->cvptr_ProjectedNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) {
            child_p->children[i] = 0;
        }

        child_p->nodeIndex = NodeIndex<D>(parent.getNodeIndex(), cIdx);
        child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

        child_p->n_coefs = this->sizeNodeCoeff;
        child_p->coefs = coefs_p;

        child_p->lockX = 0;
        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->status = 0;

        child_p->clearNorms();
        child_p->setIsLeafNode();
        child_p->setIsAllocated();
        child_p->clearHasCoefs();
        child_p->setIsEndNode();

        child_p->tree->incrementNodeCount(child_p->getScale());

        sIx++;
        child_p++;
        coefs_p += this->sizeNodeCoeff;
    }
}

template<int D>
void SerialFunctionTree<D>::allocGenChildren(MWNode<D> &parent) {
    int sIx;
    double *coefs_p;
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    GenNode<D>* child_p = this->allocGenNodes(nChildren, &sIx, &coefs_p);

    //position of first child
    parent.childSerialIx = sIx;//not used fro Gennodes?
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char**)(child_p) = this->cvptr_GenNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) {
            child_p->children[i] = 0;
        }

        child_p->nodeIndex = NodeIndex<D>(parent.getNodeIndex(), cIdx);
        child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

        child_p->n_coefs = this->sizeGenNodeCoeff;
        child_p->coefs = coefs_p;

        child_p->lockX = 0;
        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->status = 0;

        child_p->clearNorms();
        child_p->setIsLeafNode();
        child_p->setIsAllocated();
        child_p->clearHasCoefs();
        child_p->setIsGenNode();

        child_p->tree->incrementGenNodeCount();

        sIx++;
        child_p++;
        coefs_p += this->sizeGenNodeCoeff;
    }
}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* SerialFunctionTree<D>::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    *serialIx = this->nNodes;
    int chunkIx = *serialIx%(this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx+nAlloc > this->maxNodesPerChunk ) {
        //we want nodes allocated simultaneously to be allocated on the same piece.
        //possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk*((this->nNodes+nAlloc-1)/this->maxNodesPerChunk);//start of next chunk

        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk

        //careful: nodeChunks.size() is an unsigned int
        if (chunk+1 > this->nodeChunks.size()){
            //need to allocate new chunk
            double *sNodesCoeff;
            if (this->isShared()) {
                //for coefficients, take from the shared memory block
                sNodesCoeff = this->shMem->sh_end_ptr;
                this->shMem->sh_end_ptr += (this->sizeNodeCoeff*this->maxNodesPerChunk);
                //may increase size dynamically in the future
                if (this->shMem->sh_max_ptr < this->shMem->sh_end_ptr) {
                    MSG_FATAL("Shared block too small");
                }
            } else {
                sNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];
            }

            this->nodeCoeffChunks.push_back(sNodesCoeff);
            this->sNodes = (ProjectedNode<D>*) new char[this->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
            for (int i = 0; i < this->maxNodesPerChunk; i++) {
                this->sNodes[i].serialIx = -1;
                this->sNodes[i].parentSerialIx = -1;
                this->sNodes[i].childSerialIx = -1;
            }
            this->nodeChunks.push_back(this->sNodes);

            //allocate new chunk in nodeStackStatus
            int oldsize = this->nodeStackStatus.size();
            int newsize = oldsize + this->maxNodesPerChunk;
            for (int i = oldsize; i < newsize; i++) this->nodeStackStatus.push_back(0);
            this->maxNodes = newsize;

            if (chunk%100==99 and D==3) println(10, std::endl<<" number of nodes "<<this->nNodes <<",number of Nodechunks now " << this->nodeChunks.size()<<", total size coeff  (MB) "<<(this->nNodes * this->sizeNodeCoeff)/1024/128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx%(this->maxNodesPerChunk);
    }
    assert((this->nNodes+nAlloc-1)/this->maxNodesPerChunk < this->nodeChunks.size());

    ProjectedNode<D> *newNode  = this->lastNode;
    ProjectedNode<D> *newNode_cp  = newNode;

    int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk
    *coefs_p = this->nodeCoeffChunks[chunk] + chunkIx*this->sizeNodeCoeff;

    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx+i] != 0)
            println(0, *serialIx+i<<" NodeStackStatus: not available " << this->nodeStackStatus[*serialIx+i]);
        this->nodeStackStatus[*serialIx+i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode += nAlloc;

    return newNode;
}

template<int D>
void SerialFunctionTree<D>::deallocNodes(int serialIx) {
    if (this->nNodes < 0) {
        println(0, "minNodes exceeded " << this->nNodes);
        this->nNodes++;
    }
    this->nodeStackStatus[serialIx] = 0;//mark as available
    if (serialIx == this->nNodes-1) {//top of stack
        int topStack = this->nNodes;
        while (this->nodeStackStatus[topStack-1] == 0){
            topStack--;
            if (topStack < 1) break;
        }
        this->nNodes = topStack;//move top of stack
        //has to redefine lastNode
        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
    }
}
/*
void MWTree<D>::crop(double prec, double splitFac, bool absPrec) {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = getRootMWNode(i);
        root.crop(prec, splitFac, absPrec);
    }
    resetEndNodeTable();
    calcSquareNorm();
    }*/

/** Reduce the accuracy of the tree by deleting nodes
  * which have a higher precision than the requested precison.
  * By default, the relative precision of the tree is used.
template<int D>
void SerialFunctionTree<D>::crop(double prec, double splitFac, bool absPrec) {
    MWTree<D>* tree = this->tree_p;
    for (int i = 0; i < tree->rootBox.size(); i++) {
        MWNode<D> &root = tree->getRootMWNode(i);
        root.crop(prec, splitFac, absPrec);
    }
    cropChunks();
    tree->resetEndNodeTable();
    tree->calcSquareNorm();
    }*/

/** Fill all holes in the chunks with occupied nodes, then remove all empty chunks */
template<int D>
int SerialFunctionTree<D>::cropChunks() {
    if (this->maxNodesPerChunk*this->nodeChunks.size() <= this->getTree()->getNNodes()+this->maxNodesPerChunk) {
        return 0; //no chunks to remove
    }

    double *coefs_p;
    MWTree<D>* tree = this->tree_p;
    int posavail = tree->getRootBox().size();//start after root nodes
    int posocc = 0;
    int nChunksStart = this->nodeChunks.size();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==1){
        int Ix=2;
        int n_ichunk = Ix/this->maxNodesPerChunk;
        int n_inode = Ix%this->maxNodesPerChunk;
        ProjectedNode<D>* Node = this->nodeChunks[n_ichunk] + n_inode;
        //        if(rank==1)std::cout<<rank<<" START "<<Node->serialIx<<" "<<Node->parentSerialIx<<" child "<<Node->childSerialIx<<std::endl;

    }
    if(true){
    while (true) {
        //find first available spot
        while (this->nodeStackStatus[posavail]!=0 and posavail<this->nNodes) posavail++;
        if( posavail>=this->nNodes ) break;//treated all nodes
        //next allocated spot after available
        posocc = posavail;
        while (this->nodeStackStatus[posocc]==0 and posavail<this->nNodes) posocc++;
        if( posocc>=this->nNodes ) break;//treated all nodes

        int nAlloc = (1<<D);

        //move node from posocc to posavail
        int ichunk = posocc/this->maxNodesPerChunk;
        int inode = posocc%this->maxNodesPerChunk;
        ProjectedNode<D>* NodeOcc = this->nodeChunks[ichunk] + inode;

        //check that all siblings are consecutive. Should never be root node.
        for(int i = 0; i < nAlloc; i++) assert(this->nodeStackStatus[posavail+i]==0);
        for(int i = 1; i < nAlloc; i++) assert((NodeOcc+i)->parent->serialIx == (NodeOcc)->parent->serialIx);//siblings

        ichunk = posavail/this->maxNodesPerChunk;
        inode = posavail%this->maxNodesPerChunk;
        ProjectedNode<D>* NodeAvail = this->nodeChunks[ichunk] + inode;
        for(int i = 0; i < nAlloc*this->sizeNodeCoeff; i++) NodeAvail->coefs[i] = NodeOcc->coefs[i];

        //just copy everything "as is"
        for(int i = 0; i < nAlloc*sizeof(ProjectedNode<D>); i++) ((char*) NodeAvail)[i] = ((char*) NodeOcc)[i];

        //new nodes have another adress
        for(int i = 0; i < nAlloc; i++)  (NodeAvail + i)->serialIx = posavail + i;

        //new nodes have another adress. Tell parent
        NodeAvail->parent->childSerialIx = posavail;
        for(int i = 0; i < nAlloc; i++) NodeAvail->parent->children[i] = NodeAvail + i;
                                               //        if(rank==1 and NodeAvail->parent->serialIx==2)std::cout<<NodeAvail->serialIx<<" writing parent "<<NodeAvail->parent->serialIx<<" parentchild "<<NodeAvail->parent->childSerialIx<<std::endl;
        //Tell children too
        for(int i = 0; i < nAlloc; i++) {
            for(int j = 0; j < (NodeAvail+i)->getNChildren(); j++) (NodeAvail+i)->children[j]->parentSerialIx = posavail+i;
            for(int j = 0; j < (NodeAvail+i)->getNChildren(); j++) (NodeAvail+i)->children[j]->parent = NodeAvail+i;
        }

        //mark moved nodes as occupied
        for(int i = 0; i < nAlloc; i++) this->nodeStackStatus[posavail+i] = 1;
        posavail += nAlloc;

        //delete "old" nodes
        for(int i = 0; i < nAlloc; i++) this->nodeStackStatus[posocc+i] = 0;
        for(int i = 0; i < nAlloc; i++) (NodeOcc +i)->serialIx = -1;

        //check
                                               /*        if(NodeAvail->parent->childSerialIx!=NodeAvail->serialIx)std::cout<<NodeAvail->parent->childSerialIx<<" "<<NodeAvail->parent->serialIx<<" EError parent child ix "<<NodeAvail->serialIx<<" "<<NodeAvail->parentSerialIx<<std::endl;
        if(rank==1 and(NodeAvail->serialIx==56 or NodeAvail->serialIx==21 or NodeAvail->serialIx==29 or NodeAvail->parentSerialIx==21 or  NodeAvail->serialIx==104) )std::cout<<rank<<" "<<NodeAvail->serialIx<<" parent "<<NodeAvail->parent->serialIx<<" child "<<NodeAvail->childSerialIx<<" parentchild "<<NodeAvail->parent->childSerialIx<<" Occ "<<NodeOcc->serialIx<<" Occparentchild"<<NodeOcc->parent->childSerialIx<<std::endl;
                                               */
    }
    }
    /*    if(rank==1){
        int Ix=21;
        int n_ichunk = Ix/this->maxNodesPerChunk;
        int n_inode = Ix%this->maxNodesPerChunk;
        ProjectedNode<D>* Node = this->nodeChunks[n_ichunk] + n_inode;
        if(rank==1)std::cout<<rank<<" "<<Node->serialIx<<" parent "<<Node->parent->serialIx<<" "<<Node->parentSerialIx<<" child "<<Node->childSerialIx<<" parentchild "<<Node->parent->childSerialIx<<std::endl;
        Ix=56;
        n_ichunk = Ix/this->maxNodesPerChunk;
        n_inode = Ix%this->maxNodesPerChunk;
        Node = this->nodeChunks[n_ichunk] + n_inode;
        if(rank==1)std::cout<<rank<<" "<<Node->serialIx<<" parent "<<Node->parent->serialIx<<" child "<<Node->childSerialIx<<" parentchild "<<Node->parent->childSerialIx<<std::endl;
        Ix=2;
        n_ichunk = Ix/this->maxNodesPerChunk;
        n_inode = Ix%this->maxNodesPerChunk;
        Node = this->nodeChunks[n_ichunk] + n_inode;
        if(rank==1)std::cout<<rank<<" "<<Node->serialIx<<" child "<<Node->childSerialIx<<std::endl;


        }*/

    //find the last used node
    posocc = this->nNodes-1;
    while (this->nodeStackStatus[posocc]==0 and posocc>0) posocc--;
    this->nNodes = posocc + 1;
    int ichunk = posocc/this->maxNodesPerChunk;
    int inode = posocc%this->maxNodesPerChunk;
    this->lastNode = this->nodeChunks[ichunk] + inode;

    //assert(posavail >= tree.getRootBox().size());//at least root nodes must be used
    int nChunks = posocc/this->maxNodesPerChunk + 1;//number of occupied chunks
    for(int i = nChunks; i < this->nodeChunks.size(); i++) delete[] (char*)(this->nodeChunks[i]);//remove unused chunks
    if(not this->isShared())//if the data is shared, it must be freed by MPI_Win_free
        for (int i = nChunks; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    //shrink the stacks
    this->nodeChunks.resize(nChunks);
    this->nodeCoeffChunks.resize(nChunks);

    this->maxNodes = this->nodeStackStatus.size();

    this->getTree()->resetEndNodeTable();

    rewritePointers(nChunks);//clean up all valid addresses WHY NEEDED??
    return nChunksStart - nChunks;
}

//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* SerialFunctionTree<D>::allocGenNodes(int nAlloc, int *serialIx, double **coefs_p) {
    omp_set_lock(&Sfunc_tree_lock);
    *serialIx = this->nGenNodes;
    int chunkIx = *serialIx%(this->maxNodesPerChunk);

    //Not necessarily wrong, but new:
    assert(nAlloc == (1<<D));

    if(chunkIx == 0  or  chunkIx+nAlloc > this->maxNodesPerChunk ){
        //start on new chunk
        //we want nodes allocated simultaneously to be allocated on the same chunk.
        //possibly jump over the last nodes from the old chunk
        this->nGenNodes=this->maxNodesPerChunk*((this->nGenNodes+nAlloc-1)/this->maxNodesPerChunk);//start of next chunk

        int chunk = this->nGenNodes/this->maxNodesPerChunk;//find the right chunk

        //careful: nodeChunks.size() is an unsigned int
        if(chunk+1 > this->genNodeChunks.size()){
            //need to allocate new chunk
            this->sGenNodes = (GenNode<D>*) new char[this->maxNodesPerChunk*sizeof(GenNode<D>)];
            for (int i = 0; i < this->maxNodesPerChunk; i++) {
                this->sGenNodes[i].serialIx = -1;
                this->sGenNodes[i].parentSerialIx = -1;
                this->sGenNodes[i].childSerialIx = -1;
            }
            this->genNodeChunks.push_back(this->sGenNodes);
            auto *sGenNodesCoeff = new double[this->sizeGenNodeCoeff*this->maxNodesPerChunk];
            this->genNodeCoeffChunks.push_back(sGenNodesCoeff);
            //allocate new chunk in nodeStackStatus
            int oldsize = this->genNodeStackStatus.size();
            int newsize = oldsize + this->maxNodesPerChunk;
            for (int i = oldsize; i < newsize; i++) this->genNodeStackStatus.push_back(0);
            this->maxGenNodes = newsize;

            if(chunk%100==99 and D==3)println(10, "\n number of GenNodes "<<this->nGenNodes <<",number of GenNodechunks now " << this->genNodeChunks.size()<<", total size coeff  (MB) "<<(this->nGenNodes/1024) * this->sizeGenNodeCoeff/128);
        }
        this->lastGenNode = this->genNodeChunks[chunk] + this->nGenNodes%(this->maxNodesPerChunk);
        *serialIx = this->nGenNodes;
        chunkIx = *serialIx%(this->maxNodesPerChunk);
    }
    assert((this->nGenNodes+nAlloc-1)/this->maxNodesPerChunk < this->genNodeChunks.size());

    GenNode<D>* newNode  = this->lastGenNode;
    GenNode<D>* newNode_cp  = newNode;

    int chunk = this->nGenNodes/this->maxNodesPerChunk;//find the right chunk
    *coefs_p = this->genNodeCoeffChunks[chunk] + chunkIx*this->sizeGenNodeCoeff;

    for (int i = 0; i < nAlloc; i++) {
        newNode_cp->serialIx = *serialIx+i;//Until overwritten!
        if (this->genNodeStackStatus[*serialIx+i] != 0)
            println(0, *serialIx+i<<" NodeStackStatus: not available " << this->genNodeStackStatus[*serialIx+i]);
        this->genNodeStackStatus[*serialIx+i] = 1;
        newNode_cp++;
    }
    this->nGenNodes += nAlloc;
    this->lastGenNode += nAlloc;

    omp_unset_lock(&Sfunc_tree_lock);
    return newNode;
}

template<int D>
void SerialFunctionTree<D>::deallocGenNodes(int serialIx) {
    omp_set_lock(&Sfunc_tree_lock);
    if (this->nGenNodes <0) {
        println(0, "minNodes exceeded " << this->nGenNodes);
        this->nGenNodes++;
    }
    this->genNodeStackStatus[serialIx] = 0;//mark as available
    if (serialIx == this->nGenNodes-1) {//top of stack
        int topStack = this->nGenNodes;
        while (this->genNodeStackStatus[topStack-1] == 0) {
            topStack--;
            if (topStack < 1) {
                //remove all the GenNodeChunks once there are noe more genNodes
                this->deallocGenNodeChunks();
                break;
            }
        }
        this->nGenNodes = topStack;//move top of stack
        //has to redefine lastGenNode
        int chunk = this->nGenNodes/this->maxNodesPerChunk;//find the right chunk
        this->lastGenNode = this->genNodeChunks[chunk] + this->nGenNodes%(this->maxNodesPerChunk);
    }
    omp_unset_lock(&Sfunc_tree_lock);
}

template<int D>
void SerialFunctionTree<D>::deallocGenNodeChunks() {
    //if(mpiOrbRank==0 and (this->genNodeCoeffChunks.size()*2)*1024/8>10000)cout<<"deallocate genchunks MB "<<(this->genNodeCoeffChunks.size()*2)*1024/1024/8<<endl;
    for (int i = 0; i < this->genNodeCoeffChunks.size(); i++) delete[] this->genNodeCoeffChunks[i];
    for (int i = 0; i < this->genNodeChunks.size(); i++) delete[] (char*)(this->genNodeChunks[i]);
    this->genNodeCoeffChunks.clear();
    this->genNodeChunks.clear();
    this->genNodeStackStatus.clear();
}

/** Overwrite all pointers defined in the tree.
 * Necessary after sending the tree
 * could be optimized. Should reset other counters? (GenNodes...) */
template<int D>
void SerialFunctionTree<D>::rewritePointers(int nChunks){

    this->getTree()->nNodes = 0;
    this->getTree()->nodesAtDepth.clear();
    this->getTree()->squareNorm = 0.0;

    //reinitialize stacks
    int nodecount = nChunks * this->maxNodesPerChunk;
    while (nodecount > this->nodeStackStatus.size()) this->nodeStackStatus.push_back(0);
    this->maxNodes = nodecount;

    //clear all gennodes and their chunks:
    for (int i = 0; i < this->genNodeCoeffChunks.size(); i++) delete[] this->genNodeCoeffChunks[i];
    for (int i = 0; i < this->genNodeChunks.size(); i++) delete[] (char*)(this->genNodeChunks[i]);
    this->genNodeCoeffChunks.clear();
    this->genNodeStackStatus.clear();
    this->genNodeChunks.clear();
    this->nGenNodes = 0;
    this->maxGenNodes = 0;
    this->nNodes = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(int ichunk = 0 ; ichunk < nChunks; ichunk++) {
        for(int inode = 0 ; inode < this->maxNodesPerChunk; inode++) {
            ProjectedNode<D>* node = (this->nodeChunks[ichunk]) + inode;
            if (node->serialIx >= 0) {
                this->nNodes = ichunk*this->maxNodesPerChunk + (inode + 1);

                //Node is part of tree, should be processed
                this->getTree()->incrementNodeCount(node->getScale());
                if (node->isEndNode()) this->getTree()->squareNorm += node->getSquareNorm();

                //normally (intel) the virtual table does not change, but we overwrite anyway
                *(char**)(node) = this->cvptr_ProjectedNode;

                node->tree = this->getTree();

                //"adress" of coefs is the same as node, but in another array
                node->coefs = this->nodeCoeffChunks[ichunk]+ inode*this->sizeNodeCoeff;

                //adress of parent and children must be corrected
                //can be on a different chunks
                if (node->parentSerialIx >= 0) {
                    int n_ichunk = node->parentSerialIx/this->maxNodesPerChunk;
                    int n_inode = node->parentSerialIx%this->maxNodesPerChunk;
                    node->parent = this->nodeChunks[n_ichunk] + n_inode;
                    if(node->parent->serialIx!=node->parentSerialIx)std::cout<<node->parentSerialIx<<" error parent ix "<<node->parent->serialIx<<std::endl;
                    if(rank==1 and node->parent->childSerialIx!=node->serialIx-node->serialIx%8)std::cout<<" "<<node->parent->childSerialIx<<" "<<node->parent->serialIx<<" Error parent child ix "<<node->serialIx<<std::endl;
                } else {
                    node->parent = 0;
                }

                for (int i = 0; i < node->getNChildren(); i++) {
                    int n_ichunk = (node->childSerialIx+i)/this->maxNodesPerChunk;
                    int n_inode = (node->childSerialIx+i)%this->maxNodesPerChunk;
                    node->children[i] = this->nodeChunks[n_ichunk] + n_inode;
                }
                this->nodeStackStatus[node->serialIx] = 1;//occupied
            }

        }
        this->lastNode = this->nodeChunks[ichunk] + this->nNodes%(this->maxNodesPerChunk);
    }

    //update other MWTree data
    auto* Tree = static_cast<FunctionTree<D>*> (this->tree_p);

    NodeBox<D> &rBox = Tree->getRootBox();
    MWNode<D> **roots = rBox.getNodes();

    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        roots[rIdx] = (this->nodeChunks[0]) + rIdx;//adress of roots are at start of NodeChunks[0] array
    }

    this->getTree()->resetEndNodeTable();
 }

template<int D>
int SerialFunctionTree<D>::getNChunksUsed() const {
    int lastUsed = 0;
    for (int iChunk = 0; iChunk < getNChunks(); iChunk++) {
        int iShift = iChunk * this->maxNodesPerChunk;
        bool chunkUsed = false;
        for (int i = 0; i < this->maxNodesPerChunk; i++) {
            if (this->nodeStackStatus[iShift+i] == 1) {
                chunkUsed = true;
            }
        }
        if (chunkUsed) lastUsed = iChunk + 1;
    }
    return lastUsed;
}

template class SerialFunctionTree<1>;
template class SerialFunctionTree<2>;
template class SerialFunctionTree<3>;
} // namespace mrcpp
