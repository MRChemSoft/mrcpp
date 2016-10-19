#include <iostream>
#include "parallel.h"
#include "TelePrompter.h"
#include "Timer.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "../mrchem/qmfunctions/Orbital.h"
#include "SerialTree.h"
#include "InterpolatingBasis.h"
#include "MultiResolutionAnalysis.h"

#include "mrchem.h"


using namespace std;

const MultiResolutionAnalysis<3> *default_mra=0;//used to define mra when not explicitely defined

int MPI_rank = 0;
int MPI_size = 1;

/** Send or receive a serial tree using MPI
 */
void MPI_Initializations(){

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
#endif

}

MultiResolutionAnalysis<3>* initializeMRA() {
    // Constructing world box
  const int D=3;
    int scale = Input.get<int>("World.scale");
    int max_depth = Input.get<int>("max_depth");
    vector<int> corner = Input.getIntVec("World.corner");
    vector<int> boxes = Input.getIntVec("World.boxes");
    NodeIndex<D> idx(scale, corner.data());
    BoundingBox<D> world(idx, boxes.data());

    // Constructing scaling basis
    int order = Input.get<int>("order");
    InterpolatingBasis basis(order);

    // Initializing MRA
    return new MultiResolutionAnalysis<D>(world, basis, max_depth);
}

void SendRcv_Orbital(Orbital* Orb, int source, int dest, int tag){
#ifdef HAVE_MPI
  MPI_Status status;
  MPI_Comm comm=MPI_COMM_WORLD;

  struct Metadata{
    int spin;
    int occupancy;
    int NchunksReal;
    int NchunksImag;
    double error;
  };

  Metadata Orbinfo;

  if(MPI_rank==source){
    Orbinfo.spin=Orb->getSpin();
    Orbinfo.occupancy=Orb->getOccupancy();
    Orbinfo.error=Orb->getError();
    if(Orb->hasReal()){
      Orbinfo.NchunksReal = Orb->re().getSerialTree()->NodeChunks.size();//should reduce to actual number of chunks
    }else{Orbinfo.NchunksReal = 0;}
    if(Orb->hasImag()){
      Orbinfo.NchunksImag = Orb->im().getSerialTree()->NodeChunks.size();//should reduce to actual number of chunks
    }else{Orbinfo.NchunksImag = 0;}

    int count=sizeof(Metadata);
    MPI_Send(&Orbinfo, count, MPI_BYTE, dest, 0, comm);

    if(Orb->hasReal())SendRcv_SerialTree(&Orb->re(), Orbinfo.NchunksReal, source, dest, tag, comm);
    if(Orb->hasImag())SendRcv_SerialTree(&Orb->im(), Orbinfo.NchunksImag, source, dest, tag*10000, comm);

  }
  if(MPI_rank==dest){
    int count=sizeof(Metadata);
    MPI_Recv(&Orbinfo, count, MPI_BYTE, source, 0, comm, &status);
    Orb->setSpin(Orbinfo.spin);
    Orb->setOccupancy(Orbinfo.occupancy);
    Orb->setError(Orbinfo.error);

    if(Orbinfo.NchunksReal>0){
      if(not Orb->hasReal()){
	//We must have a tree defined for receiving nodes. Define one:
	if(default_mra==0){
	  default_mra = initializeMRA();
	  cout<<" defined new MRA with depth "<<default_mra->getMaxDepth()<<endl;
	}
	Orb->real = new FunctionTree<3>(*default_mra,MAXALLOCNODES);
      }
    SendRcv_SerialTree(&Orb->re(), Orbinfo.NchunksReal, source, dest, tag, comm);}

    if(Orbinfo.NchunksImag>0){
      SendRcv_SerialTree(&Orb->im(), Orbinfo.NchunksImag, source, dest, tag*10000, comm);
    }else{
      //&(Orb->im())=0;
    }
  }

#endif

}
#ifdef HAVE_MPI
template<int D>
void SendRcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int dest, int tag, MPI_Comm comm){
  MPI_Status status;
  Timer timer;
  SerialTree<D>* STree = Tree->getSerialTree();
  

  cout<<MPI_rank<<" STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<endl;
  int count = 1;
  int STreeMeta[count];

  if(MPI_rank==source){

    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
      count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
      cout<<MPI_rank<<" sending "<<STree->maxNodesPerChunk<<" ProjectedNodes to "<<dest<<endl;
      //set SNodeIx of the unused nodes to -1
      int ishift = ichunk*STree->maxNodesPerChunk;
      for (int i = 0; i < STree->maxNodesPerChunk; i++){
	if(STree->NodeStackStatus[ishift+i]!=1)(STree->NodeChunks[ichunk])[i].SNodeIx=-1;
      }
      MPI_Send(STree->NodeChunks[ichunk], count, MPI_BYTE, dest, tag+1+ichunk, comm);
      count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
      cout<<MPI_rank<<" sending "<<count<<" doubles to "<<dest<<endl;
      MPI_Send(STree->NodeCoeffChunks[ichunk], count, MPI_DOUBLE, dest, tag+ichunk*1000, comm);
    }

     timer.stop();
    cout<<" time send     " << timer<<endl;
  }
  if(MPI_rank==dest){

    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
      if(ichunk<STree->NodeChunks.size()){
	STree->SNodesCoeff = STree->NodeCoeffChunks[ichunk];
	STree->SNodes = STree->NodeChunks[ichunk];
      }else{
        STree->SNodesCoeff = new double[STree->sizeNodeCoeff*STree->maxNodesPerChunk];
        STree->NodeCoeffChunks.push_back(STree->SNodesCoeff);
	STree->SNodes = (ProjectedNode<D>*) new char[STree->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	STree->NodeChunks.push_back(STree->SNodes);
      }      
      count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
      MPI_Recv(STree->NodeChunks[ichunk], count, MPI_BYTE, source, tag+1+ichunk, comm, &status);
      cout<<MPI_rank<<" received  "<<count/1024<<" kB with ProjectedNodes from "<<source<<endl;
      count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
      MPI_Recv(STree->NodeCoeffChunks[ichunk], count, MPI_DOUBLE, source, tag+ichunk*1000, comm, &status);
      cout<<MPI_rank<<" received  "<<count<<" coefficients from "<<source<<endl;
    }
    timer.stop();
    cout<<" time receive  " << timer<<endl;
    timer.start();
    STree->RewritePointers(Nchunks);
    timer.stop();
    cout<<" time rewrite pointers  " << timer<<endl;
  }
}
#endif

/** Assign work among MPI processes
 * For a NxN loop. Does not assume symmetry
 * At each iteration data is processed, data is sent and data is received.
 * input: 
 * N, the size of one loop 
 * iter, the iteration 
 * output: 
 * which (i,j) pair to compute next
 * to whom send own data, from who receive data
 */
void Assign_NxN(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvfrom, int* rcvorb){
  int iter = 0;
  int NBlock=N/MPI_size;//number of subblocks
  int M = N%MPI_size; //number of lines and column left after subblocks are taken out
  int Maxiter=-1;
  int iter_base;
  int iBlock = 0;
  int jBlock = 0;
  int MaxiterSofar=0;
  for (int dBlock = 0; dBlock < (N+MPI_size-1)/MPI_size; dBlock++) {

  //1)  //diagonal block
    int imax=MPI_size;
    if(dBlock==NBlock)imax=M;
    int jmax=MPI_size;
    if(dBlock==NBlock)jmax=M;
    for (int i = 0; i < imax; i++) {
      for (int j = 0; j < jmax; j++) {
	iter_base = (i+j)%MPI_size * ((N+MPI_size-1)/MPI_size);
	if(dBlock==NBlock)iter_base = (i+j)%M * ((N+MPI_size-1)/MPI_size);
	iter = iter_base + dBlock*MPI_size*((N+MPI_size-1)/MPI_size);
	int i_glob = i + dBlock*MPI_size;
	int j_glob = j + dBlock*MPI_size;
	if(i==MPI_rank){
	  //this processor receive orbital j_glob from MPIrank=j, and send its own orbital i_glob to MPIrank=j.
	  if(iter>Maxiter)Maxiter=iter;
	  doi[iter]=i_glob;
	  doj[iter]=j_glob;
	  if(i==j){
	    sendto[iter]=-1;
	    sendorb[iter]=-1;
	    rcvorb[iter]=-1;
	  }else{
	    sendto[iter]=j;
	    sendorb[iter]=i_glob;
	    rcvorb[iter]=j_glob;}// contains info both which orbital and indirectly who from
	}
  //2)	//treat all rows and columns that can be generated from received data	
	if(i<=j and i==MPI_rank){
	  //do the corresponding column 
	  //for (iBlock = 0; iBlock < NBlock; iBlock++) {
	  for (int icol = i; icol < N; icol+=MPI_size) {
	    iBlock=icol/MPI_size;
	    if(iBlock!=dBlock){
	      int sh=1;
	      if(iBlock>dBlock)sh=0;//not very elegant
	      int i_glob = i + iBlock*MPI_size;
	      int j_glob = j + dBlock*MPI_size;
	      doi[iter+iBlock+sh] = i_glob;
	      doj[iter+iBlock+sh] = j_glob;
	      rcvorb[iter+iBlock+sh]=-1;//indicates "use same as before"
	      sendorb[iter+iBlock+sh]=-1;//indicates "do not send anything"
	      sendto[iter+iBlock+sh]=-1;//indicates "do not send anything"
	      if(iter+iBlock+sh>Maxiter)Maxiter=iter+iBlock+sh;
	    }
	  }
	}else if(i==MPI_rank){
	  //do the corresponding row
	  //for (jBlock = 0; jBlock < NBlock; jBlock++) {
	  for (int jrow = j; jrow < N; jrow+=MPI_size) {
	    jBlock=jrow/MPI_size;
	    if(jBlock!=dBlock){
	      int sh=1;
	      if(jBlock>dBlock)sh=0;//not very elegant
	    int i_glob = i + dBlock*MPI_size;
	    int j_glob = j + jBlock*MPI_size;
	    
	    doi[iter+jBlock+sh] = i_glob;
	    doj[iter+jBlock+sh] = j_glob;
	    rcvorb[iter+jBlock+sh]=-1;//indicates "use same as before"
	    sendorb[iter+jBlock+sh]=-1;//indicates "do not send anything"
	    sendto[iter+jBlock+sh]=-1;//indicates "do not send anything"
	    if(iter+jBlock+sh>Maxiter)Maxiter=iter+jBlock+sh;
	    }
	  }
	}
	MaxiterSofar=Maxiter;
      }
    }
  }

}

/** Define all the different MPI groups
 */
void define_groups(){

#ifdef HAVE_MPI
  int MPI_worldsize;
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_worldsize);
  int MPI_orbital_group_size = MPI_worldsize;//temporary definition
  println(0,"MPI world size: "<< MPI_worldsize ); 
  int MPI_worldrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_worldrank);
  int i = (MPI_worldrank) % MPI_orbital_group_size;

  //divide the world into groups
  //NB: each group has its own group communicator definition!
  MPI_Comm MPI_Comm_Orbital;
  MPI_Comm_split(MPI_COMM_WORLD, 0, i, &MPI_Comm_Orbital);//0 is color (same color->in same group), i is key (orders rank in the group)
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_worldsize);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_orbital_group_size);
  println(0,"orbital group size: "<< MPI_orbital_group_size); 
#endif
  
}

/** Define MPI groups that share the same compute-node i.e. the same memory and array with shared memory between
 * sh_size is size in Bytes
 * returns a pointer to the shared memory
void Share_memory(MPI_Comm ncomm, MPI_Comm ncomm_sh, int sh_size, double * d_ptr){

  //split ncomm into groups that can share memory
  MPI_Comm_split_type(ncomm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &ncomm_sh);

  int shrank;
  MPI_Comm_rank(ncomm_sh, &shrank);
  MPI_Aint size = (shrank==0) ? sh_size : 0;//rank 0 defines length of segment 
  MPI_Win win = MPI_WIN_NULL;
  int disp_unit = 16;//in order for the compiler to keep aligned. 
  //sh_size is in bytes
  MPI_Win_allocate_shared(sh_size, disp_unit, MPI_INFO_NULL, ncomm_sh, &d_ptr, &win);

  MPI_Win_fence(0, win);//wait until finished

  MPI_Aint qsize = 0;
  int qdisp = 0;
  int * qbase = NULL;
  MPI_Win_shared_query(win, 0, &qsize, &qdisp, &d_ptr);
  printf("me = %d, allocated: size=%zu, bytes at = %p\n", shrank, qsize, d_ptr);

  qbase[shrank]=shrank;
  MPI_Barrier(ncomm);
  for (int i=0; i<sh_size; i++){
    printf("shared: me=%d, i=%d, val=%d\n", shrank, i, qbase[i]);
  }
  //MPI_Win_free(&win);//to include at the end of the program
  //MPI_Comm_free(&ncomm_sh);//to include at the end of the program?

} */

#ifdef HAVE_MPI
template void SendRcv_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int source, int dest, int tag, MPI_Comm comm);
template void SendRcv_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int source, int dest, int tag, MPI_Comm comm);
template void SendRcv_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int source, int dest, int tag, MPI_Comm comm);
#endif

