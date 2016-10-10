#include <iostream>
#include "parallel.h"
#include "TelePrompter.h"
#include "Timer.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "../mrchem/qmfunctions/Orbital.h"
#include "SerialTree.h"

using namespace std;

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
#ifdef HAVE_MPI

void SendRcv_Orbital(Orbital* Orb, int source, int dest, int tag, MPI_Comm comm){
  MPI_Status status;

  //for now brute force...
  if(MPI_rank==source){
    int spin=Orb->getSpin();
    int occupancy=Orb->getOccupancy();
    double error=Orb->getError();
    MPI_Send(&spin, 1, MPI_INT, dest, 0, comm);
    MPI_Send(&occupancy, 1, MPI_INT, dest, 1, comm);
    MPI_Send(&error, 1, MPI_DOUBLE, dest, 2, comm);

    SendRcv_SerialTree(&Orb->re(), source, dest, tag, comm);
    SendRcv_SerialTree(&Orb->im(), source, dest, tag*10000, comm);
  }
  if(MPI_rank==dest){
    int spin;
    int occupancy;
    double error;
    MPI_Recv(&spin, 1, MPI_INT, source, 0, comm, &status);
    MPI_Recv(&occupancy, 1, MPI_INT, source, 1, comm, &status);
    MPI_Recv(&error, 1, MPI_DOUBLE, source, 2, comm, &status);
    Orb->setSpin(spin);
    Orb->setOccupancy(occupancy);
    Orb->setError(error);

    SendRcv_SerialTree(&Orb->re(), source, dest, tag, comm);
    SendRcv_SerialTree(&Orb->im(), source, dest, tag*10000, comm);
  }


}
template<int D>
void SendRcv_SerialTree(FunctionTree<D>* Tree, int source, int dest, int tag, MPI_Comm comm){
  MPI_Status status;
  Timer timer;
  SerialTree<D>* STree = Tree->getSerialTree();
  

  cout<<MPI_rank<<" STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<endl;
  int count = 1;
  int STreeMeta[count];

  if(MPI_rank==source){

    //send first metadata
    int Nchunks = STree->NodeChunks.size();//number of chunks to send
    STreeMeta[0] = Nchunks ;
    MPI_Send(STreeMeta, count, MPI_INT, dest, tag, comm);

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
    MPI_Recv(STreeMeta, count, MPI_INT, source, tag, comm, &status);
    int Nchunks =STreeMeta[0] ;//number of chunks to receive
    
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
template void SendRcv_SerialTree<1>(FunctionTree<1>* STree, int source, int dest, int tag, MPI_Comm comm);
template void SendRcv_SerialTree<2>(FunctionTree<2>* STree, int source, int dest, int tag, MPI_Comm comm);
template void SendRcv_SerialTree<3>(FunctionTree<3>* STree, int source, int dest, int tag, MPI_Comm comm);
#endif

