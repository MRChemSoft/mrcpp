#include <iostream>
#include "parallel.h"
#include "TelePrompter.h"
#include "Timer.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "SerialFunctionTree.h"
#include "SerialTree.h"
#include "InterpolatingBasis.h"
#include "MultiResolutionAnalysis.h"


using namespace std;

int MPI_Orb_rank = 0;
int MPI_Orb_size = 1;
int MPI_SH_rank = 0;
int MPI_SH_size = 1;
int MPI_SH_group_rank = 0;
int MPI_SH_group_size = 1;

#ifdef HAVE_MPI
MPI_Comm MPI_Comm_Orb;
MPI_Comm MPI_Comm_SH;
MPI_Comm MPI_Comm_SH_group;
#endif

/** Send or receive a serial tree using MPI
 */
void MPI_Initializations(){

#ifdef HAVE_MPI
    MPI_Init(NULL, NULL);
    define_MPI_groups();
#endif
}

#ifdef HAVE_MPI
template<int D>
void Send_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  
    println(10," STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<" sending to "<<dest);
    int count = 1;
    int STreeMeta[count];
  
    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	//set serialIx of the unused nodes to -1
	int ishift = ichunk*STree->maxNodesPerChunk;
	for (int i = 0; i < STree->maxNodesPerChunk; i++){
	    if(STree->nodeStackStatus[ishift+i]!=1)(STree->nodeChunks[ichunk])[i].setSerialIx(-1);
	}
	MPI_Send(STree->nodeChunks[ichunk], count, MPI_BYTE, dest, tag+1+ichunk, comm);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if (not STree->isShared or MPI_SH_group_rank != dest/MPI_SH_size){
	    MPI_Send(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, dest, tag+ichunk*1000, comm);}
    }
  
    timer.stop();
    println(10," time send     " << timer); 
}
#endif

#ifdef HAVE_MPI
template<int D>
void ISend_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  
    println(10,MPI_Orb_rank<<" STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<" sending to "<<dest);
    int count = 1;
    int STreeMeta[count];

    timer.start();
    for (int ichunk = 0 ; ichunk < Nchunks ; ichunk++){
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	//set serialIx of the unused nodes to -1
	int ishift = ichunk*STree->maxNodesPerChunk;
	for (int i = 0; i < STree->maxNodesPerChunk; i++){
	    if (STree->nodeStackStatus[ishift+i] != 1)(STree->nodeChunks[ichunk])[i].setSerialIx(-1);
	}
	MPI_Isend(STree->nodeChunks[ichunk], count, MPI_BYTE, dest, tag+1+ichunk, comm, &request);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if (not STree->isShared or MPI_SH_group_rank != dest/MPI_SH_size){
	    MPI_Isend(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, dest, tag+ichunk*1000, comm, &request);}
    }
    timer.stop();
    println(10, " time send     " << timer); 
}
#endif

#ifdef HAVE_MPI
template<int D>
void Rcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  
    println(10, MPI_Orb_rank<<" STree  at "<<STree<<" receiving  "<<Nchunks<<" chunks from "<<source);
    int count = 1;
    int STreeMeta[count];

    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
	if(ichunk<STree->nodeChunks.size()){
	    STree->sNodes = STree->nodeChunks[ichunk];
	}else{
	    double *sNodesCoeff = new double[STree->sizeNodeCoeff*STree->maxNodesPerChunk];
	    STree->nodeCoeffChunks.push_back(sNodesCoeff);
	    STree->sNodes = (ProjectedNode<D>*) new char[STree->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	    STree->nodeChunks.push_back(STree->sNodes);
	}      
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	MPI_Recv(STree->nodeChunks[ichunk], count, MPI_BYTE, source, tag+1+ichunk, comm, &status);
	println(10, MPI_Orb_rank<<" received  "<<count/1024<<" kB with ProjectedNodes from "<<source);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if( not STree->isShared or MPI_SH_group_rank != source/MPI_SH_size){
	    MPI_Recv(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, source, tag+ichunk*1000, comm, &status);}
	println(10, " received  "<<count<<" coefficients from "<<source);
    }
    timer.stop();
    println(10, " time receive  " << timer);
    timer.start();
    STree->rewritePointers(Nchunks);
    timer.stop();
    println(10, " time rewrite pointers  " << timer);

}
#endif

#ifdef HAVE_MPI
template<int D>
void IRcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm){
    MPI_Status status;
    Timer timer;
    SerialFunctionTree<D>* STree = Tree->getSerialFunctionTree();
  

    println(10, MPI_Orb_rank<<" STree  at "<<STree<<" number of nodes = "<<STree->nNodes<<" receiving from "<<source);
    int count = 1;
    int STreeMeta[count];

    MPI_Request request=MPI_REQUEST_NULL;

    timer.start();
    for(int ichunk = 0 ; ichunk <Nchunks ; ichunk++){
	if(ichunk<STree->nodeChunks.size()){
	    STree->sNodes = STree->nodeChunks[ichunk];
	}else{
	    double *sNodesCoeff = new double[STree->sizeNodeCoeff*STree->maxNodesPerChunk];
	    STree->nodeCoeffChunks.push_back(sNodesCoeff);
	    STree->sNodes = (ProjectedNode<D>*) new char[STree->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	    STree->nodeChunks.push_back(STree->sNodes);
	}      
	count=STree->maxNodesPerChunk*sizeof(ProjectedNode<D>);
	MPI_Irecv(STree->nodeChunks[ichunk], count, MPI_BYTE, source, tag+1+ichunk, comm, &request);
	println(10, MPI_Orb_rank<<" received  "<<count/1024<<" kB with ProjectedNodes from "<<source);
	count=STree->sizeNodeCoeff*STree->maxNodesPerChunk;
	if (not STree->isShared or MPI_SH_group_rank != source/MPI_SH_size){
	    MPI_Irecv(STree->nodeCoeffChunks[ichunk], count, MPI_DOUBLE, source, tag+ichunk*1000, comm, &request);}
	println(10, MPI_Orb_rank<<" received  "<<count<<" coefficients from "<<source);
    }


    timer.stop();
    println(10, " time receive  " << timer);
    timer.start();
    STree->rewritePointers(Nchunks);
    timer.stop();
    println(10, " time rewrite pointers  " << timer);

}
#endif


/** Assign work among MPI processes
 * For a NxN loop. Does assume symmetry
 * At each iteration data is processed, data is sent and data is received.
 * input: 
 * N, the size of one loop 
 * iter, the iteration 
 * output: 
 * which (i,j) pair to compute next
 * to whom send own data, from who receive data
 * NB: does not work for N even when N is not a multiple of MPI_Orb_size
 */
void Assign_NxN_sym(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter){
    int iter = 0;
    int NBlock=N/MPI_Orb_size;//number of subblocks
    int M = N%MPI_Orb_size; //number of lines and column left after subblocks are taken out
    int Maxiter=-1;
    int iter_base;
    int iBlock = 0;
    int jBlock = 0;

    if(N%2==0 and N%MPI_Orb_size!=0 )MSG_FATAL("Assign_NxN_sym not implemented for this case");
    for (iter = 0; iter < *MaxIter; iter++) {
	doi[iter] = -1;
	doj[iter] = -1;
	sendto[iter] = -1;
	sendorb[iter] = -1;
	rcvorb[iter] = -1;
    }

    for (int dBlock = 0; dBlock < (N+MPI_Orb_size-1)/MPI_Orb_size; dBlock++) {
	//treat diagonal block first
	int imax=MPI_Orb_size;//size of the diagonal block
	if(dBlock==NBlock)imax=M;
	int jmax=MPI_Orb_size;
	if(dBlock==NBlock)jmax=M;
	for (int i = 0; i < imax; i++) {
	    for (int j = 0; j < jmax; j++) {

		//used in a diagonal block: 1+MPI_Orb_size/2 
		//used in a block below diagonal: (1+MPI_Orb_size)*0.5
		//   if even MPI_Orb_rank: 1+MPI_Orb_size/2 if even block, and MPI_Orb_size/2 if odd block, 0.5*(1+MPI_Orb_size) on average, 
		//used in a block over diagonal: (1+MPI_Orb_size)*0.5-1
		//number of entire blocks in a column: N/MPI_Orb_size 
		//possibly used in lowest reduced block in column: M/2
		//used in all blocks in one column d (starting at 0): (1+(N/MPI_Orb_size)*(1+MPI_Orb_size))/2-d+M/2
		//used in all blocks up to column d, without d: d*(1+(N/MPI_Orb_size)*(1+MPI_Orb_size))/2-(d*(d-1))/2+d*M/2


		//(MPI_Orb_size+j-i)%MPI_Orb_size : for one block
		//((N+MPI_Orb_size-1)/MPI_Orb_size) : number of blocks in one column
		//-dBlock: The diagonal is not counted again in each block
		iter_base = (MPI_Orb_size+j-i)%MPI_Orb_size * ((N+MPI_Orb_size-1)/MPI_Orb_size)-dBlock;
		if(i==j)iter_base = 0;
		//the last block may be smaller
		if(dBlock==NBlock)iter_base = (imax+j-i)%M * ((N+MPI_Orb_size-1)/MPI_Orb_size);
		//Shift for all the columns already taken into account
		iter = iter_base + dBlock*((1+((N+MPI_Orb_size-1)/MPI_Orb_size)*(1+MPI_Orb_size))/2) - (dBlock*(dBlock-1))/2 + dBlock*M/2;

		if( (imax+j-i)%imax<=imax/2 ){
		    int i_glob = i + dBlock*MPI_Orb_size;
		    int j_glob = j + dBlock*MPI_Orb_size;

		    if(i==MPI_Orb_rank and ( imax%2!=0 or (imax+j-i)%imax<imax/2 or i<j )) {
			if(iter>Maxiter)Maxiter=iter;
			doi[iter]=i_glob;
			doj[iter]=j_glob;
			if(i!=j)rcvorb[iter]=j_glob;
			if(imax%2==0 and (imax+j-i)%imax==imax/2 and i<j){
			    sendto[iter]=j;
			    sendorb[iter]=i_glob;	    
			}
		    }else if(j==MPI_Orb_rank and ( imax%2!=0 or (imax+j-i)%imax<imax/2 or i<j )){
			sendto[iter]=i;
			sendorb[iter]=j_glob;
			if(imax%2==0 and (imax+j-i)%imax==imax/2 and i<j){
			    rcvorb[iter]=i;
			}
		    }

		    //all columns that can be generated from received data.	
		    if((i==MPI_Orb_rank and (i<=j or (MPI_Orb_size%2==0 and (imax+j-i)%imax==imax/2) ))){
			for (int icol = i; icol < N; icol+=MPI_Orb_size) {
			    iBlock=icol/MPI_Orb_size;
			    if((i!=j or iBlock>dBlock) and  (MPI_Orb_size%2!=0 or ((imax+j-i)%imax<imax/2) or ((iBlock+dBlock)%2!=0 and i>j) or ((iBlock+dBlock)%2==0 and i<j)  )){
				if(iBlock!=dBlock){
				    int sh=1;//diagonal block is the first
				    if(iBlock>dBlock)sh=0;//diagonal block accounted for
				    if(i==j)sh=-dBlock;//only below diagonal are counted
				    if(MPI_Orb_size%2==0 and (imax+j-i)%imax==imax/2)sh=-(iBlock+1)/2;//only half are present
				    if(MPI_Orb_size%2==0 and (imax+j-i)%imax==imax/2 and (dBlock%2==0 and iBlock==0))sh=1;//same iter as diagonalblock
				    if(MPI_Orb_size%2==0 and (imax+j-i)%imax==imax/2 and (dBlock%2!=0 and iBlock==1))sh=1;//iBlock=0 same iter as diagonal block, start next iter
				    int i_glob = i + iBlock*MPI_Orb_size;
				    int j_glob = j + dBlock*MPI_Orb_size;
				    if((imax+j-i)%imax<=imax/2 ){
					doi[iter+iBlock+sh] = i_glob;
					doj[iter+iBlock+sh] = j_glob;}
				    //rcvorb[iter+iBlock+sh]=-1;//indicates "use same as before"
				    //sendorb[iter+iBlock+sh]=-1;//indicates "do not send anything"
				    //sendto[iter+iBlock+sh]=-1;//indicates "do not send anything"
				    if(iter+iBlock+sh>Maxiter)Maxiter=iter+iBlock+sh;
				}
			    }
			}
		    }else if(i==MPI_Orb_rank){
			//part of blocks below blockdiagonal
			//for (jBlock = 0; jBlock < NBlock; jBlock++) {
			for (int jrow = j; jrow < N; jrow+=MPI_Orb_size) {
			    jBlock=jrow/MPI_Orb_size;
			    if(jBlock!=dBlock  and (MPI_Orb_size%2 or (imax+j-i)%imax<(imax)/2 or (jBlock%2==0))){
				int sh=1;
				if(jBlock>dBlock)sh=0;//not very elegant
				int i_glob = i + dBlock*MPI_Orb_size;
				int j_glob = j + jBlock*MPI_Orb_size;
				if((imax+j-i)%imax<(imax+1)/2){
				    doi[iter+jBlock+sh] = i_glob;
				    doj[iter+jBlock+sh] = j_glob;}
				//rcvorb[iter+jBlock+sh]=-1;//indicates "use same as before"
				//sendorb[iter+jBlock+sh]=-1;//indicates "do not send anything"
				//sendto[iter+jBlock+sh]=-1;//indicates "do not send anything"
				if(iter+jBlock+sh>Maxiter)Maxiter=iter+jBlock+sh;
			    }
			}
		    }
		}
	    }
	}
    }

    *MaxIter = Maxiter;
}


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
void Assign_NxN(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter){
    int iter = 0;
    int NBlock=N/MPI_Orb_size;//number of subblocks
    int M = N%MPI_Orb_size; //number of lines and column left after subblocks are taken out
    int Maxiter = -1;
    int iter_base;
    int iBlock = 0;
    int jBlock = 0;

    for (iter = 0; iter < *MaxIter; iter++) {
	doi[iter] = -1;
	doj[iter] = -1;
	sendto[iter] = -1;
	sendorb[iter] = -1;
	rcvorb[iter] = -1;
    }

    for (int dBlock = 0; dBlock < (N+MPI_Orb_size-1)/MPI_Orb_size; dBlock++) {

	//1)  //diagonal block
	int imax=MPI_Orb_size;
	if(dBlock==NBlock)imax=M;
	int jmax=MPI_Orb_size;
	if(dBlock==NBlock)jmax=M;
	for (int i = 0; i < imax; i++) {
	    for (int j = 0; j < jmax; j++) {
		iter_base = (i+j)%MPI_Orb_size * ((N+MPI_Orb_size-1)/MPI_Orb_size);
		if(dBlock==NBlock)iter_base = (i+j)%M * ((N+MPI_Orb_size-1)/MPI_Orb_size);
		iter = iter_base + dBlock*MPI_Orb_size*((N+MPI_Orb_size-1)/MPI_Orb_size);
		int i_glob = i + dBlock*MPI_Orb_size;
		int j_glob = j + dBlock*MPI_Orb_size;
		if(i==MPI_Orb_rank){
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
		if(i<=j and i==MPI_Orb_rank){
		    //do the corresponding column 
		    //for (iBlock = 0; iBlock < NBlock; iBlock++) {
		    for (int icol = i; icol < N; icol+=MPI_Orb_size) {
			iBlock=icol/MPI_Orb_size;
			if(iBlock!=dBlock){
			    int sh=1;
			    if(iBlock>dBlock)sh=0;//not very elegant
			    int i_glob = i + iBlock*MPI_Orb_size;
			    int j_glob = j + dBlock*MPI_Orb_size;
			    doi[iter+iBlock+sh] = i_glob;
			    doj[iter+iBlock+sh] = j_glob;
			    rcvorb[iter+iBlock+sh]=-1;//indicates "use same as before"
			    sendorb[iter+iBlock+sh]=-1;//indicates "do not send anything"
			    sendto[iter+iBlock+sh]=-1;//indicates "do not send anything"
			    if(iter+iBlock+sh>Maxiter)Maxiter=iter+iBlock+sh;
			}
		    }
		}else if(i==MPI_Orb_rank){
		    //do the corresponding row
		    //for (jBlock = 0; jBlock < NBlock; jBlock++) {
		    for (int jrow = j; jrow < N; jrow+=MPI_Orb_size) {
			jBlock=jrow/MPI_Orb_size;
			if(jBlock!=dBlock){
			    int sh=1;
			    if(jBlock>dBlock)sh=0;//not very elegant
			    int i_glob = i + dBlock*MPI_Orb_size;
			    int j_glob = j + jBlock*MPI_Orb_size;
	    
			    doi[iter+jBlock+sh] = i_glob;
			    doj[iter+jBlock+sh] = j_glob;
			    rcvorb[iter+jBlock+sh]=-1;//indicates "use same as before"
			    sendorb[iter+jBlock+sh]=-1;//indicates "do not send anything"
			    sendto[iter+jBlock+sh]=-1;//indicates "do not send anything"
			    if(iter+jBlock+sh>Maxiter)Maxiter=iter+jBlock+sh;
			}
		    }
		}
	    }
	}
    }
    int NiterMax = ((N+MPI_Orb_size-1)/MPI_Orb_size)*((N+MPI_Orb_size-1)/MPI_Orb_size)*MPI_Orb_size;
    if(Maxiter>NiterMax)cout<<"CHECK Assign_NxN "<<endl;
    *MaxIter = Maxiter;
}

/** Define all the different MPI groups
 */
void define_MPI_groups(){

#ifdef HAVE_MPI
    int MPI_worldsize;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_worldsize);
    println(10,"MPI world size: "<< MPI_worldsize ); 
    int MPI_worldrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_worldrank);

    //divide the world into groups
    //each group has its own group communicator definition
 
    //split world into groups that can share memory
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &MPI_Comm_SH);

    MPI_Comm_rank(MPI_Comm_SH, &MPI_SH_rank);//rank in the SH group
    MPI_Comm_size(MPI_Comm_SH, &MPI_SH_size);
    println(10,"MPI SH size: "<< MPI_SH_size ); 

    //define a rank of the group. NB: we assume for now that all groups have same size
    MPI_Comm MPI_Comm_SH_group;
    MPI_Comm_split(MPI_COMM_WORLD, MPI_SH_rank, MPI_worldrank, &MPI_Comm_SH_group);//MPI_SH_rank is color (same color->in same group),  MPI_worldrank is key (orders rank within the groups)

    //we define a new orbital rank, so that the orbitals within a shared memory group, have consecutive ranks.
    MPI_Comm_rank(MPI_Comm_SH_group, &MPI_SH_group_rank);//rank in the SH group
    MPI_Comm_size(MPI_Comm_SH_group, &MPI_SH_group_size);
    if(MPI_worldsize != MPI_SH_size*MPI_SH_group_size )MSG_FATAL("all MPI SH groups must be of same size");
    println(10,"MPI SH group size: "<< MPI_SH_group_size ); 
    MPI_Orb_rank = MPI_SH_rank + MPI_SH_group_rank*MPI_SH_size ;

    MPI_Comm_split(MPI_COMM_WORLD, 0, MPI_Orb_rank, &MPI_Comm_Orb);//0 is color (same color->in same group), MPI_Orb_rank is key (orders rank in the group)
    int MPI_Orb_rank;
    MPI_Comm_rank(MPI_Comm_Orb, &MPI_Orb_rank);//should be the same
    MPI_Comm_size(MPI_Comm_Orb, &MPI_Orb_size);
    println(10,"orbital group size: "<< MPI_Orb_size); 
    //    cout<<"world rank="<<MPI_worldrank<<"  Orbital rank="<<MPI_Orb_rank<<"  Shared mem rank="<<MPI_SH_rank<<"  Shared memory group rank="<<MPI_SH_group_rank<<endl;

#endif
  
}

#ifdef HAVE_MPI
/** Define MPI groups that share the same compute-node i.e. the same memory and array with shared memory between
 * sh_size is size in MBytes
 * returns a pointer to the shared memory*/
void Share_memory(int sh_size, double * &d_ptr, MPI_Win &win){

    //MPI_Aint types are used for adresses (can be larger than int).
    MPI_Aint size = (MPI_SH_rank==0) ? 1024*1024*sh_size : 0;//rank 0 defines length of segment 
    //MPI_Win win = MPI_WIN_NULL;
    int disp_unit = 16;//in order for the compiler to keep aligned. 
    //size is in bytes
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, MPI_Comm_SH, &d_ptr, &win);

    MPI_Win_fence(0, win);//wait until finished

    MPI_Aint qsize = 0;
    int qdisp = 0;
    int * qbase = NULL;
    MPI_Win_shared_query(win, 0, &qsize, &qdisp, &d_ptr);
    //printf("me = %d, allocated: size=%zu, bytes at = %p\n", shrank, qsize, d_ptr);
    MPI_Win_fence(0, win);
    //if(MPI_Orb_rank==0)d_ptr[18]=2.34;
    //printf("shared: me=%d,  val=%lf\n", shrank, d_ptr[18]);

    //MPI_Win_free(&win);//to include in destructor

} 
#endif

#ifdef HAVE_MPI
template void Rcv_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void Rcv_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void Rcv_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void Send_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int dest, int tag, MPI_Comm comm);
template void Send_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int dest, int tag, MPI_Comm comm);
template void Send_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int dest, int tag, MPI_Comm comm);
template void IRcv_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void IRcv_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void IRcv_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int source, int tag, MPI_Comm comm);
template void ISend_SerialTree<1>(FunctionTree<1>* STree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
template void ISend_SerialTree<2>(FunctionTree<2>* STree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
template void ISend_SerialTree<3>(FunctionTree<3>* STree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
#endif

