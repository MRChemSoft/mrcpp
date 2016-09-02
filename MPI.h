#ifndef MPI_H_
#define MPI_H_

#include "SerialTree.h"

extern int MPI_rank;
extern int MPI_size;

void define_groups();
void MPI_Initializations();
//template<int D>
void SendRcv_SerialTree(SerialTree<1>* STree, int source, int dest, int tag, MPI_Comm comm);
void SendRcv_SerialTree(SerialTree<2>* STree, int source, int dest, int tag, MPI_Comm comm);
void SendRcv_SerialTree(SerialTree<3>* STree, int source, int dest, int tag, MPI_Comm comm);
//void Share_memory(MPI_Comm ncomm, MPI_Comm ncomm_sh, int sh_size, double * d_ptr);

#ifdef HAVE_MPI


#else

namespace mpi {
	struct communicator {
		int rank() const { return 0; }
		int size() const { return 1; }
		void barrier() const { }
	};
	struct environment {
		environment(int argc, char **argv) { }
	};
	struct timer {
		void restart() { }
		double elapsed() { return 0.0; }
	};
	typedef int request;
}
#endif
#endif /* MPI_H_*/
