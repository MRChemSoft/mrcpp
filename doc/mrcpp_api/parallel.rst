--------
Parallel
--------

The core features of MRCPP are parallelized using a shared memory model *only*
(OpenMP). This means that there is *no intrinsic* MPI parallelization (e.i. no
data distribution across machines) *within* the library routines. However, the
code comes with a small set of features that facilitate MPI work and data
distribution in the host program, in the sense that *entire* ``FunctionTree``
objects can be located on different machines and communicated between them.
Also, a ``FunctionTree`` can be *shared* between several MPI processes
that are located on the *same* machine. This means that several processes have
read access to the same ``FunctionTree``, thus reducing the memory footprint,
as well as the need for communication.

The MPI features are available by including:

.. code-block:: cpp

    #include "MRCPP/Parallel"


The host program
----------------

In order to utilize the MPI features of MRCPP, the MPI instance must be
initialized (and finalized) by the host program, as usual:

.. code-block:: cpp

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get MPI world size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get MPI world rank

    MPI_Finalize();


For the shared memory features we must make sure that the ranks within a
communicator is actually located on the same machine. When running on
distributed architectures this can be achieved by creating separate
communicators for each physical machine, e.g. to split *MPI_COMM_WORLD*
into a new communicator group called *MPI_COMM_SHARED* that share the
same physical memory space:

.. code-block:: cpp

    // Initialize a new communicator called MPI_COMM_SHARE
    MPI_Comm MPI_COMM_SHARE;

    // Split MPI_COMM_WORLD into sub groups and assign to MPI_COMM_SHARE
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &MPI_COMM_SHARE);


Note that the main purpose of the shared memory feature of MRCPP is to avoid
memory duplication and reduce the memory footprint, it will **not**
automatically provide any work sharing parallelization for the construction of
the shared ``FunctionTree``.


Blocking communication
----------------------

A blocking send/receive means that the function call does not return until the
communication is completed. This is a simple and safe option, but can lead to
significant overhead if the communicating MPI prosesses are not synchronized.

send_tree
  Send a complete ``FunctionTree`` to a given MPI rank using blocking
  communication.

recv_tree
  Receive a complete ``FunctionTree`` from a given MPI rank using blocking
  communication.


.. code-block:: cpp

    mrcpp::FunctionTree<3> tree(MRA);

    
    // At this point tree is uninitialized on both rank 0 and 1


    // Only rank 0 projects the function
    if (rank == 0) mrcpp::project(prec, tree, func);


    // At this point tree is projected on rank 0 but still uninitialized on rank 1


    // Sending tree from rank 0 to rank 1
    int tag = 111111; // Unique tag for each communication
    int src=0, dst=1; // Source and destination ranks
    if (rank == src) mrcpp::send_tree(tree, dst, tag, MPI_COMM_WORLD);
    if (rank == dst) mrcpp::revc_tree(tree, src, tag, MPI_COMM_WORLD);


    // At this point tree is projected on both rank 0 and 1


    // Rank 0 clear the tree
    if (rank == 0) mrcpp::clear(tree);


    // At this point tree is uninitialized on rank 0 but still projected on rank 1


Non-blocking communication
--------------------------

A non-blocking send means that the sending MPI process will return from the
function call immediately and carry on computing, even if the receiving part
is not ready yet.

isend_tree
  Send a complete ``FunctionTree`` from a given MPI rank using non-blocking
  communication. Should be combined with *MPI_Wait* to verify that the data
  has been received.

There is no corresponding non-blocking receive. The reason for this is that the
communication happens in several steps, each depending on the previous. This
means that *MPI_Wait* statements would be required in between, and the purpose
of non-blocking is lost. The following example shows the importance of the wait
statement: without it the send could potentially be performed *after* the tree
has been cleared locally on rank 0.

.. code-block:: cpp

    mrcpp::FunctionTree<3> tree(MRA);


    // At this point tree is uninitialized on both rank 0 and 1


    // Only rank 0 projects the function
    if (rank == 0) mrcpp::project(prec, tree, func);


    // At this point tree is projected on rank 0 but still uninitialized on rank 1


    // Sending tree from rank 0 to rank 1
    int tag = 222222; // Unique tag for each communication
    int src=0, dst=1; // Source and destination ranks
    MPI_Request request = MPI_REQUEST_NULL;
    if (rank == src) mrcpp::isend_tree(tree, dst, tag, MPI_COMM_WORLD, &request);
    if (rank == dst) mrcpp::revc_tree(tree, src, tag, MPI_COMM_WORLD);

    // Source rank waits until the send is complete
    MPI_Wait(&request, MPI_STATUS_NULL);


    // At this point tree is projected on both rank 0 and 1


    // Rank 0 clear the tree
    if (rank == 0) mrcpp::clear(tree);


    // At this point tree is uninitialized on rank 0 but still projected on rank 1


Shared memory
-------------

The sharing of a ``FunctionTree`` happends in three steps: first a
``SharedMemory`` object is initialized with the appropriate shared memory
communicator; then this object is used in the ``FunctionTree`` constructor;
finally, *after* the ``FunctionTree`` has been properly computed, a call must
be made to the ``share_tree`` function. The reason for the last function call
is that the internal memory pointers needs to be updated *locally* on each MPI
process whenever the shared memory window has been updated.


share_tree
  Share a ``FunctionTree`` among MPI processes that share the same physical
  memory. This function should be called every time a *shared* ``FunctionTree``
  is updated, in order to update the local memory of each MPI process.

.. code-block:: cpp

    // Get rank within the shared group
    int rank;
    MPI_Comm_rank(MPI_COMM_SHARE, &rank);

    // Define master and worker ranks
    int master = 0;
    int worker = 1;

    // The tree will be shared within the given communicator
    int mem_size = 1000; //MB
    mrcpp::SharedMemory shared_mem(MPI_COMM_SHARE, mem_size);
    mrcpp::FunctionTree<3> tree(MRA, shared_mem);

    // Master rank projects the tree
    if (rank == master) mrcpp::project(prec, tree, func);

    // When a shared function is updated, it must be re-shared
    int tag = 333333; // Unique tag for each communication
    mrcpp::share_tree(tree, master, tag, MPI_COMM_SHARE); 

    // Other ranks within the shared group can update the tree
    if (rank == worker) tree.rescale(2.0);

    // When a shared function is updated, it must be re-shared
    mrcpp::share_tree(tree, worker, tag, MPI_COMM_SHARE); 


