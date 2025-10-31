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

#pragma once
/**
 * @file
 * @brief Distributed “Bank” service for sharing functions and raw data across MPI ranks.
 *
 * This header declares a minimal runtime that lets multiple MPI ranks exchange:
 * - Multiresolution functions (`CompFunction<3>`) and
 * - Raw numeric buffers (`double` / `ComplexDouble`)
 *
 * The service is organized as a central **Bank** that maintains per-client
 * accounts. Each rank interacts with the Bank through a lightweight RAII
 * client, **BankAccount**. A simple **TaskManager** piggybacks on the same
 * infrastructure to distribute integer-indexed tasks and collect “ready”
 * notifications.
 *
 * @par High-level design
 * - The **Bank** lives on one or more designated MPI ranks (see `mpi::is_bank`
 *   in the runtime). Non-bank ranks act as clients.
 * - Clients open an **account** and then `put_*` or `get_*` by integer IDs or
 *   by `NodeIndex<3>` keys. The Bank tracks sizes in kB for accounting.
 * - The **TaskManager** provides a tiny work-queue: clients request the next
 *   task, mark tasks as ready, and optionally consume ready items.
 *
 * @note The concrete message-passing, blocking semantics, and memory
 *       ownership rules are implemented in the Bank source (MPI-based).
 *       This header documents intent and call contracts at a high level.
 */

#include "CompFunction.h"
#include "parallel.h"
#include "trees/NodeIndex.h"

namespace mrcpp {

using namespace mpi;

/**
 * @brief A deposited item stored by the Bank.
 *
 * A deposit can represent either a multiresolution function (`orb`) or a
 * raw data buffer (`data`). Exactly one of them is expected to be active
 * for a given deposit.
 */
struct deposit {
    /** Pointer to a deposited function (3D component function). */
    CompFunction<3> *orb = nullptr;
    /** Pointer to a deposited plain data buffer (contiguous). */
    double *data = nullptr; // for pure data arrays
    /** True if this deposit contains a raw data buffer in @ref data. */
    bool hasdata = false;
    /** Size (number of elements) for @ref data when @ref hasdata is true. */
    int datasize = 0;
    /** Application-defined identifier used to name and retrieve this deposit. */
    int id = -1; // to identify what is deposited
    /** MPI rank that originally deposited the item. */
    int source = 0;  // mpi rank from the source of the data
};

/**
 * @brief Queue bookkeeping for task-ready notifications.
 *
 * Associates a queue identifier with a list of client ranks that registered
 * interest or contributed ready items.
 */
struct queue_struct {
    /** Queue identifier (application-defined). */
    int id = 0;
    /** Ranks that have entries or are waiting on this queue. */
    std::vector<int> clients;
};

/**
 * @brief Command codes exchanged between clients and the Bank/TaskManager.
 *
 * These enumerators are used as operation selectors in MPI messages. Listed
 * values are stable and intentionally explicit to simplify debugging.
 */
enum {
    CLOSE_BANK,              ///< 0 — Shut down the Bank service.
    CLEAR_BANK,              ///< 1 — Remove all accounts and deposits.
    NEW_ACCOUNT,             ///< 2 — Open a new client account.
    CLOSE_ACCOUNT,           ///< 3 — Close (delete) an existing account.
    GET_ORBITAL,             ///< 4 — Retrieve an orbital (internal legacy op).
    GET_FUNCTION_AND_WAIT,   ///< 5 — Blocking fetch of a function until available.
    GET_FUNCTION_AND_DELETE, ///< 6 — Fetch and erase a function.
    SAVE_ORBITAL,            ///< 7 — Store an orbital (internal legacy op).
    GET_FUNCTION,            ///< 8 — Non-blocking fetch of a function if available.
    SAVE_FUNCTION,           ///< 9 — Store a function.
    GET_DATA,                ///< 10 — Fetch a raw data buffer.
    SAVE_DATA,               ///< 11 — Store a raw data buffer.
    SAVE_NODEDATA,           ///< 12 — Store node-scoped raw data.
    GET_NODEDATA,            ///< 13 — Fetch node-scoped raw data.
    GET_NODEBLOCK,           ///< 14 — Fetch a contiguous block for a node id across IDs.
    GET_ORBBLOCK,            ///< 15 — Fetch a contiguous block for an orbital id across nodes.
    CLEAR_BLOCKS,            ///< 16 — Clear block caches/aggregations.
    GET_MAXTOTDATA,          ///< 17 — Query max total stored size (kB).
    GET_TOTDATA,             ///< 18 — Query per-account total sizes (kB).
    INIT_TASKS,              ///< 19 — Initialize TaskManager with N tasks.
    GET_NEXTTASK,            ///< 20 — Acquire the next task index.
    PUT_READYTASK,           ///< 21 — Mark (i,j) as ready.
    DEL_READYTASK,           ///< 22 — Remove a ready marker (i,j).
    GET_READYTASK,           ///< 23 — Get ready list for i (keep).
    GET_READYTASK_DEL,       ///< 24 — Get ready list for i and consume.
};

/**
 * @brief Central repository for distributed function/data sharing and task queues.
 *
 * The Bank owns per-client accounts, holds deposits, and tracks memory usage
 * in kB. Only designated Bank ranks instantiate and `open()` the service;
 * clients interact via @ref BankAccount and @ref TaskManager on worker ranks.
 *
 * @par Thread-safety
 * Bank methods are orchestrated via MPI; within a single rank, methods are not
 * inherently thread-safe unless otherwise guarded at the call site.
 */
class Bank {
public:
    /** @brief Construct an unopened Bank instance on a Bank rank. */
    Bank() = default;

    /** @brief Destructor. Ensures resources are released if not already closed. */
    ~Bank();

    /**
     * @brief Start the Bank service (receive loop, state init).
     *
     * Must be called on the Bank rank(s). After `open()`, the service listens
     * for client commands and manages account state.
     */
    void open();

    /**
     * @brief Stop the Bank service and release all resources.
     *
     * Closes all accounts and clears all deposits.
     */
    void close();

    /**
     * @brief Maximum total footprint observed so far, in kB.
     * @return Peak cumulative size (across all accounts) since `open()`.
     */
    int get_maxtotalsize();

    /**
     * @brief Current total sizes per account, in kB.
     * @return A vector of sizes aligned with internal account ordering.
     */
    std::vector<int> get_totalsize();

private:
    friend class BankAccount;
    friend class TaskManager;

    // ---- Account control (called by clients through Bank's command loop) ----

    /**
     * @brief Create a new account for client rank @p iclient.
     * @param iclient Rank creating the account (logical owner).
     * @param comm    Communicator the client uses to reach the Bank.
     * @return Integer account identifier (>0 on success).
     */
    int openAccount(int iclient, MPI_Comm comm);

    /**
     * @brief Clear and reinitialize an existing account.
     *
     * Equivalent to closing and reopening the account, preserving the account id.
     *
     * @param account Account id to clear.
     * @param iclient Requesting client rank.
     * @param comm    Client communicator.
     * @return 0 on success, negative on error.
     */
    int clearAccount(int account, int iclient, MPI_Comm comm);

    /**
     * @brief Permanently remove an account and all of its deposits.
     * @param account_id Account identifier.
     */
    void closeAccount(int account_id);

    // ---- Task manager control (internal) ----

    /**
     * @brief Initialize task bookkeeping for @p ntasks items.
     * @param ntasks  Number of tasks available (0..ntasks-1).
     * @param iclient Requesting rank.
     * @param comm    Client communicator.
     * @return Account id of the task manager context.
     */
    int openTaskManager(int ntasks, int iclient, MPI_Comm comm);

    /**
     * @brief Close and remove a TaskManager context.
     * @param account_id Associated account id.
     */
    void closeTaskManager(int account_id);

    // ---- Internal utilities ----

    /** @brief Remove all accounts and deposits (global reset). */
    void clear_bank();

    /**
     * @brief Remove a single account and its content.
     * @param account Account id to erase.
     */
    void remove_account(int account);

    // ---- Accounting & indices ----
    long long totcurrentsize = 0ll;                     ///< Sum of all account sizes (kB).
    std::vector<int> accounts;                          ///< Active account ids.

    /** Map: account id → vector of deposits. */
    std::map<int, std::vector<deposit> *> get_deposits;

    /** Map: account id → (item id → index in deposits vector). */
    std::map<int, std::map<int, int> *> get_id2ix;

    /** Map: account id → (queue id → index in queue vector). */
    std::map<int, std::map<int, int> *> get_id2qu;

    /** Map: account id → queue collection (task-ready queues). */
    std::map<int, std::vector<queue_struct> *> get_queue;

    /** Map: account id → (i → vector of j ready items). */
    std::map<int, std::map<int, std::vector<int>> *> get_readytasks;

    /** Map: account id → current size in kB (without container overhead). */
    std::map<int, long long> currentsize;

    /** Peak total size (kB) observed since last reset. */
    long long maxsize = 0;
};

/**
 * @brief RAII client-side view of a Bank account.
 *
 * A `BankAccount` encapsulates a live account and offers typed methods to
 * deposit and retrieve functions or raw buffers. Most methods are thin
 * request wrappers; the Bank performs the actual storage.
 *
 * @note By default, rank and communicator are taken from the MPI worker
 *       context (`mpi::wrk_rank`, `mpi::comm_wrk`).
 *
 * @par Ownership & lifetime
 * Returned raw pointers (e.g., from `get_orbblock`) typically reference
 * storage owned by the Bank. Callers should copy data if it must outlive
 * subsequent Bank interactions. See the implementation for exact details.
 */
class BankAccount {
public:
    /**
     * @brief Open a new account for @p iclient on communicator @p comm.
     * @param iclient Client rank that owns the account (default: @ref mpi::wrk_rank).
     * @param comm    Communicator used to contact the Bank (default: @ref mpi::comm_wrk).
     */
    BankAccount(int iclient = wrk_rank, MPI_Comm comm = comm_wrk);

    /** @brief Close the account and release any client-side resources. */
    ~BankAccount();

    /** @brief Bank-assigned account identifier (≥0 when open). */
    int account_id = -1;

    /**
     * @brief Clear and reinitialize this account.
     * @param i    Client rank issuing the request.
     * @param comm Client communicator.
     */
    void clear(int i = wrk_rank, MPI_Comm comm = comm_wrk);

    // --- Function storage/retrieval ---

    /**
     * @brief Fetch a function by @p id and delete it on the Bank.
     * @param id   Application-level identifier of the function.
     * @param orb  Output destination; resized/assigned by the Bank.
     * @return 0 on success, negative on error.
     */
    int get_func_del(int id, CompFunction<3> &orb);

    /**
     * @brief Deposit a function under identifier @p id.
     * @param id   Application-level identifier.
     * @param func Function object to store (copied/serialized by Bank).
     * @return 0 on success, negative on error.
     */
    int put_func(int id, CompFunction<3> &func);

    /**
     * @brief Fetch a function by @p id.
     * @param id   Application-level identifier.
     * @param func Output destination; resized/assigned by the Bank.
     * @param wait If nonzero, block until available; otherwise return immediately if missing.
     * @return 0 on success; negative on error; positive (e.g. 1) if not found and @p wait==0.
     */
    int get_func(int id, CompFunction<3> &func, int wait = 0);

    // --- Raw data buffers by plain id ---

    /**
     * @brief Deposit a real-valued buffer.
     * @param id   Application-level identifier.
     * @param size Number of elements in @p data.
     * @param data Pointer to contiguous buffer (copied by the Bank).
     * @return 0 on success, negative on error.
     */
    int put_data(int id, int size, double *data);

    /**
     * @brief Deposit a complex-valued buffer.
     * @copydetails put_data(int,int,double*)
     */
    int put_data(int id, int size, ComplexDouble *data);

    /**
     * @brief Retrieve a real-valued buffer by @p id.
     * @param id   Identifier previously used with @ref put_data.
     * @param size Expected number of elements; used for validation.
     * @param data Destination buffer; must have room for @p size elements.
     * @return 0 on success; negative on error; positive if not found.
     */
    int get_data(int id, int size, double *data);

    /**
     * @brief Retrieve a complex-valued buffer by @p id.
     * @copydetails get_data(int,int,double*)
     */
    int get_data(int id, int size, ComplexDouble *data);

    // --- Raw data scoped by node index (spatial addressing) ---

    /**
     * @brief Deposit real-valued data associated with a node index.
     * @param nIdx Spatial node key.
     * @param size Number of elements.
     * @param data Buffer pointer (copied by the Bank).
     * @return 0 on success, negative on error.
     */
    int put_data(NodeIndex<3> nIdx, int size, double *data);

    /**
     * @brief Deposit complex-valued data for a node index.
     * @copydetails put_data(NodeIndex<3>,int,double*)
     */
    int put_data(NodeIndex<3> nIdx, int size, ComplexDouble *data);

    /**
     * @brief Retrieve real-valued data for a node index.
     * @param nIdx Node key.
     * @param size Expected number of elements.
     * @param data Output buffer.
     * @return 0 on success; negative on error; positive if not found.
     */
    int get_data(NodeIndex<3> nIdx, int size, double *data);

    /**
     * @brief Retrieve complex-valued data for a node index.
     * @copydetails get_data(NodeIndex<3>,int,double*)
     */
    int get_data(NodeIndex<3> nIdx, int size, ComplexDouble *data);

    // --- Node-scoped data grouped under an object id (e.g., orbital id) ---

    /**
     * @brief Deposit real-valued data for a specific node @p nodeid within object @p id.
     * @param id     Object (e.g., orbital) identifier.
     * @param nodeid Node identifier within the object.
     * @param size   Number of elements.
     * @param data   Buffer pointer (copied by the Bank).
     * @return 0 on success, negative on error.
     */
    int put_nodedata(int id, int nodeid, int size, double *data);

    /**
     * @brief Deposit complex-valued data for a node within object @p id.
     * @copydetails put_nodedata(int,int,int,double*)
     */
    int put_nodedata(int id, int nodeid, int size, ComplexDouble *data);

    /**
     * @brief Retrieve real-valued data for (@p id, @p nodeid).
     * @param id      Object identifier.
     * @param nodeid  Node identifier.
     * @param size    Expected element count.
     * @param data    Output buffer.
     * @param idVec   (Out) List of object ids actually present in the block, if aggregated.
     * @return 0 on success; negative on error; positive if not found.
     */
    int get_nodedata(int id, int nodeid, int size, double *data, std::vector<int> &idVec);

    /**
     * @brief Retrieve complex-valued data for (@p id, @p nodeid).
     * @copydetails get_nodedata(int,int,int,double*,std::vector<int>&)
     */
    int get_nodedata(int id, int nodeid, int size, ComplexDouble *data, std::vector<int> &idVec);

    // --- Block retrieval helpers ---

    /**
     * @brief Retrieve a contiguous block of all real node data for @p nodeid across ids.
     * @param nodeid  Node identifier to gather.
     * @param data    (Out) Pointer to contiguous storage; copy data before next call.
     * @param idVec   (Out) List of ids participating in the block.
     * @return Number of elements in @p data on success; negative on error.
     */
    int get_nodeblock(int nodeid, double *data, std::vector<int> &idVec);

    /**
     * @brief Retrieve a contiguous block of all complex node data for @p nodeid across ids.
     * @copydetails get_nodeblock(int,double*,std::vector<int>&)
     */
    int get_nodeblock(int nodeid, ComplexDouble *data, std::vector<int> &idVec);

    /**
     * @brief Retrieve all real-valued node data for an orbital id into a single contiguous block.
     * @param orbid      Orbital (object) id.
     * @param data       (Out) Pointer reference to contiguous storage.
     * @param nodeidVec  (Out) Node ids represented in the block.
     * @param bankstart  Starting index/offset within the Bank’s internal storage.
     * @return Number of elements in the returned block; negative on error.
     *
     * @note Copy out the data if it must persist beyond this call or subsequent Bank calls.
     */
    int get_orbblock(int orbid, double *&data, std::vector<int> &nodeidVec, int bankstart);

    /**
     * @brief Retrieve all complex-valued node data for an orbital id into a contiguous block.
     * @copydetails get_orbblock(int,double*&,std::vector<int>&,int)
     */
    int get_orbblock(int orbid, ComplexDouble *&data, std::vector<int> &nodeidVec, int bankstart);
};

/**
 * @brief Minimal distributed task queue associated with a Bank account.
 *
 * The TaskManager assigns task indices in [0, @ref n_tasks). Clients can:
 * - Request the next task (`next_task()`),
 * - Mark specific items as ready (`put_readytask(i,j)` / `del_readytask(i,j)`),
 * - Retrieve ready lists (`get_readytask(i, del)`).
 *
 * The actual synchronization and distribution are performed by the Bank.
 */
class TaskManager {
public:
    /**
     * @brief Construct and initialize a task context with @p ntasks tasks.
     * @param ntasks  Total number of tasks available (0..ntasks-1).
     * @param iclient Client rank that opens the context.
     * @param comm    Communicator for Bank interaction.
     */
    TaskManager(int ntasks, int iclient = wrk_rank, MPI_Comm comm = comm_wrk);

    /** @brief Destructor; closes the TaskManager context. */
    ~TaskManager();

    /**
     * @brief Obtain the next task index to process.
     * @return Task index in [0, @ref n_tasks), or negative if none available.
     */
    int next_task();

    /**
     * @brief Mark item (@p i, @p j) as ready.
     * @param i Primary key (e.g., task group/channel).
     * @param j Secondary key (e.g., item id).
     */
    void put_readytask(int i, int j);

    /**
     * @brief Remove ready marker (@p i, @p j).
     * @param i Primary key.
     * @param j Secondary key.
     */
    void del_readytask(int i, int j);

    /**
     * @brief Retrieve the ready list for key @p i.
     * @param i   Primary key (queue id).
     * @param del If nonzero, consume (erase) the ready list; otherwise keep it.
     * @return Vector of secondary keys (j values) that are ready.
     */
    std::vector<int> get_readytask(int i, int del);

    /** @brief Bank account id associated with this task context. */
    int account_id = -1;

    /** @name Serial fallbacks
     *  These are used if the runtime is not using MPI distribution.
     *  @{ */
    int task = 0;    ///< Current task pointer (serial mode only).
    int n_tasks = 0; ///< Total tasks (serial mode only).
    /** @} */
};

/**
 * @brief Fixed size of control messages exchanged with the Bank.
 *
 * @details This constant is used by the MPI layer to size control payloads.
 */
int const message_size = 7;

} // namespace mrcpp
